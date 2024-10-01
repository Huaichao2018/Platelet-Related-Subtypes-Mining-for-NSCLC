ConsensusClusterPlus_km <- function(genelist = NULL, od = NULL, exp = NULL, clin = NULL, timecol = "time",
                                    statuscol = "status", seed = 123456, maxK = 5, cluster_character = "Cluster", plotFormat = "pdf",
                                    color_fun = color_fun1, xlab = "days", input_distance = NULL, input_clusteralg = NULL) {
    if (!dir.exists(od)) {
        dir.create(od)
    }
    # 过滤生存信息不全的样本
    colnames(clin)[colnames(clin) == timecol] <- "time"
    colnames(clin)[colnames(clin) == statuscol] <- "status"
    clin <- clin %>%
        mutate(time = as.numeric(time), status = as.numeric(status)) %>%
        filter(time > 0 & status != "" & status != "NA")
    # 过滤样本
    dat_exp <- exp[which(rownames(exp) %in% genelist), which(colnames(exp) %in% clin$sample)]
    clin <- clin[match(colnames(dat_exp), clin$sample), ]
    infor <- cbind.data.frame(clin, t(dat_exp))
    suppressPackageStartupMessages(library(ConsensusClusterPlus))
    suppressPackageStartupMessages(library(survival))
    suppressPackageStartupMessages(library(survminer))
    res_list <- list()

    if (is.null(input_distance)) {
        distances <- c("euclidean", "pearson", "spearman")
    } else {
        distances <- input_distance
    }
    if (is.null(input_clusteralg)) {
        clusteralgs <- c("km", "pam", "hc")
    } else {
        clusteralgs <- input_clusteralg
    }
    distalgs <- expand.grid(clusteralgs, distances, stringsAsFactors = FALSE)
    suppressPackageStartupMessages(library(magrittr))
    distalgs %<>% dplyr::filter(!(Var1 == "km" & Var2 != "euclidean"))
    clusteralg <- dplyr::pull(distalgs, 1)
    distance <- dplyr::pull(distalgs, 2)

    for (i in 1:length(distance)) {
        j <- i
        pic_km <- list()
        kmun <- 1
        # 无监督聚类
        conClust <- ConsensusClusterPlus(
            as.matrix(dat_exp),
            maxK = maxK,
            reps = 1000,
            pItem = 0.8,
            pFeature = 1,
            clusterAlg = clusteralg[j], # hc,pam,km
            distance = distance[i], # pearson,spearman,euclidean,binary,maximum,canberra,minkowski
            innerLinkage = "ward.D2",
            seed = seed,
            plot = plotFormat,
            title = paste0(od, "/", clusteralg[j], "_", distance[i]),
            writeTable = FALSE
        )
        # 计算2-maxK中最优分组的k值
        bestk <- getOptK(conClust = conClust, maxCls = maxK)
        # 循环输出最优分组、2-maxK其它分组的预后曲线
        for (k in c(bestk, setdiff(2:maxK, bestk))) {
            cluster <- data.frame(
                sample = as.character(names(conClust[[k]]$consensusClass)),
                Cluster = str_c(cluster_character, conClust[[k]]$consensusClass)
            )
            dat <- data.frame(
                time = as.numeric(infor$time),
                status = as.numeric(infor$status),
                Cluster = cluster$Cluster
            ) %>% arrange(`Cluster`)
            kmplot <- survfit(Surv(time, status) ~ Cluster, data = dat)
            surv_diff <- survdiff(Surv(time, status) ~ Cluster, data = dat)
            p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
            if (length(unique(dat$Cluster)) > 3) {
                fontsize <- 2
            } else {
                fontsize <- 4
            }
            p <- ggsurvplot(kmplot,
                data = dat, conf.int = FALSE, pval = TRUE, conf.int.style = "step",
                risk.table = "absolute",
                pval.size = 5, palette = color_fun,
                legend.title = cluster_character,
                legend.labs = unique(dat$Cluster),
                fontsize = fontsize,
                risk.table.y.text = FALSE, ncensor.plot = FALSE,
                xlab = xlab
            )
            if (length(unique(dat$Cluster)) > 3) {
                p <- p
            }
            pic_km[[kmun]] <- p
            kmun <- kmun + 1
            res_list[[distance[i]]][[clusteralg[j]]][[paste0("k", k)]][["p.val"]] <- p.val
            colnames(cluster)[2] <- cluster_character
            res_list[[distance[i]]][[clusteralg[j]]][[paste0("k", k)]][["cluster"]] <- cluster
            res_list[[distance[i]]][[clusteralg[j]]][[paste0("k", k)]][["best"]] <- paste0("k", bestk)
            res_list[[distance[i]]][[clusteralg[j]]][[paste0("k", k)]][["kmplot"]] <- p
        }
        pdf(file = paste0(od, "/", clusteralg[j], "_", distance[i], "/km.pdf"), width = 4.5, height = 4.5)
        print(pic_km, newPage = F)
        dev.off()
    }
    # 结果整理
    res_select <- map_dfr(names(res_list), function(distance) {
        map_dfr(names(res_list[[distance]]), function(clusteralg) {
            map_dfr(names(res_list[[distance]][[clusteralg]]), function(k) {
                tmp <- cbind.data.frame(distance = distance, clusteralg = clusteralg, k = k, p.val = res_list[[distance]][[clusteralg]][[k]][["p.val"]], bestk = res_list[[distance]][[clusteralg]][[k]][["best"]])
                return(tmp)
            })
        })
    })
    return(list(res_list = res_list, res_select = res_select))
}

getOptK <- function(conClust = NULL, minCls = 2, maxCls = 5) {
    # 最佳分类数
    Kvec <- minCls:maxCls
    x1 <- 0.1
    x2 <- 0.9 # threshold defining the intermediate sub-interval
    PAC <- rep(NA, length(Kvec))
    names(PAC) <- paste("K=", Kvec, sep = "") # from 2 to maxK
    for (i in Kvec) {
        M <- conClust[[i]][["consensusMatrix"]]
        Fn <- ecdf(M[lower.tri(M)])
        PAC[i - 1] <- Fn(x2) - Fn(x1)
    }
    optK <- Kvec[which.min(PAC)]
    return(optK)
}

signature_cox <-
    function(
        signaturelist = NULL, exp = NULL, clin = NULL, coxp = 0.05,
        bygroup = TRUE, timecol = "time", statuscol = "status", color_fun = color_fun1,
        xlab = "days", savekmplot = TRUE, best_cut = FALSE, minprop = 0.2) {
        library(survival)
        library(survminer)
        uselist <- intersect(rownames(exp), signaturelist)
        message(paste0("num of signature used for cox is ", length(uselist)))
        colnames(clin)[colnames(clin) == timecol] <- "time"
        colnames(clin)[colnames(clin) == statuscol] <- "status"
        clin <- clin %>%
            mutate(time = as.numeric(time), status = as.numeric(status)) %>%
            dplyr::filter(time > 0 & status != "" & status != "NA")
        infor <- exp[uselist, , drop = FALSE] %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column(., "sample") %>%
            merge(., clin %>%
                dplyr::select(sample, time, status), by = "sample")
        coxResult <- sapply(uselist, function(x) {
            dat <- cbind.data.frame(
                time = infor$time, status = infor$status,
                Value = infor[, x]
            )
            if (bygroup) {
                if (best_cut) {
                    sur.cut <- surv_cutpoint(
                        data = dat, time = "time",
                        event = "status", variables = "Value", minprop = minprop
                    )
                    cut <- summary(sur.cut)$cutpoint
                } else {
                    cut <- median(dat[, 3])
                }
                dat[, 3] <- factor(ifelse(dat[, 3] > cut, "High",
                    "Low"
                ), levels = c("Low", "High"))
            } else {
                cut <- NA
            }
            cox <- as.formula("Surv(time,status) ~ Value") %>% coxph(data = dat)
            pvalue <- summary(cox)$wald["pvalue"]
            HR <- summary(cox)$coef[, 2]
            lower95 <- summary(cox)$conf.int[, "lower .95"]
            upper95 <- summary(cox)$conf.int[, "upper .95"]
            res <- c(round(pvalue, 5), round(HR, 5), round(
                lower95,
                5
            ), round(upper95, 5), cut)
            return(res)
        }) %>%
            t() %>%
            as.data.frame()
        colnames(coxResult) <- c(
            "pvalue", "HR", "Low 95%CI", "High 95%CI",
            "cutoff"
        )
        sigcoxResult <- coxResult %>%
            dplyr::arrange(pvalue) %>%
            dplyr::filter(pvalue < coxp)
        cox_signature <- rownames(sigcoxResult)
        message(c("cox signature num is ", length(cox_signature)))
        if (savekmplot) {
            kmplot <- lapply(cox_signature, function(x) {
                dat <- cbind.data.frame(
                    time = infor$time, status = infor$status,
                    Group = infor[, x]
                )
                if (best_cut) {
                    sur.cut <- surv_cutpoint(
                        data = dat, time = "time",
                        event = "status", variables = "Group", minprop = minprop
                    )
                    cut <- summary(sur.cut)$cutpoint
                } else {
                    cut <- median(dat[, 3])
                }
                dat[, 3] <- factor(ifelse(dat[, 3] > cut, "High",
                    "Low"
                ), levels = c("Low", "High"))
                data_plot <- survfit(Surv(time, status) ~ Group,
                    data = dat
                )
                data_survdiff <- survdiff(Surv(time, status) ~ Group,
                    data = dat
                )
                pvalue <- 1 - pchisq(data_survdiff$chisq, length(data_survdiff$n) -
                    1)
                HR <- (data_survdiff$obs[2] / data_survdiff$exp[2]) / (data_survdiff$obs[1] / data_survdiff$exp[1])
                lower95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data_survdiff$exp[1] +
                    1 / data_survdiff$exp[2]))
                upper95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data_survdiff$exp[1] +
                    1 / data_survdiff$exp[2]))
                if (pvalue < 0.001) {
                    pval <- paste0(
                        "p < 0.001", "\n", "HR = ", round(
                            HR,
                            2
                        ), "\n", "95%CI (", round(lower95, 2), "-",
                        round(upper95, 2), ")"
                    )
                } else {
                    pval <- paste0(
                        "p = ", round(pvalue, 4), "\n",
                        "HR = ", round(HR, 2), "\n", "95%CI (", round(
                            lower95,
                            2
                        ), "-", round(upper95, 2), ")"
                    )
                }
                ggsurvplot(data_plot,
                    data = dat, conf.int = FALSE,
                    pval = pval, conf.int.style = "step", legend = "top",
                    legend.title = x, risk.table = "absolute", palette = color_fun1,
                    pval.size = 5, tables.theme = theme_cleantable(),
                    risk.table.y.text = FALSE, ncensor.plot = FALSE,
                    xlab = xlab
                ) + ggplot2::guides(color = guide_legend(nrow = 2))
            })
            names(kmplot) <- cox_signature
        } else {
            kmplot <- NULL
        }
        return(list(
            coxResult = coxResult, cox_signature = cox_signature,
            sigcoxResult = sigcoxResult, kmplot = kmplot
        ))
    }


lasso_model <-
    function(
        signaturelist = NULL, od = NULL, exp = NULL, clin = NULL,
        signaturetype = NULL, seed = 123456, timecol = "time", statuscol = "status",
        color_fun = color_fun1) {# 初始化 picture 变量为一个空的 list
      picture <<- list() 
      
        if (!dir.exists(od)) {
            dir.create(od)
        }
        message(sprintf("input signature num is %s", length(signaturelist)))
        message(sprintf("exist signature num is %s", sum(rownames(exp) %in%
            signaturelist)))
        library(glmnet)
        library(survival)
        library(ggpubr)
        colnames(clin)[colnames(clin) == timecol] <- "time"
        colnames(clin)[colnames(clin) == statuscol] <- "status"
        clin <- clin %>%
            mutate(time = as.numeric(time), status = as.numeric(status)) %>%
            filter(time > 0 & status != "" & status != "NA")
        infor <- exp %>%
            .[(rownames(exp) %in% signaturelist), ] %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column(., "sample") %>%
            merge(., clin, by = "sample")
        x <- as.matrix(infor[, signaturelist])
        rownames(x) <- infor$sample
        y <- data.matrix(Surv(infor$time, infor$status))
        set.seed(seed)
        fit <- glmnet(x, y, family = "cox")
        picture[[1]] <<- fit
        cv.fit <- cv.glmnet(x, y, family = "cox")
        picture[[2]] <<- cv.fit
        Coefficients <- coef(fit, s = cv.fit$lambda.min)
        Active.Index <- which(Coefficients != 0)
        Active.Coefficients <- Coefficients[Active.Index]
        lasso_signature <- row.names(Coefficients)[Active.Index]
        message(sprintf("after lasso signature num is %s", length(lasso_signature)))
        signature_coef <- data.frame(gene = lasso_signature, coef = Active.Coefficients)
        signature_coef$gene <- as.character(signature_coef$gene)
        signature_coef[, "group"] <- ifelse(signature_coef$coef >
            0, "pos", "neg")
        write.table(signature_coef,
            file = paste0(
                od, "/SupplementaryTable_",
                signaturetype, "_lasso_Coefficients.txt"
            ), quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = TRUE
        )
        coef_bar <- ggplot(data = signature_coef, mapping = aes(
            x = coef,
            y = factor(reorder(gene, coef)), fill = factor(group)
        )) +
            geom_bar(
                stat = "identity", color = "black", width = 0.8,
                linewidth = 0.6
            ) +
            theme_pubr() +
            guides(fill = "none") +
            ylab("") +
            xlab("Coeffecients") +
            scale_fill_manual(values = color_fun)
        p <- cowplot::plot_grid(~ plot(picture[[1]], xvar = "lambda"),
            ~ plot(picture[[2]]), coef_bar,
            labels = NA, ncol = 3,
            label_size = 24, scale = c(1, 1, 0.9)
        )
        plotout(p = p, h = 6, w = 18, od = od, num = paste0(
            "_",
            signaturetype, "_lasso"
        ))
        return(signature_coef)
    }



modelscore_km_roc <-
    function(
        signature_coef = NULL, od = NULL, exp = NULL, clin = NULL,
        best_cut = FALSE, dataset = NULL, roc_time = c(1, 3, 5),
        timecol = "time", statuscol = "status", color_fun = color_fun1,
        surtime_unit = 365, no_roc_pheatmap = FALSE, xlab = "days",
        minprop = 0.2) {
        if (!dir.exists(od)) {
            dir.create(od)
        }
        colnames(clin)[colnames(clin) == timecol] <- "time"
        colnames(clin)[colnames(clin) == statuscol] <- "status"
        clin <- clin %>%
            mutate(time = as.numeric(time), status = as.numeric(status)) %>%
            filter(time > 0 & status != "" & status != "NA")
        exists_sign <- intersect(signature_coef$signature, rownames(exp))
        message(sprintf(
            "input signature num is %s,exists signature num is %s",
            length(signature_coef$signature), length(exists_sign)
        ))
        if (length(exists_sign) > 0) {
            if (length(setdiff(signature_coef$signature, exists_sign)) >
                0) {
                exp <- exp %>% as.data.frame()
                for (i in setdiff(signature_coef$signature, exists_sign)) {
                    exp[i, ] <- 0
                }
            }
            library(survivalROC)
            library(survival)
            library(survminer)
            infor <- exp %>%
                t() %>%
                as.data.frame() %>%
                rownames_to_column(
                    .,
                    "sample"
                ) %>%
                merge(., clin, by = "sample")
            rownames(infor) <- infor$sample
            sigdata <- as.matrix(subset(infor, select = signature_coef$signature))
            Score <- sigdata %*% signature_coef$coef
            colnames(Score) <- "Score"
            write.table(Score,
                file = paste0(
                    od, "/SupplementaryTable_",
                    dataset, "_Score.txt"
                ), quote = FALSE, sep = "\t",
                col.names = TRUE, row.names = T
            )
            dat <- cbind.data.frame(
                time = infor$time, status = infor$status,
                Score = Score
            )
            if (best_cut) {
                sur.cut <- surv_cutpoint(
                    data = dat, time = "time",
                    event = "status", variables = "Score", minprop = minprop
                )
                cut <- summary(sur.cut)$cutpoint
            } else {
                cut <- median(dat[, 3])
            }
            Group <- factor(ifelse(dat[, 3] > cut, "High", "Low"),
                levels = c("Low", "High")
            ) %>% as.data.frame()
            if (length(unique(Group[, 1])) > 1 & sum(summary(Group[
                ,
                1
            ]) > 5) == length(unique(Group[, 1]))) {
                colnames(Group) <- "Group"
                rownames(Group) <- rownames(Score)
                write.table(Group,
                    file = paste0(
                        od, "/SupplementaryTable_",
                        dataset, "_Group.txt"
                    ), quote = FALSE, sep = "\t",
                    col.names = TRUE, row.names = T
                )
                dat <- cbind.data.frame(
                    time = infor$time, status = infor$status,
                    Score = Group$Group
                )
                data_plot <- survfit(Surv(time, status) ~ Score,
                    data = dat
                )
                data_survdiff <- survdiff(Surv(time, status) ~ Score,
                    data = dat
                )
                pvalue <- 1 - pchisq(data_survdiff$chisq, length(data_survdiff$n) -
                    1)
                HR <- (data_survdiff$obs[2] / data_survdiff$exp[2]) / (data_survdiff$obs[1] / data_survdiff$exp[1])
                lower95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data_survdiff$exp[1] +
                    1 / data_survdiff$exp[2]))
                upper95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data_survdiff$exp[1] +
                    1 / data_survdiff$exp[2]))
                if (pvalue < 0.001) {
                    pval <- paste0(
                        "p < 0.001", "\n", "HR = ", round(
                            HR,
                            2
                        ), "\n", "95%CI (", round(lower95, 2), "-",
                        round(upper95, 2), ")"
                    )
                } else {
                    pval <- paste0(
                        "p = ", round(pvalue, 4), "\n",
                        "HR = ", round(HR, 2), "\n", "95%CI (", round(
                            lower95,
                            2
                        ), "-", round(upper95, 2), ")"
                    )
                }
                res <- ggsurvplot(data_plot,
                    data = dat, conf.int = FALSE,
                    pval = pval, conf.int.style = "step", legend = "top",
                    legend.title = dataset, risk.table = "absolute",
                    palette = color_fun, pval.size = 5, risk.table.y.text = FALSE,
                    ncensor.plot = FALSE, xlab = xlab
                )
                p1 <- ggarrange(res$plot, res$table,
                    ncol = 1, align = "v",
                    heights = c(0.75, 0.3)
                )
                dat <- cbind.data.frame(
                    time = infor$time, status = infor$status,
                    Score = Score
                )
                if (!is.na(HR)) {
                    if (HR < 1) {
                        dat <- dat %>% mutate(status = ifelse(status >
                            0, 0, 1))
                    }
                    library(timeROC)
                    roc_model <- timeROC(
                        T = dat$time, delta = dat$status,
                        marker = dat$Score, cause = 1, weighting = "marginal",
                        times = c(roc_time[1] * surtime_unit, roc_time[2] *
                            surtime_unit, roc_time[3] * surtime_unit),
                        ROC = TRUE, iid = FALSE
                    )
                    roc_plot <- data.frame(year1x = roc_model$FP[
                        ,
                        1
                    ], year1y = roc_model$TP[, 1], year2x = roc_model$FP[
                        ,
                        2
                    ], year2y = roc_model$TP[, 2], year3x = roc_model$FP[
                        ,
                        3
                    ], year3y = roc_model$TP[, 3])
                    AUC_anno_1 <- sprintf(
                        "AUC at %s year = %s",
                        roc_time[1], sprintf("%.3f", roc_model$AUC[[1]])
                    )
                    AUC_anno_2 <- sprintf(
                        "AUC at %s year = %s",
                        roc_time[2], sprintf("%.3f", roc_model$AUC[[2]])
                    )
                    AUC_anno_3 <- sprintf(
                        "AUC at %s year = %s",
                        roc_time[3], sprintf("%.3f", roc_model$AUC[[3]])
                    )
                    p2 <- ggplot(data = roc_plot) +
                      geom_line(aes(
                        x = year1x,
                        y = year1y
                      ), linewidth = 1.2, color = "#E31A1C") +
                      geom_line(aes(x = year2x, y = year2y),
                                linewidth = 1.2,
                                color = "#377EB8"
                      ) +
                      geom_line(aes(
                        x = year3x,
                        y = year3y
                      ), linewidth = 1.2, color = "#007947") +
                      geom_abline(
                        slope = 1, intercept = 0, color = "grey",
                        linewidth = 1, linetype = 2
                      )+
                        egg::theme_article() +
                        theme(plot.background = element_rect(fill = "white")) +
                        annotate(
                            geom = "line", x = c(0.5, 0.54), y = 0.17,
                            colour = "#E31A1C", size = 1.2
                        ) +
                        annotate("text",
                            x = 0.55, y = 0.17, size = 5, label = AUC_anno_1,
                            color = "black", hjust = "left"
                        ) +
                        annotate(
                            geom = "line",
                            x = c(0.5, 0.54), y = 0.11, colour = "#377EB8",
                            size = 1.2
                        ) +
                        annotate("text",
                            x = 0.55, y = 0.11,
                            size = 5, label = AUC_anno_2, color = "black",
                            hjust = "left"
                        ) +
                        annotate(geom = "line", x = c(
                            0.5,
                            0.54
                        ), y = 0.05, colour = "#007947", size = 1.2) +
                        annotate("text",
                            x = 0.55, y = 0.05, size = 5,
                            label = AUC_anno_3, color = "black", hjust = "left"
                        ) +
                        labs(x = "1-Specificity", y = "Sensitivity") +
                        theme(axis.text.x = element_text(
                            face = "plain",
                            size = 12, color = "black"
                        ), axis.text.y = element_text(
                            face = "plain",
                            size = 12, color = "black"
                        ), axis.title.x = element_text(
                            face = "plain",
                            size = 14, color = "black"
                        ), axis.title.y = element_text(
                            face = "plain",
                            size = 14, color = "black"
                        ))
                    dat <- cbind.data.frame(
                        time = infor$time, status = infor$status,
                        Group = Group, Score = Score
                    ) %>%
                        arrange(Score) %>%
                        mutate(x = 1:nrow(dat), status = ifelse(status ==
                            0, "Alive", "Death"))
                    v <- as.numeric(table(dat[, 3])["Low"]) + 0.5
                    h <- cut
                    p3 <- ggplot(data = dat) +
                      geom_point(mapping = aes(
                        x = x,
                        y = Score, color = Group
                      )) +
                      theme_bw() +
                      theme(
                        panel.grid = element_blank(),  # 添加 panel.grid 的设置
                        legend.position = c(0.1, 0.85)  # 使用 c(x, y) 指定图例位置，而不是 "inside"
                      ) +
                      xlab("") +
                      scale_color_manual(values = color_fun, name = NULL) +
                      geom_hline(aes(yintercept = h),
                                 colour = "#000000",
                                 linetype = "dashed"
                      ) +
                      geom_vline(aes(xintercept = v),
                                 colour = "#000000", linetype = "dashed"
                      )
                    
                    p4 <- ggplot(data = dat) +
                        geom_point(mapping = aes(
                            x = x,
                            y = time, color = status
                        )) +
                        theme_bw() +
                        theme(
                            panel.grid = element_blank(),
                            legend.position = c(0.1, 0.85)
                        ) +
                        xlab("") +
                        ylab(xlab) +
                        scale_color_manual(
                            values = color_fun,
                            name = NULL
                        ) +
                        geom_vline(aes(xintercept = v),
                            colour = "#000000", linetype = "dashed"
                        )
                    library(pheatmap)
                    library(ggplotify)
                    annotation_col <- as.data.frame(Group %>% dplyr::rename(Score = Group))
                    ann_colors <- list(Score = c(
                        High = color_fun[2],
                        Low = color_fun[1]
                    ))
                    col <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
                    plotdata <- sigdata[order(Score), ] %>%
                        as.matrix() %>%
                        t()
                    rownames(plotdata) <- colnames(sigdata)
                    p <- pheatmap::pheatmap(plotdata,
                        border = F,
                        scale = "row", cluster_rows = F, cluster_cols = F,
                        show_colnames = F, fontsize = 8, legend = TRUE,
                        annotation_legend = FALSE, annotation_col = annotation_col,
                        annotation_colors = ann_colors, color = col
                    )
                    dev.off()
                    p5 <- ggplotify::as.ggplot(p)
                    g1 <- cowplot::plot_grid(p1, p2, ncol = 1, scale = c(
                        0.95,
                        0.95
                    ), labels = c("AUTO"), label_size = 20)
                    g2 <- cowplot::plot_grid(p3 + theme(plot.margin = unit(c(
                        1,
                        5, 1, 1
                    ), "lines")), p4 + theme(plot.margin = unit(c(
                        1,
                        5, 1, 1
                    ), "lines")), p5 + theme(plot.margin = unit(c(
                        1,
                        1.4, 1, 4.5
                    ), "lines")), ncol = 1, labels = c(
                        "C",
                        "D", "E"
                    ), label_size = 20)
                    p <- cowplot::plot_grid(g1, g2,
                        ncol = 2, labels = NA,
                        rel_widths = c(6, 8)
                    )
                    plotout(p = p, od = od, num = sprintf(
                        "_%s_model",
                        dataset
                    ), w = 14, h = 12)
                    if (no_roc_pheatmap) {
                        g1 <- p1
                        g2 <- cowplot::plot_grid(p3 + theme(plot.margin = unit(c(
                            1,
                            5, 1, 1
                        ), "lines")), p4 + theme(plot.margin = unit(c(
                            1,
                            5, 1, 1
                        ), "lines")), ncol = 1, labels = c(
                            "B",
                            "C"
                        ), label_size = 20)
                        p <- cowplot::plot_grid(g1, g2, ncol = 2, labels = c(
                            "A",
                            NA
                        ), scale = c(0.95, 1), label_size = 20)
                        plotout(p = p, od = od, num = sprintf(
                            "_%s_model_noroc_heat",
                            dataset
                        ), w = 12, h = 7)
                    }
                    stat_infor <- data.frame(
                        dataset = dataset, exists_sign_num = length(exists_sign),
                        km_pval = pvalue, HR = HR, auc_1 = roc_model$AUC[[1]],
                        auc_2 = roc_model$AUC[[2]], auc_3 = roc_model$AUC[[3]],
                        best_cut = best_cut, cut_value = cut
                    )
                    return(list(Score = Score, Group = Group, infor = stat_infor))
                } else {
                    return(NULL)
                }
            }
        } else {
            return(NULL)
        }
    }

immunotherapy <-
    function(
        signature_coef = NULL, od = getwd(), exp = NULL, clin = NULL,
        best_cut = FALSE, dataset = NULL, timecol = "time", statuscol = "status",
        color_fun = color_fun1, immunotherapy.responsecol = c(
            "Immunotherapy.Response",
            "Immunotherapy.Response.Group"
        ), km_xlab = "days", SinglePlotWidth = 3.5,
        SinglePlotHigh = 3.8) {
        if (!dir.exists(od)) {
            dir.create(od)
        }
        signature_coef <- as.data.frame(signature_coef)
        index_tmp <- map_lgl(1:ncol(signature_coef), function(x) {
            a <- signature_coef[, x]
            is.numeric(a)
        })
        colnames(signature_coef)[which(index_tmp)] <- "coef"
        colnames(signature_coef)[which(!index_tmp)] <- "signature"
        if (any(colnames(clin) %in% "status")) {
            colnames(clin)[colnames(clin) == "status"] <- "raw_status"
        }
        if (any(colnames(clin) %in% "time")) {
            colnames(clin)[colnames(clin) == "time"] <- "raw_time"
        }
        library(tidyverse)
        library(survival)
        library(survminer)
        library(ggpubr)
        common_sample <- intersect(clin$sample, colnames(exp))
        if (length(common_sample) > 10) {
            message(dataset)
            exp <- exp %>%
                as.data.frame() %>%
                dplyr::select(any_of(common_sample))
            clin <- clin %>% filter(sample %in% common_sample)
            exists_sign <- intersect(signature_coef$signature, rownames(exp))
            message(sprintf(
                "input signature num is %s,exists signature num is %s",
                crayon::bold(length(signature_coef$signature)), crayon::bold(length(exists_sign))
            ))
            if (length(exists_sign) > 0) {
                if (length(setdiff(signature_coef$signature, exists_sign)) >
                    0) {
                    exp <- exp %>% as.data.frame()
                    for (i in setdiff(signature_coef$signature, exists_sign)) {
                        exp[i, ] <- 0
                    }
                }
                sigdata <- as.matrix(subset(exp %>% t() %>% as.data.frame(),
                    select = signature_coef$signature
                ))
                Score <- sigdata %*% signature_coef$coef
                colnames(Score) <- "Score"
                write.table(Score,
                    file = paste0(
                        od, "/SupplementaryTable_",
                        dataset, "_Score.txt"
                    ), quote = FALSE, sep = "\t",
                    col.names = TRUE, row.names = T
                )
            }
            res <- lapply(immunotherapy.responsecol, function(x) {
                if (any(colnames(clin) %in% x)) {
                    if (length(unique(clin[, x])) > 1) {
                        use_clin <- clin %>%
                            mutate(response.col = .[
                                ,
                                x
                            ]) %>%
                            dplyr::filter(!is.na(response.col)) %>%
                            merge(., Score %>% as.data.frame() %>% rownames_to_column(
                                .,
                                "sample"
                            )) %>%
                            arrange(response.col)
                        table(use_clin$response.col)
                        p_all_sample_boxplot <- ggboxplot(
                            data = use_clin,
                            x = "response.col", y = "Score", fill = "response.col"
                        ) +
                            stat_compare_means(na.rm = TRUE) + xlab("Response") +
                            ylab("Score") + theme_bw() + theme(
                                panel.grid = element_blank(),
                                legend.title = element_blank(), legend.position = "top",
                                axis.text.x = element_text(
                                    angle = 0, hjust = 0.5,
                                    vjust = 0
                                ), axis.title.x.bottom = element_blank()
                            ) +
                            scale_fill_manual(values = color_fun[-c(
                                1,
                                2
                            )])
                        plotout(
                            od = od, name = str_glue("{dataset}_{x}_all_sample_boxplot"),
                            w = SinglePlotWidth, h = SinglePlotHigh,
                            p = p_all_sample_boxplot, plot_tif = FALSE
                        )
                        if (any(colnames(clin) %in% timecol) & any(colnames(clin) %in%
                            statuscol)) {
                            colnames(use_clin)[colnames(use_clin) ==
                                timecol] <- "time"
                            colnames(use_clin)[colnames(use_clin) ==
                                statuscol] <- "status"
                            infor <- use_clin %>%
                                mutate(
                                    time = as.numeric(time),
                                    status = as.numeric(status)
                                ) %>%
                                filter(time >
                                    0 & status != "" & status != "NA")
                            rownames(infor) <- infor$sample
                            dat <- cbind.data.frame(
                                sample = infor$sample,
                                time = infor$time, status = infor$status,
                                Score = infor$Score
                            )
                            if (best_cut) {
                                sur.cut <- surv_cutpoint(
                                    data = dat, time = "time",
                                    event = "status", variables = "Score"
                                )
                                cut <- summary(sur.cut)$cutpoint
                            } else {
                                cut <- median(dat$Score)
                            }
                            Group <- factor(ifelse(dat$Score > cut, "High",
                                "Low"
                            ), levels = c("Low", "High")) %>%
                                as.data.frame()
                            if (length(unique(Group[, 1])) > 1 & (sum(summary(Group[
                                ,
                                1
                            ]) > 5) == length(unique(Group[, 1])))) {
                                colnames(Group) <- "Group"
                                rownames(Group) <- infor$sample
                                write.table(Group,
                                    file = paste0(
                                        od, "/SupplementaryTable_",
                                        dataset, "_sample_in_", x, "_Group.txt"
                                    ),
                                    quote = FALSE, sep = "\t", col.names = TRUE,
                                    row.names = T
                                )
                                dat <- cbind.data.frame(
                                    time = infor$time,
                                    status = infor$status, Score = Group$Group
                                )
                                data_plot <- survfit(Surv(time, status) ~
                                    Score, data = dat)
                                data_survdiff <- survdiff(Surv(time, status) ~
                                    Score, data = dat)
                                pvalue <- 1 - pchisq(
                                    data_survdiff$chisq,
                                    length(data_survdiff$n) - 1
                                )
                                HR <- (data_survdiff$obs[2] / data_survdiff$exp[2]) / (data_survdiff$obs[1] / data_survdiff$exp[1])
                                lower95 <- exp(log(HR) - qnorm(0.975) *
                                    sqrt(1 / data_survdiff$exp[1] + 1 / data_survdiff$exp[2]))
                                upper95 <- exp(log(HR) + qnorm(0.975) *
                                    sqrt(1 / data_survdiff$exp[1] + 1 / data_survdiff$exp[2]))
                                if (pvalue < 0.001) {
                                    pval <- paste0(
                                        "p < 0.001", "\n", "HR = ",
                                        round(HR, 2), "\n", "95%CI (", round(
                                            lower95,
                                            2
                                        ), "-", round(upper95, 2), ")"
                                    )
                                } else {
                                    pval <- paste0(
                                        "p = ", round(
                                            pvalue,
                                            4
                                        ), "\n", "HR = ", round(HR, 2), "\n",
                                        "95%CI (", round(lower95, 2), "-",
                                        round(upper95, 2), ")"
                                    )
                                }
                                res <- ggsurvplot(data_plot,
                                    data = dat,
                                    conf.int = FALSE, pval = pval, conf.int.style = "step",
                                    legend = "top", legend.title = dataset,
                                    surv.median.line = "hv", risk.table = "absolute",
                                    palette = color_fun, pval.size = 5, risk.table.y.text = FALSE,
                                    ncensor.plot = FALSE, xlab = km_xlab
                                )
                                p_group_km <- ggarrange(res$plot + theme(axis.title.x.bottom = element_blank()),
                                    res$table,
                                    ncol = 1, align = "v", heights = c(
                                        0.7,
                                        0.3
                                    )
                                )
                                plotout(
                                    p = p_group_km, od = od, name = str_glue("{dataset}_sample_in_{x}_group_km"),
                                    w = SinglePlotWidth + 0.7, h = SinglePlotHigh +
                                        0.5, plot_tif = FALSE
                                )
                                table(infor$response.col)
                                p_boxplot <- ggboxplot(
                                    data = infor, x = "response.col",
                                    y = "Score", fill = "response.col"
                                ) +
                                    stat_compare_means(na.rm = TRUE) + xlab("Response") +
                                    ylab("Score") + theme_bw() + theme(
                                        panel.grid = element_blank(),
                                        legend.title = element_blank(), legend.position = "top",
                                        axis.text.x = element_text(
                                            angle = 0,
                                            hjust = 0.5, vjust = 0
                                        ), axis.title.x.bottom = element_blank()
                                    ) +
                                    scale_fill_manual(values = color_fun[-c(
                                        1,
                                        2
                                    )])
                                plotout(
                                    od = od, name = str_glue("{dataset}_sample_in_{statuscol}_{x}_boxplot"),
                                    w = SinglePlotWidth, h = SinglePlotHigh,
                                    p = p_boxplot, plot_tif = FALSE
                                )
                                use_data <- Group %>%
                                    rownames_to_column("sample") %>%
                                    merge(., infor) %>%
                                    select(
                                        sample, Group,
                                        response.col
                                    )
                                p.value <- use_data %>%
                                    select(Group, response.col) %>%
                                    table() %>%
                                    fisher.test(simulate.p.value = TRUE) %>%
                                    .$p.value %>%
                                    format.pval(.,
                                        digits = 2,
                                        eps = 0.001
                                    )
                                if (str_detect(string = p.value, pattern = "<")) {
                                    p.value <- paste0("Fisher.p ", p.value)
                                } else {
                                    p.value <- paste0("Fisher.p = ", p.value)
                                }
                                data1 <- Group %>%
                                    rownames_to_column("sample") %>%
                                    merge(., infor) %>%
                                    select(
                                        sample, response.col,
                                        Group
                                    ) %>%
                                    group_by(Group, response.col) %>%
                                    summarise(n = n()) %>%
                                    mutate(Percent = n / sum(n) *
                                        100) %>%
                                    mutate(txt = str_c(sprintf(
                                        "%.2f",
                                        Percent
                                    )))
                                p_sample_percentage <- ggplot(
                                    data = data1,
                                    aes(Group, Percent, fill = response.col)
                                ) +
                                    geom_col(width = 0.7) +
                                    geom_text(aes(label = txt),
                                        position = position_stack(vjust = 0.5),
                                        size = 4, color = "white"
                                    ) +
                                    scale_fill_manual(values = color_fun[-c(1:2)]) +
                                    annotate("text",
                                        label = p.value, color = "black",
                                        y = 104, x = 1.15, size = 3.7
                                    ) +
                                    theme_test(13) +
                                    theme(legend.position = "top") +
                                    labs(
                                        fill = "",
                                        y = "Percent(%)"
                                    ) +
                                    theme(
                                        axis.text.x.bottom = element_text(size = 11),
                                        axis.title.x.bottom = element_blank()
                                    )
                                plotout(
                                    od = od, name = str_glue("{dataset}_{statuscol}_{x}_sample_percentage"),
                                    w = SinglePlotWidth, h = SinglePlotHigh,
                                    p = p_sample_percentage, plot_tif = FALSE
                                )
                                data2 <- Group %>%
                                    rownames_to_column("sample") %>%
                                    merge(., infor) %>%
                                    select(
                                        sample, response.col,
                                        Group
                                    ) %>%
                                    group_by(response.col, Group) %>%
                                    summarise(n = n()) %>%
                                    mutate(Percent = n / sum(n) *
                                        100) %>%
                                    mutate(txt = str_c(sprintf(
                                        "%.2f",
                                        Percent
                                    )))
                                p_sample_percentage2 <- ggplot(
                                    data = data2,
                                    aes(response.col, Percent, fill = Group)
                                ) +
                                    geom_col(width = 0.7) +
                                    geom_text(aes(label = txt),
                                        position = position_stack(vjust = 0.5),
                                        size = 4, color = "white"
                                    ) +
                                    scale_fill_manual(values = color_fun) +
                                    annotate("text",
                                        label = p.value, color = "black",
                                        y = 104, x = 1.15, size = 3.7
                                    ) +
                                    theme_test(13) +
                                    theme(legend.position = "top") +
                                    labs(
                                        fill = "",
                                        y = "Percent(%)"
                                    ) +
                                    theme(
                                        axis.text.x.bottom = element_text(size = 11),
                                        axis.title.x.bottom = element_blank()
                                    )
                                plotout(
                                    od = od, name = str_glue("{dataset}_{statuscol}_{x}_sample_percentage2"),
                                    w = SinglePlotWidth, h = SinglePlotHigh,
                                    p = p_sample_percentage2, plot_tif = FALSE
                                )
                            } else {
                                return(NULL)
                            }
                        } else {
                            return(NULL)
                        }
                    } else {
                        return(NULL)
                    }
                }
                projects_chara <- c(
                    "p_group_km", "p_all_sample_boxplot",
                    "p_boxplot", "p_sample_percentage"
                )
                index <- map_lgl(projects_chara, function(y) {
                    exists(x = y)
                })
                projects_chara_exits <- projects_chara[which(index)]
                projects_chara_exits_project <- map(
                    projects_chara_exits,
                    function(x) {
                        get0(x)
                    }
                )
                names(projects_chara_exits_project) <- projects_chara_exits
                return(projects_chara_exits_project)
            })
            if (any(colnames(clin) %in% timecol) & any(colnames(clin) %in%
                statuscol)) {
                colnames(clin)[colnames(clin) == timecol] <- "time"
                colnames(clin)[colnames(clin) == statuscol] <- "status"
                clin <- clin %>%
                    mutate(
                        time = as.numeric(time),
                        status = as.numeric(status)
                    ) %>%
                    filter(time >
                        0 & status != "" & status != "NA")
                library(survival)
                library(survminer)
                infor <- Score %>%
                    as.data.frame() %>%
                    rownames_to_column(
                        .,
                        "sample"
                    ) %>%
                    merge(., clin, by = "sample")
                rownames(infor) <- infor$sample
                dat <- cbind.data.frame(
                    sample = infor$sample, time = infor$time,
                    status = infor$status, Score = infor$Score
                )
                if (best_cut) {
                    sur.cut <- surv_cutpoint(
                        data = dat, time = "time",
                        event = "status", variables = "Score"
                    )
                    cut <- summary(sur.cut)$cutpoint
                } else {
                    cut <- median(dat$Score)
                }
                Group <- factor(ifelse(dat$Score > cut, "High", "Low"),
                    levels = c("Low", "High")
                ) %>% as.data.frame()
                if (length(unique(Group[, 1])) > 1 & (sum(summary(Group[
                    ,
                    1
                ]) > 5) == length(unique(Group[, 1])))) {
                    colnames(Group) <- "Group"
                    rownames(Group) <- infor$sample
                    write.table(Group,
                        file = paste0(
                            od, "/SupplementaryTable_",
                            dataset, "all_sample_Group.txt"
                        ), quote = FALSE,
                        sep = "\t", col.names = TRUE, row.names = T
                    )
                    dat <- cbind.data.frame(
                        time = infor$time, status = infor$status,
                        Score = Group$Group
                    )
                    data_plot <- survfit(Surv(time, status) ~ Score,
                        data = dat
                    )
                    data_survdiff <- survdiff(Surv(time, status) ~
                        Score, data = dat)
                    pvalue <- 1 - pchisq(data_survdiff$chisq, length(data_survdiff$n) -
                        1)
                    HR <- (data_survdiff$obs[2] / data_survdiff$exp[2]) / (data_survdiff$obs[1] / data_survdiff$exp[1])
                    lower95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data_survdiff$exp[1] +
                        1 / data_survdiff$exp[2]))
                    upper95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data_survdiff$exp[1] +
                        1 / data_survdiff$exp[2]))
                    if (pvalue < 0.001) {
                        pval <- paste0(
                            "p < 0.001", "\n", "HR = ",
                            round(HR, 2), "\n", "95%CI (", round(
                                lower95,
                                2
                            ), "-", round(upper95, 2), ")"
                        )
                    } else {
                        pval <- paste0(
                            "p = ", round(pvalue, 4), "\n",
                            "HR = ", round(HR, 2), "\n", "95%CI (", round(
                                lower95,
                                2
                            ), "-", round(upper95, 2), ")"
                        )
                    }
                    res <- ggsurvplot(data_plot,
                        data = dat, conf.int = FALSE,
                        pval = pval, conf.int.style = "step", legend = "top",
                        legend.title = dataset, risk.table = "absolute",
                        palette = color_fun, pval.size = 5, surv.median.line = "hv",
                        risk.table.y.text = FALSE, ncensor.plot = FALSE,
                        xlab = km_xlab
                    )
                    p_all_sample_group_km <- ggarrange(
                        res$plot +
                            theme(axis.title.x.bottom = element_blank()),
                        res$table,
                        ncol = 1, align = "v", heights = c(
                            0.7,
                            0.3
                        )
                    )
                    plotout(
                        p = p_all_sample_group_km, od = od, name = str_glue("{dataset}_all_sample_group_km"),
                        w = SinglePlotWidth + 0.7, h = SinglePlotHigh +
                            0.5, plot_tif = FALSE
                    )
                } else {
                    return(NULL)
                }
            }
        } else {
            message(str_glue("{crayon::bold('WARNING:') } 表达谱中与临床信息表中样本交集小于10.分析跳过."))
        }
    }


limma_deg = 
function (od = NULL, DEG_exp = NULL, DEG_pdata = NULL, controlLabel = NULL, 
    caseLabel = NULL, DEG_FC = 1, DEG_P = 0.01, pvalue = NULL, 
    saveplot = FALSE, color_fun = NULL) 
{
    if (!dir.exists(od)) {
        dir.create(od)
    }
    library(limma)
    usesample <- intersect(colnames(DEG_exp), DEG_pdata[, 1])
    DEG_pdata <- DEG_pdata[match(usesample, DEG_pdata[, 1]), 
        ]
    DEG_exp <- DEG_exp[, usesample]
    design <- model.matrix(~0 + factor(DEG_pdata[, 2]))
    colnames(design) <- levels(factor(DEG_pdata[, 2]))
    rownames(design) <- colnames(DEG_exp)
    contrasts <- paste0(caseLabel, "-", controlLabel)
    contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)
    fit <- lmFit(DEG_exp, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2, 0.01)
    nrDEG <- topTable(fit2, adjust = "fdr", sort.by = "B", number = nrow(DEG_exp)) %>% 
        na.omit()
    if (is.null(pvalue)) {
        DEGs <- rownames(nrDEG)[which(abs(nrDEG$logFC) >= DEG_FC & 
            nrDEG$adj.P.Val < DEG_P)]
        up <- rownames(nrDEG)[which(nrDEG$logFC >= DEG_FC & nrDEG$adj.P.Val < 
            DEG_P)]
        down <- rownames(nrDEG)[which(nrDEG$logFC <= (-DEG_FC) & 
            nrDEG$adj.P.Val < DEG_P)]
        nrDEG[up, "Diff"] <- "up"
        nrDEG[down, "Diff"] <- "down"
        nrDEG[!(nrDEG$Diff %in% c("up", "down")), "Diff"] <- "not-significant"
        message(paste0(caseLabel, "_vs_", controlLabel))
        message(paste0("log2FC=", DEG_FC, ",", "FDR=", DEG_P))
        message("DEGs up down")
        message(paste0(length(DEGs), " ", length(up), " ", length(down)))
    }
    else {
        DEGs <- rownames(nrDEG)[which(abs(nrDEG$logFC) >= DEG_FC & 
            nrDEG$P.Value < DEG_P)]
        up <- rownames(nrDEG)[which(nrDEG$logFC >= DEG_FC & nrDEG$P.Value < 
            DEG_P)]
        down <- rownames(nrDEG)[which(nrDEG$logFC <= (-DEG_FC) & 
            nrDEG$P.Value < DEG_P)]
        nrDEG[up, "Diff"] <- "up"
        nrDEG[down, "Diff"] <- "down"
        nrDEG[!(nrDEG$Diff %in% c("up", "down")), "Diff"] <- "not-significant"
        message(paste0(caseLabel, "_vs_", controlLabel))
        message(paste0("log2FC=", DEG_FC, ",", "P.Value=", DEG_P))
        message("DEGs up down")
        message(paste0(length(DEGs), " ", length(up), " ", length(down)))
    }
    write.table(DEGs, file = paste0(od, "/SupplementaryTable_", 
        caseLabel, "_vs_", controlLabel, "_DEGs.txt"), quote = FALSE, 
        sep = "\t", col.names = FALSE, row.names = FALSE)
    write.table(nrDEG, file = paste0(od, "/SupplementaryTable_", 
        caseLabel, "_vs_", controlLabel, "_nrDEG.txt"), quote = FALSE, 
        sep = "\t", col.names = TRUE, row.names = TRUE)
    if (saveplot) {
        library(ggplot2)
        library(ggrepel)
        if (is.null(pvalue)) {
            p <- ggplot(nrDEG, aes(x = logFC, y = -log10(adj.P.Val))) + 
                geom_point(aes(color = Diff), alpha = 0.4, size = 1.5) + 
                theme_bw() + theme(panel.grid = element_blank(), 
                legend.title = element_blank(), plot.title = element_text(hjust = 0.5, 
                  size = 20)) + scale_colour_manual(limits = c("up", 
                "down", "not-significant"), values = c("#E31A1C", 
                "#007947", "gray40"), labels = c("up", "down", 
                "not-significant"))
            up <- nrDEG %>% filter(Diff == "up") %>% dplyr::arrange(adj.P.Val) %>% 
                head(10)
            down <- nrDEG %>% filter(Diff == "down") %>% dplyr::arrange(adj.P.Val) %>% 
                head(10)
            data <- rbind(up, down)
            data$gene <- rownames(data)
            p1 <- p + theme(legend.position = "right") + geom_text_repel(data = data, 
                aes(x = logFC, y = -log10(adj.P.Val), label = gene), 
                size = 3, box.padding = unit(0.5, "lines"), segment.color = "black", 
                show.legend = FALSE)
            plotout(p = p1, num = paste0("_", caseLabel, "_vs_", 
                controlLabel, "_volcano"), h = 6, w = 7, od = od)
        }
        else {
            p <- ggplot(nrDEG, aes(x = logFC, y = -log10(P.Value))) + 
                geom_point(aes(color = Diff), alpha = 0.4, size = 1.5) + 
                theme_bw() + theme(panel.grid = element_blank(), 
                legend.title = element_blank(), plot.title = element_text(hjust = 0.5, 
                  size = 20)) + scale_colour_manual(limits = c("up", 
                "down", "not-significant"), values = c("#E31A1C", 
                "#007947", "gray40"), labels = c("up", "down", 
                "not-significant"))
            up <- nrDEG %>% filter(Diff == "up") %>% dplyr::arrange(P.Value) %>% 
                head(10)
            down <- nrDEG %>% filter(Diff == "down") %>% dplyr::arrange(P.Value) %>% 
                head(10)
            data <- rbind(up, down)
            data$gene <- rownames(data)
            p1 <- p + theme(legend.position = "right") + geom_text_repel(data = data, 
                aes(x = logFC, y = -log10(P.Value), label = gene), 
                size = 3, box.padding = unit(0.5, "lines"), segment.color = "black", 
                show.legend = FALSE)
            plotout(p = p1, num = paste0("_", caseLabel, "_vs_", 
                controlLabel, "_volcano"), h = 6, w = 7, od = od)
        }
        suppressPackageStartupMessages(library(ComplexHeatmap))
        if (length(DEGs) > 1000) {
            plotlist <- nrDEG[DEGs, ] %>% arrange(desc(abs(logFC))) %>% 
                rownames_to_column(var = "gene") %>% pull(gene) %>% 
                head(1000)
        }
        else {
            plotlist <- DEGs
        }
        exp <- DEG_exp[plotlist, match(DEG_pdata[, 1], colnames(DEG_exp))]
        plotdata <- apply(as.matrix(exp), 1, function(x) scale(x, 
            scale = T, center = T)) %>% t()
        colnames(plotdata) <- colnames(exp)
        p <- Heatmap(as.matrix(plotdata), column_split = DEG_pdata[, 
            2], name = "zscore", show_row_names = F, show_column_names = F, 
            column_title_gp = gpar(fill = c(caseLabel = color_fun[2], 
                controlLabel = color_fun[1])))
        plotout(p = p, num = paste0("_", caseLabel, "_vs_", controlLabel, 
            "_heatmap"), h = 6, w = 10, od = od)
        p2 <- ggplotify::as.ggplot(p)
        p <- cowplot::plot_grid(p1, p2, rel_widths = c(7, 10), 
            ncol = 2, labels = c("AUTO"), labels = NA, scale = c(1, 
                0.95))
        plotout(p = p, num = paste0("_", caseLabel, "_vs_", controlLabel, 
            "_DEGs"), h = 6, w = 17, od = od)
    }
    return(list(DEGs = DEGs, nrDEG = nrDEG))
}

plotout = 
function (od = NULL, name = NULL, num = NULL, w = 6, h = 6, p = NULL, 
    plot_tif = FALSE) 
{
    if (!dir.exists(od)) {
        dir.create(od, recursive = TRUE)
    }
    if (!is.null(name)) {
        pdf(file = paste0(od, "/Figure_", name, ".pdf"), width = w, 
            height = h)
        invisible(print(p))
        dev.off()
        if (plot_tif) {
            tiff(file = paste0(od, "/Figure_", name, "-72.tif"), 
                width = w, height = h, units = "in", res = 72)
            print(p)
            dev.off()
            tiff(file = paste0(od, "/Figure_", name, "-300.tif"), 
                width = w, height = h, units = "in", res = 300)
            print(p)
            dev.off()
        }
    }
    else if (!is.null(num)) {
        pdf(file = paste0(od, "/Figure", num, ".pdf"), width = w, 
            height = h)
        print(p)
        dev.off()
        if (plot_tif) {
            tiff(file = paste0(od, "/Figure", num, "-72.tif"), 
                width = w, height = h, units = "in", res = 72)
            print(p)
            dev.off()
            tiff(file = paste0(od, "/Figure", num, "-300.tif"), 
                width = w, height = h, units = "in", res = 300)
            print(p)
            dev.off()
        }
    }
    else {
        message("there is no name/num for figure")
    }
}

