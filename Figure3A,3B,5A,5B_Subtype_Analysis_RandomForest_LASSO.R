
suppressMessages(library(tidyverse)); suppressMessages(library(magrittr))

c(
    list.files("/Pub/Users/cuiye/RCodes/UserCode", recursive = T, full.names = T, pattern = "\\.R$")
) %>%
    walk(source)

# 设置结果路径
out_home <- "D:/r_workplace/OUT"
setwd(out_home)
load("D:/r_workplace/BEST3/out/1.molecular_characteristics/xgene.RData")
load("D:/r_workplace/BEST3/data/train_data.RData", verbose = T)
od <- "out/xgene_subtypes.2"
dir.create(od, recursive = TRUE)

color_fun1  <- color_fun3 <- color_fun3 <- RColorBrewer::brewer.pal(9,'Set1')

xgene_k <- "k3"
xgene_distance <- "pearson"
xgene_clusteralg <- "pam"
# 亚型识别
xgene_cc_res <- ConsensusClusterPlus_km(
    genelist = xgene, od = od,
    exp = train_data$tumor_exprs, clin = train_data$data_clinical,
    input_distance = xgene_distance, input_clusteralg = xgene_clusteralg,
    seed = 123456, maxK = 5
)
save(xgene_cc_res, file = str_glue("{od}/xgene_cc_res.RData"))

# 亚型间差异基因
sample_cluster <- xgene_cc_res[["res_list"]][[xgene_distance]][[xgene_clusteralg]][[xgene_k]][["cluster"]]
unique_cluster <- unique(sample_cluster$Cluster) %>% sort()


combine_deg_res <- map_dfr(unique_cluster, function(cluster) {
    sample_pdata <- sample_cluster %>% mutate(Cluster = ifelse(Cluster == cluster, cluster, "Rest"))
    deg_res <- limma_deg(
        od = glue("{od}/deg/"),
        DEG_exp = train_data$tumor_exprs, DEG_pdata = sample_pdata,
        controlLabel = "Rest", caseLabel = cluster,
        DEG_FC = 0, DEG_P = 1, pvalue = NULL, saveplot = FALSE, color_fun = color_fun1
    )
    return(deg_res$nrDEG %>% rownames_to_column("gene") %>% mutate(Cluster = cluster))
})
print(deg_res) 
print(limma_deg)


deg_gene <- combine_deg_res %>%
    dplyr::filter(logFC >= 1, adj.P.Val < 0.01) %>%
    pull(gene) %>%
    unique()
length(deg_gene) # 3454
save(deg_gene, file = str_glue("{od}/deg_gene.RData"))

load("D:/r_workplace/BEST3/out/2.xgene_subtypes/deg_gene.RData")
# 样本分布
dat <- merge(train_data$data_clinical, sample_cluster) %>% column_to_rownames("sample")
table(dat %>% select(Cluster, Histological.Type))

# 单因素分析识别预后基因
cox_res <- signature_cox(
    signaturelist = deg_gene, 
    exp = train_data$tumor_exprs, 
    clin = train_data$data_clinical,
    coxp = 0.01, 
    bygroup = TRUE, 
    savekmplot = FALSE
)

cox_res %>% names()
cox_res$sigcoxResult %>% rownames()
write_tsv(
    x = cox_res$coxResult %>% rownames_to_column("gene") %>% as_tibble(),
    file = file.path(od, "单因素cox全部结果.tsv")
)

library(survival)
library(randomForestSRC)
infor <- train_data$tumor_exprs[cox_res$cox_signature, ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    merge(., train_data$data_clinical %>% select(sample, time, status)) %>%
    select(-sample)
set.seed(123456)
rfs_obj <- rfsrc(Surv(time, status) ~ ., data = infor, nodesize = 15, importance = TRUE)
# rfs_vs <- var.select(object = rfs_obj, method = "vh", nrep = 10, nstep = 5)
rfs_vs <- var.select(object = rfs_obj, conservative = "high")
top_vars <- rfs_vs$topvars
length(top_vars)

save(list = c("rfs_obj","rfs_vs","top_vars"),file = file.path(od,"randomforest_res.RData"))

write_tsv(x = tibble("top_vars" = top_vars),file.path(od,'randomForestSRC_topvars_gene.tsv'))


library(Matrix)
# LASSO-COX选择关键基因
lasso_res <- lasso_model(
    signaturelist = cox_res$cox_signature,
    od = str_glue("{od}/lasso"),
    seed = 123456,
    exp = train_data$tumor_exprs, clin = train_data$data_clinical,
    signaturetype = "cox_gene"
)
length(lasso_res[[1]])

# 两种方法取交集
genelist <- Reduce(intersect, list(top_vars, lasso_res[[1]]))
length(genelist)

saveRDS(genelist, file = "genelist.rds")


a <- ggvenn::ggvenn(
  data = list('RandomForest get Genes' = top_vars,'LASSO get Genes' = lasso_res[[1]]),
  digits = 1,
  fill_color = ggsci::pal_jco()(2),
  fill_alpha = 0.6,
  text_size = 4
)

ggsave(plot = a,filename = file.path(od,'Figure_forest_lasso_gene_venn.pdf'),width = 4,height = 4)


get_function(modelscore_km_roc) %>% write_lines("/Pub/Users/wangyk/project/Poroject/F210823005_非小细胞肺癌血小板/aftersale3/src/1.TMP",append = T)


# 多因素cox系数构建模型
library(survival)
library(survminer)
multi_cox_model <- as.formula(paste0("Surv(time, status) ~", str_c("`", genelist, "`", collapse = "+"))) %>% coxph(data = infor)
signature_coef <- summary(multi_cox_model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("signature")

model_res <- modelscore_km_roc(
    signature_coef = signature_coef,
    od = str_glue("{od}/train/"),
    no_roc_pheatmap = TRUE,
    exp = train_data$tumor_exprs, clin = train_data$data_clinical,
    # time = roc_time,
    dataset = "TCGA"
)
save(model_res,file = file.path(od,"model_res.rdata"))

write_tsv(x = signature_coef[,1:2],file = file.path(od,'coef_file.txt'))

res <- uni_multi_cox(
    od = str_glue("{od}/train/"),
    infor = model_res$Group %>% rownames_to_column("sample") %>%
        merge(., train_data$data_clinical) %>% rename(Score = Group),
    factors = c("Histological.Type", "Stage", "Age", "Sex", "Score"),
    dataset = "TCGA", coxp = 1, w = 10
)

load("D:/r_workplace/BEST3/GSE50081_GPL570/clinical.RData")
load("D:/r_workplace/BEST3/GSE50081_GPL570/expression.RData")

vali_res <- modelscore_km_roc(
    signature_coef = signature_coef,
    od = str_glue("{od}/vali/"),
    no_roc_pheatmap = TRUE, best_cut = T,
    exp = expression, clin = clinical %>% rename(sample = 1,time = OS.Time,status = OS.Status),
    # time = roc_time,
    dataset = "GSE50081_GPL570"
)

res <- uni_multi_cox(
    od = str_glue("{od}/vali/"),
    infor = vali_res$Group %>% rownames_to_column("sample") %>%
        inner_join(clinical %>% dplyr::rename(sample = 1) %>%
            mutate(Age = ifelse(Age < 60, "<60", ">=60"))) %>% 
            rename(Score = Group,time = OS.Time,status = OS.Status),
    factors = c("Score", "Smoking", "Age", "Sex"),
    dataset = "GSE50081_GPL570", coxp = 1, w = 10
)

