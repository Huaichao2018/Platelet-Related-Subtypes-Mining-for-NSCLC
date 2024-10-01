skcm_expr <- readRDS('D:/r_workplace/BEST3/data/ICB_TCGA_SKCM_exprs.rds')
skcm_clinical <- readRDS('D:/r_workplace/BEST3/data/ICB_TCGA_SKCM_clinical.rds')
signature_coef = read.delim("D:/r_workplace/BEST3/out/2.xgene_subtypes/coef_file.txt")
  od <- "D:/r_workplace/OUT/imm"


library(dplyr)
library(tibble)
library(glue)

skcm_clinical %<>%
    select(1, OS_Time, OS_Status, Response_infor) %>%
    rename(sample = 1, time = 2, status = 3, response = 4) %>% 
    mutate(response = ifelse(response == 1,'Y','N')) 
color_fun1  <- color_fun3 <- color_fun3 <- RColorBrewer::brewer.pal(9,'Set1')
skcm_res <- modelscore_km_roc(
    signature_coef = signature_coef,
    od = glue("{od}"),
    no_roc_pheatmap = TRUE, best_cut = T,
    exp = skcm_expr, clin =skcm_clinical,
    # time = roc_time,
    dataset = "TCGA_SKCM"
)

skcm_res$Score[rownames(skcm_res$Group),,drop = F] %>% bind_cols(skcm_res$Group) -> df
df %<>% rownames_to_column('sample') %>% inner_join(skcm_clinical)

library(ggpubr)

p2 = ggboxplot(df %>% filter(!is.na(response)),'response','Score',fill = 'response')+stat_compare_means()
p2
ggsave(p2,filename = "D:/r_workplace/OUT/imm/p2.pdf",width = 3,height = 4)
remove.packages("rlang")
install.packages("rlang")

library(rlang)
library(dplyr)
library(ggplot2)
p = df %>%
    filter(!is.na(response)) %>%
    ggplot(aes(x = Group, fill = response)) +
    geom_bar(position = position_fill())+
    scale_y_continuous(labels = scales::label_percent())+
    ggstatsplot::theme_ggstatsplot()+
    theme(legend.position = 'top')+
    labs(y = NULL)
p
ggsave(p,filename = "out/immunotherapy.3/p3.pdf",width = 3,height = 4)
