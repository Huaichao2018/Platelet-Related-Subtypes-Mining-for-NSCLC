 
library(tidyverse)

# Set result path
out_home <- "~/project/Poroject/F210823005_NSCLC_platelet/aftersale3/"
setwd(out_home)
od <- str_glue("{out_home}/1.molecular_characteristics")
load("data/train_data.RData",
     verbose = T
)
dir.create(od)

# Differential expression analysis between cancer and normal samples
cancer_deg_res <- limma_deg(
  od = str_glue("{od}/exprs/"), DEG_exp = train_data$data_exprs,
  DEG_pdata = train_data$data_pdata %>% mutate(Group = ifelse(`_sample_type` == "Primary Tumor", "Tumor", "Normal")) %>% select(sample, Group) %>% as.data.frame(),
  controlLabel = "Normal", caseLabel = "Tumor",
  DEG_FC = 1, DEG_P = 0.01, pvalue = NULL, saveplot = FALSE,
  color_fun = color_fun1
)
save(cancer_deg_res, file = str_glue("{od}/cancer_deg_res.RData"))


load("data/xgene_genelist.RData")
# Intersect with platelet genes
xgene <- Reduce(intersect, list(xgene_genelist, cancer_deg_res$DEGs))
length(xgene) # 
save(xgene, file = str_glue("{od}/xgene.RData"))


a <- ggvenn::ggvenn(
  data = list('Platelet Gene' = xgene_genelist,'NSCLC DEGs\n(Tumor vs. Normal)' = cancer_deg_res$DEGs),
  digits = 1,
  fill_color = ggsci::pal_rickandmorty()(12)[3:4],
  fill_alpha = 0.6,
  text_size = 5
)

ggsave(plot = a,filename = file.path(od,'Figure_DE_PLATELET_gene_venn.pdf'),width = 6,height = 6)

