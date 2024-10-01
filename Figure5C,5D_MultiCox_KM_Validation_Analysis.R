
load("D:/r_workplace/BEST3/out/2.xgene_subtypes/model_res.rdata")
load("D:/r_workplace/BEST3/data/train_data.RData", verbose = T)
od <- "D:/r_workplace/OUT"
dir.create(od)
library(glue)
res <- uni_multi_cox(
    od = glue("{od}/train/"),
    infor = model_res$Group %>% rownames_to_column("sample") %>%
        merge(., clinical) %>% rename(Score = Group),
    factors = c("Sample.Type", "Stage", "Age", "Sex", "Score"),
    dataset = "TCGA", coxp = 1, w = 8
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

library(survival)
library(dplyr)
# 检查 clinical 数据框的列名
colnames(clinical)
# 检查 vali_res$Group 的结构
head(vali_res$Group)

# 确保 vali_res$Group 的行名为样本名，并将其转换为数据框
group_df <- vali_res$Group %>% rownames_to_column("Sample")

# 合并数据框
clinical_with_group <- clinical %>%
  inner_join(group_df, by = "Sample") %>%
  rename(Group = Group) # 确保重命名以便后续使用

# 计算 Kaplan-Meier 生存模型
km_fit <- survfit(Surv(OS.Time, OS.Status) ~ Group, data = clinical_with_group)

# 检查 km_fit 的结构
summary(km_fit)
library(survminer)

p=ggsurv <- ggsurvplot(
  km_fit,
  data = clinical_with_group,
  risk.table = TRUE, # 显示风险表
  pval = TRUE, # 显示p值
  conf.int = TRUE, # 显示置信区间
  xlab = "Time (days)", # x轴标签
  ylab = "Survival Probability" # y轴标签
)
p
# 保存生存曲线图
ggsave(filename = str_glue("{od}/vali/survival_plot.png"), plot = ggsurv$plot, width = 8, height = 6)


###########################GES37745验证#########################################

load("D:/r_workplace/BEST3/output_data.RData")
library(dplyr)
unique(new_clinical$OS.Status)

new_clinical <- new_clinical %>%
  filter(!is.na(OS.Status) & !is.na(OS.Time) & 
           OS.Status != "not known" & OS.Time != "not known" & 
           OS.Status != "" & OS.Time != "") %>%
  mutate(OS.Status = ifelse(OS.Status == "yes", 0, 1),  # 转换为 1 和 0
         OS.Time = as.numeric(OS.Time)) %>%
  filter(!is.na(OS.Time))  # 过滤出 OS.Time 非 NA 的行


vali_res1 <- modelscore_km_roc(
  signature_coef = signature_coef,
  od = str_glue("{od}/vali/"),
  no_roc_pheatmap = TRUE, best_cut = T,
  exp = exp, clin = new_clinical %>% rename(sample = 1,time = OS.Time,status = OS.Status),
  # time = roc_time,
  dataset = "output_data"
)

res1 <- uni_multi_cox(
  od = str_glue("{od}/vali/"),
  infor = vali_res1$Group %>% rownames_to_column("sample") %>%
    inner_join(new_clinical %>% dplyr::rename(sample = 1) %>%
                 mutate(Age = ifelse(Age < 60, "<60", ">=60"))) %>% 
    rename(Score = Group,time = OS.Time,status = OS.Status),
  factors = c("Score", "Smoking", "Age", "Sex"),
  dataset = "output_data", coxp = 1, w = 10
)



group_df <- vali_res1$Group %>% rownames_to_column("Sample")

# 合并数据框
clinical_with_group <- new_clinical %>%
  inner_join(group_df, by = "Sample") %>%
  rename(Group = Group) # 确保重命名以便后续使用


# 计算 Kaplan-Meier 生存模型
km_fit <- survfit(Surv(OS.Time, OS.Status) ~ Group, data = clinical_with_group)

# 检查 km_fit 的结构
summary(km_fit)
library(survminer)

p=ggsurv <- ggsurvplot(
  km_fit,
  data = clinical_with_group,
  risk.table = TRUE, # 显示风险表
  pval = TRUE, # 显示p值
  conf.int = TRUE, # 显示置信区间
  xlab = "Time (days)", # x轴标签
  ylab = "Survival Probability" # y轴标签
)
p
# 保存生存曲线图
ggsave(filename = str_glue("{od}/vali/survival_plot1.png"), plot = ggsurv$plot, width = 10, height = 6)


##################GSE132133########################################


load("D:/r_workplace/OUT/GSE13213.RData")
library(dplyr)
unique(clinical_2$OS.Status)

#clinical_2 <- clinical_2 %>%
  #filter(!is.na(OS.Status) & !is.na(OS.Time) & 
          # OS.Status != "not known" & OS.Time != "not known" & 
          # OS.Status != "" & OS.Time != "") %>%
  #mutate(OS.Status = ifelse(OS.Status == "yes", 0, 1),  # 转换为 1 和 0
        # OS.Time = as.numeric(OS.Time)) %>%
  #filter(!is.na(OS.Time))  # 过滤出 OS.Time 非 NA 的行
#summary(clinical_2$OS.Time)
#summary(clinical_2$OS.Status)


vali_res2 <- modelscore_km_roc(
  signature_coef = signature_coef,
  od = str_glue("{od}/vali/"),
  no_roc_pheatmap = TRUE, best_cut = T,
  exp = exp1, clin = clinical_2%>% rename(sample = 1,time = OS.Time,status = OS.Status),
  # time = roc_time,
  dataset = "GSE13213.RData"
)

res2 <- uni_multi_cox(
  od = str_glue("{od}/vali/"),
  infor = vali_res2$Group%>% rownames_to_column("sample") %>%
    inner_join(clinical_2%>% dplyr::rename(sample = 1) %>%
                 mutate(Age = ifelse(Age < 60, "<60", ">=60"))) %>% 
    rename(Score = Group,time = OS.Time,status = OS.Status),
  factors = c("Score", "Smoking", "Age", "Sex"),
  dataset = "GSE13213.RData", coxp = 1, w = 10
)



group_df <- vali_res2$Group %>% rownames_to_column("Sample")

# 合并数据框
clinical_with_group2 <- clinical_2 %>%
  inner_join(group_df, by = "Sample") %>%
  rename(Group = Group) # 确保重命名以便后续使用

clinical_with_group2 <- clinical_with_group2 %>%
  mutate(OS.Time = as.numeric(OS.Time)) %>%
  mutate(OS.Status = as.numeric(OS.Status))

# 计算 Kaplan-Meier 生存模型
km_fit2 <- survfit(Surv(OS.Time, OS.Status) ~ Group, data = clinical_with_group2)

# 检查 km_fit 的结构
summary(km_fit2)
library(survminer)

p=ggsurv <- ggsurvplot(
  km_fit2,
  data = clinical_with_group2,
  risk.table = TRUE, # 显示风险表
  pval = TRUE, # 显示p值
  conf.int = TRUE, # 显示置信区间
  xlab = "Time (days)", # x轴标签
  ylab = "Survival Probability" # y轴标签
)
p
# 保存生存曲线图
ggsave(filename = str_glue("{od}/vali/survival_plot2.png"), plot = ggsurv$plot, width = 10, height = 6)

#################################GSE11969########################

load("D:/r_workplace/OUT/GSE11969.RData")
library(dplyr)
unique(clinical_3$OS.Status)

#clinical_2 <- clinical_2 %>%
#filter(!is.na(OS.Status) & !is.na(OS.Time) & 
# OS.Status != "not known" & OS.Time != "not known" & 
# OS.Status != "" & OS.Time != "") %>%
#mutate(OS.Status = ifelse(OS.Status == "yes", 0, 1),  # 转换为 1 和 0
# OS.Time = as.numeric(OS.Time)) %>%
#filter(!is.na(OS.Time))  # 过滤出 OS.Time 非 NA 的行
#summary(clinical_2$OS.Time)
#summary(clinical_2$OS.Status)


vali_res3 <- modelscore_km_roc(
  signature_coef = signature_coef,
  od = str_glue("{od}/vali/"),
  no_roc_pheatmap = TRUE, best_cut = T,
  exp = exp2, clin = clinical_3%>% rename(sample = 1,time = OS.Time,status = OS.Status),
  # time = roc_time,
  dataset = "GSE11969.RData"
)

res3 <- uni_multi_cox(
  od = str_glue("{od}/vali/"),
  infor = vali_res3$Group%>% rownames_to_column("sample") %>%
    inner_join(clinical_3%>% dplyr::rename(sample = 1) %>%
                 mutate(Age = ifelse(Age < 60, "<60", ">=60"))) %>% 
    rename(Score = Group,time = OS.Time,status = OS.Status),
  factors = c("Score", "Smoking", "Age", "Sex"),
  dataset = "GSE11969.RData", coxp = 1, w = 10
)



group_df <- vali_res3$Group %>% rownames_to_column("Sample")

# 合并数据框
clinical_with_group3 <- clinical_3 %>%
  inner_join(group_df, by = "Sample") %>%
  rename(Group = Group) # 确保重命名以便后续使用

clinical_with_group3<- clinical_with_group3 %>%
  mutate(OS.Time = as.numeric(OS.Time)) %>%
  mutate(OS.Status = as.numeric(OS.Status))

# 计算 Kaplan-Meier 生存模型
km_fit3<- survfit(Surv(OS.Time, OS.Status) ~ Group, data = clinical_with_group3)

# 检查 km_fit 的结构
summary(km_fit3)
library(survminer)

p=ggsurv <- ggsurvplot(
  km_fit2,
  data = clinical_with_group3,
  risk.table = TRUE, # 显示风险表
  pval = TRUE, # 显示p值
  conf.int = TRUE, # 显示置信区间
  xlab = "Time (days)", # x轴标签
  ylab = "Survival Probability" # y轴标签
)
p
#################################GSE42127#########################


load("D:/r_workplace/OUT/GSE42127.RData")
library(dplyr)
unique(clinical_3$OS.Status)

#clinical_2 <- clinical_2 %>%
#filter(!is.na(OS.Status) & !is.na(OS.Time) & 
# OS.Status != "not known" & OS.Time != "not known" & 
# OS.Status != "" & OS.Time != "") %>%
#mutate(OS.Status = ifelse(OS.Status == "yes", 0, 1),  # 转换为 1 和 0
# OS.Time = as.numeric(OS.Time)) %>%
#filter(!is.na(OS.Time))  # 过滤出 OS.Time 非 NA 的行
#summary(clinical_2$OS.Time)
#summary(clinical_2$OS.Status)


vali_res4 <- modelscore_km_roc(
  signature_coef = signature_coef,
  od = str_glue("{od}/vali/"),
  no_roc_pheatmap = TRUE, best_cut = T,
  exp = exp3, clin = clinical_3%>% rename(sample = 1,time = OS.Time,status = OS.Status),
  # time = roc_time,
  dataset = "GSE42127.RData"
)

res4 <- uni_multi_cox(
  od = str_glue("{od}/vali/"),
  infor = vali_res4$Group%>% rownames_to_column("sample") %>%
    inner_join(clinical_3%>% dplyr::rename(sample = 1) %>%
                 mutate(Age = ifelse(Age < 60, "<60", ">=60"))) %>% 
    rename(Score = Group,time = OS.Time,status = OS.Status),
  factors = c("Score", "Smoking", "Age", "Sex"),
  dataset = "GSE42127.RData", coxp = 1, w = 10
)



group_df <- vali_res4$Group %>% rownames_to_column("Sample")

# 合并数据框
clinical_with_group4 <- clinical_3 %>%
  inner_join(group_df, by = "Sample") %>%
  rename(Group = Group) # 确保重命名以便后续使用

clinical_with_group4<- clinical_with_group4 %>%
  mutate(OS.Time = as.numeric(OS.Time)) %>%
  mutate(OS.Status = as.numeric(OS.Status))

# 计算 Kaplan-Meier 生存模型
km_fit4<- survfit(Surv(OS.Time, OS.Status) ~ Group, data = clinical_with_group4)

# 检查 km_fit 的结构
summary(km_fit4)
library(survminer)

p=ggsurv <- ggsurvplot(
  km_fit4,
  data = clinical_with_group4,
  risk.table = TRUE, # 显示风险表
  pval = TRUE, # 显示p值
  conf.int = TRUE, # 显示置信区间
  xlab = "Time (days)", # x轴标签
  ylab = "Survival Probability" # y轴标签
)
p
##################GSE30219#####################

load("D:/r_workplace/OUT/GSE30219.RData")
library(dplyr)
unique(cli$OS.Status)

#clinical_2 <- clinical_2 %>%
#filter(!is.na(OS.Status) & !is.na(OS.Time) & 
# OS.Status != "not known" & OS.Time != "not known" & 
# OS.Status != "" & OS.Time != "") %>%
#mutate(OS.Status = ifelse(OS.Status == "yes", 0, 1),  # 转换为 1 和 0
# OS.Time = as.numeric(OS.Time)) %>%
#filter(!is.na(OS.Time))  # 过滤出 OS.Time 非 NA 的行
#summary(clinical_2$OS.Time)
#summary(clinical_2$OS.Status)


vali_res5 <- modelscore_km_roc(
  signature_coef = signature_coef,
  od = str_glue("{od}/vali/"),
  no_roc_pheatmap = TRUE, best_cut = T,
  exp = exp, clin = cli%>% rename(sample = 1,time = OS.Time,status = OS.Status),
  # time = roc_time,
  dataset = "GSE30219.RData"
)

res5 <- uni_multi_cox(
  od = str_glue("{od}/vali/"),
  infor = vali_res5$Group%>% rownames_to_column("sample") %>%
    inner_join(cli%>% dplyr::rename(sample = 1) %>%
                 mutate(Age = ifelse(Age < 60, "<60", ">=60"))) %>% 
    rename(Score = Group,time = OS.Time,status = OS.Status),
  factors = c("Score", "Smoking", "Age", "Sex"),
  dataset = "GSE30219.RData", coxp = 1, w = 10
)



group_df <- vali_res5$Group %>% rownames_to_column("Sample")

# 合并数据框
clinical_with_group5<- cli%>%
  inner_join(group_df, by = "Sample") %>%
  rename(Group = Group) # 确保重命名以便后续使用

clinical_with_group5<- clinical_with_group5 %>%
  mutate(OS.Time = as.numeric(OS.Time)) %>%
  mutate(OS.Status = as.numeric(OS.Status))

# 计算 Kaplan-Meier 生存模型
km_fit5<- survfit(Surv(OS.Time, OS.Status) ~ Group, data = clinical_with_group5)

# 检查 km_fit 的结构
summary(km_fit5)
library(survminer)

p=ggsurv <- ggsurvplot(
  km_fit5,
  data = clinical_with_group5,
  risk.table = TRUE, # 显示风险表
  pval = TRUE, # 显示p值
  conf.int = TRUE, # 显示置信区间
  xlab = "Time (days)", # x轴标签
  ylab = "Survival Probability" # y轴标签
)
p
##################################





























