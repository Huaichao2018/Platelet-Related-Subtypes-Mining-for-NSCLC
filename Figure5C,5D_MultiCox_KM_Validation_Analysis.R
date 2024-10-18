# Load required data
load("D:/r_workplace/BEST3/out/2.xgene_subtypes/model_res.rdata")
load("D:/r_workplace/BEST3/data/train_data.RData", verbose = TRUE)

# Set output directory and create it if it does not exist
od <- "D:/r_workplace/OUT"
dir.create(od)

# Load required libraries
library(glue)

# Perform univariate and multivariate Cox regression on the training dataset
res <- uni_multi_cox(
  od = glue("{od}/train/"),  # Output directory for results
  infor = model_res$Group %>% 
    rownames_to_column("sample") %>%
    merge(., clinical) %>%
    rename(Score = Group),  # Rename 'Group' to 'Score'
  factors = c("Sample.Type", "Stage", "Age", "Sex", "Score"),  # Covariates for Cox model
  dataset = "TCGA", 
  coxp = 1, 
  w = 8  # Width of output plot
)

# Load clinical and expression data for validation
load("D:/r_workplace/BEST3/GSE50081_GPL570/clinical.RData")
load("D:/r_workplace/BEST3/GSE50081_GPL570/expression.RData")

# Perform model validation with Kaplan-Meier and ROC analysis
vali_res <- modelscore_km_roc(
  signature_coef = signature_coef,
  od = str_glue("{od}/vali/"),  # Validation output directory
  no_roc_pheatmap = TRUE,  # Disable ROC heatmap
  best_cut = TRUE,  # Use the best cutoff value
  exp = expression,  # Expression data
  clin = clinical %>% 
    rename(sample = 1, time = OS.Time, status = OS.Status),  # Rename clinical columns
  dataset = "GSE50081_GPL570"
)

# Perform univariate and multivariate Cox regression on the validation dataset
res <- uni_multi_cox(
  od = str_glue("{od}/vali/"),
  infor = vali_res$Group %>% 
    rownames_to_column("sample") %>%
    inner_join(
      clinical %>%
        dplyr::rename(sample = 1) %>%
        mutate(Age = ifelse(Age < 60, "<60", ">=60"))  # Categorize age
    ) %>% 
    rename(Score = Group, time = OS.Time, status = OS.Status),
  factors = c("Score", "Smoking", "Age", "Sex"),
  dataset = "GSE50081_GPL570", 
  coxp = 1, 
  w = 10
)

# Prepare data for Kaplan-Meier analysis
group_df <- vali_res$Group %>% rownames_to_column("Sample")

# Merge clinical data with group information
clinical_with_group <- clinical %>%
  inner_join(group_df, by = "Sample") %>%
  rename(Group = Group)

# Compute Kaplan-Meier survival curves
km_fit <- survfit(Surv(OS.Time, OS.Status) ~ Group, data = clinical_with_group)

# Load survival plotting library
library(survminer)

# Generate Kaplan-Meier plot with risk table, p-value, and confidence intervals
p <- ggsurv <- ggsurvplot(
  km_fit,
  data = clinical_with_group,
  risk.table = TRUE,  # Display risk table
  pval = TRUE,  # Display p-value
  conf.int = TRUE,  # Display confidence interval
  xlab = "Time (days)",  # X-axis label
  ylab = "Survival Probability"  # Y-axis label
)

# Display plot
p

# Save survival curve plot
ggsave(
  filename = str_glue("{od}/vali/survival_plot.png"),
  plot = ggsurv$plot, 
  width = 8, height = 6
)


















