library(tidyverse)
library(ggthemes)
library(scales)
setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd3525082/TCGA_rnaseq")

tibble_expr_sample_clin <- read_tsv("expr_data_select_clin_annot_FPKM-UQ.txt")

target_gene <- "MAGEA9B"
y_title <- paste(target_gene, " gene expression level", "(normalized FPKM-UQ values)")
my.labels <- c("Adrenocortical carcinoma (n=79)","B-cell lymphoma (DLBCL) (n=48)","Bladder carcinoma (n=414)","Breast carcinoma (n=1108)","Cervical carcinoma (n=306)","Cholangiocarcinoma (n=36)","Colon adenocarcinoma (n=477)","Esophageal carcinoma (n=160)","Glioblastoma (n=169)","Head&neck carcinoma (n=502)","Kidney chromophobe tumor (n=65)","Kidney renal cell carcinoma (n=539)","Kidney papillary carcinoma (n=289)","Liver hepatocellular carcinoma (n=368)","Low grade glioma (n=1055)","Lung adenocarcinoma (n=535)","Lung squamous carcinoma (n=502)","Mesothelioma (n=86)","Ovarian carcinoma (n=379)","Pancreas carcinoma (n=178)","Pheochromocytoma paraganglioma (n=183)","Prostate adenocarcinoma (n=496)","Rectal adenocarcinoma (n=167)","Sarcoma (n=263)","Skin melanoma (n=471)","Stomach adenocarcinoma (n=375)","Testicular germ cell tumor (n=156)","Thyroid carcinoma (n=509)","Thymoma (n=119)","Uterine endometrial carcinoma (n=550)","Uterine carcinosarcoma (n=56)","Uveal melanoma (n=78)")
plot_table  <- dplyr::filter(tibble_expr_sample_clin, tissue_type == "primary_tumor" | tissue_type ==  "metastasis")

ggplot(plot_table, aes(x = tumor_type, y = get(target_gene))) +
  geom_boxplot(outlier.alpha = 0, color = "red") +
  geom_jitter(aes(color = tissue_type), size = 0.8, width = 0.3, alpha = 0.3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(face = "bold", size = 12, color = "black"), axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(size = 10)) +
#  coord_cartesian(ylim = c(0, 10000000)) +
  labs(x ="Tumor type (number of samples)", y = y_title) + 
  scale_x_discrete(labels = my.labels) +
  scale_y_continuous(labels = comma) +
  theme(panel.grid.major.x = element_blank()) +
  geom_vline(xintercept=seq(0.5,33,1), lwd = 0.2, colour = "grey")

