library(tidyverse)
library(stringr)

setwd("/media/nicolas/ce71f585-e4d4-43a3-adea-f042dd3525082/TCGA_rnaseq")                              #set working directory as appropriate (it must contain the source files) 

#Step1: generate a data frame (tibble_expr) with genes (identified by Ensembl gene identifier) as rows and samples (identified by their GDC manifest file identifier) as columns.
file_names <- list.files(path = "TCGA_rnaseq_FPKM-UQ", pattern = ".FPKM-UQ.txt$", all.files = FALSE, full.names = TRUE)     #generate a vector listing the names of the files
gene_list <- read_tsv(file_names[1], col_names = FALSE)                                                                     #create a new tibble, copy the content of the first data file
colnames(gene_list) <- c("ensembl_gene_id", "x")                                                                            #name the 1st column with the gene identifier, 2nd column x, which will be deleted soon)
gene_list <- dplyr::select(gene_list, ensembl_gene_id) %>%                                                                  #keep only 1st column; add dplyr:: to the select function, which has another action in base R
arrange(ensembl_gene_id)                                                                                                    #sort by gene ID

tibble_expr <- gene_list                                                                                    #create a new tibble to host the data, copy the gene identifiers
counter <- 0                                                                                                #the counter that will allow to follow the progress of the loop is set to zero
for (sample_name in file_names) {                                                                           #start of the loop: sample files are processed 1 by 1
  tibble_temp <- read_tsv(sample_name, col_names = FALSE)                                                   #import sample file
  buffer_name <- str_replace(sample_name, "TCGA_rnaseq_FPKM-UQ/", "")                                       #generate a sample name without the prefix
  buffer_name <- str_replace(buffer_name, ".FPKM-UQ.txt", "")                                               #remove the suffix from this sample name
  colnames(tibble_temp) <- c("ensembl_gene_id", buffer_name)                                                #rename the columns
  tibble_expr <- full_join(tibble_expr, tibble_temp, by = "ensembl_gene_id")                                #join this table to the growing table with the gene names and the sample values already added
  counter <- counter + 1                                                                                    #increase counter value by 1
  print(counter)                                                                                            #print counter value on screen
}                                                                                                           #move to next sample and process similarly
tibble_expr$ensembl_gene_id <- str_sub(tibble_expr$ensembl_gene_id,1,15)                                    #remove suffix (version) from the gene identifiers
tibble_expr <- arrange(tibble_expr, ensembl_gene_id)                                                        #sort by gene ID

#write.table(tibble_expr, file = "tibble_expr.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)      #if necessary (to avoid repeated lengthy runs) save data as new table

#tibble_expr <- read_tsv("tibble_expr.txt")                                                                              #load table data (start here if previous step was done and data saved)

#Step2: replace the ensembl gene identifier by the HGNC symbol retrieved from BioMart, select only protein-coding genes
ensembl <- biomaRt::useMart("ensembl")
ensembl = biomaRt::useDataset ("hsapiens_gene_ensembl", mart=ensembl)
#select the Biomart gene attributes HGNC symbol and type, matching the Ensembl genes from the gene list
mart_annot <- as.tibble(biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"), filters = "ensembl_gene_id", values = tibble_expr$ensembl_gene_id, mart=ensembl)) %>%
  dplyr::filter(gene_biotype == "protein_coding" | hgnc_symbol == "XIST") %>%     #select only the protein-coding genes, select also XIST (control for female tissues)
  dplyr::filter(!is.na(hgnc_symbol)) %>%                                          #remove genes without HGNC symbol (labeled "NA")
  dplyr::filter(hgnc_symbol!="") %>%                                              #remove genes without HGNC symbol (empty field)
  dplyr::select(ensembl_gene_id, hgnc_symbol) %>%                                 #select the columns named gene_id and HGNC_symbol
  arrange(hgnc_symbol)                                                            #sort by HGNC symbol
tibble_hgnc <- left_join(mart_annot, tibble_expr, by = "ensembl_gene_id") %>%     #join this table with the gene annotation table by matching gene identifier
  dplyr::select(2:11094) %>%                                                      #select all the columns except the gene ID column (i.e. delete that column)
  arrange(hgnc_symbol) %>%                                                        #sort by gene symbol
# two genes are represented more than once, for obscure reasons (PRAMEF7 twice, PINX many times); the next step replaces each duplicated gene rows by a single row and sums up the expression values.
  group_by(hgnc_symbol) %>%                                                       #take rows with same gene symbol (i.e. duplicates) as one group
  summarize_all(sum) %>%                                                          #replace that group by a single row with summed expression values
  ungroup()                                                                       #remove grouping step

#Step3: transpose the tibble (each row is a sample, each column is a gene)
tibble_hgnc <- gather(tibble_hgnc, key = manifest_file_id, value = value, 2:ncol(tibble_hgnc))        #transpose the data frame (1st step: transform in tidy data)
tibble_hgnc <- spread(tibble_hgnc, key = hgnc_symbol, value = 'value')                                #transpose the data frame (2nd step: untidy with gene symbols as column headers)

manifest <- read.table("Downloaded_from_GDC/gdc_manifest.2018-09-06_FPKM-UQ.txt", header = TRUE)    #import manifest file obtained from GDC website
file_uuids <- manifest$id                                                                   #retrieve the file identifiers
#the following loop was copied from a website (https://www.biostars.org/p/306400/), I do not master its instructions, but it works
library(GenomicDataCommons)
TCGAtranslateID = function(file_ids, legacy = FALSE) {
  info <- files(legacy = legacy) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  id_list <- lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  barcodes_per_file <- sapply(id_list,length)
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}
corresp_sampleUUID_barcode <- as_tibble(TCGAtranslateID(file_uuids)) %>%            #change the data frame into a tibble
  arrange(file_id)                                                                  #sort by file_id identifier

corresp_sampleUUID_barcode <- corresp_sampleUUID_barcode %>%
  dplyr::filter(!is.na(submitter_id)) %>%                               #remove rows with NA values in the submitter_id column
  mutate(bcr_patient_barcode = str_sub(submitter_id, 1, 12)) %>%        #add a new row with the patient barcode identifier (= first 12 characters of sample barcode ID)
  mutate(tissue_type = as.factor(str_sub(submitter_id, 14, 15))) %>%    #add a new row with the tissue type identifier (= last 2 but 1 characters of sample barcode ID)
  mutate_all(as.character)                                              #change all values in character type

corresp_sampleUUID_barcode$tissue_type <- recode_factor(corresp_sampleUUID_barcode$tissue_type, "01" = "primary_tumor","02" = "primary_tumor","03" = "primary_tumor","05" = "primary_tumor","06"="metastasis","07"="metastasis","11"="adjacent_normal")

manifest <- rename(manifest, file_id = id) %>%
  rename(manifest_file_id = filename)
manifest$manifest_file_id <- str_replace(manifest$manifest_file_id,".FPKM-UQ.txt.gz", "")
corresp_manifestUUID_barcode <- full_join(manifest, corresp_sampleUUID_barcode, by = "file_id") %>%
  dplyr::select(manifest_file_id, submitter_id, bcr_patient_barcode, tissue_type)

sample_annot <- read_tsv("Downloaded_from_GDC/slide.tsv") %>%                                          #import sample annotations (obtained from)
  mutate(tumor_code = factor(str_sub(project_id, 6, 9))) %>%                                    #add column with 4-letter tumor code (derived from truncated TCGA-XXXX code)
  rename(submitter_id = sample_submitter_id) %>%                                                       #change column name
  dplyr::filter(section_location == "TOP")                                                             #select only the samples labeled TOP (to avoid doublet rows for individual samples)
full_sample_annot <- left_join(corresp_manifestUUID_barcode, sample_annot, by = "submitter_id") %>%    #join both sample annotation tables
  group_by(submitter_id) %>%                                                                           #group rowwise by submitter_id identifier
  slice(1) %>%                                                                                         #remove first line
  ungroup() %>%
  arrange(manifest_file_id)

#create vector with tumor codes in defined order
#library(TCGAbiolinks)
#list_tumor_types <- c("TCGA-LAML","TCGA-ACC","TCGA-BLCA","TCGA-LGG","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-DLBC","TCGA-MESO","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THYM","TCGA-THCA","TCGA-UCS","TCGA-UCEC","TCGA-UVM")
#for each tumor type, generate a table with tumor-specific clinical annotations
#for (tumor_type in list_tumor_types) {
#  query <- GDCquery(project = tumor_type, data.category = "Clinical", file.type = "xml")
#  GDCdownload(query)
#  clinical <- GDCprepare_clinic(query, clinical.info = "patient")
#tumor_type <- paste("TCGA_xml_clinical_data/", tumor_type, sep = "")
#write.table(clinical, file = tumor_type, sep = "\t", quote = FALSE, row.names = FALSE)
#}

tibble_clinannot <- tibble::tibble()  
list_tumor_types <- list.files(path = "TCGA_xml_clinical_data", pattern = "TCGA-", all.files = FALSE, full.names = TRUE)  #
for (tumor_type in list_tumor_types) {
print(tumor_type)
tibble_buffer <- read_tsv(tumor_type, quote = "") %>%
  mutate_all(as.character)
tibble_clinannot <- bind_rows(tibble_clinannot, tibble_buffer, .id = NULL)
}

tibble_clinannot_all <- full_join(corresp_manifestUUID_barcode, tibble_clinannot, by ="bcr_patient_barcode")
tibble_clinannot_all <- dplyr::select(tibble_clinannot_all, -2,-3,-4)                                                                #removes duplicate columns

tibble_expr_sample <- left_join(tibble_hgnc, full_sample_annot, by = "manifest_file_id")
tibble_expr_sample_clin <- left_join(tibble_expr_sample, tibble_clinannot_all, by = "manifest_file_id")
tibble_expr_sample_clin <-tibble_expr_sample_clin %>%
  mutate(patient_id = as.character(tumor_code)) %>%
  dplyr::filter(!is.na(bcr_patient_barcode)) %>%
  dplyr::filter(!is.na(tumor_code)) %>%
  arrange(tissue_type, bcr_patient_barcode) %>%
  dplyr::select(bcr_patient_barcode, tissue_type, A1BG:ZZEF1,tumor_code,days_to_birth,gender,race_list,ethnicity,other_dx,history_of_neoadjuvant_treatment,person_neoplasm_cancer_status,vital_status,days_to_last_followup,days_to_death,stage_event_pathologic_stage,stage_event_tnm_categories,days_to_last_known_alive,tumor_tissue_site,histological_type,diagnosis_subtype,neoplasm_histologic_grade,anatomic_neoplasm_subdivision,days_to_initial_pathologic_diagnosis,age_at_initial_pathologic_diagnosis,year_of_initial_pathologic_diagnosis,metastatic_site_list,other_metastatic_site,breast_carcinoma_progesterone_receptor_status,breast_carcinoma_estrogen_receptor_status,human_papillomavirus_types,patient_death_reason,death_cause_text,treatment,microsatellite_instability,kras_mutation_found,braf_gene_analysis_result,loss_expression_of_mismatch_repair_proteins_by_ihc_results,h_pylori_infection,ldh1_mutation_found,p53_gene_analysis,egfr_amplication_status,hpv_status_by_p16_testing,hpv_status_by_ish_testing,days_from_date_of_initial_pathologic_diagnosis_to_date_of_birth,tumor_type,leukemia_french_american_british_morphology_code,cytogenetic_abnormalities,child_pugh_classification_grade,viral_hepatitis_serologies,kras_mutation_result,egfr_mutation_result,eml4_alk_translocation_result,breslow_depth_value,melanoma_ulceration_indicator,malignant_neoplasm_mitotic_count_rate,distant_metastasis_anatomic_site,percent_tumor_nuclei,percent_monocyte_infiltration,percent_normal_cells,percent_eosinophil_infiltration,percent_lymphocyte_infiltration,percent_neutrophil_infiltration,percent_necrosis,percent_granulocyte_infiltration,number_proliferating_cells,percent_stromal_cells,percent_inflam_infiltration,percent_tumor_cells)

lookup_table <- c("ACC"="Adrenocortical carcinoma (n=79)","BLCA"="Bladder carcinoma (n=414)","BRCA"="Breast carcinoma (n=1108)","CESC"="Cervical carcinoma (n=306)","CHOL"="Cholangiocarcinoma (n=36)","COAD"="Colon adenocarcinoma (n=477)","DLBC"="B-cell lymphoma (DLBCL) (n=48)","ESCA"="Esophageal carcinoma (n=160)","GBM"="Glioblastoma (n=169)","HNSC"="Head&neck carcinoma (n=502)","KICH"="Kidney chromophobe tumor (n=65)","KIRC"="Kidney renal cell carcinoma (n=539)","KIRP"="Kidney papillary carcinoma (n=289)","LGG"="Low grade glioma (n=1055)","LIHC"="Liver hepatocellular carcinoma (n=368)","LUAD"="Lung adenocarcinoma (n=535)","LUSC"="Lung squamous carcinoma (n=502)","MESO"="Mesothelioma (n=86)","OV"="Ovarian carcinoma (n=379)","PAAD"="Pancreas carcinoma (n=178)","PCPG"="Pheochromocytoma paraganglioma (n=183)","PRAD"="Prostate adenocarcinoma (n=496)","READ"="Rectal adenocarcinoma (n=167)","SARC"="Sarcoma (n=263)","SKCM"="Skin melanoma (n=471)","STAD"="Stomach adenocarcinoma (n=375)","TGCT"="Testicular germ cell tumor (n=156)","THCA"="Thyroid carcinoma (n=509)","THYM"="Thymoma (n=119)","UCEC"="Uterine endometrial carcinoma (n=550)","UCS"="Uterine carcinosarcoma (n=56)","UVM"="Uveal melanoma (n=78)")
tibble_expr_sample_clin$tumor_type <- lookup_table[tibble_expr_sample_clin$tumor_code]

write.table(tibble_expr_sample_clin, file = "expr_data_select_clin_annot_FPKM-UQ.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)