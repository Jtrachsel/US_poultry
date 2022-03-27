library(tidyverse)
library(phyloseq)



tax <- read_tsv('./processed_data/ASV_taxonomy.tsv') %>%
  column_to_rownames(var='ASV_seq') %>%
  as.matrix()

ASV_counts <- read_tsv('processed_data/ASV_counts.tsv') %>%
  mutate(sample_ID=gsub('-','_', sample_ID)) %>%
  filter(sample_ID != 'Undetermined') %>%
  column_to_rownames(var = 'sample_ID') %>%
  as.matrix()


samdat <-
  read_csv('./US_poultry2021_metadata.csv') %>%
  mutate(sample_ID=gsub('-','_',SampleID)) %>%
  select(sample_ID, everything(), -SampleID) %>%
  mutate(
    sample_type=case_when(
                 grepl('mock_plate', sample_ID)   ~ 'mock',
                 grepl('neg_', sample_ID)   ~ 'neg',
                 TRUE ~ 'exp'),
    Challenge=ifelse(sample_type == 'exp', Challenge, 'none'),
    Vaccine  =ifelse(sample_type == 'exp', Vaccine, 'none'),
    Vaccine_Booster_Route=ifelse(sample_type == 'exp', Vaccine_Booster_Route, 'none'),
    is_neg=ifelse(sample_type == 'neg', TRUE, FALSE)) %>%
  sample_data()

sample_names(samdat) <- samdat$sample_ID


all(rownames(ASV_counts) %in% samdat$sample_ID)

exp_obj <-
  phyloseq(sample_data(samdat),
           otu_table(ASV_counts, taxa_are_rows = F),
           tax_table(tax))


# phyloseq::sample_sums(exp_obj)

sample_data(exp_obj)$Read_depth <- sample_sums(exp_obj)

saveRDS(exp_obj, file = './processed_data/phyloseq_raw.rds')





