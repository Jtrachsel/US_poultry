library(tidyverse)
library(phyloseq)



tax <-
  read_tsv('./processed_data/ASV_taxonomy.tsv') %>%
  mutate(ASV=paste0('ASV_', seq_along(ASV_seq)))


fix_domain <- function(tax){

  tax %>%
    mutate(domain = case_when(
      is.na(domain) ~ 'unclassified',
      TRUE ~ domain))

}


fix_level <- function(tax, level){
  # browser()
  position <- which(colnames(tax) == level)

  tax[[level]] <- ifelse(is.na(tax[[level]]), paste0(tax[[position-1]], '_unclassified'),tax[[level]])
  tax[[level]] <- sub('unclassified_unclassified', 'unclassified',  tax[[level]])
  return(tax)
}


tax <-
  tax %>%
  fix_domain() %>%
  fix_level( 'phylum') %>%
  fix_level('class') %>%
  fix_level('order') %>%
  fix_level('family') %>%
  fix_level('genus') %>%
  fix_level('species')






swip_swapper <- tax$ASV
names(swip_swapper) <- tax$ASV_seq

tax <-
  tax %>%
  column_to_rownames(var='ASV') %>%
  as.matrix()

rownames(tax)

ASV_counts <- read_tsv('processed_data/ASV_counts.tsv') %>%
  mutate(sample_ID=gsub('-','_', sample_ID)) %>%
  filter(sample_ID != 'Undetermined') %>%
  column_to_rownames(var = 'sample_ID') %>%
  as.matrix()

colnames(ASV_counts) <- swip_swapper[colnames(ASV_counts)]


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
    is_neg=ifelse(sample_type == 'neg', TRUE, FALSE),
    dpc=day_post_challenge) %>%
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





