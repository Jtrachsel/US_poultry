library(tidyverse)
library(phyloseq)
library(cowplot)

#function to :
# 1) subset phyloseq object
# 2) make ordination (NMDS)
# 3) PERMANOVAs
# 4) calc diffabund OTUs
#  a) DESeq2
#  b) ANCOM-BC
# 5) calc diffs in alpha diversity

#########

ps <- read_rds('processed_data/phyloseq_decontam.rds')

tax <-
  ps@tax_table %>%
  as(., Class = 'matrix') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV')



#
# fix_domain <- function(tax){
#
#   tax %>%
#     mutate(domain = case_when(
#       is.na(domain) ~ 'unclassified',
#       TRUE ~ domain))
#
# }
#
#
# fix_level <- function(tax, level){
#   # browser()
#   position <- which(colnames(tax) == level)
#
#   tax[[level]] <- ifelse(is.na(tax[[level]]), paste0(tax[[position-1]], '_unclassified'),tax[[level]])
#   tax[[level]] <- sub('unclassified_unclassified', 'unclassified',  tax[[level]])
#   return(tax)
# }
#
#
# tax <-
#   tax %>%
#   fix_domain() %>%
#   fix_level( 'phylum') %>%
#   fix_level('class') %>%
#   fix_level('order') %>%
#   fix_level('family') %>%
#   fix_level('genus') %>%
#   fix_level('species')




sample_data(ps)$Vaccine <- factor(sample_data(ps)$Vaccine, levels = c('Mock', 'BBS866', 'AVIPRO'))


sample_data(ps)$Challenge <- factor(sample_data(ps)$Vaccine, levels = c('Mock', 'BBS866', 'AVIPRO'))



# A: vaccination effects on Turkey microbiome
# BBS866 vaccine
A_comp1 <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge == -1 &
                 ps@sam_data$Vaccine %in% c('Mock', 'BBS866')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)


NMDS_from_phyloseq <-
  function(ps, grouping_set){
    library(funfuns)
    ps <- ps %>% rarefy_even_depth()
    otu_mat <- otu_table(ps) %>% as(., Class='matrix')
    sam_dat <- sample_data(ps) %>% as(., Class = 'data.frame')

    if(length(grouping_set > 1)){

    }


    NMDS <- NMDS_ellipse(sam_dat, otu_mat, grouping_set = grouping_set)

    return(NMDS)
  }

# tst <- A_comp1 %>% NMDS_from_phyloseq('Vaccine')
#
# tst[[1]] %>%
#   ggplot(aes(x=MDS1, y=MDS2, fill=Vaccine)) +
#   geom_point(shape=21, color='white', size=4) +
#   geom_segment(aes(xend=centroidX, yend=centroidY, color=Vaccine), alpha=.5)+
#   theme_bw()
#
#



A_comp2 <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge == -1 &
                   ps@sam_data$Vaccine %in% c('Mock', 'AVIPRO')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)


#
# tst <- A_comp2 %>% NMDS_from_phyloseq('Vaccine')
#
# tst[[1]] %>%
#   ggplot(aes(x=MDS1, y=MDS2, fill=Vaccine)) +
#   geom_point(shape=21, color='white', size=4) +
#   geom_segment(aes(xend=centroidX, yend=centroidY, color=Vaccine), alpha=.5)+
#   theme_bw()
#

#### OK FINE IVE GOT TO SEE IF THE EFFECT IS THE SAME ###
A_comp3 <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge == -1) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)

##### THIS ONE

tst <- A_comp3 %>% NMDS_from_phyloseq('Vaccine')

tst[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=Vaccine)) +
  geom_point(shape=21, color='white', size=4) +
  geom_segment(aes(xend=centroidX, yend=centroidY, color=Vaccine), alpha=.5)+
  theme_bw()


A_mat <- A_comp3@otu_table %>% as(Class='matrix')
A_dat <- A_comp3@sam_data %>% as(Class='data.frame')

adonis(A_mat ~ Vaccine, data = A_dat)

pairwise.adonis(A_mat, factors = A_dat[['Vaccine']]) %>%
  ggplot(aes(y=pairs, x=F.Model)) + geom_col() +
  theme_half_open()


A_3_DE <- phyloseq_to_deseq2(A_comp3, design = ~ Vaccine)

A_3_DE <- A_3_DE %>% DESeq2::DESeq()

library(DESeq2)

resultsNames(A_3_DE)
# BBS866_diffabund <- results(A_3_DE, name ='Vaccine_BBS866_vs_Mock')
# AVIPRO_diffabund <- results(A_3_DE, name ='Vaccine_AVIPRO_vs_Mock')

BBS866_diffabund <-
  lfcShrink(A_3_DE, coef ='Vaccine_BBS866_vs_Mock') %>%
  as(Class = 'data.frame') %>%
  rownames_to_column(var='ASV') %>%
  rownames_to_column() %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > .5) %>%
  left_join(tax) %>%
  mutate(vaccine='BBS866')


AVIPRO_diffabund <-
  lfcShrink(A_3_DE, coef ='Vaccine_AVIPRO_vs_Mock') %>%
  as(Class = 'data.frame') %>%
  rownames_to_column(var='ASV') %>%
  rownames_to_column() %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > .5) %>%
  left_join(tax) %>%
  mutate(vaccine='AVIPRO')


both_diff <- bind_rows(AVIPRO_diffabund, BBS866_diffabund)



both_diff %>% filter(ASV %in% in_both) %>%
  ggplot(aes(x=log2FoldChange, y=genus, fill=vaccine)) + geom_point(shape=21) +
  geom_vline(xintercept = 0)



AVIPRO_diffabund %>% ggplot(aes(x=log2FoldChange, y=genus)) + geom_point()

BBS866_diffabund %>% ggplot(aes(x=log2FoldChange, y=genus)) + geom_point()
#

in_both <- intersect(AVIPRO_diffabund$ASV, BBS866_diffabund$ASV)


# B. Infantis Challenge effect on turkey microbiome
# really only 1 timepoint here even though its 10 and 14 in the metadata

B <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge %in% c(10, 14) &
                   ps@sam_data$Vaccine %in% c('Mock') &
                   ps@sam_data$Challenge %in% c('Infantis', 'Mock')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)


###



# C. BBS 866 vaccination effect on Infantis-challenged turkey microbiome

#
# 13 mock-vaccinated/Infantis-challenged turkeys at 10dpi:    Minfantis_D10 #1-13
# VERSUS
# 14 BBS866-vaccinated/Infantis-challenged turkeys at 10dpi:    BBS866_OO_D10 #1-14
#

C <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge %in% c(10) &
                   ps@sam_data$Vaccine %in% c('Mock', 'BBS866') &
                   ps@sam_data$Challenge %in% c('Infantis')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)


###
#
# D. AviPro vaccination effect on Infantis-challenged turkey microbiome
# Samples: Mock-vaccinated/Infantis-challenged and AviPro-vaccinated/Infantis-challenged
# Question: Did AviPro vaccination modify the Infantis-induced alterations of the microbial composition of the turkey cecum (cecal contents)
#
# 13 mock-vaccinated/Infantis-challenged turkeys at 10dpi:    Minfantis_D10 #1-13
# VERSUS
# 14 -AviPro-vaccinated/Infantis-challenged turkeys at 10dpi:    Avipro_OO_D10 #1-14
#


D <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge %in% c(10) &
                   ps@sam_data$Vaccine %in% c('Mock', 'AVIPRO') &
                   ps@sam_data$Challenge %in% c('Infantis')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)




######

#
# E. Reading Challenge effect on turkey microbiome
# Samples: Mock-vaccinated/mock-challenged and mock-vaccinated/Reading-challenged @ 2,7, 14, 21 dpi
# Question: did challenge with Salmonella Reading alter the microbial composition of the turkey cecum (cecal contents)
#
# Comparison #1 @ 2dpi:
# 12 Mock-vaccinated/mock-challenged turkeys:        MM_D2 #1-12
# VERSUS
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D2 #1-12
#
# Comparison #2 @ 7dpi:
# 12 Mock-vaccinated/mock-challenged turkeys:        MM_D7 #1-12
# VERSUS
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D7 #1-12
#
# Comparison #3 @ 14dpi:
# 12 Mock-vaccinated/mock-challenged turkeys:        MM_D14 #1-12
# VERSUS
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D14 #1-12
#
# Comparison #4 @ 21dpi:
# 16 Mock-vaccinated/mock-challenged turkeys:        MM_D21 #1-16
# VERSUS
# 16 mock-vaccinated/Reading-challenged turkeys:    Mreading_D21 #1-16
#



E <-
  ps %>%
  subset_samples(ps@sam_data$Vaccine %in% c('Mock') &
                 ps@sam_data$Challenge %in% c('Mock', 'Reading')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)


#####


F_ <-
  ps %>%
  subset_samples(ps@sam_data$Vaccine %in% c('BBS866', 'Mock') &
                   ps@sam_data$Challenge %in% c('Reading')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)




######
G <-
  ps %>%
  subset_samples(ps@sam_data$Vaccine %in% c('AVIPRO', 'Mock') &
                   ps@sam_data$Challenge %in% c('Reading')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)



##########

# can use this to look at differences between the 2 vaccinations and the mocks
A_comp3 <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge == -1) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)

