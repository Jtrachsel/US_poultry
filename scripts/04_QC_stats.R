library(tidyverse)
library(phyloseq)

track <- read_tsv('./processed_data/QC_stats_reads_per_sample.tsv') %>%
  filter(sample_ID != 'Undetermined') %>%
  mutate(sample_type=case_when(
    grepl('mock', sample_ID) ~ 'mock',
    grepl('neg', sample_ID)  ~ 'neg',
    TRUE                     ~ 'exp'
  ))

sort(track$nonchim)
track$sample_type

### THIS ONE ###

P_reads_histogram <-
  track %>%
  ggplot(aes(nonchim,fill=sample_type)) +
  geom_histogram(color='black', bins = 50) +
  ylab('number of samples') +
  xlab('number of reads passing QC') +
  xlim(0,150000) +
  theme_bw()
P_reads_histogram
# track %>% filter(nonchim)
ggsave('figures/01_QC_hist_prefilt.jpeg', bg='white')

# tax <- read_tsv('./processed_data/ASV_taxonomy.tsv')

# ASV_counts <- read_tsv('processed_data/ASV_counts.tsv')

ps <- readRDS('./processed_data/phyloseq_decontam.rds')

ASV_long <-
  ps %>%
  # rarefy_even_depth() %>%
  # transform_sample_counts(fun = function(x) x / sum(x)) %>%
  psmelt()

# ASV_long2 <-
#   ASV_counts %>%
#   pivot_longer(cols = -sample_ID,names_to = 'ASV_seq', values_to = 'count') %>%
#   left_join(tax) %>%
#   mutate(sample_type=case_when(
#     grepl('mock', sample_ID) ~ 'mock',
#     grepl('neg', sample_ID)  ~ 'neg',
#     TRUE                     ~ 'exp'
#   ))  %>%
#   group_by(sample_ID) %>%
#   mutate(percent_abund=count/sum(count)*100) %>%
#   ungroup()

# mean_proportion_per_sample
# for each ASV within each sample:
  # calculate what proportion of the sampletotal that ASV occupies
  # ASV_count / sample_total_count
# ASV_long2 %>%
#   filter(percent_abund > 1) %>%
#   ggplot(aes(x=sample_type, y=percent_abund)) + geom_point()



# ALL_ASV_summary <-
#   ASV_long %>%
#   group_by(OTU,sample_type) %>%
#   summarise(count_in_type=sum(Abundance),
#             num_samples_in_type=sum(Abundance > 0),
#             av_abund_in_type=mean(percent_abund),
#             med_abund_in_type = median(percent_abund),
#             .groups = 'drop') %>%
#   arrange(desc(med_abund*num_samples)) %>%
#   mutate(ASV_ID=paste0('ASV_', 1:n()),
#          overall_rank= 1:n()) %>%
#   group_by(sample_type) %>%
#   mutate(within_type_rank=1:n())
#
# look <-
#   ALL_ASV_summary %>% arrange((within_type_rank)) %>%
#   filter(sample_type == 'neg') %>%
#   slice_head(n = 100)





# one for the asvs that occur in each sample type
# ASV_long
# ASV_long %>% arrange((count))
# ASV_long$Abundance
# ASV_long %>%
#   filter(Abundance>0) %>%
#   group_by(sample_ID, Challenge) %>%
#   summarise(num_ASVs=n()) %>%
#   ggplot(aes(x=Challenge, y=num_ASVs)) + geom_boxplot() +
#   # scale_y_log10() +
#   geom_jitter()

#
# ASV_long %>%
#   filter(Abundance>0) %>%
#   group_by(sample_ID, Challenge) %>%
#   summarise(num_ASVs=n()) %>%
#   ggplot(aes(x=Challenge, y=num_ASVs)) + geom_boxplot() +
#   # scale_y_log10() +
#   geom_jitter()
#
#

#
# mock_ASVs <- ASV_long %>% filter(sample_type == 'mock' & Abundance > 0) %>% pull(OTU) %>% unique()
# neg_ASVs <- ASV_long %>% filter(sample_type == 'neg' & Abundance > 0) %>% pull(OTU) %>% unique()
# exp_ASVs <- ASV_long %>% filter(sample_type == 'exp' & Abundance > 0) %>% pull(OTU) %>% unique()


#
# ALL_ASV_summary %>%
#   filter(ASV_seq  %in% neg_ASVs) %>%
#   filter(tot_count > 0) %>%
#   filter(num_samples >1) %>%
#   ggplot(aes(x=ASV_rank, y=(med_abund))) +
#   geom_point(aes(color=sample_type))
#
#
#



# problems
#
# MOCK_ASV_summary <-
#   ASV_long %>%
#   filter(OTU %in% mock_ASVs) %>%
#   group_by(OTU,sample_type) %>%
#   summarise(tot_count=sum(Abundance),
#             num_samples=sum(Abundance > 0),
#             av_abund=mean(percent_abund),
#             med_abund = median(percent_abund))# %>%
# # arrange(desc(tot_count))



# PRINT THIS SUMMARY
# unclass_ASV_summary <-
#   ASV_long %>%
#   filter(is.na(domain)) %>%
#   filter(Abundance != 0) %>%
#   group_by(OTU, ASV_seq) %>%
#   summarise(tot_counts=sum(Abundance),
#             num_samples=n(),
#             av_per_sample=mean(Abundance),
#             med_per_sample=median(Abundance)) %>%
#   arrange(desc(num_samples))
#
#
# unclass_ASV_summary %>% pull(ASV_seq)

# unclass_ASV_summary %>% slice_head() %>% pull(ASV)

# bad_unclassified_ASVs <- unclass_ASV_summary %>% filter(num_samples < 20) %>% pull(OTU)


# ps <- prune_taxa(!(taxa_names(ps) %in% bad_unclassified_ASVs), x = ps)
########


#
# mock_melt <-
#   ps_clean %>%
#   prune_samples(x = ., samples = sample_data(.)$sample_type == 'mock') %>%
#   prune_taxa(x=., taxa = taxa_sums(.) > 0) %>%
#   transform_sample_counts(function(x){x/sum(x)}) %>%
#   psmelt()


#
# mock_melt %>% ggplot(aes(x=Sample, y=(Abundance))) + geom_col(aes(fill=phylum), color='black')
#

# mock_seqs <-
#   mock_melt %>%
#   group_by(OTU, phylum, class, order, family, genus) %>%
#   summarise(tot_abund=sum(Abundance)) %>%
#   arrange(desc(tot_abund)) %>%
#   pull(OTU)


# mock_seqs[2]

########### pre-post decontam compare


ps <- read_rds('./processed_data/phyloseq_raw.rds')


otu_mat_pre <- otu_table(ps) %>% as(., Class='matrix')
sam_dat_pre <- sample_data(ps) %>% as(., Class = 'data.frame')

library(funfuns)

NMDSpre <- NMDS_ellipse(sam_dat_pre, otu_mat_pre, grouping_set = 'sample_type', MDS_trymax = 50)

# plot(NMDS[[3]])

# THIS ONE LOW READ SAMPLES MORE SIMILAR TO NTCS
NMDSpre[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(fill=sample_type, size=Read_depth), color='black', shape=21)

ggsave('figures/02_NMDS_pre_filter.jpeg', bg='white')

########

ps_post <- read_rds('./processed_data/phyloseq_decontam.rds')


otu_mat_post <- otu_table(ps_post) %>% as(., Class='matrix')
sam_dat_post <- sample_data(ps_post) %>% as(., Class = 'data.frame')


NMDSpost <- NMDS_ellipse(sam_dat_post, otu_mat_post, grouping_set = 'sample_type', MDS_trymax = 50)


NMDSpost[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(fill=sample_type, size=Read_depth), color='white', shape=21)


phyloseq::ntaxa(ps)
phyloseq::ntaxa(ps_post)

min(sample_sums(ps_post))
###### mocks removed #####




ps_post <- ps_post %>% prune_samples(samples = ps_post@sam_data$sample_type != 'mock') %>%
  rarefy_even_depth(rngseed = 7)





otu_mat_post <- otu_table(ps_post) %>% as(., Class='matrix')
sam_dat_post <- sample_data(ps_post) %>% as(., Class = 'data.frame')


NMDSpost <- NMDS_ellipse(sam_dat_post, otu_mat_post, grouping_set = 'sample_type', MDS_trymax = 200)

### THIS ONE
NMDSpost[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(fill=sample_type, size=Read_depth), color='black', shape=21)

ggsave('figures/03_NMDS_post_filter.jpeg')





######

# merge colonization data
col_dat <-  read_tsv('processed_data/Infantis_colonization_clean.tsv')


ps <- read_rds('./processed_data/phyloseq_decontam.rds')



SAMDAT <-
  ps@sam_data %>%
  as(object = ., Class = 'matrix') %>%
  as.data.frame() %>%
  as_tibble()

SAMDAT <-
  SAMDAT %>%
  mutate(turkey_ID=sub('.*_([0-9]+$)','\\1',sample_ID)) %>%
  mutate(ID=paste0(turkey_ID, '_', Vaccine, '_', Challenge)) %>%
  select(-starts_with('CFU'))

col_dat <-
  col_dat %>%
  mutate(Vaccine=sub('mock', 'Mock', Vaccine),
         Vaccine=sub('AviPro', 'AVIPRO', Vaccine),
         Vaccine=sub('866', 'BBS866', Vaccine),
         ID=paste0(turkey_ID, '_', Vaccine, '_', 'Infantis')) %>%
  select(ID, CFU, log10_CFU)

SAMDAT <- SAMDAT %>% left_join(col_dat) %>% sample_data()
sample_names(SAMDAT) <- SAMDAT$sample_ID
sample_data(ps) <- SAMDAT

write_rds(ps, 'processed_data/phyloseq_final.rds')
#####
# NMDSpost[[1]] %>%
#   ggplot(aes(x=MDS1, y=MDS2)) +
#   geom_point(aes(fill=Vaccine), color='white', shape=21, size=3)

#
# NMDSpost[[1]] %>%
#   ggplot(aes(x=MDS1, y=MDS2)) +
#   geom_point(aes(fill=Challenge), color='white', shape=21, size=3) +
#   facet_wrap(~day_post_challenge)


#
# AUTHORS   Gao,K.
# TITLE     Time course response of ileal and faecal microbiota and metabolic
# profiles to in-feed antibiotics in growing pigs
#

# REFERENCE   1  (bases 1 to 253)
# AUTHORS   Zhou,Q., Zhu,Y., Zhang,Y., Wang,L. and Dong,Z.
# TITLE     Studies on Bacterial Diversities in Hot Springs at Quzhuomu in
# Tibet
# JOURNAL   Unpublished

# REFERENCE   1  (bases 1 to 253)
# AUTHORS   Zhou,Q., Zhu,Y., Zhang,Y., Wang,L. and Dong,Z.
# TITLE     Studies on Bacterial Diversities in Hot Springs at Quzhuomu in
# Tibet
# JOURNAL   Unpublished

# REFERENCE   1  (bases 1 to 407)
# AUTHORS   Gao,K.
# TITLE     Time course response of ileal and faecal microbiota and metabolic
# profiles to in-feed antibiotics in growing pigs
# JOURNAL   Unpublished

# REFERENCE   1  (bases 1 to 409)
# AUTHORS   Lowe,P.P., Gyongyosi,B., Satishchandran,A., Iracheta-Vellve,A.,
# Ambade,A., Kodys,K., Catalano,D., Ward,D.V. and Szabo,G.
# TITLE     Alcohol-Related Changes in the Intestinal Microbiome Influence
# Neutrophil Infiltration, Inflammation and Steatosis in Early
# Alcoholic Hepatitis in Mice
# JOURNAL   Unpublished

# REFERENCE   2  (bases 1 to 418)
# AUTHORS   Yang,Y.
# TITLE     Direct Submission
# JOURNAL   Submitted (12-FEB-2015) Animal Science and Technology, Nanjing
# Agricultural University, Weigang 1#, Nanjing, Jiangsu 210095, China
# COMMENT     Sequences were screened for chimeras by the submitter using mothur.
#
# REFERENCE   1  (bases 1 to 407)
# AUTHORS   An,C.
# TITLE     Direct Submission
# JOURNAL   Submitted (02-AUG-2017) Department of Plague, Brucellosis and
# Etiologic Biological Vector Control, Shaanxi Center for Disease
# Control and Prevention, Jiandong, xi'an, Shaanxi 710054, China

# REFERENCE   1  (bases 1 to 251)
# AUTHORS   Walsh,P.
# TITLE     Direct Submission
# JOURNAL   Submitted (24-JUL-2019) School Of Medicine, Trinity College Dublin,
# TTMI, St. James Campus, Dublin 8, Ireland
