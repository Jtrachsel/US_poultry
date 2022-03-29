library(tidyverse)
library(decontam)
library(funfuns)
library(phyloseq)

ps <- read_rds('./processed_data/phyloseq_raw.rds')



otu_mat_pre <- otu_table(ps) %>% as(., Class='matrix')
sam_dat_pre <- sample_data(ps) %>% as(., Class = 'data.frame')

library(funfuns)

NMDSpre <- NMDS_ellipse(sam_dat_pre, otu_mat_pre, grouping_set = 'sample_type', MDS_trymax = 100)

# plot(NMDS[[3]])

# THIS ONE LOW READ SAMPLES MORE SIMILAR TO NTCS
NMDSpre[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2)) +
  geom_point(aes(fill=sample_type, size=Read_depth), color='white', shape=21) +
  theme_bw()



# these 4 samples are evident in the nmds with the mocks and negs
# they have the lowest number of reads, so makes sense that they would be most contaminated
# probably pretty much equivalent to negative controls?
# next closest probably valid sample has 7327 reads
bad_samples <- c('Avipro_OO_D14_6',
                 'BBS866_OO_D21_10',
                 'Mreading_D21_3',
                 'MM_D21_6')





# remove bad samples, remove singletons
ps <- ps %>%
  prune_samples(x=., samples = !(sample_names(.) %in% bad_samples)) %>%
  prune_taxa(x=., taxa = taxa_sums(.) > 1)

#######

# sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is_neg")
table(contamdf.prev$contaminant)

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is_neg", threshold=0.5)
table(contamdf.prev05$contaminant)

#plots
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_type %in% c("exp", "mock"), ps.pa)
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_type == "neg", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

good_ASVs <-
  contamdf.prev %>%
  filter(!contaminant) %>%
  rownames_to_column(var='ASV') %>%
  pull(ASV)





### look for closest negatives for each sample
# vegdist(otu_mat_pre) %>%
#   as.matrix() %>%
#   as.data.frame() %>%
#   rownames_to_column(var='from') %>%
#   pivot_longer(cols=-from, names_to = 'to', values_to = 'bray') %>%
#   filter(grepl('neg',from)) %>%
#   group_by(to) %>%
#   summarise(best_neg=from[which.min(bray)],
#             bray=bray[which.min(bray)]) %>% arrange(bray)
#





ps_clean <-
  prune_taxa(x = ps, taxa = (taxa_names(ps) %in% good_ASVs)) %>%
  prune_samples(x=., samples = !(sample_names(.) %in% bad_samples)) %>%
  prune_samples(x=., samples = sample_data(.)$sample_type != 'neg') %>%
  prune_samples(x=., samples = sample_data(.)$sample_type != 'mock') %>%
  prune_taxa(x=., taxa = taxa_sums(.) > 1)


sample_sums(ps_clean) %>% sort()



ASV_long <-
  ps %>%
  # rarefy_even_depth() %>%
  # transform_sample_counts(fun = function(x) x / sum(x)) %>%
  psmelt()

ASV_long %>%
  filter(Abundance>0) %>%
  group_by(sample_ID, Challenge) %>%
  summarise(num_ASVs=n()) %>%
  ggplot(aes(x=Challenge, y=num_ASVs)) + geom_boxplot() +
  # scale_y_log10() +
  geom_jitter()



ASV_long %>%
  filter(Abundance>0) %>%
  group_by(sample_ID, Challenge) %>%
  summarise(num_ASVs=n()) %>%
  ggplot(aes(x=Challenge, y=num_ASVs)) + geom_boxplot() +
  # scale_y_log10() +
  geom_jitter()



unclass_ASV_summary <-
  ASV_long %>%
  filter(is.na(domain)) %>%
  filter(Abundance != 0) %>%
  group_by(OTU) %>%
  summarise(tot_counts=sum(Abundance),
            num_samples=n(),
            av_per_sample=mean(Abundance),
            med_per_sample=median(Abundance)) %>%
  arrange(desc(num_samples))


bad_unclassified_ASVs <- unclass_ASV_summary %>% filter(num_samples < 20) %>% pull(OTU)


write_rds(ps_clean, './processed_data/phyloseq_decontam.rds')


