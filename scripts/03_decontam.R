library(decontam)
library(funfuns)


ps <- read_rds('./processed_data/phyloseq_object.rds')

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
  contamdf.prev05 %>%
  filter(!contaminant) %>%
  rownames_to_column(var='ASV_seq') %>%
  pull(ASV_seq)





### look for closest negatives for each sample
vegdist(otu_mat_pre) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column(var='from') %>%
  pivot_longer(cols=-from, names_to = 'to', values_to = 'bray') %>%
  filter(grepl('neg',from)) %>%
  group_by(to) %>%
  summarise(best_neg=from[which.min(bray)],
            bray=bray[which.min(bray)]) %>% arrange(bray)



ps_clean <-
  prune_taxa(x = ps, taxa = (taxa_names(ps) %in% good_ASVs)) %>%
  prune_samples(x=., samples = !(sample_names(.) %in% bad_samples)) %>%
  prune_samples(x=., samples = sample_data(.)$sample_type != 'neg') %>%
  prune_taxa(x=., taxa = taxa_sums(.) > 1)


sample_sums(ps_clean) %>% sort()

write_rds(ps_clean, './processed_data/phyloseq_decontam.rds')
