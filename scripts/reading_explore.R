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


E@sam_data$dpc <- E@sam_data$day_post_challenge
E@sam_data$set <- paste(E@sam_data$day_post_challenge, E@sam_data$Challenge, sep = '_')

tmp <- NMDS_from_phyloseq(E, 'set')

tmp[[2]] <-
  tmp[[2]] %>%
  mutate(day_post_challenge=sub('(-?[0-9]+)_([A-Za-z]+)','\\1',group),
         Challenge=sub('(-?[0-9]+)_([A-Za-z]+)','\\2',group))

tmp[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, color=Challenge)) +
  geom_segment(aes(xend=centroidX, yend=centroidY))+
  geom_point(aes(color=Challenge)) +
  facet_wrap(~factor(day_post_challenge))

tmp[[1]] %>%
  select(centroidX, centroidY, Challenge, day_post_challenge) %>%
  unique() %>%
  arrange(day_post_challenge) %>%
  ggplot(aes(x=centroidX, y=centroidY, color=Challenge)) +
  geom_path(aes(group=Challenge))+
  geom_point(size=6) +
  geom_text(aes(label=day_post_challenge), color='black')

# Should try lrt tests here
E@sam_data$dpc <- factor(E@sam_data$dpc, levels=c('-1', '2', '7', '14', '21'))

dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + dpc + Challenge:dpc)

dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Challenge + dpc)


# for a LRT, lfc values depend on the contrast/coef but the p values do not?
clustering_sig_genes <-
  dds_lrt_time %>%
  results() %>%
  as.data.frame() %>%
  rownames_to_column(var='ASV') %>%
  filter(padj < 0.01) %>%
  left_join(tax)


# look for ASVs with similar patterns of abundance

# rlog take a while...
if (!file.exists('E_rlog.rds')){
  rld <- rlog(dds, blind=FALSE)
  write_rds(rld, 'E_rlog.rds')
} else {
  rld <- read_rds('E_rlog.rds')
}

library(DEGreport)


cluster_rlog <- rld[clustering_sig_genes$ASV, ]
rlog_mat <- as.matrix(cluster_rlog@assays@data[[1]])

clusters <- degPatterns(rlog_mat, metadata = data.frame(E@sam_data),minc=5, time = "dpc", col=NULL)

clusters_col <- degPatterns(rlog_mat, metadata = data.frame(E@sam_data),minc=5, time = "dpc", col='Challenge')



look <- clusters_col[['raw']]
look$Challenge %>% unique()


look %>%
  filter(cluster == 1) %>%
  ggplot(aes(x=dpc, y=value, color=Challenge))+
  geom_point()
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
