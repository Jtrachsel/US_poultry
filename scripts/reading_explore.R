library(tidyverse)
library(phyloseq)
library(cowplot)
library(DESeq2)
library(ANCOMBC)
library(funfuns)

source('scripts/functions.R')

ps <- read_rds('processed_data/phyloseq_final.rds')

tax <-
  ps@tax_table %>%
  as(., Class = 'matrix') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV')


# set reference levels
sample_data(ps)$Vaccine <- factor(sample_data(ps)$Vaccine, levels = c('Mock', 'BBS866', 'AVIPRO'))
sample_data(ps)$Challenge <- factor(sample_data(ps)$Challenge, levels = c('Mock', 'Infantis', 'Reading'))

##### tax glom
# ps@tax_table <- ps@tax_table[,-1] # drop ASVs so glomming works
# ps_genus <- tax_glom(ps, taxrank = rank_names(ps)[6]) # genus
# ps_fam <- tax_glom(ps, taxrank = rank_names(ps)[5]) # family

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
  # ps_genus %>%
  ps %>%
  subset_samples(ps@sam_data$Vaccine %in% c('Mock') &
                   ps@sam_data$Challenge %in% c('Mock', 'Reading')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  prune_taxa(taxa=colSums(.@otu_table > 0) > 5) # taxa detected in more than 5 samples


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
  geom_point(data=tmp[[1]], aes(x=MDS1, y=MDS2, color=Challenge), alpha=.5) +
  geom_segment(data=tmp[[1]], aes(x=MDS1, y=MDS2,xend=centroidX, yend=centroidY, color=Challenge), alpha=.25) +
  geom_path(aes(group=Challenge))+
  geom_point(size=6) +
  geom_text(aes(label=day_post_challenge), color='black') +
  ggtitle('NMDS of Mock-Mock vs Mock-Reading') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock = 'green', Reading = 'skyblue')) +
  theme_bw()

system('mkdir -p Reading_figures')
ggsave('Reading_figures/01_E_NMDS.jpeg', width = 7, height = 4, units = 'in', bg='white')
# hard to say what the 'challenge' effect is with this subset, because over the course of the
# challenge both groups change quite a bit.
# mock mock changes a lot
# mock Reading changes a lot
# Time effect?





# Should try lrt tests here

# THESE OTUS ARE ONES THAT DONT CHANGE IN THE SAME WAY OVER TIME
# OTUS THAT CHANGE OVER TIME ARE IGNORED == TIME EFFECT
# OTUS THAT ARE ONLYDIFFERENT BETWEEN THE GROUPS ARE IGNORED == ROOM EFFECT +  CHALLENGE EFFECT



E@sam_data$dpc <- factor(E@sam_data$dpc, levels=c('-1', '2', '7', '14', '21'))

# THIS IDS OTUS THAT BEHAVE DIFFERENTLY OVER TIME BETWEEN THE TWO GROUPS
dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge + Challenge:day_post_challenge)

# dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + dpc + Challenge:dpc)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Challenge + dpc)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Challenge + day_post_challenge)


# # this IDs OTUS that change over time
# dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Challenge)
# #
#
# # this IDS OTUS that differ between groups
# dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ day_post_challenge)
#


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

# clusters <- degPatterns(rlog_mat, metadata = data.frame(E@sam_data),minc=5, time = "dpc", col=NULL)

clusters_col <- degPatterns(rlog_mat, metadata = data.frame(E@sam_data),minc=5, time = "dpc", col='Challenge')

clusters_col$plot + scale_color_manual(values=c(Mock = 'green', Reading = 'skyblue')) +
  theme_bw()
ggsave('Reading_figures/02_E_LRT_OTUs.jpeg', width = 7, height = 4, units = 'in', bg='white')


E_interesting <-
  clusters_col$df %>%
  filter(cluster %in% c(6)) %>%
  transmute(ASV=genes, cluster = cluster) %>%
  left_join(tax)


E_LRT_ASVs <-
  E_interesting %>%
  arrange(cluster) %>%
  select(ASV, cluster, phylum, class,family)

#####
#
# F. BBS 866 vaccination effect on Reading-challenged turkey microbiome
# Samples: Mock-vaccinated/Reading-challenged and BBS866-vaccinated/Reading-challenged
# Question: Did BBS 866 vaccination modify the Reading-induced alterations of the microbial composition of the turkey cecum (cecal contents)
#
# Comparison #1 @ 2dpi:
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D2 #1-12
# VERSUS
# 12 BBS866-vaccinated/Reading-challenged turkeys:    BBS866_OO_D2 #1-12
#
# Comparison #2 @ 7dpi:
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D7 #1-12
# VERSUS
# 12 BBS866-vaccinated/Reading-challenged turkeys:    BBS866_OO_D7 #1-12
#
# Comparison #3 @ 14dpi:
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D14 #1-12
# VERSUS
# 12 BBS866-vaccinated/Reading-challenged turkeys:    BBS866_OO_D14 #1-12
#
# Comparison #4 @ 21dpi:
# 16 mock-vaccinated/Reading-challenged turkeys:    Mreading_D21 #1-16
# VERSUS
# 16 BBS866-vaccinated/Reading-challenged turkeys:    BBS866_OO_D21 #1-16


F_ <-
  ps %>%
  subset_samples(ps@sam_data$Vaccine %in% c('BBS866', 'Mock') &
                   ps@sam_data$Challenge %in% c('Reading')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  prune_taxa(taxa=colSums(.@otu_table > 0) > 5) # taxa detected in more than 5 samples



F_@sam_data$dpc <- F_@sam_data$day_post_challenge
F_@sam_data$set <- paste(F_@sam_data$day_post_challenge, F_@sam_data$Vaccine, sep = '_')

tmp <- NMDS_from_phyloseq(F_, 'set')

tmp[[2]] <-
  tmp[[2]] %>%
  mutate(day_post_challenge=sub('(-?[0-9]+)_([A-Za-z]+)','\\1',group),
         Vaccine=sub('(-?[0-9]+)_([A-Za-z]+)','\\2',group))

tmp[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, color=Vaccine)) +
  geom_segment(aes(xend=centroidX, yend=centroidY))+
  geom_point(aes(color=Vaccine)) +
  facet_wrap(~factor(day_post_challenge))

tmp[[1]] %>%
  select(centroidX, centroidY, Vaccine, day_post_challenge) %>%
  unique() %>%
  arrange(day_post_challenge) %>%
  ggplot(aes(x=centroidX, y=centroidY, color=Vaccine)) +
  geom_point(data=tmp[[1]], aes(x=MDS1, y=MDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=tmp[[1]], aes(x=MDS1, y=MDS2,xend=centroidX, yend=centroidY, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine))+
  geom_point(size=6) +
  geom_text(aes(label=day_post_challenge), color='black') +
  ggtitle('NMDS of Mock-Reading vs BBS866-Reading') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock = 'skyblue', BBS866 = 'orange')) +
  theme_bw()


ggsave('Reading_figures/03_F_NMDS.jpeg', width = 7, height = 4, units = 'in', bg='white')


#### INSERT LRT HERE ####
# Should try lrt tests here

# THESE OTUS ARE ONES THAT DONT CHANGE IN THE SAME WAY OVER TIME
# OTUS THAT CHANGE OVER TIME ARE IGNORED == TIME EFFECT
# OTUS THAT ARE ONLYDIFFERENT BETWEEN THE GROUPS ARE IGNORED == ROOM EFFECT +  CHALLENGE EFFECT



F_@sam_data$dpc <- factor(F_@sam_data$dpc, levels=c('-1', '2', '7', '14', '21'))

# THIS IDS OTUS THAT BEHAVE DIFFERENTLY OVER TIME BETWEEN THE TWO GROUPS
dds <- phyloseq_to_deseq2(physeq = F_, design = ~ Vaccine + day_post_challenge + Vaccine:day_post_challenge)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Vaccine + day_post_challenge)


# # this IDs OTUS that change over time
# dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Challenge)
# #
#
# # this IDS OTUS that differ between groups
# dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ day_post_challenge)
#


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
# a bit less than 10 mins
if (!file.exists('F_rlog.rds')){
  rld <- rlog(dds, blind=FALSE)
  write_rds(rld, 'F_rlog.rds')
} else {
  rld <- read_rds('F_rlog.rds')
}

library(DEGreport)


cluster_rlog <- rld[clustering_sig_genes$ASV, ]
rlog_mat <- as.matrix(cluster_rlog@assays@data[[1]])

# clusters <- degPatterns(rlog_mat, metadata = data.frame(E@sam_data),minc=5, time = "dpc", col=NULL)

clusters_col <- degPatterns(rlog_mat, metadata = data.frame(F_@sam_data),minc=5, time = "dpc", col='Vaccine')
clusters_col$plot +scale_color_manual(values=c(Mock = 'skyblue', BBS866 = 'orange')) +
  theme_bw()

ggsave('Reading_figures/04_F_LRT_OTUs.jpeg', width = 8, height = 4, units = 'in', bg='white')


F_interesting <-
  clusters_col$df %>%
  filter(cluster %in% c(5)) %>%
  transmute(ASV=genes, cluster = cluster) %>%
  left_join(tax)


F_LRT_ASVs <-
  F_interesting %>%
  arrange(cluster) %>%
  select(ASV, cluster, phylum, class,family)


# clusters_col$df %>%
# filter(cluster == 2) %>%
# transmute(ASV=genes, cluster = cluster) %>%
# left_join(tax)

#
# look <- clusters_col[['raw']]
# look$Challenge %>% unique()

#
# look %>%
#   filter(cluster == 1) %>%
#   ggplot(aes(x=dpc, y=value, color=Challenge))+
#   geom_point()


#### END LRT INSERT ####

######
# G. AviPro vaccination effect on Reading-challenged turkey microbiome
# Samples: Mock-vaccinated/Reading-challenged and AviPro-vaccinated/Reading-challenged
# Question: Did AviPro vaccination modify the Reading-induced alterations of the microbial composition of the turkey cecum (cecal contents)

# Comparison #1 @ 2dpi:
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D2 #1-12
# VERSUS
# 12 AviPro -vaccinated/Reading-challenged turkeys:    Avipro_OO_D2 #1-12
#
# Comparison #2 @ 7dpi:
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D7 #1-12
# VERSUS
# 12 AviPro -vaccinated/Reading-challenged turkeys:    Avipro_OO_D7 #1-12
#
# Comparison #3 @ 14dpi:
# 12 mock-vaccinated/Reading-challenged turkeys:    Mreading_D14 #1-12
# VERSUS
# 12 AviPro -vaccinated/Reading-challenged turkeys:    Avipro_OO_D14 #1-12
#
# Comparison #4 @ 21dpi:
# 16 mock-vaccinated/Reading-challenged turkeys:    Mreading_D21 #1-16
# VERSUS
# 13 AviPro -vaccinated/Reading-challenged turkeys:    Avipro_OO_D21 #1-13
#
# ***I don’t plan to necessarily compare the vaccines to one another, but instead create 2 manuscripts (if all the data allow). Although, it might be interesting to know, if both vaccines alter the microbiota, do they do it similarly; i.e.:
#   1. Compare in A: comparison #1 & #2
# 2. Compare F to G on the corresponding dpi





G <-
  ps %>%
  subset_samples(ps@sam_data$Vaccine %in% c('AVIPRO', 'Mock') &
                   ps@sam_data$Challenge %in% c('Reading')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  prune_taxa(taxa=colSums(.@otu_table > 0) > 5) # taxa detected in more than 5 samples



G@sam_data$dpc <- G@sam_data$day_post_challenge
G@sam_data$set <- paste(G@sam_data$day_post_challenge, G@sam_data$Vaccine, sep = '_')

tmp <- NMDS_from_phyloseq(G, 'set')

tmp[[2]] <-
  tmp[[2]] %>%
  mutate(day_post_challenge=sub('(-?[0-9]+)_([A-Za-z]+)','\\1',group),
         Vaccine=sub('(-?[0-9]+)_([A-Za-z]+)','\\2',group))

tmp[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, color=Vaccine)) +
  geom_segment(aes(xend=centroidX, yend=centroidY))+
  geom_point(aes(color=Vaccine)) +
  facet_wrap(~factor(day_post_challenge))

# tmp[[1]] %>%
#   select(centroidX, centroidY, Vaccine, day_post_challenge) %>%
#   unique() %>%
#   arrange(day_post_challenge) %>%
#   ggplot(aes(x=centroidX, y=centroidY, color=Vaccine)) +
#   geom_path(aes(group=Vaccine))+
#   geom_point(size=6) +
#   geom_text(aes(label=day_post_challenge), color='black')

tmp[[1]] %>%
  select(centroidX, centroidY, Vaccine, day_post_challenge) %>%
  unique() %>%
  arrange(day_post_challenge) %>%
  ggplot(aes(x=centroidX, y=centroidY, color=Vaccine)) +
  geom_point(data=tmp[[1]], aes(x=MDS1, y=MDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=tmp[[1]], aes(x=MDS1, y=MDS2,xend=centroidX, yend=centroidY, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine))+
  geom_point(size=6) +
  geom_text(aes(label=day_post_challenge), color='black') +
  ggtitle('NMDS of Mock-Reading vs AVIPRO-Reading') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock = 'skyblue', AVIPRO = 'purple')) +
  theme_bw()




ggsave('Reading_figures/05_G_NMDS.jpeg', width = 7, height = 4, units = 'in', bg='white')



#### INSERT LRT HERE ####
# Should try lrt tests here

# THESE OTUS ARE ONES THAT DONT CHANGE IN THE SAME WAY OVER TIME
# OTUS THAT CHANGE OVER TIME ARE IGNORED == TIME EFFECT
# OTUS THAT ARE ONLYDIFFERENT BETWEEN THE GROUPS ARE IGNORED == ROOM EFFECT +  CHALLENGE EFFECT



G@sam_data$dpc <- factor(G@sam_data$dpc, levels=c('-1', '2', '7', '14', '21'))

# THIS IDS OTUS THAT BEHAVE DIFFERENTLY OVER TIME BETWEEN THE TWO GROUPS
dds <- phyloseq_to_deseq2(physeq = G, design = ~ Vaccine + day_post_challenge + Vaccine:day_post_challenge)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Vaccine + day_post_challenge)


# # this IDs OTUS that change over time
# dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Challenge)
# #
#
# # this IDS OTUS that differ between groups
# dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ day_post_challenge)
#


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
# 10ish mins?
if (!file.exists('G_rlog.rds')){
  rld <- rlog(dds, blind=FALSE)
  write_rds(rld, 'G_rlog.rds')
} else {
  rld <- read_rds('G_rlog.rds')
}

library(DEGreport)


cluster_rlog <- rld[clustering_sig_genes$ASV, ]
rlog_mat <- as.matrix(cluster_rlog@assays@data[[1]])

# clusters <- degPatterns(rlog_mat, metadata = data.frame(E@sam_data),minc=5, time = "dpc", col=NULL)

clusters_col <- degPatterns(rlog_mat, metadata = data.frame(G@sam_data),minc=5, time = "dpc", col='Vaccine')

clusters_col$plot +scale_color_manual(values=c(Mock = 'skyblue', AVIPRO = 'purple')) +
  theme_bw()


ggsave('Reading_figures/06_G_LRT_OTUs.jpeg', width = 8, height = 4, units = 'in', bg='white')




G_interesting <-
  clusters_col$df %>%
  filter(cluster %in% c(10)) %>%
  transmute(ASV=genes, cluster = cluster) %>%
  left_join(tax)



G_LRT_ASVs <-
  G_interesting %>%
  arrange(cluster) %>%
  select(ASV, cluster, phylum, class,family)


#
# clusters_col$df %>%
#   filter(cluster == 2) %>%
#   transmute(ASV=genes, cluster = cluster) %>%
#   left_join(tax)
#
#
# look <- clusters_col[['raw']]
# look$Challenge %>% unique()
#
#
# look %>%
#   filter(cluster == 1) %>%
#   ggplot(aes(x=dpc, y=value, color=Challenge))+
#   geom_point()

#### END LRT INSERT HERE ####

### GOTTA DO THEM ALL TOGETHER


ps@sam_data$Vaccine %>% unique()
ps@sam_data$Challenge %>% unique() # c('Mock', 'Reading')
ps@sam_data$day_post_challenge

G@sam_data$day_post_challenge %>% unique()

EFG <-
  ps %>%
  subset_samples(ps@sam_data$Vaccine %in% c('BBS866','AVIPRO', 'Mock') &
                   ps@sam_data$Challenge %in% c('Mock','Reading')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  prune_taxa(taxa=colSums(.@otu_table > 0) > 5) # taxa detected in more than 5 samples


EFG




EFG@sam_data$dpc <- EFG@sam_data$day_post_challenge
EFG@sam_data$set <- paste(EFG@sam_data$day_post_challenge, EFG@sam_data$Vaccine,EFG@sam_data$Challenge, sep = '_')
EFG@sam_data$GROUP <- paste(EFG@sam_data$Vaccine,EFG@sam_data$Challenge, sep = '_')


tmp <- NMDS_from_phyloseq(EFG, 'set')

tmp[[2]] <-
  tmp[[2]] %>%
  mutate(day_post_challenge=sub('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+)','\\1',group),
         Vaccine=sub('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+)','\\2',group),
         Challenge=sub('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+)','\\3',group),
         GROUP=paste(Vaccine, Challenge, sep='_'))


### THIS ONE ###

tmp[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, color=GROUP)) +
  geom_segment(aes(xend=centroidX, yend=centroidY))+
  geom_point(aes(color=GROUP)) +
  facet_wrap(~factor(day_post_challenge))

tmp[[1]] %>%
  select(centroidX, centroidY, GROUP, day_post_challenge) %>%
  unique() %>%
  arrange(day_post_challenge) %>%
  ggplot(aes(x=centroidX, y=centroidY, color=GROUP)) +
  geom_path(aes(group=GROUP))+
  geom_point(size=6) +
  geom_text(aes(label=day_post_challenge), color='black') +
  ggtitle('NMDS of all groups over time') +
  ylab('MDS2')+
  xlab('MDS1') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_Reading = 'skyblue',
                              AVIPRO_Reading = 'purple',
                              BBS866_Reading = 'orange')) +
  theme_bw()



ggsave('Reading_figures/07_EFG_NMDS.jpeg', width = 7, height = 4, units = 'in', bg='white')
### PERMANOVA SECTION ##

library(pairwiseAdonis)


min(sample_sums(EFG))

# 10ish mins?
if (!file.exists('avg_bray.rds')){
  AVG_BRAY <- avgdist(otu_table(EFG)@.Data, min(sample_sums(EFG)))
  write_rds(AVG_BRAY, 'avg_bray.rds')
} else {
  AVG_BRAY <- read_rds('avg_bray.rds')
}


# rownames(sample_data(EFG))  == names(AVG_BRAY)


# 10ish mins?
if (!file.exists('all_pwadon.rds')){
  PW_ADON <- pairwise.adonis(x = AVG_BRAY, factors = sample_data(EFG)$set,perm =9999)
  write_rds(PW_ADON, 'all_pwadon.rds')
} else {
  PW_ADON <- read_rds('all_pwadon.rds')
}


PW_ADON %>%
  mutate(FROM=sub('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+) vs (-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+)','',pairs),
         TO=sub('','',pairs))

# WITHIN DAY COMPARISONS
within_day <- grepl('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+) vs (\\1)_([A-Za-z0-9]+)_([A-Za-z]+)', PW_ADON$pairs)

# WITHIN GROUP COMPARISONS
within_room <- grepl('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+) vs (-?[0-9]+)_(\\2)_(\\3)', PW_ADON$pairs)


## SET REFERENCE GROUP HERE, NOW ITS MOCK MOCK
pw_day_Mock_Mock <-
  PW_ADON[within_day,] %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',pairs),
         G2=sub('(.*) vs (.*)','\\2',pairs)) %>%
  filter(grepl('Mock_Mock',G1) |grepl('Mock_Mock',G2)) %>%
  mutate(REF=ifelse(grepl('Mock_Mock',G1), G1,G2),
         COMP=ifelse(grepl('Mock_Mock',G1), G2,G1),
         day=sub('(-?[0-9]+)_.*','\\1',REF)%>% as.numeric(),
         room=sub('(-?[0-9]+)_(.*)','\\2',COMP),
         fdr=p.adjust(p.value, method = 'fdr'))

# this is with Mock_Reading as reference
pw_day_Mock_Reading <-
  PW_ADON[within_day,] %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',pairs),
         G2=sub('(.*) vs (.*)','\\2',pairs)) %>%
  filter(grepl('Mock_Reading',G1) |grepl('Mock_Reading',G2)) %>%
  mutate(REF=ifelse(grepl('Mock_Reading',G1), G1,G2),
         COMP=ifelse(grepl('Mock_Reading',G1), G2,G1),
         day=sub('(-?[0-9]+)_.*','\\1',REF) %>% as.numeric(),
         room=sub('(-?[0-9]+)_(.*)','\\2',COMP),
         fdr=p.adjust(p.value, method = 'fdr'))
#
### GROUP COMPARISONS AT EACH DAY
pw_day_Mock_Reading %>%
  ggplot(aes(x=day, y=F.Model, color=room))+
  geom_point() +
  geom_line() +
  ylim(0,25) +
  ggtitle('Difference from Mock_Reading at each day', 'All comparisons within each day are significant')+
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_Reading = 'skyblue',
                              AVIPRO_Reading = 'purple',
                              BBS866_Reading = 'orange')) +
  theme_bw()
ggsave('Reading_figures/12_PERMANOVA_vs_MOCK_Reading.jpeg', width=7,height = 5, units = 'in', bg='white')


pw_day_Mock_Mock %>%
  ggplot(aes(x=day, y=F.Model, color=room)) +
  geom_point() +
  geom_line() +
  ylim(0,25) +
  ggtitle('Difference from Mock_Mock at each day', 'All comparisons within day are significant') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_Reading = 'skyblue',
                              AVIPRO_Reading = 'purple',
                              BBS866_Reading = 'orange')) +
  theme_bw()

ggsave('Reading_figures/13_PERMANOVA_vs_MOCK_MOCK.jpeg', width=7,height = 5, units = 'in', bg='white')


pw_day_all <- PW_ADON[within_day,] %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',pairs),
         G2=sub('(.*) vs (.*)','\\2',pairs)) %>%
  # filter(grepl('Mock_Reading',G1) |grepl('Mock_Reading',G2)) %>%
  mutate(REF=ifelse(grepl('Mock_Mock',G1), G1,G2),
         COMP=ifelse(grepl('Mock_Mock',G1), G2,G1),
         day=sub('(-?[0-9]+)_.*','\\1',G1) %>% as.numeric(),
         room1=sub('(-?[0-9]+)_(.*)','\\2',REF),
         room2=sub('(-?[0-9]+)_(.*)','\\2',COMP),
         fdr=p.adjust(p.value, method = 'fdr'),
         group=paste(room1, room2, sep=' vs '),
         group=factor(group, levels = c("Mock_Mock vs Mock_Reading",
                                        "Mock_Mock vs BBS866_Reading",
                                        "Mock_Mock vs AVIPRO_Reading",
                                        "Mock_Reading vs BBS866_Reading",
                                        "Mock_Reading vs AVIPRO_Reading",
                                        "BBS866_Reading vs AVIPRO_Reading")))




# pw_day_all %>%
#   ggplot(aes(x=F.Model, y=fct_rev(group))) +
#   geom_col() +
#   facet_wrap(~day, nrow = 1)
#


pw_day_all %>%
  ggplot(aes(x=F.Model, y=fct_rev(group))) +
  geom_col() +
  facet_wrap(~day, ncol = 1)

ggsave('Reading_figures/08_PERMANOVA_groups_at_each_day.jpeg',width=5, height=6, units='in', bg='white')



pw_day_all %>%
  filter(day == -1) %>%
  ggplot(aes(x=F.Model, y=fct_rev(group))) +
  geom_col() +
  geom_text(aes(label=paste0('P=',signif(fdr,2))),hjust=1 ) +
  facet_wrap(~day, ncol = 1)


#
pw_room <- PW_ADON[within_room,] %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',pairs),
         G2=sub('(.*) vs (.*)','\\2',pairs)) %>%
  filter(grepl('-1',G1) |grepl('-1',G2)) %>%
  mutate(REF=ifelse(grepl('-1',G1), G1,G2),
         COMP=ifelse(grepl('-1',G1), G2,G1),
         day=sub('(-?[0-9]+)_.*','\\1',COMP) %>% as.numeric(),
         room=sub('(-?[0-9]+)_(.*)','\\2',REF),
         fdr_pval=p.adjust(p.value, method = 'fdr'))

# pw_room %>% ggplot(aes(x=day, y=R2, color=room)) + geom_point() + geom_line()

d_m1tib <-
  tibble(room= pw_room %>% pull(room) %>% unique(),
         day=-1,
         fdr_pval=1,
         F.Model=1)


pw_room %>%
  bind_rows(d_m1tib) %>%
  mutate(sig=ifelse(fdr_pval < 0.05, TRUE, NA)) %>%
  ggplot(aes(x=day, y=F.Model, color=room))+
  geom_point(aes(shape=sig)) +
  geom_line() +
  xlim(-1,25) +
  ggtitle('Within group changes over time',
          'All groups compared back to themselves at D -1') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_Reading = 'skyblue',
                              AVIPRO_Reading = 'purple',
                              BBS866_Reading = 'orange')) +
  theme_bw()


ggsave('Reading_figures/09_within_group_changes_time.jpeg',width=7, height=5, units='in', bg='white')

library(vegan)


### CUMULATIVE DIFF FROM PREVIOUS DAY

pw_room_cumulative <- PW_ADON[within_room,] %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',pairs),
         G2=sub('(.*) vs (.*)','\\2',pairs)) %>%
  # filter(grepl('-1',G1) |grepl('-1',G2)) %>%
  mutate(REF=ifelse(grepl('-1',G1), G1,G2),
         COMP=ifelse(grepl('-1',G1), G2,G1),
         d1=sub('(-?[0-9]+)_.*','\\1',REF) %>% as.numeric(),
         d2=sub('(-?[0-9]+)_.*','\\1',COMP) %>% as.numeric(),
         room=sub('(-?[0-9]+)_(.*)','\\2',REF),
         day1=ifelse(d1<d2, d1, d2),
         day2=ifelse(d1<d2, d2, d1),
         fdr_pval=p.adjust(p.value, method = 'fdr')) %>%
  filter(day1 == -1 & day2 == 2  |
           day1 == 2  & day2 == 7  |
           day1 == 7  & day2 == 14 |
           day1 == 14 & day2 == 21) %>%
  arrange(day1) %>%
  group_by(room) %>%
  mutate(cumulative_F.model=cumsum(F.Model)) %>%
  select(pairs, F.Model, R2, p.value, room, day1, day2, fdr_pval, cumulative_F.model) %>%
  mutate(period=paste0(day1, ' vs ', day2),
         period=factor(period, levels = c('-1 vs 2','2 vs 7','7 vs 14','14 vs 21')))



pw_room_cumulative %>%
  ggplot(aes(x=period, y=cumulative_F.model, color=room, group=room)) +
  geom_line() + geom_point() +
  ggtitle('Cumulative change in communities within a room') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_Reading = 'skyblue',
                              AVIPRO_Reading = 'purple',
                              BBS866_Reading = 'orange')) +
  theme_bw()


ggsave('Reading_figures/10_cumulative_change_within_room.jpeg',width=7, height=5, units='in', bg='white')

### Write out PERMANOVA tables ###


### INCLUDE THIS

PERMANOVA_WITHIN_DAY <-
  pw_day_all %>%
  select(pairs, F.Model, R2, p.value, fdr) %>%
  mutate(across(where(is.numeric), signif, digits=2)) %>%
  write_tsv('output/Reading_within_day_PERMANOVAs.tsv')




### INCLUDE THIS

PERMANOVA_WITHIN_ROOM <- pw_room_tab <-
  pw_room %>%
  select(pairs, F.Model, R2, p.value, fdr_pval) %>%
  mutate(across(where(is.numeric), signif, digits=2)) %>%
  write_tsv('output/Reading_within_room_PERMANOVAs.tsv')




# EFG ALPHA #


EFG_rare <- rarefy_even_depth(EFG)

EFG_rare@sam_data$shannon <- diversity(EFG_rare@otu_table)

EFG_rare@sam_data %>%
  ggplot(aes(x=GROUP, y=shannon, fill=GROUP, group=set)) +
  geom_boxplot() +
  facet_wrap(~day_post_challenge, nrow=1)

mod_dat <-
  EFG_rare@sam_data %>%
  data.frame() %>%
  mutate(GROUP=factor(GROUP, levels=c('Mock_Mock', 'Mock_Reading', 'BBS866_Reading', 'AVIPRO_Reading')),
         day_post_challenge = factor(day_post_challenge, levels = c('-1', '2', '7','14', '21')))

colnames(EFG_rare@sam_data)

shannon_mod <- lm(data = mod_dat, formula= shannon ~ day_post_challenge * GROUP)

shannon_mod %>% summary()
library(emmeans)


EMMEANS <- emmeans(shannon_mod, ~ GROUP | day_post_challenge)

shan_contrasts <- emmeans::contrast(EMMEANS, method='pairwise', adjust='fdr')

contrasts_table <- shan_contrasts %>% as.data.frame()

plot_dat <- EMMEANS %>% as.data.frame()


## THIS ONE
plot_dat %>%
  ggplot(aes(x=day_post_challenge,
             y=emmean,ymin=lower.CL,
             ymax=upper.CL,
             group=GROUP, color=GROUP)) +
  geom_pointrange(position = position_dodge(width = .2)) +
  geom_line() + ggtitle('Alpha diversity over time') +
  ylab('Shannon') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_Reading = 'skyblue',
                              AVIPRO_Reading = 'purple',
                              BBS866_Reading = 'orange')) +
  theme_bw()

ggsave('Reading_figures/11_alpha_div.jpeg',width=7, height=5, units='in', bg='white')

# INCLUDE THIS
alpha_contrasts_table <-
  contrasts_table %>%
  mutate(across(where(is.numeric), signif, digits=2),
         sig=ifelse(p.value < 0.05, TRUE, FALSE)) %>%
  write_tsv('output/Reading_alpha_tests.tsv')


###

### Beta dispersion section ###


BETA_DISPER <- betadisper(AVG_BRAY, group = sample_data(EFG)$set)

EFG_rare@sam_data$group_disper <- BETA_DISPER$distances


EFG_rare@sam_data %>% ggplot(aes(x=GROUP, y=group_disper)) + geom_boxplot() +
  facet_wrap(~day_post_challenge, ncol=1)

### LRT Section ###

# EFG <- prune_samples(EFG@sam_data$GROUP != 'Mock_Mock', x = EFG)
#
# EFG@sam_data$dpc <- factor(EFG@sam_data$dpc, levels=c('-1', '2', '7', '14', '21'))
#
# # THIS IDS OTUS THAT BEHAVE DIFFERENTLY OVER TIME BETWEEN THE TWO GROUPS
# dds <- phyloseq_to_deseq2(physeq = EFG, design = ~ Vaccine + day_post_challenge + Vaccine:day_post_challenge)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Vaccine + day_post_challenge)
#
#
# # # this IDs OTUS that change over time
# # dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge)
# # dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Challenge)
# # #
# #
# # # this IDS OTUS that differ between groups
# # dds <- phyloseq_to_deseq2(physeq = E, design = ~ Challenge + day_post_challenge)
# # dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ day_post_challenge)
# #
#
#
# # for a LRT, lfc values depend on the contrast/coef but the p values do not?
# clustering_sig_genes <-
#   dds_lrt_time %>%
#   results() %>%
#   as.data.frame() %>%
#   rownames_to_column(var='ASV') %>%
#   filter(padj < 0.01) %>%
#   left_join(tax)
#
#
# # look for ASVs with similar patterns of abundance
#
# # rlog take a while...
# # 10ish mins?
# if (!file.exists('EFG_rlog.rds')){
#   rld <- rlog(dds, blind=FALSE)
#   write_rds(rld, 'EFG_rlog.rds')
# } else {
#   rld <- read_rds('EFG_rlog.rds')
# }
#
# library(DEGreport)
#
#
# cluster_rlog <- rld[clustering_sig_genes$ASV, ]
# rlog_mat <- as.matrix(cluster_rlog@assays@data[[1]])
#
# # clusters <- degPatterns(rlog_mat, metadata = data.frame(E@sam_data),minc=5, time = "dpc", col=NULL)
#
# clusters_col <- degPatterns(rlog_mat, metadata = data.frame(EFG@sam_data),minc=5, time = "dpc", col='Vaccine')
# #
# # EFG_interesting <-
# #   clusters_col$df %>%
# #   filter(cluster %in% c(10)) %>%
# #   transmute(ASV=genes, cluster = cluster) %>%
# #   left_join(tax)

bind_rows(E_LRT_ASVs %>% mutate(comp='E'),
          F_LRT_ASVs%>% mutate(comp='F'),
          G_LRT_ASVs%>% mutate(comp='G')) %>%
  group_by(family) %>% tally() %>%
  arrange(desc(n))


bind_rows(E_LRT_ASVs %>% mutate(comp='E'),
          F_LRT_ASVs%>% mutate(comp='F'),
          G_LRT_ASVs%>% mutate(comp='G')) %>%
  group_by(ASV) %>% tally() %>%
  arrange(desc(n))



