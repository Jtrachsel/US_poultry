library(funfuns)
library(tidyverse)
library(phyloseq)
library(cowplot)
library(DESeq2)
library(ANCOMBC)






#function to :
# 1) subset phyloseq object
# 2) make ordination (NMDS)
# 3) PERMANOVAs
# 4) calc diffabund OTUs
#  a) DESeq2
#  b) ANCOM-BC
# 5) calc diffs in alpha diversity

DESeq_single_contrast <-
  function(DESeq_obj,
           COEF,
           padj_cut=0.05,
           l2fc_cut=.5,
           tax_tab,
           rowname_var='ASV'){
    # browser()
    res_names <- sub('_vs_','_', COEF) %>%
      strsplit(split = '_') %>% unlist()

    COL=res_names[1]

    res <-
      lfcShrink(DESeq_obj, coef=COEF) %>%
      as(Class = 'data.frame') %>%
      rownames_to_column(var=rowname_var) %>%
      filter(padj < 0.05) %>%
      filter(abs(log2FoldChange) > .5) %>%
      left_join(tax_tab) %>%
      mutate({{COL}}:=ifelse(log2FoldChange > 0, res_names[2], res_names[3]))

    return(res)
  }



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

make_difabund_plot <- function(diffabund_res, variable){
  # browser()
  diffabund_res %>%
    # filter(ASV %in% in_both) %>%
    ggplot(aes(x=log2FoldChange, y=genus, fill=!!sym(variable))) +
    geom_point(size=2.5, shape=21) +
    geom_vline(xintercept = 0)
}

explore <- function(ps, variable){
  # TODO
  # add adonis2 results object as part of return list
  # browser()

  tax <-
    ps@tax_table %>%
    as(., Class = 'matrix') %>%
    as.data.frame() %>%
    rownames_to_column(var = 'ASV')


  tst <- ps %>% NMDS_from_phyloseq(variable)

  p1 <- tst[[1]] %>%
    ggplot(aes(x=MDS1, y=MDS2, fill=!!sym(variable))) +
    geom_point(shape=21, color='white', size=4) +
    geom_segment(aes(xend=centroidX, yend=centroidY, color=!!sym(variable)), alpha=.5)+
    theme_bw()


  A_mat <- ps@otu_table %>% as(Class='matrix')
  A_dat <- ps@sam_data %>% as(Class='data.frame')


  # ADON <- adonis2(A_mat ~ A_dat[[variable]], data = A_dat)

  p2 <- pairwise.adonis(A_mat, factors = A_dat[[variable]]) %>%
    ggplot(aes(y=pairs, x=F.Model)) + geom_col() +
    theme_half_open()



  ### THIS SHOULD BE A FUNCTION




  A_3_DE <- phyloseq_to_deseq2(ps, design = as.formula(paste0('~', variable)))

  A_3_DE <- A_3_DE %>% DESeq2::DESeq()

  res_tibble <-
    tibble(comparison=resultsNames(A_3_DE)[-1]) %>%
    mutate(res=map(.x = comparison,
                   .f=~DESeq_single_contrast(DESeq_obj = A_3_DE, COEF = .x, tax_tab = tax)),
           difabund_plot=
             map(.x=res, .f=~make_difabund_plot(.x, variable)))
  # res_tibble$difabund_plot[[1]]


  # this could have more than 2, need a way to deal with that
  # resultsNames(A_3_DE)
  # ref_level <- levels(ps@sam_data[[variable]])[1]
  # other_level <- levels(ps@sam_data[[variable]])[2]
  # #
  # # BBS866_diffabund <- results(A_3_DE, name ='Vaccine_BBS866_vs_Mock')
  # # AVIPRO_diffabund <- results(A_3_DE, name ='Vaccine_AVIPRO_vs_Mock')
  #
  # diffabund <-
  #   lfcShrink(A_3_DE, coef =2) %>%
  #   as(Class = 'data.frame') %>%
  #   rownames_to_column(var='ASV') %>%
  #   filter(padj < 0.05) %>%
  #   filter(abs(log2FoldChange) > .5) %>%
  #   left_join(tax) %>%
  #   mutate({{variable}}:=ifelse(log2FoldChange > 0, other_level, ref_level))


  # p3 <- diffabund %>%
  #   # filter(ASV %in% in_both) %>%
  #   ggplot(aes(x=log2FoldChange, y=genus, fill=!!sym(variable))) +
  #   geom_point(size=2.5, shape=21) +
  #   geom_vline(xintercept = 0)
  # end function
  return(list(p1, p2, res_tibble))
}


#########

ps <- read_rds('processed_data/phyloseq_decontam.rds')

tax <-
  ps@tax_table %>%
  as(., Class = 'matrix') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV')


# set reference levels
sample_data(ps)$Vaccine <- factor(sample_data(ps)$Vaccine, levels = c('Mock', 'BBS866', 'AVIPRO'))
sample_data(ps)$Challenge <- factor(sample_data(ps)$Challenge, levels = c('Mock', 'Infantis', 'Reading'))




# taking stock of the available samples
ps@sam_data %>%
  group_by(day_post_challenge, Challenge, Vaccine) %>%
  tally() %>%
  pivot_wider(names_from = Challenge, values_from = n)



# A: vaccination effects on Turkey microbiome
# BBS866 vaccine
A_comp1 <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge == -1 &
                 ps@sam_data$Vaccine %in% c('Mock', 'BBS866')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)



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

### ANCOMBC explore ##

ANCOM_RES <- ancombc(A_comp3, formula = "Vaccine", group = 'Vaccine', struc_zero = TRUE, global = TRUE)

ANCOM_RES$res$beta$VaccineBBS866 %>% hist()
ANCOM_RES$res$beta$VaccineAVIPRO %>% hist()

ANCOM_RES$res$lfc


##### THIS ONE

tst <- A_comp3 %>% NMDS_from_phyloseq('Vaccine')

tst[[1]] %>%
  ggplot(aes(x=MDS1, y=MDS2, fill=Vaccine)) +
  geom_point(shape=21, color='white', size=4) +
  geom_segment(aes(xend=centroidX, yend=centroidY, color=Vaccine), alpha=.5) +
  theme_bw()


A_mat <- A_comp3@otu_table %>% as(Class='matrix')
A_dat <- A_comp3@sam_data %>% as(Class='data.frame')

look <- adonis2(A_mat ~ Vaccine, data = A_dat)

pairwise.adonis(A_mat, factors = A_dat[['Vaccine']]) %>%
  ggplot(aes(y=pairs, x=F.Model)) + geom_col() +
  theme_half_open()


A_3_DE <- phyloseq_to_deseq2(A_comp3, design = ~ Vaccine)

A_3_DE <- A_3_DE %>% DESeq2::DESeq()


num_res <- length(resultsNames(A_3_DE))

res_tibble <- tibble(comparison=resultsNames(A_3_DE)[-1])

res_tibble %>%
  mutate(res=map(.x = comparison,
                 .f=~DESeq_single_contrast(DESeq_obj = A_3_DE, COEF = .x, tax_tab = tax)))


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

in_both <- intersect(AVIPRO_diffabund$ASV, BBS866_diffabund$ASV)


# difabund in either vaccine vs mock
both_diff %>%
  # filter(ASV %in% in_both) %>%
  ggplot(aes(x=log2FoldChange, y=genus, fill=vaccine)) +
  geom_point(size=2.5, shape=21) +
  geom_vline(xintercept = 0)

# diffabund in both vaccines vs mock
both_diff %>%
  filter(ASV %in% in_both) %>%
  ggplot(aes(x=log2FoldChange, y=genus, fill=vaccine)) +
  geom_point(size=2.5, shape=21) +
  geom_vline(xintercept = 0)



# AVIPRO_diffabund %>% ggplot(aes(x=log2FoldChange, y=genus)) + geom_point()

# BBS866_diffabund %>% ggplot(aes(x=log2FoldChange, y=genus)) + geom_point()
#


###########

# B. Infantis Challenge effect on turkey microbiome
# really only 1 timepoint here even though its 10 and 14 in the metadata

B <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge %in% c(10, 14) &
                   ps@sam_data$Vaccine %in% c('Mock') &
                   ps@sam_data$Challenge %in% c('Infantis', 'Mock')) %>%
  prune_taxa(taxa = taxa_sums(.) > 10)


B_res <- explore(ps = B, variable = 'Challenge')

B_res[[1]]
B_res[[2]]
B_res[[3]]$difabund_plot[[1]]
B_res[[3]]$res

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

C_res <- explore(C, 'Vaccine')
###
#

#
# D. AviPro vaccination effect on Infantis-challenged turkey microbiome
# Samples: Mock-vaccinated/Infantis-challenged and AviPro-vaccinated/Infantis-challenged
# Question: Did AviPro vaccination modify the Infantis-induced alterations of the microbial composition of the turkey cecum (cecal contents)
#
# 13 mock-vaccinated/Infantis-challenged turkeys at 10dpi:    Minfantis_D10 #1-13
# VERSUS
# 14 -AviPro-vaccinated/Infantis-challenged turkeys at 10dpi:    Avipro_OO_D10 #1-14
#

# This one is a little tricky because we really need to know what the 'infantis effect'
# is
# If we had mock infected turkeys to compare all infected groups to, this would
#  be more effective.



D <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge %in% c(10) &
                   ps@sam_data$Vaccine %in% c('Mock', 'AVIPRO') &
                   ps@sam_data$Challenge %in% c('Infantis')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)


D_res <- explore(D, 'Vaccine')




# 10 days into an infantis infection, how different do the groups look?

CD_comb <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge %in% c(10) &
                   ps@sam_data$Vaccine %in% c('Mock', 'AVIPRO', 'BBS866') &
                   ps@sam_data$Challenge %in% c('Infantis')) %>%
  prune_taxa(taxa = taxa_sums(.) > 10)


CD_res <- explore(CD_comb, 'Vaccine')

CD_res[[1]]
CD_res[[2]]
CD_res[[3]]$difabund_plot[[1]]
CD_res[[3]]$difabund_plot[[2]]

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

