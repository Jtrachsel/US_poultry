library(tidyverse)
library(phyloseq)
library(cowplot)
library(DESeq2)
library(ANCOMBC)
library(funfuns)


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
      # filter(padj < 0.05) %>%
      # filter(abs(log2FoldChange) > .5) %>%
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


#function to :
# 1) subset phyloseq object
# 2) make ordination (NMDS)
# 3) PERMANOVAs
# 4) calc diffabund OTUs
#  a) DESeq2
#  b) ANCOM-BC
# 5) calc diffs in alpha diversity



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


  PWADON_res <- pairwise.adonis(A_mat, factors = A_dat[[variable]])
  p2 <-  PWADON_res %>%
    ggplot(aes(y=pairs, x=F.Model)) +
    geom_col() +
    geom_text(aes(label=p.adjusted))+
    theme_half_open()



  ### THIS SHOULD BE A FUNCTION
  # Alpha Diversity


  ALPHA <-
    ps %>%
    rarefy_even_depth() %>%
    estimate_richness() %>%
    rownames_to_column(var = 'sample_ID') %>%
    left_join(ps@sam_data %>% as(Class='matrix') %>% as.data.frame())

  # ALPHA %>%
  #   ggplot(aes(x=!!sym(variable)), y=Shannon, fill=!!sym(variable)) +
  #   geom_boxplot()

  # lm(data = ALPHA, formula = Shannon~!!sym(variable))) %>% broom::tidy()


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
  return(list(NMDS=p1, PWADON=p2, DIFFABUND=res_tibble, ALPHA_dat=ALPHA))
}


#########


ps <- read_rds('processed_data/phyloseq_final.rds')

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


ps@sam_data %>%
  group_by(day_post_vaccination,day_post_challenge, Challenge, Vaccine) %>%
  tally() %>%
  pivot_wider(names_from = Challenge, values_from = n)



# A: vaccination effects on Turkey microbiome
# BBS866 vaccine
A_comp1 <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge == -1 &
                 ps@sam_data$Vaccine %in% c('Mock', 'BBS866')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0)







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
#
# ANCOM_RES <- ancombc(A_comp3, formula = "Vaccine", group = 'Vaccine', struc_zero = TRUE, global = TRUE)
#
# ANCOM_RES$res$beta$VaccineBBS866 %>% hist()
# ANCOM_RES$res$beta$VaccineAVIPRO %>% hist()
#
# ANCOM_RES$res$lfc


##### THIS ONE

# alpha Diversity

A_COMP_3_RES <- explore(A_comp3, variable = 'Vaccine')
A_COMP_3_RES$NMDS
A_COMP_3_RES$PWADON
A_COMP_3_RES$DIFFABUND

A_COMP_3_RES$DIFFABUND %>% select(-difabund_plot) %>%
  mutate(filename=paste0('output/', comparison,'.tsv'),
         WRITE=map2(.x=filename, .y=res, .f=~write_tsv(x = .y, file = .x)))


A_COMP_3_RES$ALPHA_dat <-
  A_COMP_3_RES$ALPHA_dat %>%
  mutate(Vaccine=factor(Vaccine, levels=c('Mock', 'AVIPRO', 'BBS866')))

lm(data = A_COMP_3_RES$ALPHA_dat, formula = Shannon~Vaccine) %>% broom::tidy()

A_COMP_3_RES$ALPHA_dat %>%
  ggplot(aes(x=Vaccine, y=Shannon, fill=Vaccine)) + geom_boxplot() +
  theme_half_open() +
  ggtitle('', 'AVIPRO vs Mock P=0.0002, BBS866 vs Mock P=0.20')



ggsave('figures/05_vaccine_effect_alpha.jpeg', bg='white', height=4, width=7, units = 'in')

A_COMP_3_RES$NMDS

# tst <- A_comp3 %>% NMDS_from_phyloseq('Vaccine')
#
# tst[[1]] %>%
#   ggplot(aes(x=MDS1, y=MDS2, fill=Vaccine)) +
#   geom_point(shape=21, color='white', size=4) +
#   geom_segment(aes(xend=centroidX, yend=centroidY, color=Vaccine), alpha=.5) +
#   theme_bw()
ggsave('figures/06_vaccine_effect_NMDS.jpeg', bg='white', height=5, width=7, units = 'in')


A_COMP_3_RES$PWADON
ggsave('figures/07_vaccine_effect_PWadon.jpeg', bg='white', height=3, width = 7, units = 'in')


#
# A_3_DE <- phyloseq_to_deseq2(A_comp3, design = ~ Vaccine)
#
# A_3_DE <- A_3_DE %>% DESeq2::DESeq()
#
#
# num_res <- length(resultsNames(A_3_DE))
#
# res_tibble <- tibble(comparison=resultsNames(A_3_DE)[-1])
#
# res_tibble %>%
#   mutate(res=map(.x = comparison,
#                  .f=~DESeq_single_contrast(DESeq_obj = A_3_DE, COEF = .x, tax_tab = tax)))
#

# BBS866_diffabund <- results(A_3_DE, name ='Vaccine_BBS866_vs_Mock')
# AVIPRO_diffabund <- results(A_3_DE, name ='Vaccine_AVIPRO_vs_Mock')

BBS866_diffabund <- A_COMP_3_RES$DIFFABUND$res[[1]]%>% mutate(vaccine='BBS866')

# BBS866_diffabund <-
#   lfcShrink(A_3_DE, coef ='Vaccine_BBS866_vs_Mock') %>%
#   as(Class = 'data.frame') %>%
#   rownames_to_column(var='ASV') %>%
#   rownames_to_column() %>%
#   filter(padj < 0.05) %>%
#   filter(abs(log2FoldChange) > .5) %>%
#   left_join(tax) %>%
#   mutate(vaccine='BBS866')
AVIPRO_diffabund <- A_COMP_3_RES$DIFFABUND$res[[2]] %>% mutate(vaccine='AVIPRO')

#
# AVIPRO_diffabund <-
#   lfcShrink(A_3_DE, coef ='Vaccine_AVIPRO_vs_Mock') %>%
#   as(Class = 'data.frame') %>%
#   rownames_to_column(var='ASV') %>%
#   rownames_to_column() %>%
#   filter(padj < 0.05) %>%
#   filter(abs(log2FoldChange) > .5) %>%
#   left_join(tax) %>%
#   mutate(vaccine='AVIPRO')


both_diff <- bind_rows(AVIPRO_diffabund, BBS866_diffabund)

in_both <- intersect(AVIPRO_diffabund$ASV, BBS866_diffabund$ASV)

library(cowplot)
# difabund in either vaccine vs mock
both_diff %>%
  # filter(ASV %in% in_both) %>%
  ggplot(aes(x=log2FoldChange, y=genus, fill=vaccine)) +
  geom_point(size=2.5, shape=21) +
  geom_vline(xintercept = 0) +
  theme_half_open() +
  annotate(geom='label', x=-5, y=-1, label='Mock')+
  annotate(geom='text', x=-5, y=-3, label='')+
  annotate(geom='label', x=5, y=-1, label='Vaccine')+
  theme(panel.grid.major.y = element_line(color='grey'))

ggsave('./figures/08_vaccine_diffabund_all.jpeg', bg='white', width=9, height=9, units = 'in')

# diffabund in both vaccines vs mock
both_diff %>%
  filter(ASV %in% in_both) %>%
  ggplot(aes(x=log2FoldChange, y=genus, fill=vaccine)) +
  geom_jitter(size=2.5, shape=21, height = .1) +
  geom_vline(xintercept = 0)+
  theme_half_open() +
  annotate(geom='label', x=-5, y=-.5, label='Mock')+
  annotate(geom='text', x=-5, y=-1, label='')+
  annotate(geom='label', x=5, y=-.5, label='Vaccine')+
  theme(panel.grid.major.y = element_line(color='grey'))

ggsave('./figures/09_vaccine_diffabund_BOTH.jpeg', bg='white', width=9, height=5, units = 'in')



#### TO DO!!! FROM THE OTUS THAT WENT UP IN INFANTIS INFECTED GROUPS LABEL THE
# ONES THAT SHOW UP IN THE VACCINE ASSOCIATED GROUPS


# AVIPRO_diffabund %>% ggplot(aes(x=log2FoldChange, y=genus)) + geom_point()

# BBS866_diffabund %>% ggplot(aes(x=log2FoldChange, y=genus)) + geom_point()
#


###########

# B. Infantis Challenge effect on turkey microbiome
# really only 1 timepoint here even though its 10 and 14 in the metadata

# ONLY MOCK VACCINATED
# MOCK AND INFANTIS INFECTED

B <-
  ps %>%
  subset_samples(ps@sam_data$day_post_challenge %in% c(10, 14) &
                   ps@sam_data$Vaccine %in% c('Mock') &
                   ps@sam_data$Challenge %in% c('Infantis', 'Mock')) %>%
  prune_taxa(taxa = taxa_sums(.) > 10)


B_res <- explore(ps = B, variable = 'Challenge')

B_res$NMDS + ggtitle('Mock Vaccinated Only')
ggsave('./figures/10_Mock_vaccine_Infantis_challenge_NMDS.jpeg', bg='white', height=5, width=7, units = 'in')
B_res$PWADON

B_res$DIFFABUND$difabund_plot[[1]]
ggsave('./figures/11_Mock_vaccine_Infantis_challenge_difabund.jpeg', bg='white', width=9, height=8, units = 'in')

B_res$ALPHA_dat %>%
  ggplot(aes(x=Challenge, y=Shannon)) +
  geom_boxplot(aes(fill=Challenge)) +
  ggtitle('','Infantis vs Mock P=0.00008')
ggsave('./figures/12_Mock_vaccine_Infantis_challenge_alpha.jpeg', bg='white', width=7, height=4)

lm(data=B_res$ALPHA_dat, formula = Shannon~Challenge) %>% broom::tidy()

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


C_res$NMDS

ggsave('figures/13_BBSvaccine_effect_on_infection_NMDS.jpeg', bg='white', height=5, width=7, units = 'in')
#
C_res$PWADON
#
C_res$DIFFABUND$difabund_plot
ggsave('figures/14_BBSvaccine_effect_on_infection_difabund.jpeg', bg='white', height=8, width=7, units = 'in')
#


lm(data=C_res$ALPHA_dat, formula= Shannon~Vaccine) %>% broom::tidy()
C_res$ALPHA_dat %>%
  ggplot(aes(x=Vaccine, y=Shannon, fill=Vaccine)) +
  geom_boxplot()+
  ggtitle('','P=0.002')
#
ggsave('figures/15_BBSvaccine_effect_on_infection_alpha.jpeg', bg='white', height=4, width=7, units = 'in')



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


D_res$NMDS

ggsave('figures/16_AVIPROvaccine_effect_on_infection_NMDS.jpeg', bg='white', height=5, width=7, units = 'in')
#
D_res$PWADON
#
D_res$DIFFABUND$difabund_plot
ggsave('figures/17_AVIPROvaccine_effect_on_infection_difabund.jpeg', bg='white', height=9, width=7, units = 'in')
#

lm(data=D_res$ALPHA_dat, formula = Shannon~Vaccine) %>% broom::tidy()
D_res$ALPHA_dat %>%
  ggplot(aes(x=Vaccine, y=Shannon, fill=Vaccine)) +
  geom_boxplot() +
  ggtitle('','P=0.069')
#
ggsave('figures/18_AVIPROvaccine_effect_on_infection_alpha.jpeg', bg='white', height=4, width=7, units = 'in')

# # 10 days into an infantis infection, how different do the groups look?
#
# CD_comb <-
#   ps %>%
#   subset_samples(ps@sam_data$day_post_challenge %in% c(10) &
#                    ps@sam_data$Vaccine %in% c('Mock', 'AVIPRO', 'BBS866') &
#                    ps@sam_data$Challenge %in% c('Infantis')) %>%
#   prune_taxa(taxa = taxa_sums(.) > 10)
#
#
# CD_res <- explore(CD_comb, 'Vaccine')
#
# CD_res$NMDS
# CD_res$PWADON
# CD_res$DIFFABUND$difabund_plot[[1]]
# CD_res$DIFFABUND$difabund_plot[[2]]

MOCK_MOCKD7D14 <-
  prune_samples(
    samples = ps@sam_data$Vaccine == 'Mock' &
      ps@sam_data$Challenge == 'Mock' &
      ps@sam_data$day_post_challenge %in% c(14), x = ps)
#ps@sam_data$day_post_challenge %in% c(7), x = ps)
#ps@sam_data$day_post_challenge %in% c(14), x = ps)

INFANTIS_INFECT <-
  prune_samples(
    samples = ps@sam_data$Challenge == 'Infantis',
    x=ps)


INFANTIS <- merge_phyloseq(MOCK_MOCKD7D14, INFANTIS_INFECT)

INFANTIS@sam_data$Vac.Inf <- paste0(INFANTIS@sam_data$Vaccine,'.', INFANTIS@sam_data$Challenge)

INFANTIS@sam_data$Vac.Inf <-
  factor(INFANTIS@sam_data$Vac.Inf,
         levels = c('Mock.Mock', 'Mock.Infantis','AVIPRO.Infantis', 'BBS866.Infantis'))

INFANTIS2 <- INFANTIS

INFANTIS2@sam_data$Vac.Inf <-
  factor(INFANTIS@sam_data$Vac.Inf,
         levels = c('Mock.Infantis','Mock.Mock', 'AVIPRO.Infantis', 'BBS866.Infantis'))




library(cowplot)

library(DESeq2)

tst <- explore(ps = INFANTIS, variable = 'Vac.Inf')

tst2 <- explore(ps = INFANTIS2, variable = 'Vac.Inf')

tst$DIFFABUND
tst2$DIFFABUND %>%
  mutate(filename=paste0('output/', comparison,'.tsv'),
         WRITE=map2(.x=filename, .y=res, .f=~write_tsv(x = .y, file = .x)))


bind_rows(tst$DIFFABUND, tst2$DIFFABUND) %>%
  filter(comparison != 'Vac.Inf_Mock.Mock_vs_Mock.Infantis') %>%
  select(-difabund_plot) %>% unnest(res) %>%
  mutate(padj2=p.adjust(pvalue, method = 'fdr')) %>%
  filter(padj2 < 0.05) %>%
  group_by(comparison) %>%
  nest() %>%
  mutate(filename=paste0('output/', comparison,'.tsv'),
         WRITE=map2(.x=filename, .y=data, .f=~write_tsv(x = .y, file = .x)))







tst$NMDS
ggsave('figures/19_Vaccine_infection_NMDS.jpeg', bg='white', height=5, width=7, units = 'in')



tst$PWADON
ggsave('figures/20_Vaccine_infection_PERMANOVA.jpeg', bg='white', height=2.5, width=6, units = 'in')

p <- tst$ALPHA_dat %>% ggplot(aes(x=Vac.Inf, y=Shannon, fill=Vac.Inf)) +
  geom_boxplot()
p
ggsave('figures/21_Vaccine_infection_ALPHA.jpeg', bg='white', height=4, width=7, units = 'in')



MOD <- lm(data=tst$ALPHA_dat, formula= Shannon ~ Vac.Inf)
EMM <- emmeans::emmeans(MOD, specs = c('Vac.Inf'))

CONTRASTS <-
  contrast(EMM, method = 'pairwise', adjust = 'BH') %>%
  broom::tidy() #%>%
  # arrange(adj.p.value)


CONTRASTS %>% select(contrast, estimate, std.error, adj.p.value) %>%
  write_tsv('output/All_infantis_alpha_tests.tsv')
# VECTOR <- CONTRASTS$adj.p.value
# names(VECTOR) <- CONTRASTS$contrast
#
# LETS <- multcompLetters(VECTOR)
# LETS <- tibble(NAMES=names(LETS$Letters),LETS=LETS$monospacedLetters) %>%
  # mutate(NAMES=gsub(' ','',NAMES))



# ANNO_DAT <-
  # tst$ALPHA_dat %>%
  # group_by(Vac.Inf) %>%
  # summarise(NAMES=unique(Vac.Inf), Shannon=max(Shannon)) %>%
  # ungroup() %>%
  # left_join(LETS)

# p+geom_text(data=ANNO_DAT, aes(label=LETS), nudge_y = .2)



# tst$DIFFABUND$res[[1]]
# tst$DIFFABUND %>% select(-difabund_plot) %>% unnest(res)
# tst$DIFFABUND$res[[1]] %>% as_tibble() %>% arrange(desc(padj))



# VAC_CHAL_DESEQ <- INFANTIS %>% phyloseq_to_deseq2(design = ~ Vaccine + Challenge) %>%
#   DESeq()
#
#
# resultsNames(VAC_CHAL_DESEQ)
# infantis_effect_diffabund <-
#   DESeq2::lfcShrink(VAC_CHAL_DESEQ, coef ='Challenge_Infantis_vs_Mock' ) %>%
#   as.data.frame() %>%
#   rownames_to_column(var='ASV') %>%
#   left_join(tax) %>%
#   filter(padj <0.05) %>%
#   filter(abs(log2FoldChange) > .5) %>%
#   mutate(condition=ifelse(log2FoldChange > 0, 'Infantis', 'Mock-challenged'))
#
#
# infantis_effect_diffabund %>% ggplot(aes(x=log2FoldChange, y=genus, fill=condition)) + geom_point(shape=21)
#
#
# UP_IN_SALMONELLA <- infantis_effect_diffabund %>% filter(log2FoldChange > 0) %>% pull(ASV)
# UP_IN_MOCKS <- infantis_effect_diffabund %>% filter(log2FoldChange < 0) %>% pull(ASV)
#

#Very interesting... Vaccines seem to have made Salmonella colonization
# worse.
# Significantly higher Sal Shedding in both AviPro and BBS866 vaccinated groups
# compared to Mock Vaccinated Infantis challenged group.

# Large differences in bacterial communities from cecal contents
# when comparing infantis infected groups back to mock-infected controls, the
# mock vaccinated are the most similar to the mock-challenged. and  both
# vaccine groups are quite different from the mock vaccine and each other.

# seems like each vaccine enhanced the 'salmonella effect' on microbial communities
#

#in summary, not only is vaccination associated with worse colonization, it
# is also associated with greater difference in microbial communities after
# salmonella Infantis infection.


# SPECULATION ZONE:
# could it be that being vaccinated enhanced the inflammatory immune response
# so when salmonella was encountered again, increased gut inflammation allowed more
# extensive colonization?

# Alternatively, could vaccination have effected microbial communities? then
# these changed microbial communities were less effective at excluding Salmonella
# upon infection?



### COULD THROW THESE IN ###


PRE_INFECTION <- prune_samples(ps, samples=ps@sam_data$day_post_challenge == -1)



INFANTIS <- merge_phyloseq(MOCK_MOCKD7D14, INFANTIS_INFECT, PRE_INFECTION)
INFANTIS@sam_data$TIME <- ifelse(INFANTIS@sam_data$day_post_challenge > 0, 'PREINFECTION', 'POSTINFECTION')
INFANTIS@sam_data$Vac.Inf.time <- paste0(INFANTIS@sam_data$Vaccine,'.', INFANTIS@sam_data$Challenge, '.',INFANTIS@sam_data$TIME)
#
# INFANTIS@sam_data$Vac.Inf <-
#   factor(INFANTIS@sam_data$Vac.Inf,
#          levels = c('Mock.Mock', 'Mock.Infantis','AVIPRO.Infantis', 'BBS866.Infantis'))
#
# INFANTIS@sam_data$Challenge
library(cowplot)

library(DESeq2)

tst <- explore(ps = INFANTIS, variable = 'Vac.Inf.time')

tst$NMDS



ps@sam_data %>%
  as.data.frame() %>%
  filter(day_post_challenge == -1) %>%
  group_by(Vaccine) %>%
  tally()



##### linear assoc
# when you account for vaccine/room only two ASVs are linearly associated with salmonella
# colonization in the cecum
# this maybe becasue communities are so different within each room,


AVIvac <- prune_samples(INFANTIS_INFECT, samples = INFANTIS_INFECT@sam_data$Vaccine == 'AVIPRO')
BBSvac <- prune_samples(INFANTIS_INFECT, samples = INFANTIS_INFECT@sam_data$Vaccine == 'BBS866')
MOCKvac <- prune_samples(INFANTIS_INFECT, samples = INFANTIS_INFECT@sam_data$Vaccine == 'Mock')
#


# INFANTIS_glom <- tax_glom(INFANTIS_INFECT, taxrank = 'genus')
INFANTIS_DESEQ <- phyloseq_to_deseq2(MOCKvac, ~ log10_CFU) %>% DESeq()
resultsNames(INFANTIS_DESEQ)

lfcShrink(INFANTIS_DESEQ, coef = 'log10_CFU') %>%
  as.data.frame() %>%
  rownames_to_column(var='ASV') %>%
  left_join(tax) %>%
  filter(padj < 0.05)


scale(INFANTIS_INFECT@sam_data$log10_CFU)



