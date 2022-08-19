library(phyloseq)
library(tidyverse)

ps <- read_rds('processed_data/phyloseq_final.rds')

MOCK_MOCKD7D14 <-
  prune_samples(
    samples = ps@sam_data$Vaccine == 'Mock' &
              ps@sam_data$Challenge == 'Mock' &
              ps@sam_data$day_post_challenge %in% c(7,14), x = ps)
             #ps@sam_data$day_post_challenge %in% c(7), x = ps)
             #ps@sam_data$day_post_challenge %in% c(14), x = ps)

INFANTIS_INFECT <-
  prune_samples(
    samples = ps@sam_data$Challenge == 'Infantis',
    x=ps)


INFANTIS <- merge_phyloseq(MOCK_MOCKD7D14, INFANTIS_INFECT)

INFANTIS@sam_data$Vac_Inf <- paste0(INFANTIS@sam_data$Vaccine,'_', INFANTIS@sam_data$Challenge)
library(cowplot)

library(DESeq2)

tst <- explore(ps = INFANTIS, variable = 'Vac_Inf')

tst[[3]]$difabund_plot[[1]]



# Infantis vaccine trial

### PULL IN MOCK INFECTED MOCK VACCINATED FROM D7 and D14
# MOCK MOCK D7
# MOCK MOCK D14
# MOCK INF D10
# AVIPRO INF D10
# BBS866 INF D10

# first answer SHawns Qs as written


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










#### now try a Jules Remix ###
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
