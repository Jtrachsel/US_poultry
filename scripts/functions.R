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
