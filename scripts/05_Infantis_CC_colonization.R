library(phyloseq)
library(tidyverse)

set.seed(10)

# read in shedding data

colonization_dat <-
  read_tsv('processed_data/US_Poultry_Infantis_CC_colonization.txt')

col_dat <-
  colonization_dat %>%
  transmute(turkey_ID=`Turkey ID#`,
            Room,
            Vaccine=sub(' \\(SB 358\\)','',Vaccine),
            Day=`Day (T= )`,
            dillution=`Dilution (x-fold)`, CFU=`CFU/g`) %>%
  mutate(CFU=map_dbl(.x=CFU, .f=~ifelse(is.na(.x),sample(1:20, 1), .x )),
         log10_CFU=log10(CFU),
         Vaccine=factor(Vaccine, levels = c('mock', 'AviPro', '866')))



bootstrap_colonization_analysis <- function(seed, colonization_dat){
  set.seed(seed)
  col_dat <-
    colonization_dat %>%
    transmute(turkey_ID=`Turkey ID#`,
              Room,
              Vaccine=sub(' \\(SB 358\\)','',Vaccine),
              Day=`Day (T= )`,
              dillution=`Dilution (x-fold)`, CFU=`CFU/g`) %>%
    mutate(CFU=map_dbl(.x=CFU, .f=~ifelse(is.na(.x),sample(1:20, 1), .x )),
           log10_CFU=log10(CFU),
           Vaccine=factor(Vaccine, levels = c('mock', 'AviPro', '866')))

  # col_dat %>% ggplot(aes(x=log10_CFU, fill=Vaccine)) + geom_histogram()

  p <- col_dat %>% ggplot(aes(x=Vaccine, y=log10_CFU, fill=Vaccine)) + geom_boxplot()
  lm_summary <-
    broom::tidy(lm(data = col_dat, formula = log10_CFU ~ Vaccine), conf.int=TRUE) %>%
    mutate(SEED=seed)

  results <- list(PLOT=list(p),LM=lm_summary)

  return(results)
}


BS_res <-
  tibble(seed=1:100,
         res=map(.x=seed, .f=~bootstrap_colonization_analysis(.x, colonization_dat))) %>%
  unnest_wider(res) %>%
  unnest(LM)



BS_res %>%
  filter(term != '(Intercept)') %>%
  ggplot(aes(x=p.value, fill=term)) + geom_histogram()

ggsave('figures/SHED_BS_pval_hist.jpeg', bg='white', height=2.5, width=6.5)

BS_res %>%
  filter(term != '(Intercept)') %>%
  ggplot(aes(x=estimate, y=SEED, xmin=conf.low, xmax=conf.high, color=term)) +
  geom_pointrange(position = position_dodge(width = .5)) +
  geom_point(aes(fill=term), shape=21, color='black', size=2)+
  geom_vline(xintercept = 0)

ggsave('figures/SHED_BS_conf_int.jpeg', bg='white', height=9.5, width=6.5)

BS_res %>% filter(term == 'Vaccine866')


col_dat %>% write_tsv('processed_data/Infantis_colonization_clean.tsv')


library(emmeans)
library(cowplot)

col_dat %>%
  ggplot(aes(x=Vaccine, y=log10_CFU, fill=Vaccine)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .15,size=2.5, shape=21, color='black',aes(fill=Vaccine))+
  ylab('log10 CFU/g Cecal Contents') +
  theme_half_open() +
  ylim(0,6)

ggsave('figures/04_Infantis_CC_Vaccine.jpeg', bg='white', width = 7, height=4)

LM <- lm(data=col_dat, formula=log10_CFU~Vaccine)
EMMEANS <- emmeans::emmeans(LM, specs = 'Vaccine')
CONTRAST <- contrast(EMMEANS, method = 'pairwise', adjust = 'BH')
broom::tidy(CONTRAST) %>% write_tsv('output/Infantis_CC_Vaccine_tests.tsv')
