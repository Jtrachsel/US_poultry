library(usethis)


use_directory('raw_data')
use_directory('processed_data')
use_directory('reference_data')
use_directory('scripts')
use_readme_rmd()

options(timeout=6000)
getOption('timeout')
download.file(url = 'http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData',
              destfile = './reference_data/SILVA_SSU_r138_2019.RData' )

