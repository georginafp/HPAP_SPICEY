##--------------------------------------------------------##
## LIFTOVER SUPERENHANCERS, ENHANCERS-HUBS, ENHANCER-STRETCH
##--------------------------------------------------------##

library(magrittr)
library(GenomicRanges)
source("R/utils.R")


# load enhancers datasets (hg19)
files <- list.files("/homes/users/gfuentes/scratch/projects/spicey/test/data/TEST_ENHANCERS",
                  full.names = T,
                  pattern = "hg19")

# read and import them as GRanges
enh_datasets <- GRangesList(
  setNames(lapply(files, function(file) {
    load(file)
    get(ls()[1])
  }),
  tools::file_path_sans_ext(basename(files)))
)

# Liftover from hg19 to hg38
enh_datasets_hg38 <- enh_datasets |>
  lapply(liftOver_hg19_hg38) |>
  setNames(gsub("hg19", "hg38", names(enh_datasets))) |>
  GRangesList()
saveRDS(enh_datasets_hg38, "/homes/users/gfuentes/scratch/projects/spicey/test/data/TEST_ENHANCERS/islet_enhancers_datasets_hg38.rds")

