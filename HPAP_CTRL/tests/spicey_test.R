library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(Seurat)
library(tidyr)
library(regioneR)
library(purrr)
library(Signac)
library(purrr)
library(ggside)
library(viridis)
library(tibble)
library(ggsignif)
library(RColorBrewer)
library(ggplotify)
library(pheatmap)
library(patchwork)
library(scales)
library(data.table)
library(here)
library(InteractionSet)
library(stringr)
library(tools)


source("R/getsi.R")
source("R/retsi.R")
source("R/utils.R")
source("R/link_spicey.R")
source("R/run_spicey.R")




# Co-accessibility mode
retsi <- compute_retsi_only("data/ATAC.rds")
getsi <- compute_getsi_only("data/RNA.rds")

spicey_linked_coaccess <- join_spicey_coaccessibility(retsi, getsi, "data/LINKS.rds", )
saveRDS(result_coacc, "/homes/users/gfuentes/scratch/projects/spicey_old/HPAP_CTRL/data/HPAP_CTRL_SPICEY_COACC.rds")

spicey_linked_nearest <- join_spicey_nearest(retsi, getsi)
saveRDS(result_nearest, "/homes/users/gfuentes/scratch/projects/spicey_old/HPAP_CTRL/data/HPAP_CTRL_SPICEY_NEAREST.rds")


