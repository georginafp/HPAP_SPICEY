library(SPICEY)

# Co-accessibility mode
result_coacc <- run_spicey(
  atac_path = "data/FINAL_ATAC.rds",
  rna_path = "data/FINAL_RNA.rds",
  links_path = "data/COACC_LINKS.rds",
  linking_method = "coaccessibility"
)

# Nearest-gene mode
result_nearest <- run_spicey(
  atac_path = "data/FINAL_ATAC.rds",
  rna_path = "data/FINAL_RNA.rds",
  linking_method = "nearest"
)
