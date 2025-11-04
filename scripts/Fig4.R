

library(here)
library(ggtree)
library(ggpubr)
source(here("src/setup.R"))
source(here("src/landscape_plots.R"))
theme_set(theme_cowplot(font_size = 8, line_size = 0.25))
library(glue)
library(grid)
library(ComplexHeatmap)
rm(list = ls())

#setwd("/Users/yvesgreatti/github/brca-driver-estimation")

config <- read_yaml(here("config.yaml"))
cnarm <- fread(file.path(config$basedir,"zenododirectory/normal_breast/cn_arm.csv.gz"))
metrics <- fread(file.path(config$basedir, "zenododirectory/normal_breast/metrics.csv.gz"))
cndat_md <- distinct(metrics, cell_id,cell_type1, sample, genotype, organoid)
cn_arm <- left_join(cn_arm, cndat_md)

## Identify extreme aneuploid cells 

### Find cutoff
df_cum <- metrics %>% filter(keep_cell == TRUE)

total_cells <- nrow(df_cum %>% filter(n_aneuploid_arm > 0))

df_cum %>% 
  filter(n_aneuploid_arm > 0) %>% 
  arrange(n_aneuploid_arm) %>% 
  mutate(x = 1) %>% 
  mutate(x = cumsum(x) / total_cells) %>% 
  #filter(n_aneuploid_arm > 2) %>% 
  ggplot(aes(y = x, x = n_aneuploid_arm)) + 
  geom_line()

df_cum %>% 
  filter(n_aneuploid_arm > 0) %>% 
  arrange(n_aneuploid_arm) %>% 
  mutate(x = 1) %>% 
  mutate(x = cumsum(x) / total_cells) %>% 
  filter(x > 0.95)

metrics %>% 
  #filter(cohort1 != "Bulk") %>% 
  filter(keep_cell == TRUE) %>% 
  replace_na(list(n_aneuploid_arm =0)) %>% 
  filter(n_aneuploid_arm > 0) %>% 
  group_by(n_aneuploid_arm) %>% 
  summarise(n = n()) %>% 
  mutate(n = n / sum(n)) %>% 
  mutate(cumsum(n))

#set extreme cutoff as. > 95%
extreme_cutoff <- 6 

cancer <- fread(paste0(config$cancerdir, "nikzainal_bins.csv.gz")) %>% 
  mutate(sample = str_remove(sample, "a")) %>% 
  mutate(sample = str_remove(sample, "_2"))

cancer_metadata <- fread(paste0(config$cancerdir, "nikzainal_metadata.csv.gz"))

#only include samples with purity > 0.5
samples_to_keep <- cancer_metadata %>% 
  filter(Aberrant.cell.fraction > 0.5) %>% 
  pull(sample)
cancer <- filter(cancer, sample %in% samples_to_keep)
cancer <- left_join(cancer %>% select(-ploidy), cancer_metadata %>% select(sample, ploidy))
cancer$ploidy <- round(cancer$ploidy)

cancer_arm <- signals::per_chrarm_cn(cancer %>% rename(cell_id = sample) %>% mutate(copy = state))
cancer_arm <- cancer_arm %>% 
  left_join(distinct(cancer, sample, ploidy) %>% rename(cell_id = sample)) %>% 
  mutate(ploidy = round(ploidy))

cn_arm_format <- metrics %>% mutate(narm = n_aneuploid_arm)



mycells <- metrics %>% filter(keep_cell == TRUE) %>% 
  filter(n_aneuploid_arm > extreme_cutoff | 
           cell_id %in% c("SA1150-A95735B-R34-C35", 
                          "SA1150-A95735B-R44-C60",
                          "SA1150-A95735B-R42-C51",
                          "SA1150-A95735B-R43-C30"))
extreme_cells <- mycells

pl_cancer <- cancer %>%
  mutate(ploidy = round(ploidy)) %>% 
  as.data.table() %>%
  .[, list(Gain = sum(state > ploidy, na.rm = TRUE) / .N,
           Loss = sum(state < ploidy, na.rm = TRUE) / .N,
           n = .N), by = .(chr, start, end)] %>%
  plottinglist_()

cn <- fread(file.path(config$basedir, "zenododirectory/normal_breast/cn.csv.gz"))

md <- metrics[keep_cell == TRUE] %>% 
  add_count(sample) %>% 
  filter(n >= config$cutoff)
cell_list <- extreme_cells %>% 
  filter(cell_id %in% md$cell_id) %>% 
  add_count(sample)

chroms <- c(paste0(1:11), "13", "15", "17", "20", "X")
glist <- list()
for (mysample in unique(cell_list$sample)){
  print(mysample)
  cells_to_keep <- cell_list %>% filter(sample == mysample)
  print(nrow(cells_to_keep))
  cell_md <- df %>% filter(cell_id %in% cells_to_keep$cell_id)
  cn <- dat$cn %>% filter(cell_id %in% cells_to_keep$cell_id)
  # if (nrow(cells_to_keep) < 2){
  #   next
  # }
  if (nrow(cells_to_keep) < 4){
    cl <- list(tree = NULL)
    cl$clusters <- cells_to_keep %>% mutate(clone_id = "0")
  } else{
    cl <- umap_clustering(cn, minPts = 3,n_neighbors = 20, min_dist = 0, field = "state", umapmetric = "euclidean")
  }
  ptid <- strsplit(mysample, "-")[[1]][2]
  my_title = paste0("Patient ", ptid, ", ", cell_md$genotype[1], "+/-", ", ", nrow(cell_md), "/", nrow(md[sample == mysample]), " cells")
  pdf(here(paste0("Figures/Heatmaps/", mysample, "_heatmap_extremecells.pdf")), 
      w = 89 * 0.039,
      h = 2)
  p <- plotHeatmap(cn,
                   column_title = my_title,
                   clone_pal = extreme_cols,
                   column_title_gp = gpar(fontsize = 8),
                   linkheight = 2,
                   chrlabels = chroms,
                   show_heatmap_legend = F,
                   plotfrequency = F, 
                   frequency_height = 0.5,
                   anno_width = 0.1,
                   annofontsize = 7,
                   show_legend = F,
                   show_clone_text = F,
                   show_library_label = F,
                   #tree = cl$tree,
                   plottree = F,
                   reorderclusters = T,
                   clusters = cell_md)
  print(p)
  glist[[mysample]] <- grid.grabExpr(draw(p), width = 89 * 0.039, height = 2)
  dev.off()
}

