## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = TRUE,
  fig.width = 11,
  tidy.opts = list(width.cutoff = 95),
  message = FALSE,
  warning = FALSE,
  time_it = TRUE
)
all_times <- list()  # store the time for each chunk
knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res <- difftime(Sys.time(), now, units = "secs")
      all_times[[options$label]] <<- res
    }
  }
}))

## ----eval=FALSE---------------------------------------------------------------
#  githubURL <- "https://github.com/feiyoung/PRECAST/blob/main/vignettes_data/bc2.rda?raw=true"
#  download.file(githubURL,"bc2.rda",mode='wb')

## ----eval = FALSE-------------------------------------------------------------
#  load("bc2.rda")

## ----eval=FALSE---------------------------------------------------------------
#  # library(DR.SC)
#  # dir.file <- "Section" ## the folders Section1 and Section2, and each includes two folders spatial and filtered_feature_bc_matrix
#  # seuList <- list()
#  # for (r in 1:2) {
#  #   message("r = ", r)
#  #   seuList[[r]] <- read10XVisium(paste0(dir.file, r))
#  # }
#  # bc2 <- seuList

## ----eval = FALSE-------------------------------------------------------------
#  library(PRECAST)
#  library(Seurat)

## ----eval = FALSE-------------------------------------------------------------
#  bc2 ## a list including two Seurat object

## ----eval = FALSE-------------------------------------------------------------
#  head(bc2[[1]])

## ----eval = FALSE-------------------------------------------------------------
#  ## Create PRECASTObject.
#  set.seed(2022)
#  PRECASTObj <- CreatePRECASTObject(bc2)

## ----eval = FALSE-------------------------------------------------------------
#  ## check the number of genes/features after filtering step
#  PRECASTObj@seulist
#  
#  ## Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
#  PRECASTObj <-  AddAdjList(PRECASTObj, platform = "Visium")
#  
#  ## Add a model setting in advance for a PRECASTObj object. verbose =TRUE helps outputing the information in the algorithm.
#  PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal=FALSE, verbose=TRUE, int.model=NULL)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ### Given K
#  PRECASTObj <- PRECAST(PRECASTObj, K=14)
#  

## ----eval = FALSE-------------------------------------------------------------
#  ## backup the fitting results in resList
#  resList <- PRECASTObj@resList
#  PRECASTObj <- selectModel(PRECASTObj)
#  

## ----eval = FALSE-------------------------------------------------------------
#  seuInt <- IntegrateSpaData(PRECASTObj, species='Human')
#  seuInt
#  ## The low-dimensional embeddings obtained by PRECAST are saved in PRECAST reduction slot.

## ----eval = FALSE, fig.height = 4, fig.width=7.5------------------------------
#  cols_cluster <-  c( "#FD7446" ,"#709AE1", "#31A354","#9EDAE5",
#   "#DE9ED6" ,"#BCBD22", "#CE6DBD" ,"#DADAEB" ,
#   "yellow", "#FF9896","#91D1C2", "#C7E9C0" ,"#6B6ECF", "#7B4173" )
#  p12 <- SpaPlot(seuInt, batch=NULL,point_size=1, cols=cols_cluster, combine=TRUE)
#  p12
#  # users can plot each sample by setting combine=FALSE

## ----eval = FALSE, fig.height = 3, fig.width=5--------------------------------
#  seuInt <- AddUMAP(seuInt)
#  p13 <- SpaPlot(seuInt, batch=NULL,item='RGB_UMAP',point_size=2, combine=TRUE, text_size=15)
#  p13
#  #seuInt <- AddTSNE(seuInt)
#  #SpaPlot(seuInt, batch=NULL,item='RGB_TSNE',point_size=2, combine=T, text_size=15)

## ----eval = FALSE-------------------------------------------------------------
#  seuInt <- AddTSNE(seuInt, n_comp = 2)
#  library(patchwork)
#  p1 <- dimPlot(seuInt, item='cluster', point_size = 0.5, font_family='serif', cols=cols_cluster) # Times New Roman
#  p2 <- dimPlot(seuInt, item='batch', point_size = 0.5,  font_family='serif')
#  p1 + p2

## ----eval = FALSE-------------------------------------------------------------
#  library(Seurat)
#  dat_deg <- FindAllMarkers(seuInt)
#  library(dplyr)
#  n <- 10
#  dat_deg %>%
#    group_by(cluster) %>%
#    top_n(n = n, wt = avg_log2FC) -> top10
#  
#  seuInt <- ScaleData(seuInt)
#  seus <- subset(seuInt, downsample = 400)
#  
#  

## ----eval = FALSE, fig.height = 6, fig.width=9--------------------------------
#  color_id <- as.numeric(levels(Idents(seus)))
#  library(ggplot2)
#  ## HeatMap
#  p1 <- doHeatmap(seus, features = top10$gene, cell_label= "Domain",
#                  grp_label = F, grp_color = cols_cluster[color_id],
#                  pt_size=6,slot = 'scale.data') +
#    theme(legend.text = element_text(size=10),
#          legend.title = element_text(size=13, face='bold'),
#          axis.text.y = element_text(size=5, face= "italic", family='serif'))
#  p1

## ----eval = FALSE-------------------------------------------------------------
#  sessionInfo()
