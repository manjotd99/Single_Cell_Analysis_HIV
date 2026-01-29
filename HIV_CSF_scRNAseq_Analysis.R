library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

# merged Seurat object from GSE202410
# metadata includes disease_status and tissue
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

# filter low-quality cells
combined_filtrd <- subset(
  combined,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    nCount_RNA > 2000 &
    nCount_RNA < 15000 &
    percent.mt < 10

#normalize the dataset 
combined_norm <- NormalizeData(combined_filtrd)
combined_norm <- FindVariableFeatures(
  combined_norm,
  selection.method = "vst",
  nfeatures = 2000
)

# scale and PCA
variable.genes <- VariableFeatures(combined_norm)
combined_norm <- ScaleData(combined_norm, features = variable.genes)
combined_norm <- RunPCA(combined_norm, features = variable.genes)


# clustering then choose PCs based on elbow plot
combined_norm <- FindNeighbors(combined_norm, dims = 1:20)
combined_norm <- FindClusters(combined_norm, resolution = 0.5)
combined_norm <- RunUMAP(combined_norm, dims = 1:20)


#identifying cluster identities 
Idents(combined_norm) <- "seurat_clusters"
cd14_clusters <- c("6", "10", "14")
cd16_clusters <- c("11", "7")
mg_clusters   <- c("17")

# assign cell type labels to clusters based on marker expression  
combined_norm$celltype <- "Other"
combined_norm$celltype[combined_norm$seurat_clusters %in% cd14_clusters] <- "CD14"
combined_norm$celltype[combined_norm$seurat_clusters %in% cd16_clusters] <- "CD16"
combined_norm$celltype[combined_norm$seurat_clusters %in% mg_clusters]   <- "Microglia Cells"


#subset the CSF 
cd14_csf <- subset(combined_norm, subset = celltype == "CD14" & tissue == "CSF")
cd16_csf <- subset(combined_norm, subset = celltype == "CD16" & tissue == "CSF")
mg_csf   <- subset(combined_norm, subset = celltype == "Microglia Cells" & tissue == "CSF")


#proportion/tests
prop_test <- function(obj, gene = "GLS") {
  expr <- FetchData(obj, vars = gene)[,1]
  status <- obj$disease_status
  clec_pos <- expr > 0
  tbl <- table(status, clec_pos)
  print(tbl)
  
  if(any(tbl < 5)) {
    fisher.test(tbl)
  } else {
    chisq.test(tbl)
  }
}

prop_test(cd14_csf)
prop_test(cd16_csf)
prop_test(mg_csf)

ttest_gene <- function(obj, gene = "GLS") {
  expr <- FetchData(obj, vars = gene)
  expr$group <- obj$disease_status
  t.test(expr[, gene] ~ expr$group)
}

ttest_gene(cd14_csf, "GLS")
ttest_gene(cd16_csf, "GLS")
ttest_gene(mg_csf,   "GLS")


#function to create UMAP figures 
make_umap_panel <- function(obj, title_text) {
  umap <- Embeddings(obj, "umap")
  expr_vals <- FetchData(obj, vars = "GLS")[,1]
  
  df <- data.frame(
    UMAP_1 = umap[,1],
    UMAP_2 = umap[,2],
    GLS = expr_vals
  )
  
  n_expr <- sum(df$GLS > 0)
  n_total <- nrow(df)
  avg_expr <- round(mean(df$GLS), 3)
  
  label_text <- paste0(
    title_text,
    "\nCells (gene expression > 0): ", n_expr, "/", n_total,
    "\nAverage expression: ", avg_expr
  )
  
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = GLS)) +
    geom_point(size = 0.7) +
    scale_color_viridis_c(option = "plasma") +
    ggtitle(label_text) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

make_panel <- function(obj, celltype_name) {
  obj_ctrl <- subset(obj, subset = disease_status == "Uninfected control")
  obj_hiv  <- subset(obj, subset = disease_status == "HIV")
  
  p_ctrl <- make_umap_panel(obj_ctrl, paste0(celltype_name, " – Control"))
  p_hiv  <- make_umap_panel(obj_hiv,  paste0(celltype_name, " – HIV"))
  
  p_ctrl + p_hiv
}

  #create UMAP panels for CSF CD14 and CD16 monocytes stratified by disease status 
fig_cd14_csf <- make_panel(cd14_csf, "CD14 CSF")
fig_cd16_csf <- make_panel(cd16_csf, "CD16 CSF")
