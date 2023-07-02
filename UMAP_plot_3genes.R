library(Seurat)
library(ggplot2)
library(scales)
library(ggnewscale)

data.fibro <- readRDS("/path/to/Seurat_merged.rds")

# example of the FeaturePlot "blend" function
FeaturePlot(data.fibro, features = c("PLN", "CNN1"), blend = TRUE)


# extract a dataframe with UMAP coords and gene expression for the genes: PLN, CNN1, CCDC102B
# fetch data
umapCoord <- as.data.frame(Embeddings(object = data.fibro[["umap"]]))
gene1 <- FetchData(data.fibro, vars = "PLN")
gene2 <- FetchData(data.fibro, vars = "CNN1")
gene3 <- FetchData(data.fibro, vars = "CCDC102B")

# 
umapCoord$PLN <- gene1
umapCoord$CNN1 <- gene2
umapCoord$CCDC102B <- gene3

# Create DF with UMAP coords and gene expression columns
df.fibro <- data.frame(
  UMAP_1 = umapCoord$UMAP_1,
  UMAP_2 = umapCoord$UMAP_2,
  PLN = FetchData(data.fibro, vars = "PLN"),
  CNN1 = FetchData(data.fibro, vars = "CNN1"),
  CCDC102B = FetchData(data.fibro, vars = "CCDC102B"))


# Re-scale gene expression for coloring the UMAP
df.fibro2 <- cbind.data.frame(
  df.fibro,
  colour = rgb(
    rescale(df.fibro$PLN),
    rescale(df.fibro$CNN1),
    rescale(df.fibro$CCDC102B)
  )
)

# Plot the data: Approach 1
ggplot(df.fibro2, aes(UMAP_1, UMAP_2, colour = colour)) +
  geom_point() +
  scale_colour_identity() +
  new_scale_colour() +
  # shape = NA --> invisible layers
  geom_point(aes(colour = PLN), shape = NA) +
  scale_colour_gradient(low = "black", high = "red") +
  new_scale_colour() +
  geom_point(aes(colour = CNN1), shape = NA) +
  scale_colour_gradient(low = "black", high = "green") +
  new_scale_colour() +
  geom_point(aes(colour = CCDC102B), shape = NA) +
  scale_colour_gradient(low = "black", high = "blue")

# Plot the data: Approach 2
####### do not run (only for installing the library)
# library(devtools)
# devtools::install_github("teunbrand/ggchromatic")
####### 
library(ggchromatic)

# returns: Error in if (0 <= angle & angle < 90) { : the condition has length > 1 (needs debugging)
ggplot(df.fibro, aes(UMAP_1, UMAP_2, colour = rgb_spec(PLN, CNN1, CCDC102B))) +
  geom_point()


