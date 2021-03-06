# load libraries 
# SeuratDisk is a library to convert h5ad files into h5seurat files (Seurat readable format)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(pheatmap)


# The given h5ad files were created using an older annData (h5ad file) format
# to convert them to the latest annData version (h5ad file) use the python scanpy code uploaded to this git repo.
# the files I have sent you both are already updated and converted so you do not need to do the next block of code!

# # convert h5ad into h5seurat readable file (*note working directory for file will be different than mine)
# Convert("/home/matthew/datatransfer/mercier_hsc/droplet-marrow.h5ad", dest = "h5seurat", overwrite = FALSE)
# Convert("/home/matthew/datatransfer/mercier_hsc/facs-marrow.h5ad", dest = "h5seurat", overwrite = FALSE)
# 
# # import the now converted data files into Seurat
marrow.droplet <- LoadH5Seurat("/home/matthew/datatransfer/mercier_hsc/droplet-marrow.h5seurat")
marrow.facs <- LoadH5Seurat("/home/matthew/datatransfer/mercier_hsc/facs-marrow.h5seurat")

# the idents function allows you to choose which cell identity to work with
# for example, this data object has the clustering information from both leiden and louvain
# try changing them and running the dimension plot function and see what happens
Idents(marrow.droplet) <- "leiden"
Idents(marrow.facs) <- "leiden"

# dimension plot funciton, will plot pca, tsne, or umap...umap is default.
DimPlot(marrow.droplet, label = TRUE) +ggtitle("Droplet")
DimPlot(marrow.facs, label = TRUE) + ggtitle("Facs")

Idents(marrow.droplet) <- "louvain"
Idents(marrow.facs) <- "louvain"
DimPlot(marrow.droplet, label = TRUE) +ggtitle("Droplet")
DimPlot(marrow.facs, label = TRUE) + ggtitle("Facs")


# after running both dimension plots using leiden and louvain you will see a difference in number of clusters
# we can run further analysis on both identities (leiden and louvain) but for simplicity lets just pick one.

Idents(marrow.droplet) <- "louvain"
Idents(marrow.facs) <- "louvain"

# the clustering we see is all the timepoints combined. lets now subset the data into the timepoints we are interested in
# the subset fucntion allows us to seperate the data object on whatever metadata we like, in this case 3 months and 24 months.

droplet.3months <- subset(x = marrow.droplet, subset = age == "3m")
droplet.24months <-subset(x = marrow.droplet, subset = age == "24m")

facs.3months <- subset(x = marrow.facs, subset = age == "3m")
facs.24months <-subset(x = marrow.facs, subset = age == "24m")


# now we can look at the dimension plot of just 3 months or 24 months if we wanted to.
DimPlot(droplet.3months, label = TRUE) + ggtitle("droplet 3 months")
DimPlot(droplet.24months, label = TRUE) + ggtitle("droplet 24 months")



# to do differential expression between age groups we don't need to split the data up
# we can use the built in group.by function in seurat
# in our case, we will use group.by age and the first identity is "3 months" compared to the second identity which is "24 months"
# the settings for DE are somewhat default in order to just do this function fast...we may want to change the parameters for proper analysis
# this will give the DE across all genes found in each cluster, if we want to look at a specific cluster we can add subset.ident ="CLUSTER OF INTEREST" 
# to the function.

marrow.droplet.markers.age.test <- FindMarkers(marrow.droplet, ident.1 = "3m", ident.2="24m",group.by="age",
                         only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)

marrow.facs.markers.age <- FindMarkers(marrow.facs, ident.1 = "3m", ident.2="24m",group.by="age",
                    only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)



# once this is done we will have two data frame with all the DE genes. you will notice the genes are the index of the dataframe and not a column
# to fix this just run this code
# what this code does is takes all the genes from the index and makes it a column and resets the index to be numerical from 1 to x.
# the last part of the code just rearranges the order of the column so gene marker is the last column.
# the important columns besides gene marker is pct.1 and pct.2
# in our case pct.1 is age 3 months and pct.2 is age 24 months
# for example, if we have one gene that is pct.1 90% and pct.2 10% we can infer that the specific gene is more expressed in pct.1 which is timepoint 3 months


# for droplet markers
marrow.droplet.markers.age <- cbind(gene = rownames(marrow.droplet.markers.age), marrow.droplet.markers.age)
rownames(marrow.droplet.markers.age) <- 1:nrow(marrow.droplet.markers.age)
marrow.droplet.markers.age<- marrow.droplet.markers.age[, c(2,3,4,5,6,1)]

# for facs markers
marrow.facs.markers.age <- cbind(gene = rownames(marrow.facs.markers.age), marrow.facs.markers.age)
rownames(marrow.facs.markers.age) <- 1:nrow(marrow.facs.markers.age)
marrow.facs.markers.age<- marrow.facs.markers.age[, c(2,3,4,5,6,1)]


# if we wanted to now find ALL markers (per cluster) we can use the find all markers function
# this is where we would want to use the subset data by age
# this will give us another way to confirm results and if we want to look at clusters with interesting changes that we are unaware of 

# for droplet and facs markers at 3 months age.
marrow.droplets.markers.3m <- FindAllMarkers(droplet.3months, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
marrow.facs.markers.3m <- FindAllMarkers(facs.3months, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)

# for droplet and facs markers at 24 months age
marrow.droplets.markers.24m <- FindAllMarkers(droplet.24months, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
marrow.facs.markers.24m <- FindAllMarkers(facs.24months, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)


# if we want to see the level of expression of a gene on the dimension plot we can use this function

# if we use the marrow.droplet data object, it contains all the time points
FeaturePlot(marrow.droplet, features = c("Fos","Junb"),cols = c("grey", "red"))

# if we wanted to compare across timepoints we can use the split.by function

FeaturePlot(marrow.droplet, features = c("Fos"),cols = c("grey", "red"),split.by = "age")

# but as you can see, it makes the plot a little hard to read since there are 6 timepoints 
# we can use featureplot on our subset data to solve this

# if we save the feature plots as variables than we can use the + to create a facet plot

p1 <-FeaturePlot(droplet.3months, features = c("Fos"),cols = c("grey", "red")) + ggtitle("3 months")
p2<-FeaturePlot(droplet.24months, features = c("Fos"),cols = c("grey", "red")) + ggtitle("24 months")

#facet plot having both feature plots next to each other
p1+p2

# if we want to see what # cluster we can add label = TRUE

FeaturePlot(droplet.3months, features = c("Fos"),cols = c("grey", "red"), label=TRUE) + ggtitle("3 months")
FeaturePlot(droplet.24months, features = c("Fos"),cols = c("grey", "red"), label=TRUE) + ggtitle("24 months")

# now if we looked at our marrow.droplet.markers.age dataframe where we used the FindMarkers function comparing timepoints 3 and 24 we can see the decrease
# in fos which is shown on this plot

# If we wanted to compare a particular cluster between both time points we can use this function
# it is the same FindMarkers function, except we will group by timepoints (3m and 24m), and the subset.ident is the cluster we want to compare
# for this example, i want to compare cluster 14, so pct.1 will be cluster 14 at 3m and pct.2 will be cluster 14 at 24m

marrow.droplet.markers.clusterCompare <- FindMarkers(marrow.droplet, ident.1 = "3m", ident.2="24m",group.by="age",
                            only.pos = FALSE, min.pct = -Inf,logfc.threshold = 0,subset.ident = "14")

# our little hacky fix to make the gene names of the index a column of genes
marrow.droplet.markers.clusterCompare <- cbind(gene = rownames(marrow.droplet.markers.clusterCompare), marrow.droplet.markers.clusterCompare)
rownames(marrow.droplet.markers.clusterCompare) <- 1:nrow(marrow.droplet.markers.clusterCompare)
marrow.droplet.markers.clusterCompare<- marrow.droplet.markers.clusterCompare[, c(2,3,4,5,6,1)]

# if we wanted to overlay the dimplot between age groups 3 and 24 and use the group.by feature we have to subset the Seurat object again
# this time we want to subset the data by age 3m or 24m (this will give us both 3m and 24m)


droplet.3months.24months <- subset(x = marrow.droplet, subset = age == "3m" | age == "24m")

# now we can use the dimplot group.by function
DimPlot(droplet.3months.24months, group.by ="age")

# how to label clusters
# first we want to save the original cluster numbers just in case, so we will copy it to a new column in the meta data
# I named it old.ident.louvain just because it is the old cluster identity from the louvain algorithm 
# you have options, you can do this for the marrow droplet and facs objects before subsetting so when you subset the labels are already there
# if you dont want to redo that you can just do this for each object and just change the object variable name
# keep in mind facs will have different clustering labels!

marrow.droplet[["old.ident.louvain"]] <- Idents(object = marrow.droplet)

# now we make a vector of the labels, going from 0 to 25
# some repeating labels need to have an additional label so they remain unique, you can't have the same label repeating twice
# please note I just picked the first label for each cluster from the metadata, you might need to go through this carefully and figure out what the label should really be

new.cluster.ids <- c("granulocyte_1", 'granulocyte_2','granulocyte_3','erythroblast','promonocyte','granulocytopoietic cell_1',
                     'granulocytopoietic cell_2','granulocytopoietic cell_3','naive T cell','hematopoietic precursor cell_1','monocyte_1',
                     'monocyte_2','macrophage','granulocytopoietic cell_4','proerythroblast','granulocytopoietic cell_5','hematopoietic precursor cell_2',
                     'immature B cell','granulocyte_3','granulocytopoietic cell_6','precursor B cell','granulocyte_4','megakaryocyte-erythroid progenitor cell',
                     'NK cell','granulocytopoietic cell_5','basophil')

# this renames and sets the new identity
names(new.cluster.ids) <- levels(marrow.droplet)
marrow.droplet <- RenameIdents(marrow.droplet, new.cluster.ids)

# keep note of these new arguments, repel and label size, this allows control of the cluster labels
# repeal makes it such that it won't overlap the over labels, and label.size changes the font size.
# noLegend() removes the side legend which you do not need anymore if you label directly on the clusters
DimPlot(marrow.droplet, label = TRUE, repel=TRUE, label.size = 3) + NoLegend()

# dimplot zooming in on the clusters 9 and 16
# we already have a seurat object that subset the droplet marrow object into just 3m and 24m age groups
# now we can subset this again just to have clusters 9 and 16 (on my seurat object i did not rename anything yet so they are still 9 and 16)
# your seurat object has been renamed to hematopoietic precursor cell_1 and hematopoietic precursor cell_2

c.9.16 <- subset(droplet.3months.24months, idents = c("9","16") )

# shows which cells belong to which age group
DimPlot(c.9.16,group.by ="age")

# both age group clustering side by side
DimPlot(c.9.16,split.by ="age",label=TRUE)

# both age groups together
DimPlot(c.9.16,label=TRUE)


# subclustering cluster 9 from ages 3/24m, idents is 9 for cluster 9 (my object doesn't have the new labels so make sure to change this!)
cluster9.droplet <- subset(droplet.3months.24months, idents = c("9") )

# plot to make sure the cluster is fine
DimPlot(cluster9.droplet,split.by ="age",label=TRUE)

# now we need to re-do the analyze from the beginning, skipping the QC filtering steps since this object is already filtered.

# this function finds high cell-to-cell variation
# the parameters used here are the default from seurat
cluster9.droplet <- FindVariableFeatures(cluster9.droplet, selection.method = "vst", nfeatures = 2000)

# scale the data
all.genes <- rownames(cluster9.droplet)
cluster9.droplet <- ScaleData(cluster9.droplet, features = all.genes)

# perform PCA
cluster9.droplet <- RunPCA(cluster9.droplet, features = VariableFeatures(object = cluster9.droplet))

# look at the PCA dimplot now
DimPlot(cluster9.droplet, reduction = "pca")

# now we do a elbow plot to determine how many components we should include in our analysis 
# Where the bend (elbow) occurs gives us a good idea of which components hold the most variability (the important data!)
ElbowPlot(cluster9.droplet)

# we choose the dims based on the elbow plot, the resolution paramter sometimes requires tuning to give the best clustering results
# this part is always tricky and sometimes you may need to change the dims a few times, but from looking at the elbow plot lets go with 12 dims


# default algorithm is louvain (which is what we were using in the previous analysis!)
# also by looking at their results in the params section of the reduction table within the seurat object, they used a resolution of 1
# we will do the same and see how it looks
cluster9.droplet <- FindNeighbors(cluster9.droplet, dims = 1:12)
cluster9.droplet <- FindClusters(cluster9.droplet, resolution = 1)

# now we can run a umap and tsne, keep the dims parameter the same value from FindNeighbours
cluster9.droplet <- RunUMAP(cluster9.droplet, dims = 1:12)
cluster9.droplet <- RunUMAP(cluster9.droplet, dims = 1:12)

# look at the UMAP dimplot now
DimPlot(cluster9.droplet, label= TRUE)


# when looking at the metadata becareful! louvain and leiden don't mean anything anymore...the important column is seurat_clusters!

# dimplot split by age
DimPlot(cluster9.droplet, label= TRUE, split.by ="age")

# now we can fun the findallmarkers function and see if we can label the clusters (note this is using data from 3m/24m age groups!)
# you might have to adjust these parameters
cluster9.droplet.markers <- FindAllMarkers(cluster9.droplet, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.25)

# you can run the findMarkers function next to look between clusters and look between age groups

# heatmap

# we will use the droplet.3months.24months object since we want to compare those months and use subset in the heatmap function to get cluster 9
# remember its a different label for you!

# In order to use a heatmap function we need to scale the data / find variable features
# this is not required for DEG, but we scale the data and use the FindVariablefeatures when we want to cluster/heatmap stuff.

droplet.3months.24months <- FindVariableFeatures(droplet.3months.24months, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(droplet.3months.24months)
droplet.3months.24months <- ScaleData(droplet.3months.24months, features = all.genes)


# now we take a list of genes we want to use as features, since they are not common seperated we need to do some manipulation
# first we put it into a vector which will keep it seperated by new lines

gene_list <- c("Fosb
Dusp1
Eif2s3y
Neat1
Fos
Klf6
Jun
Ddx3y
Tsix
Egr1
Dhx9
Ppp1r10
Slc7a5
Cct6a
2610018G03Rik
Srp54a
Ddx23
Supt6h
Aff4
Atf7ip
Ttc37
Erbb2ip
Elf2
Ccl3
Casp3
Myl10
Ctdspl2
Cd69
Brca2
A630007B06Rik
Tdg
Tnfaip3
Kif13b
D14Abb1e
Lcorl
2810008D09Rik
AI314180
Ier3
Thbs1
Tob2
Hsph1
Dio2
Cntln
Tec
Mcph1
Gm13034
Mreg
D830031N03Rik
Dhx29
Kdm5c
BC018507
Bmp2k
Nup98
Zbed6
Tnfrsf13b
Vps13a
Serpina3g
Hist2h2ac
Rn45s
Btbd7
Zbtb38
Hmga1
Tmem88
Rev3l
Pola1
Nupr1
Kdm5d
Vegfa
Mfap3l")

# now we split the gene_list by "\n" which is the escape character letting R and the compiler to know it is a new line
gene_list <-strsplit(gene_list, "\n")
# we use unlsit to break down the list into just a vector of characters! exactly what we want!
gene_list<-unlist(gene_list)

# for some reason the DoHeatMap function is weird and doesn't ignore the empty age group variables (3m,18m,etc) so we need to remove them from the 
# factor levels...R has a function called droplevels for this exact issue

droplet.3months.24months@meta.data$age<-droplevels(droplet.3months.24months@meta.data$age) 

# now we can use the seurat doheatmap function, inside we have our subset function to select cluster 9 and we can group.by age
# now you will notice it looks a bit hard to read...we can use the downsample function of seurat's subset to select a sample of cells
# they recommend this for heatmaps so its easier to read...not sure exactly how their downsample function is performed but we can look into that!
# I gave a few examples to compare no downsample to different levels of downsampling.

DoHeatmap(subset(droplet.3months.24months,idents = c("9")) ,features=gene_list, group.by ="age")
DoHeatmap(subset(droplet.3months.24months,idents = c("9"), downsample=300) ,features=gene_list, group.by ="age")
DoHeatmap(subset(droplet.3months.24months,idents = c("9"), downsample=100) ,features=gene_list, group.by ="age")
DoHeatmap(subset(droplet.3months.24months,idents = c("9"), downsample=50) ,features=gene_list, group.by ="age")

# this block of code is to look at the average gene expression to double check that DEG's make sense

# basic feature plot split by age
FeaturePlot(marrow.droplet, features = c("Fosb","Jun"),cols = c("grey", "red"),split.by = "age") 
FeaturePlot(droplet.24months, features = c("Fosb","Jun"),cols = c("grey", "red"))

# calculating the average gene expression per cell in each cluster 
# looking at 3m and 24 seurat objects
avg.droplet.3months<- log1p(AverageExpression(droplet.3months, verbose = FALSE)$RNA)
avg.droplet.3months$gene <- rownames(avg.droplet.3months)

avg.droplet.24months<- log1p(AverageExpression(droplet.24months, verbose = FALSE)$RNA)
avg.droplet.24months$gene <- rownames(avg.droplet.24months)

# summing the rows of the gene of interest to see if there is a big difference in total average gene expression 
sum(avg.droplet.3months["Fosb",1:26])
sum(avg.droplet.24months["Fosb",1:26])

sum(avg.droplet.3months["Jun",1:26])
sum(avg.droplet.24months["Jun",1:26])


# subset facs dataset clusters 1,3,9 to recluster together
# first we need to subset the facs dataset to just have 3m/24m age groups
facs.3months.24months <- subset(x = marrow.facs, subset = age == "3m" | age == "24m")
# now we subset again but only grabbing the clusters we want
facs.c.1.3.9 <- subset(facs.3months.24months, idents = c("1","3","9") )

# now we wil repeat the exact same steps we did when we re-clustered cluster 9 from the droplet dataset
# plot to make sure the cluster is fine
DimPlot(facs.c.1.3.9,split.by ="age",label=TRUE)

# now we need to re-do the analyze from the beginning, skipping the QC filtering steps since this object is already filtered.


facs.c.1.3.9 <- FindVariableFeatures(facs.c.1.3.9, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(facs.c.1.3.9)
facs.c.1.3.9 <- ScaleData(facs.c.1.3.9, features = all.genes)

facs.c.1.3.9 <- RunPCA(facs.c.1.3.9, features = VariableFeatures(object = facs.c.1.3.9))

DimPlot(facs.c.1.3.9, reduction = "pca")

# take a look at the elbow plot and see how many dimensions you want to keep (i think 10-14 is okay for this one)
ElbowPlot(facs.c.1.3.9)

# adjust the resolution to see how it changes the clustering...at 0.8 it seems okay but clusters 10 and 11 seem a bit strange!
facs.c.1.3.9 <- FindNeighbors(facs.c.1.3.9, dims = 1:12)
facs.c.1.3.9 <- FindClusters(facs.c.1.3.9, resolution = 0.8)

# now we can run a umap and tsne, keep the dims parameter the same value from FindNeighbours
facs.c.1.3.9 <- RunUMAP(facs.c.1.3.9, dims = 1:12)

# look at the UMAP dimplot now
DimPlot(facs.c.1.3.9, label= TRUE)


DimPlot(facs.c.1.3.9, label= TRUE, split.by ="age")


facs.c.1.3.9.markers <- FindAllMarkers(facs.c.1.3.9, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.25)

