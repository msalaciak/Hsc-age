# load libraries 
# SeuratDisk is a library to convert h5ad files into h5seurat files (Seurat readable format)
library(Seurat)
library(SeuratDisk)
library(ggplot2)


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

marrow.droplet.markers.age <- FindMarkers(marrow.droplet, ident.1 = "3m", ident.2="24m",group.by="age",
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