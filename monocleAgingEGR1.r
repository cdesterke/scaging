library(monocle)
mat<-as.matrix(data)
colonnes<-as.data.frame(colnames(mat))
colonnes$group<-meta$classe

colnames(colonnes)<-c("id","group")
row.names(colonnes)<-colonnes$id
pd <- new("AnnotatedDataFrame", data = colonnes)


genes<-as.data.frame(row.names(data))
colnames(genes)<-"gene_short_name"
row.names(genes)<-genes$gene_short_name
fd <- new("AnnotatedDataFrame", data = genes)
cds <- newCellDataSet(mat, phenoData = pd,featureData = fd,expressionFamily=tobit())
cth <- newCellTypeHierarchy()

EGR1_id <- row.names(subset(fData(cds), gene_short_name == "EGR1"))


cth <- addCellType(cth, "EGR1_low", classify_func = function(x) { x[EGR1_id,] <= 3 })
cth <- addCellType(cth, "EGR1_medium", classify_func = function(x) { x[EGR1_id,] < 8 & x[EGR1_id,] > 3} )
cth <- addCellType(cth, "EGR1_high", classify_func = function(x) { x[EGR1_id,] >= 8 })

cds <- classifyCells(cds, cth, 0.1)

table(pData(cds)$CellType)


my_feat <- fData(cds)
my_feat$id<-my_feat$gene_short_name
head(my_feat)

plot_pc_variance_explained(cds, return_all = FALSE)
cds <- reduceDimension(cds, max_components = 2, num_dim = 30, reduction_method = 'tSNE', verbose = TRUE)


cds <- clusterCells(cds)
table(pData(cds)$CellType)
#ggplot graph cell type
library(ggplot2)

pie <- ggplot(pData(cds),
aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x = element_blank(), axis.title.y = element_blank())


#cluster cells with hiearchy
cds <- clusterCells(cds,cth) 

plot_cell_clusters(cds, 1, 2, color = "Cluster") +  facet_wrap(~CellType)
plot_cell_clusters(cds, 1, 2, color = "CellType",cell_size=2 )

save(cds,file="monoclehiearchy.rda")

#differential expressed genes against cell type
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
print(head(pData(cds)))

expressed_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~CellType")

sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))

my_ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
cds2 <- setOrderingFilter(cds, ordering_genes = my_ordering_genes)

#dimensional reduction and gene selection
cds2 <- reduceDimension(cds, method = 'DDRTree')
gene_to_cluster <- row.names(diff_test_res)[order(diff_test_res$qval)][1:100] 
conca<-c("EGR1",gene_to_cluster)
cds2 <- orderCells(cds2)


#build graphs on trajectory 
plot_cell_clusters(cds2, color_by = 'as.factor(CellType)',markers="EGR1",cell_size=2 )+  facet_wrap(~CellType)
plot_cell_clusters(cds2, color_by = 'as.factor(CellType)',markers="EGR1",cell_size=2 )+  facet_wrap(~group)

plot_cell_trajectory(cds2, color_by = "group",cell_size=2)
plot_cell_trajectory(cds2, color_by = "group",markers="EGR1",markers_linear = TRUE,show_branch_points=FALSE)
plot_cell_trajectory(cds2, color_by = "Pseudotime")
plot_cell_trajectory(cds2, color_by = "CellType") +  facet_wrap(~group)
my_pseudotime_cluster <- plot_pseudotime_heatmap(cds2[gene_to_cluster,],cores = 8,
show_rownames = TRUE,return_heatmap = TRUE,cluster_rows = TRUE)
plot_genes_in_pseudotime(cds2[c("CD9","LY6D","LAMA3","ANTXR1","DPT","CTGF","PLPP3","JAM2","SERPINE1"),],
cell_size = 2, color_by = "group",ncol = 2)
#save data
write.table(pData(cds),file="phenotype.txt")
write.table(diff_test_res,file="difexpEGR1trajectory.txt")
save(cds,file="monoclecds.rda")
save(cds2,file="monoclecdsfiltre.rda")
### end of the code





