library(monocle)
mat<-as.matrix(data)
colonnes<-as.data.frame(colnames(mat))
colonnes$group<-meta$group

colnames(colonnes)<-c("id","group")
row.names(colonnes)<-colonnes$id
pd <- new("AnnotatedDataFrame", data = colonnes)


genes<-as.data.frame(row.names(data))
colnames(genes)<-"gene_short_name"
row.names(genes)<-genes$gene_short_name
fd <- new("AnnotatedDataFrame", data = genes)
cds <- newCellDataSet(mat, phenoData = pd,featureData = fd)