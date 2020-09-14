data<-read.table("matrix.txt",h=T,row.names=1)

meta<-read.table("meta.txt",h=T)

head(data[1:5,1:5])
head(meta)


hsc <- CreateSeuratObject(counts = data, project = "aging", min.features = 200,min.cell=0)

hsc <- AddMetaData(object = hsc,  metadata = meta)

x = hsc[[]]


## cell cycle analysis

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


##scaling without normalization because TPM
hsc <- FindVariableFeatures(hsc, selection.method = "vst")
hsc <- ScaleData(hsc, features = rownames(hsc),check.for.norm = F)


hsc <- RunPCA(hsc, features = VariableFeatures(hsc), ndims.print = 1:10, nfeatures.print = 10)

DimHeatmap(hsc, dims = c(1, 2))

head(hsc[[]],n=30)

RidgePlot(hsc, features = c("CCNA2", "CCNB2", "CCNB1", "CCND2"), ncol = 2)


hsc <- RunPCA(hsc, features = VariableFeatures(hsc),npcs = 50, verbose = FALSE)

ElbowPlot(hsc,n=50)
DimPlot(hsc, reduction = "pca", group.by="group",pt.size=2)

hsc <- RunUMAP(hsc, reduction = "pca", dims = 1:20)
DimPlot(hsc, reduction = "umap", group.by="classe",pt.size=2)

hsc <- RunTSNE(hsc, reduction = "pca", dims = 1:10)
DimPlot(hsc, reduction = "tsne", group.by="group",pt.size=2.5)
DimPlot(hsc, reduction = "tsne", group.by="group",pt.size=2.5)
DimPlot(hsc, reduction = "umap", group.by="group",pt.size=2.5)

FeaturePlot(hsc, features = c("EGR1"),min.cutoff = "q9",cols=c("#CCFFFF","darkblue"),split.by= "orig.ident",reduction = "tsne",pt.size=2)
FeaturePlot(hsc, features = c("CDK6"),min.cutoff = "q9",cols=c("#CCFFFF","darkblue"),split.by= "orig.ident",reduction = "umap",pt.size=2)

#set identities with meta column
Idents(object = hsc) <- 'group'
Idents(object = hsc)

## marrow regression on age
marrow <- ScaleData(hsc, vars.to.regress = "age", features = rownames(hsc))
marrow <- RunPCA(marrow, features = VariableFeatures(marrow))
DimPlot(marrow,group.by="group",pt.size=2.5,reduction = "tsne")
DimPlot(hsc,group.by="group",pt.size=2.5,reduction = "tsne")
DimPlot(hsc,group.by="classe",pt.size=2.5,reduction = "tsne")


DimPlot(marrow,group.by="group",pt.size=2.5,reduction = "pca")

## differential expressed genes
#set identities with meta column
Idents(object = hsc) <- 'classe'
Idents(object = hsc)

markers <- FindMarkers(hsc, ident.1 = "old", ident.2 = "young")

VlnPlot(hsc, features = c("EGR1"), slot = "counts", log = TRUE,pt.size=1,split.by="group",col=c("darkred","red","darkgreen","green"))