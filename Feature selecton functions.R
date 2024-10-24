scmarker.fs <- function(data){
    data <- log(as.matrix(data)+1)
    res=ModalFilter(data=data,geneK=10,cellK=10)
    res=GeneFilter(obj=res)
    res=getMarker(obj=res)
    return(res$marker)
}

# weighted
m3drop.fs <- function(data){
    norm <- M3DropConvertData(data, is.counts=TRUE)
    M3Drop_genes <- M3DropFeatureSelection(norm,suppress.plot=TRUE)
    return(M3Drop_genes$Gene)
}

# weighted
nbdrop.fs <- function(data){
    count_mat <- NBumiConvertData(data, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(count_mat)
    NBDropFS <- NBumiFeatureSelectionCombinedDrop(DANB_fit, qval.thresh = 0.05, suppress.plot=TRUE)
    return(NBDropFS$Gene)
}

hvg.fs <- function(data){
    data = CreateSeuratObject(counts=data)
    data <- NormalizeData(object = data)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    return(VariableFeatures(data))
}

sct.fs <- function(data){
    data = CreateSeuratObject(counts=data)
    data <- SCTransform(data, verbose = FALSE, variable.features.n = 3000)
    return(VariableFeatures(data))
}

genebasisr.fs=function(data){
  data = as.data.frame(data)
  sce=SingleCellExperiment(assays = list(counts = data, logcounts = log2(data+1)))
  sce = retain_informative_genes(sce)
  return(gene_search(sce, n_genes_total = 50, verbose = FALSE)$gene)
}

#weighted
feast.fs=function(data, k, n_feature=2000){
  # n_feature=4000
  Y = process_Y(data)
  # ixs = FEAST(Y, k=k)
  ixs  = FEAST_fast(Y, k=k)
  return(rownames(Y)[ixs][1:n_feature])
}

hrg.fs=function(data){
  data = CreateSeuratObject(counts=data)
  data <- NormalizeData(data,verbose = FALSE)
  all.genes = rownames(data)
  data <- ScaleData(data, features = all.genes, verbose = FALSE)
  if(ncol(data)<50){
    data <- RunPCA(data, features = all.genes,verbose = FALSE,npcs=20)
  }else{
    data <- RunPCA(data, features = all.genes,verbose = FALSE)
  }
  data=FindRegionalGenes(data)
  gene_num = HRG_elbowplot(data)
  return(RegionalGenes(data,nfeatures = gene_num))
}

lpr.fs <- function(data){
    df <- ave_poi_ori(counts)
    result <- loess.reg(df)
    return(result)
}
