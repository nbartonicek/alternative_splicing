library(DEXSeq)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(BiocParallel)

setwd("/researchers/nenad.bartonicek/projects/PRMT1/scripts")
txdb<-makeTxDbFromGFF("../annotation/Mus_musculus.GRCm39.110.gtf.gz")

flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
names(flattenedAnnotation) =
  sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)

bamFiles = list.files("../results/star/", pattern="sortedByCoord.out.bam$",full.names=T)
bamFiles_df<-data.frame(files=bamFiles,Code=gsub("_.*","",basename(bamFiles)))

metadata<-read.csv("../annotation/metadata.csv") %>% 
  mutate(treatment=paste0(Genotype,"_",Condition)) %>%
  filter(treatment %in% c("WT_untreated","Prmt1ko_untreated")) %>%
  left_join(.,bamFiles_df)

bamFiles = BamFileList( bamFiles )

if(!file.exists("se.Rdata")){
se = summarizeOverlaps(
  flattenedAnnotation, BamFileList(bamFiles), singleEnd=FALSE,
  fragments=TRUE, ignore.strand=TRUE)
save(se,file="se.Rdata")
} else {
  load("se.Rdata")
}

colData(se)$condition = factor(c("WT_untreated", "Prmt1ko_untreated","WT_untreated","Prmt1ko_untreated","Prmt1ko_untreated","WT_untreated"))

options(MulticoreParam=MulticoreParam(workers=6))
BPPARAM = MulticoreParam(6)

dxd = DEXSeqDataSetFromSE( se, design= ~ sample + exon + condition:exon )
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM,fitType='local')
pdf("dispersions.pdf",width=6,height=5)
plotDispEsts( dxd )
dev.off()
dxd = testForDEU( dxd, BPPARAM=BPPARAM)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

save(dxd,file="dxd.Rdata")
plotDispEsts( dxd )



dxr1 = DEXSeqResults( dxd )

symbols<-AnnotationDbi::select(org.Mm.eg.db, keys=dxr1$groupID, columns="SYMBOL", keytype="ENSEMBL") %>%
  distinct() 

filtered<-dxr1 %>% 
  as.data.frame() %>% 
  filter(padj<0.1) %>%
  left_join(.,symbols,by=c("groupID"="ENSEMBL"))

IFN_response_activators<-c("Ticam1","Bud13","Zfr","Pml")

for(activator in IFN_response_activators){

  pdf(paste0(activator,".pdf"),width=6,height=5)
  id<-symbols$ENSEMBL[symbols$SYMBOL %in% activator]
  plotDEXSeq( dxr1, id, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
  dev.off()

}


pdf("Ccnd2.pdf",width=6,height=5)
plotDEXSeq( dxr1, "ENSMUSG00000000184", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

                               