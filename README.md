# ChAIR-Viewer
Visualization for single-cell multi-omic data


ChAIR-Viewer is developped for visualization single-cell multi-omic data. For ensemble data: bedgraph file (coverage signal), bedpe file (loop), bed file (peak/chromhmm/gene annotation) and more. For single cell data: fragment file. ChAIR-Viewer provides a variety of ways to sort cell, including: sorted by RNA UMI count, sorted by ATAC fragment count, sorted by contact/PET count, sorted by maximum contact distance, ordered of hierarchical clustering (in RNA/ATAC/PET data) and more.

ChAIR-Viewer is developped by R, some package are needed.
```
RColorBrewer
pheatmap
plotrix
dplyr
GenomicRanges
ggpubr
ggplot2
patchwork
```

# USAGE
#### 1. Only the two longest and two exon most transcripts are retained (Convenient output of gene annotation track) (Optional)
```
python getMaxlengthRNA.py ~/cellranger-7.1.0/genome/refdata-gex-mm10-2020-A/genes/genes.gtf | grep -v "^#" | awk '$3=="exon" || $3=="UTR"'|perl -lane '$id="#";$id=$1 if(/transcript_id "(.+?)"/);$name="#"; $name=$1 if(/gene_name "(.+?)"/); print "$F[0]\t".($F[3]-1)."\t$F[4]\t$id\t$name\t0\t$F[6]\t$F[2]"' > mm10.genome.max.bed
```
else
```
cat ~/cellranger-7.1.0/genome/refdata-gex-mm10-2020-A/genes/genes.gtf | grep -v "^#" | awk '$3=="exon" || $3=="UTR"'|perl -lane '$id="#";$id=$1 if(/transcript_id "(.+?)"/);$name="#"; $name=$1 if(/gene_name "(.+?)"/); print "$F[0]\t".($F[3]-1)."\t$F[4]\t$id\t$name\t0\t$F[6]\t$F[2]"' > mm10.genome.max.bed
```