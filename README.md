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

#### 2. Prepare all data (Convenient for multiple calls)
The readInfo.R script is used to organize the input files, and needs to be modified by the user according to your actual situation.
```
mkdir alldata
Rscript readInfo.R mm10.genome.max.bed mouseBrain mouseBrain.chromstat.hmm.bed
```
The program searches the current directory for the following files and stores all results as an RData object.

file|description
:---:|:---:
mouseBrain.clean.ATAC.mouseBrain.fragments.unsorted.tsv|scATAC fragment file
mouseBrain.clean.ATAC.mouseBrain.q001_peaks.final.bed|peak file (3 columns)
mouseBrain.clean.ATAC.mouseBrain.bedgraph|scATAC coverage signal
mouseBrain.clean.ATAC.mouseBrain.RPKM.bedgraph|scATAC coverage signal by using RPKM normalization
mouseBrain.clean.PET.bedpe.mouseBrain.bedpe|scPET file (11 columns: chrom1, start1, end1, chrom2, start2, end2, barcode, readID, score, strand1, strand2)
mouseBrain.clean.PET.bedpe.mouseBrain.clusters.bedpe|loop file (no given anchor model) (7 columns: chrom1, start1, end1, chrom2, start2, end2, contact count)
mouseBrain.clean.PET.bedpe.mouseBrain.clusters.twopeaks.bedpe|loop file (no given anchor model) (overlap with peak) (7 columns: chrom1, start1, end1, chrom2, start2, end2, contact count)
mouseBrain.clean.PET.bedpe.mouseBrain.bedpe.ipet.loops|loop file ( given anchor model) (10 columns: chrom1, start1, end1, chrom2, start2, end2, contact count, anchor1 coverage, anchor2 coverage, VC_sqrt score)
mouseBrain.clean.GEX.mouseBrain.bedgraph|scRNA coverage signal
mouseBrain.clean.GEX.mouseBrain.RPKM.bedgraph|scRNA coverage signal by using RPKM normalization
mouseBrain.clean.dup.GEX.mouseBrain.bedgraph|scRNA (contain mapped reads from intergenic) coverage signal. (If there is no mouseBrain.clean.dup.GEX.mouseBrain.bedgraph file, you can copy the mouseBrain.clean.GEX.mouseBrain.bedgraph file)
mouseBrain.clean.dup.GEX.mouseBrain.RPKM.bedgraph|scRNA (contain mapped reads from intergenic) coverage signal by using RPKM normalization. (If there is no mouseBrain.clean.dup.GEX.mouseBrain.RPKM.bedgraph file, you can copy the mouseBrain.clean.GEX.mouseBrain.RPKM.bedgraph file)
mouseBrain.clean.GEX.mouseBrain.fragments.unsorted.tsv|scRNA fragment file
mouseBrain.clean.dup.GEX.mouseBrain.fragments.unsorted.tsv|scRNA (contain mapped reads from intergenic) fragment file. (If there is no mouseBrain.clean.dup.GEX.mouseBrain.fragments.unsorted.tsv file, you can copy the mouseBrain.clean.GEX.mouseBrain.fragments.unsorted.tsv file)

#### 3. plot
```
Rscript plotData.R
```

#### 4. The ChAIR-viewer version of shiny is will be finished soon. 

# CONTACT
黄加祥 (Jiaxiang Huang, 12207142@zju.edu.cn)<br/>
黄星宇 (Xingyu Huang, xingyu.huang@zju.edu.cn/huang182@live.cn)<br/>