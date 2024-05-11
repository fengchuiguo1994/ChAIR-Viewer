args=commandArgs(T)

readfile = function(infile) {
  if (length(grep(".gz$",infile)) || length(grep(".GZ$",infile)) || length(grep(".gzip$",infile)) || length(grep(".GZIP$",infile))) {
    gf = gzfile(infile,'rt')
  } else {
    gf = infile
  }
  return(gf)
}


# args=c("mm10.genome.max.bed","P11_Astroependymal")
## gene annotation
genefile = args[1] # "mm10.genome.max.bed"
gene = read.table(readfile(genefile),header=F,stringsAsFactors=F)
names(gene)<-c("chrom","start","end","transcription","gene","score","strand","type")
gene$ID = paste(gene$gene,gene$transcription,sep="-")
genelist <- split(gene, gene$chrom)


## ATAC
scATACfragmentfile = paste0("mouseBrain.clean.ATAC.",args[2],".fragments.unsorted.tsv") # mouseBrain.clean.ATAC.P11_Astroependymal.fragments.unsorted.tsv
scATACfragment = read.table(readfile(scATACfragmentfile),header=F,stringsAsFactors=F)
names(scATACfragment)<-c("chrom","start","end","barcode","count")
scATACfragmentlist <- split(scATACfragment, scATACfragment$chrom)

scATACpeakfile = paste0("mouseBrain.clean.ATAC.",args[2],".q001_peaks.final.bed") # mouseBrain.clean.ATAC.P11_Astroependymal.q001_peaks.final.bed
scATACpeak = read.table(readfile(scATACpeakfile),header=F,stringsAsFactors=F)
names(scATACpeak)<-c("chrom","start","end")
scATACpeaklist <- split(scATACpeak, scATACpeak$chrom)

scATACARCCOVfile = paste0("mouseBrain.clean.ATAC.",args[2],".bedgraph") # mouseBrain.clean.ATAC.P11_Astroependymal.bedgraph
scATACARCCOV = read.table(readfile(scATACARCCOVfile),header=F,stringsAsFactors=F)
names(scATACARCCOV)<-c("chrom","start","end","value")
scATACARCCOVlist <- split(scATACARCCOV, scATACARCCOV$chrom)

scATACARCCOVrpkmfile = paste0("mouseBrain.clean.ATAC.",args[2],".RPKM.bedgraph") # mouseBrain.clean.ATAC.P11_Astroependymal.RPKM.bedgraph
scATACARCCOVrpkm = read.table(readfile(scATACARCCOVrpkmfile),header=F,stringsAsFactors=F)
names(scATACARCCOVrpkm)<-c("chrom","start","end","value")
scATACARCCOVrpkmlist <- split(scATACARCCOVrpkm, scATACARCCOVrpkm$chrom)


## PET
scATACPETfile = paste0("mouseBrain.clean.PET.bedpe.",args[2],".bedpe") # mouseBrain.clean.PET.bedpe.P11_Astroependymal.bedpe
scATACPET = read.table(readfile(scATACPETfile),header=F,stringsAsFactors=F)
names(scATACPET)<-c("chrom1","start1","end1","chrom2","start2","end2","barcode","name","score","strand1","strand2")
scATACPET = scATACPET[scATACPET$chrom1 == scATACPET$chrom2 & scATACPET$end2 - scATACPET$start1 >= 1000,]
scATACPETlist <- split(scATACPET, scATACPET$chrom1)

scATACloopfile = paste0("mouseBrain.clean.PET.bedpe.",args[2],".clusters.bedpe") # mouseBrain.clean.PET.bedpe.P11_Astroependymal.clusters.bedpe
scATACloop = read.table(readfile(scATACloopfile),header=F,stringsAsFactors=F)
names(scATACloop)<-c("chrom1","start1","end1","chrom2","start2","end2","count")
scATACloop = scATACloop[scATACloop$chrom1 == scATACloop$chrom2 & scATACloop$count > 0,]
scATAClooplist <- split(scATACloop, scATACloop$chrom1)

scATAClooptwopeakfile = paste0("mouseBrain.clean.PET.bedpe.",args[2],".clusters.twopeaks.bedpe") # mouseBrain.clean.PET.bedpe.P11_Astroependymal.clusters.twopeaks.bedpe
scATAClooptwopeak = read.table(readfile(scATAClooptwopeakfile),header=F,stringsAsFactors=F)
names(scATAClooptwopeak)<-c("chrom1","start1","end1","chrom2","start2","end2","count")
scATAClooptwopeak = scATAClooptwopeak[scATAClooptwopeak$chrom1 == scATAClooptwopeak$chrom2 & scATAClooptwopeak$count > 0,]
scATAClooptwopeaklist <- split(scATAClooptwopeak, scATAClooptwopeak$chrom1)

scATACloopgivenanchorfile = paste0("mouseBrain.clean.PET.bedpe.",args[2],".bedpe.ipet.loops") # mouseBrain.clean.PET.bedpe.P11_Astroependymal.bedpe.ipet.loops
scATACloopgivenanchor = read.table(readfile(scATACloopgivenanchorfile),header=F,stringsAsFactors=F)
names(scATACloopgivenanchor)<-c("chrom1","start1","end1","chrom2","start2","end2","count","cov1","cov2","score")
scATACloopgivenanchor = scATACloopgivenanchor[scATACloopgivenanchor$chrom1 == scATACloopgivenanchor$chrom2 & scATACloopgivenanchor$count > 0,]
scATACloopgivenanchorlist <- split(scATACloopgivenanchor, scATACloopgivenanchor$chrom1)


## RNA
scATACGEXCOVfile = paste0("mouseBrain.clean.GEX.",args[2],".bedgraph") # mouseBrain.clean.GEX.P11_Astroependymal.bedgraph
scATACGEXCOV = read.table(readfile(scATACGEXCOVfile),header=F,stringsAsFactors=F)
names(scATACGEXCOV)<-c("chrom","start","end","value")
scATACGEXCOVlist <- split(scATACGEXCOV, scATACGEXCOV$chrom)

scATACGEXCOVrpkmfile = paste0("mouseBrain.clean.GEX.",args[2],".RPKM.bedgraph") # mouseBrain.clean.GEX.P11_Astroependymal.RPKM.bedgraph
scATACGEXCOVrpkm = read.table(readfile(scATACGEXCOVrpkmfile),header=F,stringsAsFactors=F)
names(scATACGEXCOVrpkm)<-c("chrom","start","end","value")
scATACGEXCOVrpkmlist <- split(scATACGEXCOVrpkm, scATACGEXCOVrpkm$chrom)

scATACGEXdupCOVfile = paste0("mouseBrain.clean.dup.GEX.",args[2],".bedgraph") # mouseBrain.clean.dup.GEX.P11_Astroependymal.bedgraph
scATACGEXdupCOV = read.table(readfile(scATACGEXdupCOVfile),header=F,stringsAsFactors=F)
names(scATACGEXdupCOV)<-c("chrom","start","end","value")
scATACGEXdupCOVlist <- split(scATACGEXdupCOV, scATACGEXdupCOV$chrom)

scATACGEXdupCOVrpkmfile = paste0("mouseBrain.clean.dup.GEX.",args[2],".RPKM.bedgraph") # mouseBrain.clean.dup.GEX.P11_Astroependymal.RPKM.bedgraph
scATACGEXdupCOVrpkm = read.table(readfile(scATACGEXdupCOVrpkmfile),header=F,stringsAsFactors=F)
names(scATACGEXdupCOVrpkm)<-c("chrom","start","end","value")
scATACGEXdupCOVrpkmlist <- split(scATACGEXdupCOVrpkm, scATACGEXdupCOVrpkm$chrom)

scRNAfragmentfile = paste0("mouseBrain.clean.GEX.",args[2],".fragments.unsorted.tsv") # mouseBrain.clean.GEX.P11_Astroependymal.fragments.unsorted.tsv
scRNAfragment = read.table(readfile(scRNAfragmentfile),header=F,stringsAsFactors=F)
names(scRNAfragment)<-c("chrom","start","end","barcode","count")
scRNAfragmentlist <- split(scRNAfragment, scRNAfragment$chrom)

scRNAdupfragmentfile = paste0("mouseBrain.clean.dup.GEX.",args[2],".fragments.unsorted.tsv") # mouseBrain.clean.dup.GEX.P11_Astroependymal.fragments.unsorted.tsv
scRNAdupfragment = read.table(readfile(scRNAdupfragmentfile),header=F,stringsAsFactors=F)
names(scRNAdupfragment)<-c("chrom","start","end","barcode","count")
scRNAdupfragmentlist <- split(scRNAdupfragment, scRNAdupfragment$chrom)


## chromstat
chromstatfile = args[3] # /data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrain/mouseBrain.chromstat.hmm.bed
chromstat = read.table(readfile(chromstatfile),header=F,stringsAsFactors=F)
names(chromstat)<-c("chrom","start","end","name","score","strand","thickStart","thickEnd","itemRgb")
chromstatlist <- split(chromstat, chromstat$chrom)


## save
save(genelist,scATACfragmentlist,scATACpeaklist,scATACARCCOVlist,scATACARCCOVrpkmlist,scATACPETlist,scATAClooplist,scATAClooptwopeaklist,scATACloopgivenanchorlist,scATACGEXCOVlist,scATACGEXCOVrpkmlist,scATACGEXdupCOVlist,scATACGEXdupCOVrpkmlist,scRNAfragmentlist,scRNAdupfragmentlist,chromstatlist,file=paste0("mouseBrain.",args[2],".all.data.RData")) # mouseBrain.P11_Astroependymal.all.data.RData

chromref = c(paste0("chr",1:100),"chrX")
chromlistraw = unique(c(names(genelist), names(scATACfragmentlist), names(scATACpeaklist), names(scATACARCCOVlist), names(scATACARCCOVrpkmlist), names(scATACPETlist), names(scATAClooplist), names(scATAClooptwopeaklist), names(scATACloopgivenanchorlist), names(scATACGEXCOVlist), names(scATACGEXCOVrpkmlist), names(scATACGEXdupCOVlist), names(scATACGEXdupCOVrpkmlist), names(scRNAfragmentlist), names(scRNAdupfragmentlist), names(chromstatlist)))
chromlist = intersect(chromref, chromlistraw)

genelist2 = genelist
scATACfragmentlist2 = scATACfragmentlist
scATACpeaklist2 = scATACpeaklist
scATACARCCOVlist2 = scATACARCCOVlist
scATACARCCOVrpkmlist2 = scATACARCCOVrpkmlist
scATACPETlist2 = scATACPETlist
scATAClooplist2 = scATAClooplist
scATAClooptwopeaklist2 = scATAClooptwopeaklist
scATACloopgivenanchorlist2 = scATACloopgivenanchorlist
scATACGEXCOVlist2 = scATACGEXCOVlist
scATACGEXCOVrpkmlist2 = scATACGEXCOVrpkmlist
scATACGEXdupCOVlist2 = scATACGEXdupCOVlist
scATACGEXdupCOVrpkmlist2 = scATACGEXdupCOVrpkmlist
scRNAfragmentlist2 = scRNAfragmentlist
scRNAdupfragmentlist2 = scRNAdupfragmentlist
chromstatlist2 = chromstatlist

for (k in chromlist){
  genelist = list()
  genelist[[k]] = genelist2[[k]]
  scATACfragmentlist = list()
  scATACfragmentlist[[k]] = scATACfragmentlist2[[k]]
  scATACpeaklist = list()
  scATACpeaklist[[k]] = scATACpeaklist2[[k]]
  scATACARCCOVlist = list()
  scATACARCCOVlist[[k]] = scATACARCCOVlist2[[k]]
  scATACARCCOVrpkmlist = list()
  scATACARCCOVrpkmlist[[k]] = scATACARCCOVrpkmlist2[[k]]
  scATACPETlist = list()
  scATACPETlist[[k]] = scATACPETlist2[[k]]
  scATAClooplist = list()
  scATAClooplist[[k]] = scATAClooplist2[[k]]
  scATAClooptwopeaklist = list()
  scATAClooptwopeaklist[[k]] = scATAClooptwopeaklist2[[k]]
  scATACloopgivenanchorlist = list()
  scATACloopgivenanchorlist[[k]] = scATACloopgivenanchorlist2[[k]]
  scATACGEXCOVlist = list()
  scATACGEXCOVlist[[k]] = scATACGEXCOVlist2[[k]]
  scATACGEXCOVrpkmlist = list()
  scATACGEXCOVrpkmlist[[k]] = scATACGEXCOVrpkmlist2[[k]]
  scATACGEXdupCOVlist = list()
  scATACGEXdupCOVlist[[k]] = scATACGEXdupCOVlist2[[k]]
  scATACGEXdupCOVrpkmlist = list()
  scATACGEXdupCOVrpkmlist[[k]] = scATACGEXdupCOVrpkmlist2[[k]]
  scRNAfragmentlist = list()
  scRNAfragmentlist[[k]] = scRNAfragmentlist2[[k]]
  scRNAdupfragmentlist = list()
  scRNAdupfragmentlist[[k]] = scRNAdupfragmentlist2[[k]]
  chromstatlist = list()
  chromstatlist[[k]] = chromstatlist2[[k]]

  save(genelist,scATACfragmentlist,scATACpeaklist,scATACARCCOVlist,scATACARCCOVrpkmlist,scATACPETlist,scATAClooplist,scATAClooptwopeaklist,scATACloopgivenanchorlist,scATACGEXCOVlist,scATACGEXCOVrpkmlist,scATACGEXdupCOVlist,scATACGEXdupCOVrpkmlist,scRNAfragmentlist,scRNAdupfragmentlist,chromstatlist, file=paste0("alldata/mouseBrain.",args[2],".",k,".data.RData"))
}

