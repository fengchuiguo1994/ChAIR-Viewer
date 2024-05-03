library(RColorBrewer)
library(pheatmap)
library(plotrix)
library(dplyr)
library(GenomicRanges)
# library(data.table)

run <- function(region,file,cut_region1="",cut_region2="",cut= F,gene="",direction = T,showAllPET=F,cell_number=500,type = "PET",sameCell = F,PETcount = 2,VCscore=0.001,condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,readOtherData="",outfmt="png",output="",showBarcode=T,showGex=T,showArc=T,showLoop=T,adjustBin=F){

  region = gsub("-|:"," ",region)
  regionarr = strsplit(trimws(region)," |\t")
  chrom = regionarr[[1]][1]
  start = as.numeric(regionarr[[1]][2])
  end = as.numeric(regionarr[[1]][3])
  parlist = c(0,0.15,0.3,0.7325,0.8125,0.875,0.9375,1) # RNAf ATACf PETf gene RNA ATAC/peak loop # kesheng
  twosup=T # kesheng
  givenAnchor=T
  dis = 8000 # kesheng

  #### read cell cycle data
  aa_list <- list()
  for(i in file){
      if(clean==T){
      aa_list[[i]] = filterData(data_list[[i]][["genelist"]],data_list[[i]][["scATACfragmentlist"]],data_list[[i]][["scATACpeaklist"]],data_list[[i]][["scATACARCCOVrpkmlist"]],data_list[[i]][["scATACPETlist"]],data_list[[i]][["scATAClooplist"]],data_list[[i]][["scATAClooptwopeaklist"]],data_list[[i]][["scATACloopgivenanchorlist"]],data_list[[i]][["scATACGEXdupCOVrpkmlist"]],data_list[[i]][["scRNAfragmentlist"]],data_list[[i]][["chromstatlist"]],chrom,start,end,dis=dis,PETcount = PETcount,VCscore=VCscore,condition=condition,sample_num = sample_num,cleanPET=clean,cut_region1=cut_region1,cut_region2=cut_region2,cut=cut)
      }else{
      aa_list[[i]] = filterData(data_list[[i]][["genelist"]],data_list[[i]][["scATACfragmentlist"]],data_list[[i]][["scATACpeaklist"]],data_list[[i]][["scATACARCCOVrpkmlist"]],data_list[[i]][["scATACPETlist"]],data_list[[i]][["scATAClooplist"]],data_list[[i]][["scATAClooptwopeaklist"]],data_list[[i]][["scATACloopgivenanchorlist"]],data_list[[i]][["scATACGEXCOVrpkmlist"]],data_list[[i]][["scRNAdupfragmentlist"]],data_list[[i]][["chromstatlist"]],chrom,start,end,dis=dis,PETcount = PETcount,VCscore=VCscore,condition=condition,sample_num = sample_num,cleanPET=clean,cut_region1=cut_region1,cut_region2=cut_region2,cut=cut)        
      }

      if(adjustBin && (end-start >= 20000)){
        aa_list[[i]][["scATACARCCOVflt"]] <- adjustBin(aa_list[[i]][["scATACARCCOVflt"]],chrom,start,end)
        aa_list[[i]][["scATACGEXCOVflt"]] <- adjustBin(aa_list[[i]][["scATACGEXCOVflt"]],chrom,start,end)
      }

  }
  #######################
  if(direction && gene != ""){
        geneinfo <- getGene(aa_list[[1]][["transcriptinfo"]],chrom=chrom,gene)
        # start end strand
  }
  total_pic = length(aa_list)
  yaxis_list <- yaxis_cor(aa_list,twosup=twosup,givenAnchor=givenAnchor,showPETnum=showPETnum)
  cell_number_vec <- cell_number_cal(aa_list,cell_number=cell_number,type=type,direction=direction,geneinfo=geneinfo)

  if(output=="" && outfmt=="pdf"){
  # output=paste("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/",chrom,start,end,file,"gene",gene,"direction",direction,"type",type,"PETcount",PETcount,"VCscore",VCscore,"sample_num",sample_num,"clean",clean,"showPETnum",showPETnum,".pdf",sep="-")
    if(cut){
      output=paste("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/",chrom,start,end,file,"cut",cut,"cut_region",cut_region1[1],"direction",direction,"type",type,"clean",clean,"sample_num",sample_num,".pdf",sep="-")
    }else{
      output=paste("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/",chrom,start,end,file,"gene",gene,"direction",direction,"type",type,"clean",clean,"sample_num",sample_num,".pdf",sep="-")
    }
  } else if(output=="" && outfmt=="png"){
    if(cut){
      output=paste("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/",chrom,start,end,file,"cut",cut,"cut_region",cut_region1[1],"direction",direction,"type",type,"PETcount",PETcount,"VCscore",VCscore,"sample_num",sample_num,".png",sep="-")
    }else{
      output=paste("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/",chrom,start,end,file,"gene",gene,"direction",direction,"type",type,"PETcount",PETcount,"VCscore",VCscore,"sample_num",sample_num,".png",sep="-")     
    }

  }

  if(outfmt=="png"){
    png(output,height=2*480*2,width=480*total_pic*2,res=72*3)
  } else if(outfmt=="pdf"){
    pdf(output,height=14,width=7*total_pic)
  }

  for(pic in 1:total_pic){
      aa <- aa_list[[pic]]
      summaryData(aa)
      ### cal barcode order
      if(sameCell){
          bb_list_list <- getTopSameData(aa_list,cell_number=cell_number,type=type,direction=direction,geneinfo=geneinfo)
          bb_list <- bb_list_list[[1]]
          
          orderBarcode <- bb_list_list[[2]]
          #record_i <- bb_list_list[[3]] ### calucate the pic of having most cell number with PET
          #bb <- bb_list[[record_i]]
          if(direction){ 
              list_tmp <- sortedByDirection(bb_list[[1]][["scATACPETflt"]],geneinfo,showAllPET=showAllPET)
              orderPET_dataframe1 <- list_tmp[[1]]
              orderPETBarcode1 = orderPET_dataframe1$orderBarcode
              bb_list[[1]][["scATACPETflt"]] <- list_tmp[[2]]

              list_tmp <- sortedByDirection(bb_list[[2]][["scATACPETflt"]],geneinfo,showAllPET=showAllPET)
              orderPET_dataframe2 <- list_tmp[[1]]
              orderPETBarcode2 = orderPET_dataframe2$orderBarcode
              bb_list[[2]][["scATACPETflt"]] <- list_tmp[[2]]

              color_class = orderPET_dataframe1$color_class
              color_class = c(rep(0,length(setdiff(orderBarcode,orderPETBarcode1))),color_class)
              color_class = color_class + 1

              orderPETBarcode = c(setdiff(orderPETBarcode2,intersect(orderPETBarcode1,orderPETBarcode2)),intersect(orderPETBarcode1,orderPETBarcode2),setdiff(orderPETBarcode1,intersect(orderPETBarcode1,orderPETBarcode2)))
              #if(is.null(orderPETBarcode1) && pic==1) orderPETBarcode <- c()
              #if(is.null(orderPETBarcode2) && pic==2) orderPETBarcode <- c()
              #cell_number_vec <- length(orderPETBarcode)
              #if(is.null(orderPETBarcode1) && pic==1) cell_number_vec <- 0
              #if(is.null(orderPETBarcode2) && pic==2) cell_number_vec <- 0
              
          }else{
              orderPET_dataframe = clusterPET(bb[["scATACPETflt"]], start, end,nbin=50,nclass=20,direction=F)
          }
        
          bb <- bb_list[[pic]]
          #if(direction) bb[["scATACPETflt"]] <- list_tmp[[2]]
      }else{
          bb_list = getTopData(aa,cell_number=cell_number,type=type,direction=direction,geneinfo=geneinfo)
          bb <- bb_list[[1]]

          orderBarcode <- bb_list[[2]]

          if(dim(bb[["scATACPETflt"]])[1]==0){
                orderPET_dataframe <- data.frame(orderBarcode=c(),color_class=c())
          }else{    
            if(direction){
                  list_tmp <- sortedByDirection(bb[["scATACPETflt"]],geneinfo,showAllPET=showAllPET)
                  orderPET_dataframe <- list_tmp[[1]]
                  bb[["scATACPETflt"]] <- list_tmp[[2]]

                  write.table(orderPET_dataframe$orderBarcode,paste0("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/Dcc/excit_scATACPETflt_barcode_",pic,".table"),sep="\t",quote=F)
                  write.table(bb[["scATACPETflt"]],paste0("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/Dcc/excit_scATACPETflt_",pic,".table"),sep="\t",quote=F)      
            }else{
                orderPET_dataframe = clusterPET(bb[["scATACPETflt"]], start, end,nbin=50,nclass=20,direction=F)
            }        
          }
          orderPETBarcode = orderPET_dataframe$orderBarcode
          color_class = orderPET_dataframe$color_class
          color_class = c(rep(0,length(setdiff(orderBarcode,orderPETBarcode))),color_class)
          color_class = color_class + 1
          orderPETBarcode = c(setdiff(orderBarcode,orderPETBarcode),orderPETBarcode)
          orderATACBarcode = sortedByNumber(bb[["scATACfragmentflt"]])
          orderATACBarcode = c(setdiff(orderBarcode,orderATACBarcode),orderATACBarcode)
          orderRNABarcode = sortedByNumber(bb[["scRNAfragmentflt"]])
          orderRNABarcode = c(setdiff(orderBarcode,orderRNABarcode),orderRNABarcode)
      }

      plotData3(bb,color_class,orderPETBarcode,orderPETBarcode,orderPETBarcode,chrom,start,end,yaxis_list,cell_number_vec,autoscale=F,parlist,twosup=twosup,givenAnchor=givenAnchor,pic=pic,total_pic=total_pic,direction=direction,showPETnum=showPETnum,sameCell=sameCell,readOtherData=readOtherData,showBarcode=showBarcode,showGex=showGex,showArc=showArc,showLoop=showLoop) 

      if(sameCell){
        barcode_vec <- c(paste0("PlotSameCell"),orderPETBarcode)
      }else{
        if(!exists("barcode_vec")){
          barcode_vec <- c(paste0("Plot",pic),orderPETBarcode)
        }else{
          barcode_vec <- c(barcode_vec,paste0("Plot",pic),orderPETBarcode)
        }        
      }

      
      
  }

  if(outfmt=="pdf"){
    write.table(as.data.frame(barcode_vec),file=gsub("pdf","csv",output[1]),quote=F,col.names=F,row.names=F,sep=",")
  } 

  dev.off()
}

loadData <- function(file=c("patski.allG1.data.RData","patski.allS.data.RData","patski.allG2M.data.RData"),path="/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineCellCycleFinal/",addLoad="",readOtherData=""){
  data_list <- list()
  for(i in file){
      chromstatlist <- list()
      load(paste0(path,i))
      if(addLoad!=""){
        load(addLoad)
      }
      
      if(readOtherData == ""){
        data_list[[i]] <- list(genelist,scATACfragmentlist,scATACpeaklist,scATACARCCOVlist,scATACARCCOVrpkmlist,scATACPETlist,scATAClooplist,scATAClooptwopeaklist,scATACloopgivenanchorlist,scATACGEXCOVlist,scATACGEXCOVrpkmlist,scATACGEXdupCOVlist,scATACGEXdupCOVrpkmlist,scRNAfragmentlist,scRNAdupfragmentlist,chromstatlist)
      }else if(readOtherData == "bulkChIATAC"){
        data_list[[i]] <- list(genelist,scATACfragmentlist,scATACpeaklist,bulkChIATACCOVlist,bulkChIATACCOVrpkmlist,scATACPETlist,bulkChIATAClooplist,bulkChIATAClooptwopeaklist,bulkChIATACloopgivenanchorlist,scATACGEXCOVlist,scATACGEXCOVrpkmlist,scATACGEXdupCOVlist,scATACGEXdupCOVrpkmlist,scRNAfragmentlist,scRNAdupfragmentlist,chromstatlist)
      }else if(readOtherData == "bulkATAC"){
        data_list[[i]] <- list(genelist,scATACfragmentlist,scATACpeaklist,bulkATACCOVlist,bulkATACCOVrpkmlist,scATACPETlist,scATAClooplist,scATAClooptwopeaklist,scATACloopgivenanchorlist,scATACGEXCOVlist,scATACGEXCOVrpkmlist,scATACGEXdupCOVlist,scATACGEXdupCOVrpkmlist,scRNAfragmentlist,scRNAdupfragmentlist,chromstatlist)
      }else if(readOtherData == "bulkRNA"){
        data_list[[i]] <- list(genelist,scATACfragmentlist,scATACpeaklist,scATACARCCOVlist,scATACARCCOVrpkmlist,scATACPETlist,scATAClooplist,scATAClooptwopeaklist,scATACloopgivenanchorlist,bulkRNACOVlist,bulkRNACOVrpkmlist,scATACGEXdupCOVlist,scATACGEXdupCOVrpkmlist,scRNAfragmentlist,scRNAdupfragmentlist,chromstatlist)       
      }

      names(data_list[[i]]) <- c("genelist","scATACfragmentlist","scATACpeaklist","scATACARCCOVlist","scATACARCCOVrpkmlist","scATACPETlist","scATAClooplist","scATAClooptwopeaklist","scATACloopgivenanchorlist","scATACGEXCOVlist","scATACGEXCOVrpkmlist","scATACGEXdupCOVlist","scATACGEXdupCOVrpkmlist","scRNAfragmentlist","scRNAdupfragmentlist","chromstatlist")
  }
  return(data_list)
}

summaryData = function(mylist) {
  # geneflt,transcriptinfo,scATACfragmentflt,scATACpeakflt,scATACARCCOVflt,scATACPETflt,scATACloopflt,scATAClooptwopeakflt,scATACGEXCOVflt,scRNAfragmentflt
  cat(sprintf("There are %d RNA reads, and there are %d cells\n", nrow(mylist[["scRNAfragmentflt"]]),length(unique(mylist[["scRNAfragmentflt"]]$barcode))))
  cat(sprintf("There are %d ATAC fragments, and there are %d cells\n", nrow(mylist[["scATACfragmentflt"]]),length(unique(mylist[["scATACfragmentflt"]]$barcode))))
  cat(sprintf("There are %d interactions, and there are %d cells\n", nrow(mylist[["scATACPETflt"]]),length(unique(mylist[["scATACPETflt"]]$barcode))))
  tmptmp = unique(mylist[["geneflt"]]$gene)
  cat(sprintf("There are %d genes, these genes are: %s\n", length(tmptmp),paste(tmptmp,collapse = ";")))
  cat(sprintf("There are %d peaks\n", nrow(mylist[["scATACpeakflt"]])))
  cat(sprintf("There are %d loops\n", nrow(mylist[["scATACloopflt"]])))
  cat(sprintf("There are %d loops (two support)\n", nrow(mylist[["scATAClooptwopeakflt"]])))
  cat(sprintf("There are %d loops (given anchor)\n", nrow(mylist[["scATACloopgivenanchorflt"]])))
}

filterData = function(genelist,scATACfragmentlist,scATACpeaklist,scATACARCCOVlist,scATACPETlist,scATAClooplist,scATAClooptwopeaklist,scATACloopgivenanchorlist,scATACGEXCOVlist,scRNAfragmentlist,chromstatlist,chrom,start,end, dis = 8000, PETcount = 2,VCscore=0.001,condition="PETDNARNA",sample_num = -1,cleanPET=F,cut_region1="",cut_region2="",cut=F) {

  scATACfragment = scATACfragmentlist[[chrom]]
  scATACfragmentflt = scATACfragment[scATACfragment$chrom==chrom & scATACfragment$end > start & scATACfragment$start < end, ]

  chromstat = chromstatlist[[chrom]]
  chromstatflt = chromstat[chromstat$chrom==chrom & chromstat$end > start & chromstat$start < end, ]

  scATACPET = scATACPETlist[[chrom]]
  scATACPETflt = scATACPET[scATACPET$chrom1==chrom & scATACPET$end1 > start & scATACPET$start1 < end & scATACPET$chrom2==chrom & scATACPET$end2 > start & scATACPET$start2 < end, ]
  scATACPETflt = scATACPETflt[scATACPET$end2-scATACPET$start1 >= dis, ]
  
  scRNAfragment = scRNAfragmentlist[[chrom]]
  scRNAfragmentflt = scRNAfragment[scRNAfragment$chrom==chrom & scRNAfragment$end > start & scRNAfragment$start < end, ]
  
  if(sample_num>0){
  set.seed(1234)
  barcode <- sample(unique(data.frame(do.call(rbind,scATACfragmentlist))$barcode),size=sample_num)
  scRNAfragmentflt = scRNAfragmentflt[scRNAfragmentflt$barcode %in% barcode,]
  scATACfragmentflt = scATACfragmentflt[scATACfragmentflt$barcode %in% barcode,]
  scATACPETflt = scATACPETflt[scATACPETflt$barcode %in% barcode,]
  }

  #print("Barcode number")
  #print(length(unique(data.frame(do.call(rbind,scATACfragmentlist))$barcode)))
  
  if(condition == "PETDNARNA"){
    barcode = intersect(scATACfragmentflt$barcode, intersect(scATACPETflt$barcode, scRNAfragmentflt$barcode))
    scATACfragmentflt = scATACfragmentflt[scATACfragmentflt$barcode %in% barcode,]
    scATACPETflt = scATACPETflt[scATACPETflt$barcode %in% barcode,]
    scRNAfragmentflt = scRNAfragmentflt[scRNAfragmentflt$barcode %in% barcode,]
  }else if(condition == "DNARNA"){
    barcode = intersect(scATACfragmentflt$barcode, scRNAfragmentflt$barcode)
    scATACfragmentflt = scATACfragmentflt[scATACfragmentflt$barcode %in% barcode,]
    scATACPETflt = scATACPETflt[scATACPETflt$barcode %in% barcode,]
    scRNAfragmentflt = scRNAfragmentflt[scRNAfragmentflt$barcode %in% barcode,]
  }else if(condition == "RNA"){
    barcode = scRNAfragmentflt$barcode
    scATACfragmentflt = scATACfragmentflt[scATACfragmentflt$barcode %in% barcode,]
    scATACPETflt = scATACPETflt[scATACPETflt$barcode %in% barcode,]
    scRNAfragmentflt = scRNAfragmentflt[scRNAfragmentflt$barcode %in% barcode,]
  }

  gene = genelist[[chrom]]
  geneflt = gene[gene$end > start & gene$start < end & gene$type == "exon", ]
  genelist = unique(geneflt$gene)

  scATACpeak = scATACpeaklist[[chrom]]
  scATACpeakflt = scATACpeak[scATACpeak$chrom==chrom & scATACpeak$end > start & scATACpeak$start < end, ]

  scATACARCCOV = scATACARCCOVlist[[chrom]]
  scATACARCCOVflt = scATACARCCOV[scATACARCCOV$chrom==chrom & scATACARCCOV$end > start & scATACARCCOV$start < end, ]

  scATACloop = scATAClooplist[[chrom]]
  scATACloop = scATACloop[scATACloop$chrom1 == scATACloop$chrom2 & scATACloop$count >= PETcount,]
  scATACloopflt = scATACloop[scATACloop$chrom1==chrom & scATACloop$end1 > start & scATACloop$start1 < end & scATACloop$chrom2==chrom & scATACloop$end2 > start & scATACloop$start2 < end, ]

  scATAClooptwopeak = scATAClooptwopeaklist[[chrom]]
  scATAClooptwopeak = scATAClooptwopeak[scATAClooptwopeak$chrom1 == scATAClooptwopeak$chrom2 & scATAClooptwopeak$count >= PETcount,]  
  scATAClooptwopeakflt = scATAClooptwopeak[scATAClooptwopeak$chrom1==chrom & scATAClooptwopeak$end1 > start & scATAClooptwopeak$start1 < end & scATAClooptwopeak$chrom2==chrom & scATAClooptwopeak$end2 > start & scATAClooptwopeak$start2 < end, ]

  scATACloopgivenanchor = scATACloopgivenanchorlist[[chrom]]
  ## tmpadd
  # print(head(scATACloopgivenanchor))
  # names(scATACloopgivenanchor) <- c("chrom1","start1","end1","chrom2","start2","end2","count","cov1","cov2","score")

  scATACloopgivenanchor = scATACloopgivenanchor[scATACloopgivenanchor$chrom1 == scATACloopgivenanchor$chrom2 & scATACloopgivenanchor$count >= PETcount & scATACloopgivenanchor$score >= VCscore,]  
  scATACloopgivenanchorflt = scATACloopgivenanchor[scATACloopgivenanchor$chrom1==chrom & scATACloopgivenanchor$end1 > start & scATACloopgivenanchor$start1 < end & scATACloopgivenanchor$chrom2==chrom & scATACloopgivenanchor$end2 > start & scATACloopgivenanchor$start2 < end, ]

  scATACGEXCOV = scATACGEXCOVlist[[chrom]]
  scATACGEXCOVflt = scATACGEXCOV[scATACGEXCOV$chrom==chrom & scATACGEXCOV$end > start & scATACGEXCOV$start < end, ]

  tmpdata = data.frame(chrom=c(),start=c(),end=c())
  numberofgeneinfo = length(names(table(geneflt$ID)))
  namesgeneinfo = names(table(geneflt$ID))
  starts = c()
  stops  = c()
  sizes  = c()
  strands = c()
  for (i in (1:numberofgeneinfo)) {
    subgeneinfo  = geneflt[which(geneflt$ID == namesgeneinfo[i]),]
    starts = c(starts,min(subgeneinfo$start,subgeneinfo$end))
    stops  = c(stops, max(subgeneinfo$start,subgeneinfo$end))
    sizes  = c(sizes, stops[i] - starts[i])
    strands = c(strands,subgeneinfo$strand[1])
  }
  transcriptinfo = data.frame(names=namesgeneinfo,starts=starts,stops=stops,sizes=sizes,strand=strands)
  transcriptinfo = transcriptinfo[order(transcriptinfo$names),]
  typeshow = "squished"
  if (typeshow == "collapsed") {
    transcriptinfo$plotrow = 1
  } else if(typeshow == "squished") {
    if (nrow(transcriptinfo) > 0) {
      maxn = 1
      transcriptinfo$plotrow = 0
      for (i in (1:nrow(transcriptinfo))){
        overlapflag = T
        for (j in (1:maxn)) {
          transcriptinfotmp = transcriptinfo[transcriptinfo$plotrow==j &
                                               transcriptinfo$starts[i] < transcriptinfo$stops & transcriptinfo$stops[i] > transcriptinfo$starts,]
          if (nrow(transcriptinfotmp) == 0) {
            transcriptinfo$plotrow[i] = j
            overlapflag = F
            break
          }
        }
        if (overlapflag) {
          maxn = maxn + 1
          transcriptinfo$plotrow[i] = maxn
        }
      }
    }
  } else if (typeshow == "flatten") {
    transcriptinfo$plotrow = seq(1:nrow(transcriptinfo))
  }
  # genelist,scATACfragmentlist,scATACpeaklist,scATACARCCOVlist,scATACPETlist,scATAClooplist,scATAClooptwopeaklist,scATACGEXCOVlist,scRNAfragmentlist,chrom,start,end
  
  #print(scATACpeakflt)

  if(cleanPET){
    tmp_transcriptinfo <- transcriptinfo
    tmp_transcriptinfo = tmp_transcriptinfo[order(tmp_transcriptinfo$names),]
    for(i in 1:nrow(tmp_transcriptinfo)){
        tmp_transcriptinfo$names_new[i] <- strsplit(tmp_transcriptinfo[i,1],"-")[[1]][1]
    }
    tmp_transcriptinfo <- tmp_transcriptinfo[order(tmp_transcriptinfo$names_new,tmp_transcriptinfo$sizes),]
    tmp_gene = ""
    for(i in nrow(tmp_transcriptinfo):1){
      if(tmp_gene == tmp_transcriptinfo$names_new[i]){
        tmp_transcriptinfo <- tmp_transcriptinfo[-i,]
      }else{
        tmp_gene <- tmp_transcriptinfo$names_new[i]
      }
    }
    print(tmp_transcriptinfo$names_new)
    tmp <- c()
    for(i in 1:nrow(tmp_transcriptinfo)){
      if(tmp_transcriptinfo$strand[i] == "+"){
        gene_pos <- tmp_transcriptinfo$starts[i]
      }else{
        gene_pos <- tmp_transcriptinfo$stops[i]
      }     
      gene_start = gene_pos - 5000
      gene_end = gene_pos + 5000
      add <- scATACPETflt[(scATACPETflt$start1 > gene_start & scATACPETflt$end1 < gene_end) | ( scATACPETflt$start2 > gene_start & scATACPETflt$end2 < gene_end),]
      tmp <- rbind(tmp,add)
    }
  scATACPETflt <- unique(tmp)
  }

  tmp_scATACPETflt <- data.frame()
  if(cut){
    cut_region1 = gsub("-|:"," ",cut_region1)
    cut_region1arr = strsplit(trimws(cut_region1)," |\t")
    cut_chrom1 = cut_region1arr[[1]][1]
    cut_start1 = as.numeric(cut_region1arr[[1]][2])
    cut_end1 = as.numeric(cut_region1arr[[1]][3])

    for(k in 1:length(cut_region2)){

    cut_region = gsub("-|:"," ",cut_region2[k])
    cut_regionarr = strsplit(trimws(cut_region)," |\t")
    cut_chrom = cut_regionarr[[1]][1]
    cut_start = as.numeric(cut_regionarr[[1]][2])
    cut_end = as.numeric(cut_regionarr[[1]][3])

    tmp_scATACPETflt <- rbind(tmp_scATACPETflt,scATACPETflt[(scATACPETflt$start1 > cut_start1 & scATACPETflt$end1 < cut_end1 & scATACPETflt$start2 > cut_start & scATACPETflt$end2 < cut_end) | (scATACPETflt$start1 > cut_start & scATACPETflt$end1 < cut_end & scATACPETflt$start2 > cut_start1 & scATACPETflt$end2 < cut_end1),])
    }    

    scATACPETflt <- tmp_scATACPETflt
  }
  list <- list(geneflt,transcriptinfo,scATACfragmentflt,scATACpeakflt,scATACARCCOVflt,scATACPETflt,scATACloopflt,scATAClooptwopeakflt,scATACloopgivenanchorflt,scATACGEXCOVflt,scRNAfragmentflt,chromstatflt)
  names(list) <- c("geneflt","transcriptinfo","scATACfragmentflt","scATACpeakflt","scATACARCCOVflt","scATACPETflt","scATACloopflt","scATAClooptwopeakflt","scATACloopgivenanchorflt","scATACGEXCOVflt","scRNAfragmentflt","chromstatflt")
  return (list)
}

getGene <- function(transcriptinfo,chrom,gene){

  for(i in 1:nrow(transcriptinfo)){
    transcriptinfo$names_new[i] <- strsplit(transcriptinfo[i,1],"-")[[1]][1]
  }
  transcriptinfo_gene <- transcriptinfo[transcriptinfo$names_new %in% gene,]
  transcriptinfo_gene <- transcriptinfo_gene[order(-transcriptinfo_gene$sizes),]
  if(length(transcriptinfo_gene$strand)==0){
    stop(sprintf("No number in %s \n",gene))
  }
  
  return(list(transcriptinfo_gene$starts[1],transcriptinfo_gene$stops[1],transcriptinfo_gene$strand[1]))
  
}

sortedByDirection <- function(scATACPETflt,geneinfo,showAllPET=F){
  genestart = geneinfo[[1]]
  geneend = geneinfo[[2]]
  genestrand = geneinfo[[3]]
  if(genestrand == "-"){
    cut_pos = geneend
  }else{
    cut_pos = genestart
  }
  cut_start = cut_pos - 5000
  cut_end = cut_pos + 5000
 
  if(nrow(scATACPETflt)!=0){

    for(i in 1:nrow(scATACPETflt)){

      if((scATACPETflt$start1[i] > cut_start && scATACPETflt$end1[i] < cut_end) || (scATACPETflt$start2[i] > cut_start && scATACPETflt$end2[i] < cut_end)){
        dis1 <- (cut_start+cut_end) - (scATACPETflt$start1[i] + scATACPETflt$end1[i])
        dis2 <- (scATACPETflt$start2[i] + scATACPETflt$end2[i]) - (cut_start+cut_end)
        if(dis1 > dis2){
          scATACPETflt$cor[i] <- scATACPETflt$start1[i] + scATACPETflt$end1[i]
          scATACPETflt$anchor1[i] <- 2
          scATACPETflt$anchor2[i] <- 1
        }else {
          scATACPETflt$cor[i] <- scATACPETflt$start2[i] + scATACPETflt$end2[i]
          scATACPETflt$anchor1[i] <- 1
          scATACPETflt$anchor2[i] <- 2
        }
        if(scATACPETflt$start1[i] > cut_start && scATACPETflt$end1[i] < cut_end){
          scATACPETflt$anchor1[i] <- 1
        }
        if(scATACPETflt$start2[i] > cut_start && scATACPETflt$end2[i] < cut_end){
          scATACPETflt$anchor2[i] <- 1
        }    
      }else{
         scATACPETflt$cor[i] <- 0
         scATACPETflt$anchor1[i] <- 3
         scATACPETflt$anchor2[i] <- 3
      }
    }
    
    tmp <- scATACPETflt
    scATACPETflt <- scATACPETflt[(scATACPETflt$start1 > cut_start & scATACPETflt$end1 < cut_end)|(scATACPETflt$start2 > cut_start & scATACPETflt$end2 < cut_end),]

    scATACPETflt <- scATACPETflt[order(scATACPETflt$cor),]

    tmp1 <- scATACPETflt[(scATACPETflt$cor) < (cut_start+cut_end),][order(scATACPETflt$cor),]$barcode
    tmp1 <- unique(tmp1)
    tmp2 <- scATACPETflt[(scATACPETflt$cor) > (cut_start+cut_end),][order(-scATACPETflt$cor),]$barcode
    tmp2 <- unique(tmp2)
    tmp2 <- rev(tmp2)

    final_barcode <- c(na.omit(tmp1),na.omit(tmp2))
    out_order <- data.frame(orderBarcode=unique(final_barcode),color_class=rep(6,length(unique(final_barcode))))  

    if(showAllPET){
      return(list(out_order,tmp))
    }else{
      return(list(out_order,scATACPETflt))
    }
  }else{
    return(list(out_order<-data.frame(),scATACPETflt))
  }
}

cell_number_cal <- function(aa_list,cell_number,type,direction=F,geneinfo=""){

  if(type == "All") {
    print("No need to calculate cell_number")
    return(0)
  }

  cell_number_vec <- c()
  for(i in 1:length(aa_list)){
    aa <- aa_list[[i]]
    ### divide list into dataframe
    scATACfragmentflt = aa[["scATACfragmentflt"]] # DNA
    scATACPETflt = aa[["scATACPETflt"]] # PET
    scRNAfragmentflt = aa[["scRNAfragmentflt"]] # RNA

    scATACfragmentflt$barcode <- as.character(scATACfragmentflt$barcode)
    scATACPETflt$barcode <- as.character(scATACPETflt$barcode)
    scRNAfragmentflt$barcode <- as.character(scRNAfragmentflt$barcode)
    if(type=="PET"){
      obj = scATACPETflt
    } else if(type=="DNA"){
      obj = scATACfragmentflt
    } else if(type=="RNA"){
      obj = scRNAfragmentflt
    }

    if(direction){
      genestart = geneinfo[[1]]
      geneend = geneinfo[[2]]
      genestrand = geneinfo[[3]]
      if(genestrand == "-"){
        cut_pos = geneend
      }else{
        cut_pos = genestart
      }
      cut_start = cut_pos - 5000
      cut_end = cut_pos + 5000

      obj_PET <- scATACPETflt
      obj_DNA <- scATACfragmentflt
      obj_RNA <- scRNAfragmentflt

      PET_barcode <- obj_PET[((obj_PET$start1 > cut_start & obj_PET$end1 < cut_end)|(obj_PET$start2 > cut_start & obj_PET$end2 < cut_end)) ,]$barcode
      DNA_barcode <- obj_DNA[(obj_DNA$start > cut_start & obj_DNA$end < cut_end),]$barcode
      RNA_barcode <- obj_RNA[(obj_RNA$start > genestart & obj_RNA$end < geneend),]$barcode

    if(type=="PETDNARNA"){
      obj = scATACPETflt
      obj <- obj[obj$barcode %in% PET_barcode,]
      obj <- obj[obj$barcode %in% DNA_barcode,]
      obj <- obj[obj$barcode %in% RNA_barcode,]
      }else if(type=="PET"){
      obj = scATACPETflt
      obj <- obj[obj$barcode %in% PET_barcode,]
    }

    }
    cell_number_vec[i] <- 0
    if(length(table(obj$barcode)) <= cell_number){
      cell_number_vec[i] <- length(table(obj$barcode))
    }else{
      cell_number_vec[i] <- cell_number
    }
    

    cat(sprintf("List of %s \n", names(aa_list)[i]))
    print(cell_number_vec[i])
  }
  #if(cell_number > min(cell_number_vec)) cell_number <- min(cell_number_vec)
  ### this code just to show cell number
  return(cell_number_vec)
}

getTopData <- function(aa,cell_number=500,type="RNA",direction=F,geneinfo=c("")){
  ### divide list into dataframe
  geneflt = aa[["geneflt"]]
  transcriptinfo = aa[["transcriptinfo"]]
  scATACfragmentflt = aa[["scATACfragmentflt"]] # DNA
  scATACpeakflt = aa[["scATACpeakflt"]]
  scATACARCCOVflt = aa[["scATACARCCOVflt"]]
  scATACPETflt = aa[["scATACPETflt"]] # PET
  scATACloopflt = aa[["scATACloopflt"]]
  scATAClooptwopeakflt = aa[["scATAClooptwopeakflt"]]
  scATACloopgivenanchorflt = aa[["scATACloopgivenanchorflt"]]
  scATACGEXCOVflt = aa[["scATACGEXCOVflt"]]
  scRNAfragmentflt = aa[["scRNAfragmentflt"]] # RNA
  chromstatflt = aa[["chromstatflt"]]

  scATACfragmentflt$barcode <- as.character(scATACfragmentflt$barcode)
  scATACPETflt$barcode <- as.character(scATACPETflt$barcode)
  scRNAfragmentflt$barcode <- as.character(scRNAfragmentflt$barcode)

  if(type=="PET"){
      obj = scATACPETflt
    } else if(type=="DNA"){
      obj = scATACfragmentflt
    } else if(type=="RNA"){
      obj = scRNAfragmentflt
  } 

  if(direction){
    genestart = geneinfo[[1]]
    geneend = geneinfo[[2]]
    genestrand = geneinfo[[3]]
    if(genestrand == "-"){
      cut_pos = geneend
    }else{
      cut_pos = genestart
    }
    cut_start = cut_pos - 5000
    cut_end = cut_pos + 5000

    obj_PET <- scATACPETflt
    obj_DNA <- scATACfragmentflt
    obj_RNA <- scRNAfragmentflt

    PET_barcode <- obj_PET[((obj_PET$start1 > cut_start & obj_PET$end1 < cut_end)|(obj_PET$start2 > cut_start & obj_PET$end2 < cut_end)) ,]$barcode
    DNA_barcode <- obj_DNA[(obj_DNA$start > cut_start & obj_DNA$end < cut_end),]$barcode
    RNA_barcode <- obj_RNA[(obj_RNA$start > genestart & obj_RNA$end < geneend),]$barcode

    if(type=="PETDNARNA"){
      obj = scATACPETflt
      obj <- obj[obj$barcode %in% PET_barcode,]
      obj <- obj[obj$barcode %in% DNA_barcode,]
      obj <- obj[obj$barcode %in% RNA_barcode,]
      }else if(type=="PET"){
      obj = scATACPETflt
      obj <- obj[obj$barcode %in% PET_barcode,]
    } 

    # scATACPETflt <- scATACPETflt[(scATACPETflt$start1 > cut_start & scATACPETflt$end1 < cut_end)|(scATACPETflt$start2 > cut_start & scATACPETflt$end2 < cut_end),]
  }

  order <- data.frame(barcode=names(table(obj$barcode)),barcodecount=as.vector(table(obj$barcode)))
  if(!is.null(names(table(obj$barcode)))) {
    order <- order[order(-order$barcodecount),]
  }

  if(cell_number == -1){
    orderBarcode <- order$barcode
  }else if(nrow(order) <= cell_number){
    orderBarcode <- order$barcode   
  }
  else{
    orderBarcode <- order[1:cell_number,]$barcode
  }

  print(length(orderBarcode))

  bb <- list(geneflt,transcriptinfo,scATACfragmentflt[scATACfragmentflt$barcode %in% orderBarcode,],scATACpeakflt,scATACARCCOVflt,scATACPETflt[scATACPETflt$barcode %in% orderBarcode,],scATACloopflt,scATAClooptwopeakflt,scATACloopgivenanchorflt,scATACGEXCOVflt,scRNAfragmentflt[scRNAfragmentflt$barcode %in% orderBarcode,],chromstatflt)
  names(bb) <- c("geneflt","transcriptinfo","scATACfragmentflt","scATACpeakflt","scATACARCCOVflt","scATACPETflt","scATACloopflt","scATAClooptwopeakflt","scATACloopgivenanchorflt","scATACGEXCOVflt","scRNAfragmentflt","chromstatflt")

  return(list(bb,orderBarcode))
}

getTopSameData <- function(aa_list,cell_number=500,type="RNA",direction=F,geneinfo=geneinfo){
  ### if type="All", it means we will get all barcode and ignore cell_number
  ### input

  tmp_barcode <- c()
  for(i in 1:length(aa_list)){
    aa <- aa_list[[i]]
    ### divide list into dataframe
    geneflt = aa[["geneflt"]]
    transcriptinfo = aa[["transcriptinfo"]]
    scATACfragmentflt = aa[["scATACfragmentflt"]] # DNA
    scATACpeakflt = aa[["scATACpeakflt"]]
    scATACARCCOVflt = aa[["scATACARCCOVflt"]]
    scATACPETflt = aa[["scATACPETflt"]] # PET
    scATACloopflt = aa[["scATACloopflt"]]
    scATAClooptwopeakflt = aa[["scATAClooptwopeakflt"]]
    scATACloopgivenanchorflt = aa[["scATACloopgivenanchorflt"]]
    scATACGEXCOVflt = aa[["scATACGEXCOVflt"]]
    scRNAfragmentflt = aa[["scRNAfragmentflt"]] # RNA
    chromstatflt = aa[["chromstatflt"]]


    scATACfragmentflt$barcode <- as.character(scATACfragmentflt$barcode)
    scATACPETflt$barcode <- as.character(scATACPETflt$barcode)
    scRNAfragmentflt$barcode <- as.character(scRNAfragmentflt$barcode)
     
    if(direction){
      genestart = geneinfo[[1]]
      geneend = geneinfo[[2]]
      genestrand = geneinfo[[3]]
      if(genestrand == "-"){
        cut_pos = geneend
      }else{
        cut_pos = genestart
      }
      cut_start = cut_pos - 5000
      cut_end = cut_pos + 5000

      obj_PET <- scATACPETflt
      obj_DNA <- scATACfragmentflt
      obj_RNA <- scRNAfragmentflt

      PET_barcode <- obj_PET[((obj_PET$start1 > cut_start & obj_PET$end1 < cut_end)|(obj_PET$start2 > cut_start & obj_PET$end2 < cut_end)) ,]$barcode
      DNA_barcode <- obj_DNA[(obj_DNA$start > cut_start & obj_DNA$end < cut_end),]$barcode
      RNA_barcode <- obj_RNA[(obj_RNA$start > genestart & obj_RNA$end < geneend),]$barcode

    if(type=="PETDNARNA"){
      obj = scATACPETflt
      obj <- obj[obj$barcode %in% PET_barcode,]
      obj <- obj[obj$barcode %in% DNA_barcode,]
      obj <- obj[obj$barcode %in% RNA_barcode,]
    }else if(type=="PET"){
      obj = scATACPETflt
      obj <- obj[obj$barcode %in% PET_barcode,]
    }

      cat(sprintf("direction T type: There are %d intersect barcode, %d former, %d latter\n", length(intersect(tmp_barcode,obj$barcode)),length(unique(tmp_barcode)),length(unique(obj$barcode))))
      tmp_barcode <- c(tmp_barcode,obj$barcode)
      order <- data.frame(barcode=names(table(tmp_barcode)),barcodecount=as.vector(table(tmp_barcode)))
    }

    order <- order[order(-order$barcodecount),]

    if(cell_number == -1){
            orderBarcode <- order$barcode
    }else{
      if(!is.integer(order)){
        if(nrow(order) <= cell_number){
          orderBarcode <- order$barcode   
        }else{
            orderBarcode <- order[1:cell_number,]$barcode
        }
      }else{
        orderBarcode <- c()
      }

    }      
  }
  
  
  print(length(orderBarcode))
  
  ### output
  bb_list <- list()
  record_PET <- 0
  record_i <- 0
  for(i in 1:length(aa_list)){
    aa <- aa_list[[i]]
    ### divide list into dataframe
    geneflt = aa[["geneflt"]]
    transcriptinfo = aa[["transcriptinfo"]]
    scATACfragmentflt = aa[["scATACfragmentflt"]] # DNA
    scATACpeakflt = aa[["scATACpeakflt"]]
    scATACARCCOVflt = aa[["scATACARCCOVflt"]]
    scATACPETflt = aa[["scATACPETflt"]] # PET
    scATACloopflt = aa[["scATACloopflt"]]
    scATAClooptwopeakflt = aa[["scATAClooptwopeakflt"]]
    scATACloopgivenanchorflt = aa[["scATACloopgivenanchorflt"]]
    scATACGEXCOVflt = aa[["scATACGEXCOVflt"]]
    scRNAfragmentflt = aa[["scRNAfragmentflt"]] # RNA
    chromstatflt = aa[["chromstatflt"]]

    scATACfragmentflt$barcode <- as.character(scATACfragmentflt$barcode)
    scATACPETflt$barcode <- as.character(scATACPETflt$barcode)
    scRNAfragmentflt$barcode <- as.character(scRNAfragmentflt$barcode)
    
    #  print(tail(scATACfragmentflt))
    #  print(tail(scATACfragmentflt[scATACPETflt$barcode %in% orderBarcode,]))
    bb_list[[i]] <- list(geneflt,transcriptinfo,scATACfragmentflt[scATACfragmentflt$barcode %in% orderBarcode,],scATACpeakflt,scATACARCCOVflt,scATACPETflt[scATACPETflt$barcode %in% orderBarcode,],scATACloopflt,scATAClooptwopeakflt,scATACloopgivenanchorflt,scATACGEXCOVflt,scRNAfragmentflt[scRNAfragmentflt$barcode %in% orderBarcode,],chromstatflt)
    names(bb_list[[i]]) <- c("geneflt","transcriptinfo","scATACfragmentflt","scATACpeakflt","scATACARCCOVflt","scATACPETflt","scATACloopflt","scATAClooptwopeakflt","scATACloopgivenanchorflt","scATACGEXCOVflt","scRNAfragmentflt","chromstatflt")

    if(record_PET < length(unique(scATACPETflt[scATACPETflt$barcode %in% orderBarcode,]$barcode))){
      record_PET <- length(unique(scATACPETflt[scATACPETflt$barcode %in% orderBarcode,]$barcode))
      record_i <- i
    }
  }
  return(list(bb_list,orderBarcode,record_i))
}

clusterPET = function(scATACPETflt, start, end, nbin=50, nclass = 10,direction = F){
  begin = start
  last = end
  bin = (last - begin)/nbin
  binid = paste("bin", seq(1,nbin),sep="")
  tmpdat <- data.frame(matrix(0, nrow = length(unique(scATACPETflt$barcode)), ncol = length(binid), dimnames = list(unique(scATACPETflt$barcode), binid)))

  if(length(unique(scATACPETflt$barcode))>=2){

        if(length(unique(scATACPETflt$barcode)) <= nclass ){
      nclass = length(unique(scATACPETflt$barcode))
    }

    if(direction){
      for (i in (1:nrow(scATACPETflt))){
        s = floor(((scATACPETflt$start1[i] + scATACPETflt$end1[i])/2-begin)/bin) + 1
        e = floor(((scATACPETflt$start2[i] + scATACPETflt$end2[i])/2-begin)/bin) + 1
        if(s < 1){
          s = 1
        } 
        if(e > nbin){
          e = nbin
        } 
        tmpdat[scATACPETflt$barcode[i],paste("bin",s:e,sep="")] = 1
      }
    }else{
      for (i in (1:nrow(scATACPETflt))){
      s = floor(((scATACPETflt$start1[i] + scATACPETflt$end1[i])/2-begin)/bin) + 1
      e = floor(((scATACPETflt$start2[i] + scATACPETflt$end2[i])/2-begin)/bin) + 1
      if(s < 1){
        s = 1
      } 
      if(e > nbin){
        e = nbin
      } 
      tmpdat[scATACPETflt$barcode[i],paste("bin",e,sep="")] = 1
      tmpdat[scATACPETflt$barcode[i],paste("bin",s,sep="")] = 1
      }
    }
    

    PET_order <- data.frame(barcode=names(table(scATACPETflt$barcode)),barcodecount=as.vector(table(scATACPETflt$barcode)))
    PET_order <- PET_order[order(-PET_order$barcodecount),]

    PET_order_barcode <- PET_order$barcode

    #print(length(PET_order_barcode))
    tmpdat <- tmpdat[PET_order_barcode, ]
    PH<-pheatmap(tmpdat, legend = FALSE, cluster_rows=TRUE, cluster_cols=FALSE, cutree_rows = nclass, 
                border_color = NA,
                treeheight_row = 50,
                treeheight_col = 0,
                annotation_legend = FALSE,
                annotation_names_row = FALSE,
                show_rownames=FALSE,
                show_colnames = FALSE,
                silent = TRUE)
    res_clu <- cutree(PH$tree_row, k = nclass)


    tmp = as.data.frame(res_clu)
    tmp$barcode = rownames(tmp)
    tmp = tmp[order(tmp$res_clu),]

    ### count each barcode length
    barcode_length_of_each_class <- data.frame(class = c(1:nclass),length = rep(0,nclass),num = rep(0,nclass))

    Total_PET = scATACPETflt

    ## mean
    # for(i in 1:(nrow(tmp))){
    #   barcode_length_of_each_class$length[tmp$res_clu[i]] = barcode_length_of_each_class$length[tmp$res_clu[i]] + (sum(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$end1)+sum(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$end2))/2 - (sum(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$start1)+sum(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$start2))/2
    #   barcode_length_of_each_class$num[tmp$res_clu[i]] = barcode_length_of_each_class$num[tmp$res_clu[i]] + 1
    # }
    # barcode_length_of_each_class$length_aver = barcode_length_of_each_class$length / barcode_length_of_each_class$num
    # barcode_length_of_each_class <- barcode_length_of_each_class[order(barcode_length_of_each_class$length_aver),]

    ## max
    for(i in 1:(nrow(tmp))){
      record = max(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$end2) - min(Total_PET[Total_PET$barcode %in% tmp$barcode[i],]$start1)
      if(barcode_length_of_each_class$length[tmp$res_clu[i]] < record) {
          barcode_length_of_each_class$length[tmp$res_clu[i]] = record
      }
    }
    barcode_length_of_each_class <- barcode_length_of_each_class[order(barcode_length_of_each_class$length),]

    ### smallest is the first
    barcode_length_of_each_class$length_order <- c(1:nclass)
    barcode_length_of_each_class <- barcode_length_of_each_class[order(barcode_length_of_each_class$class),]  
    for(i in 1:nrow(tmp)){
      tmp$length_order[i] <- barcode_length_of_each_class$length_order[tmp$res_clu[i]]
    }
    tmp <- tmp[order(tmp$length_order),]

    out_order <- data.frame(orderBarcode=tmp$barcode,color_class=tmp$length_order)
  }else{
    out_order <- data.frame(orderBarcode=scATACPETflt$barcode,color_class=c(1))
  }
  
  return(out_order)
}

sortedByNumber = function(scATACPETflt){
  ATAC_order <- data.frame(barcode=names(table(scATACPETflt$barcode)),barcodecount=as.vector(table(scATACPETflt$barcode)))
  if(nrow(ATAC_order)!=0){
    ATAC_order <- ATAC_order[order(ATAC_order$barcodecount),]
  }
  return(ATAC_order$barcode)
}

yaxis_cor <- function(mylist_list,twosup=T,givenAnchor=T,showPETnum=F){

  mylist_list_yaxis <- list()
  for(i in 1:length(mylist_list)){
    mylist <- mylist_list[[i]]
    gexbedgraphflt = mylist[["scATACGEXCOVflt"]]
    gex = max(gexbedgraphflt$value) * 1.05
  
    arcbedgraphflt = mylist[["scATACARCCOVflt"]]
    arc = max(arcbedgraphflt$value) * 1.05

    if (twosup) {
      if(givenAnchor){
        loopflt = mylist[["scATACloopgivenanchorflt"]]

        if(showPETnum){
          loopflt$count = log10(loopflt$count + 1)
          loop = max(loopflt$count) + log10(2)
        }else{
          if(!(is.null(loopflt) || nrow(loopflt) == 0)){
            loop = max(loopflt$score) * 1.05
          }else{
            loop=0
          }
          
        } 
      }else{
        loopflt = mylist[["scATAClooptwopeakflt"]]
        loopflt$count = log10(loopflt$count + 1)
        loop = max(loopflt$count) + log10(2)
      }
    } else {
    loopflt = mylist[["scATACloopflt"]]
    loopflt$count = log10(loopflt$count + 1)
    loop = max(loopflt$count) + log10(2)
    }

    # gex=2000
    mylist_list_yaxis[[i]] <- c(gex,arc,loop)
  }

  return(mylist_list_yaxis)
}

plotData3 = function(mylist,color_class,orderPETBarcode,orderATACBarcode,orderRNABarcode,chrom,start,end,yaxis_list,cell_number_vec,autoscale=T,parlist=c(0,0.15,0.3,0.7,0.775,0.85,0.925,1),twosup=T,givenAnchor=T,pic=1,total_pic=1,direction=F,showPETnum=F,sameCell=F,readOtherData="",showBarcode=T,showGex=T,showArc=T,showLoop=T) {
  # geneflt,transcriptinfo,scATACfragmentflt,scATACpeakflt,scATACARCCOVflt,scATACPETflt,scATACloopflt,scATAClooptwopeakflt,scATACGEXCOVflt,scRNAfragmentflt
  parlist <- (parlist + 0.04)/1.04
  
  if(readOtherData == ""){
    print("")
  }else if(readOtherData == "bulkChIATAC"){
    showBarcode = F
    showGex = F
    showArc = T
    showLoop = T
  }else if(readOtherData == "bulkATAC"){
    showBarcode = F
    showGex = F
    showArc = T
    showLoop = F
  }else if(readOtherData == "bulkRNA"){
    showBarcode = F
    showGex = T
    showArc = F
    showLoop = F   
  }

  if(!sameCell){
    adjust_pos <- max(cell_number_vec) - cell_number_vec[pic]
  }else{
    adjust_pos <-  0
  }

  # if(length(orderPETBarcode) >= 150){
  #   point_ratio = 2
  # }else if(length(orderPETBarcode) >=50){
  #   point_ratio <- 150 / length(orderPETBarcode) * 2 
  # }else{
  #   point_ratio = 6
  # }

  if(max(cell_number_vec) >= 150){
    point_ratio = 2
  }else if(max(cell_number_vec) >=50){
    point_ratio <- 150 / max(cell_number_vec) * 2 
  }else{
    point_ratio = 6
  }

  if(direction){
    color_DNA = "#2F7E4F"
  }else{
    color_DNA = "brown"
  }

  if(autoscale){
    yaxis <- yaxis_list[[pic]]
  }else{
    yaxis <- c(0,0,0)
    for(m in 1:3){
      for(i in 1:total_pic){
      if(yaxis[m] < yaxis_list[[i]][m]) yaxis[m] <- yaxis_list[[i]][m]
      }
    }
  }

  x_pic1 <- (pic-1)/total_pic
  x_pic2 <- pic/total_pic

  if(length(orderRNABarcode)==0){
    upvalue = 1
  }else{
    tmpdata = mylist[["scRNAfragmentflt"]][mylist[["scRNAfragmentflt"]]$barcode %in% orderRNABarcode,]

    upvalue = floor(quantile(tmpdata$count, 0.99))    

  }
  
  if(nrow(tmpdata)!=0){
      color = colorRampPalette(brewer.pal(9,"Blues"))(upvalue+1)
  }
  
  #   if (upvalue == 1){
  #     upvalue = 3
  #   }
  #   if (upvalue == 2){
  #     upvalue = 3
  #   }  
  #   if (upvalue > 9){
  #     color = colorRampPalette(brewer.pal(9,"Greens"))(upvalue)
  #   } else {
  #     color = brewer.pal(upvalue,"Greens")
  #   }
  
  

  par(cex=0.5, mai=c(0.1,0.2,0.1,0.1))
  if(pic==1){
    par(fig=c(x_pic1,x_pic2,parlist[1],parlist[2]))
  }else{
    par(fig=c(x_pic1,x_pic2,parlist[1],parlist[2]),new=T)
  }
  #plot(c(1,1),type='n',xlim=c(start,end), ylim=c(0,length(orderRNABarcode)+1),xaxt="n",xaxs="i")
  plot(c(1,1),type='n',xlim=c(start,end), ylim=c(0,max(cell_number_vec)+1),xaxt="n",xaxs="i")  

  mtext(length(orderRNABarcode), side = 1,cex=2,line=3)
  legend("topright", legend = c(0,"","","",upvalue), fill = brewer.pal(9,"Blues")[c(1,3,5,7,9)], bty = "n",border = brewer.pal(9,"Blues")[c(1,3,5,7,9)],y.intersp = 0.5)

  if(length(orderRNABarcode)!=0 && showBarcode){
    for (i in (1:length(orderRNABarcode))) {
      tmpdata = mylist[["scRNAfragmentflt"]][mylist[["scRNAfragmentflt"]]$barcode==orderRNABarcode[i],]
      if (nrow(tmpdata) != 0){
        for (row in (1:nrow(tmpdata))) {

          if (tmpdata$end[row] - tmpdata$start[row] < 150) {
            xleft = tmpdata$start[row]
            ybottom = 1*i + 0.25 + adjust_pos
            xright = tmpdata$end[row]
            ytop = 1*i + 0.75 + adjust_pos
            if(tmpdata$count[row]+4 > upvalue+1){
                color_count=upvalue+1
            }
            else{
                color_count=tmpdata$count[row]+4
            }
            #rect(xleft, ybottom, xright, ytop,col=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),border=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),lwd=NA)

            draw.circle((xleft+xright)/2,(ybottom+ytop)/2,0.45,lwd = point_ratio * 0.15 / 0.4,col=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),border=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8))
          } else {
            xleft = tmpdata$start[row]
            ybottom = 1*i + 0.25 + adjust_pos
            xright = tmpdata$end[row]
            ytop = 1*i + 0.75 + adjust_pos
            if(tmpdata$count[row]+4 > upvalue+1){
                color_count=upvalue+1
            }
            else{
                color_count=tmpdata$count[row]+4
            }
            #rect(c(xleft,xright-80), ybottom, c(xleft+80,xright), ytop,col=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.9),border=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.9),lwd=NA)

            draw.circle((xleft+xleft+80)/2,(ybottom+ytop)/2,0.45,lwd = point_ratio * 0.15 / 0.4,col=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),border=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8))

            draw.circle((xright+xright+80)/2,(ybottom+ytop)/2,0.45,lwd = point_ratio * 0.15 / 0.4,col=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8),border=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.8))
          }

          #rect(xleft, ybottom, xright, ytop,col=color[tmpdata$count[row]],border=rgb(col2rgb(color[color_count])[1]/255,col2rgb(color[color_count])[2]/255,col2rgb(color[color_count])[3]/255,alpha = 0.9),lwd=NA)

          #rect(xleft, ybottom, xright, ytop,col=color[tmpdata$count[row]],border=color[color_count],lwd=NA)
          # draw.circle((xleft+xright)/2,(ybottom+ytop)/2,0.4,lwd = 3,col=color,border=color)
        }
      }
      
    }
  }

  par(cex=0.5, mai=c(0,0.2,0.1,0.1))
  par(fig=c(x_pic1,x_pic2,parlist[2],parlist[3]), new=T)
  #plot(c(1,1),type='n',xlim=c(start,end), ylim=c(0,length(orderATACBarcode)+1),xaxt="n",xaxs="i")
  plot(c(1,1),type='n',xlim=c(start,end), ylim=c(0,max(cell_number_vec)+1),xaxt="n",xaxs="i")

  legend("topright", legend = c(1,2),pch = c(16,16), col = c(color_DNA,"black"), bty = "n",border = c(color_DNA,"black"),y.intersp = 0.7)
  if(length(orderATACBarcode)!=0 && showBarcode){
    for (i in (1:length(orderATACBarcode))) {
      tmpdata = mylist[["scATACfragmentflt"]][mylist[["scATACfragmentflt"]]$barcode==orderATACBarcode[i],]
      if (nrow(tmpdata) != 0){
        for (row in (1:nrow(tmpdata))) {
          xleft = tmpdata$start[row]
          ybottom = 1*i - 0 + adjust_pos
          xright = tmpdata$end[row]
          ytop = 1*i + 0.05 + adjust_pos
          color = "black"
          #color = brewer.pal(9,"Reds")[9]
          if (tmpdata$count[row] == 1){
            color = color_DNA
            #color = brewer.pal(9,"Reds")[5]
          }
          #rect(xleft, ybottom, xright, ytop,col=color,border=color,lwd=3)
          draw.circle((xleft+xright)/2,(ybottom+ytop)/2,0.45,lwd = point_ratio * 0.15 / 0.4,col=color[tmpdata$count[row]],border=color[tmpdata$count[row]])
        }
      }
    }
  }
  #### copy
  #CC=c(brewer.pal(8, "Dark2"),brewer.pal(9, "Set1"),brewer.pal(12, "Set3"),brewer.pal(8, "Accent"),brewer.pal(9, "YlOrRd"),brewer.pal(9, "YlOrBr"),brewer.pal(9, "YlGnBu"),brewer.pal(9, "YlGn"),brewer.pal(9, "Reds"),brewer.pal(9, "RdPu"),brewer.pal(9, "Purples"),brewer.pal(9, "PuRd"),brewer.pal(9, "PuBuGn"),brewer.pal(9, "PuBu"),brewer.pal(9, "Oranges"),brewer.pal(9, "Greys"),brewer.pal(9, "Greens"),brewer.pal(9, "Blues"),brewer.pal(11, "BrBG"),brewer.pal(11, "RdGy"))
  #CC=c(brewer.pal(9, "YlOrRd")[7],brewer.pal(9, "YlOrBr")[7],brewer.pal(9, "YlGnBu")[7],brewer.pal(9, "YlGn")[7],brewer.pal(9, "Reds")[7],brewer.pal(9, "RdPu")[7],brewer.pal(9, "Purples")[7],brewer.pal(9, "PuRd")[7],brewer.pal(9, "PuBuGn")[7],brewer.pal(9, "PuBu")[7],brewer.pal(9, "Oranges")[7],brewer.pal(9, "Greens")[7],brewer.pal(9, "Blues")[7],brewer.pal(11, "BrBG")[9])
  #CC=c(brewer.pal(9, "YlOrBr")[7],brewer.pal(9, "Reds")[7],brewer.pal(9, "YlGn")[7],brewer.pal(9, "RdPu")[7],brewer.pal(9, "YlGnBu")[7],brewer.pal(9, "PuRd")[7],brewer.pal(9, "PuBuGn")[7],brewer.pal(9, "Purples")[7])
  ### browser color
  CC=c("#FF0000","#1E1EFF","#007A00","#800080","#FFBD46","#38D1D8","#89410E")

  par(fig=c(x_pic1,x_pic2,parlist[3],parlist[4]), new=TRUE)
  #plot(c(1,1),type='n',xlim=c(start,end), ylim=c(0,length(orderPETBarcode)+1),xaxt="n",xaxs="i")
  plot(c(1,1),type='n',xlim=c(start,end), ylim=c(0,max(cell_number_vec)+1),xaxt="n",xaxs="i")

  if(length(orderPETBarcode)!=0 && showBarcode){
    for (i in (1:length(orderPETBarcode))) {
      tmpdata = mylist[["scATACPETflt"]][mylist[["scATACPETflt"]]$barcode==orderPETBarcode[i],]
      if (nrow(tmpdata) != 0){
        for (row in (1:nrow(tmpdata))){
          x1 = (tmpdata$start1[row] + tmpdata$end1[row])/2
          x2 = (tmpdata$start2[row] + tmpdata$end2[row])/2
          height = 0.7
          vert_y  = height
          vert_x= (x1 + x2) / 2
          a = vert_y / ((vert_x - x1) * (vert_x - x2))
          if(direction){
            color_list <- c("#25CC69","#2F7E4F","#8080804A")
            color_list_curve <- c("#25CC69","#25CC69","#8080804A")
            ypos <- 1*i + 0.2 + adjust_pos
            curve(ypos + 0 * a * (x - x1) * (x - x2), from=x1, to=x2, add=TRUE, col=color_list_curve[tmpdata$anchor1[row]],lwd=0.25)
            draw.circle(x1,ypos,0.45,lwd = point_ratio,col=color_list[tmpdata$anchor1[row]],border=color_list[tmpdata$anchor1[row]])
            draw.circle(x2,ypos,0.45,lwd = point_ratio,col=color_list[tmpdata$anchor2[row]],border=color_list[tmpdata$anchor2[row]])
          }else{
            #color = CC[i %% length(CC) + 1]
            color = CC[color_class[i] %% length(CC) + 1]
            ypos <- 1*i + 0.2 + adjust_pos
            # print(color)
            curve(ypos + 0 * a * (x - x1) * (x - x2), from=x1, to=x2, add=TRUE, col=color,lwd=0.25)
            draw.circle(x1,ypos,0.45,lwd = point_ratio,col=color,border=color)
            draw.circle(x2,ypos,0.45,lwd = point_ratio,col=color,border=color)
          }

        }
      }
    }
  }
  

  color = "red"
  #rect(start, length(orderPETBarcode)+3, end, length(orderPETBarcode)+3,col=color,border=color)
  par(cex=0.5, mai=c(0.2,0.2,0.1,0.1))
  par(fig=c(x_pic1,x_pic2,parlist[4],parlist[5]), new=TRUE)
  transcriptinfo = mylist[["transcriptinfo"]]

  if(is.null(mylist[["chromstatflt"]])){
    plot(c(1,1),type='n',xlim=c(start,end), ylim=c(1,max(transcriptinfo$plotrow)+1),xaxs="i")
  }else{
    plot(c(1,1),type='n',xlim=c(start,end), ylim=c(1,max(transcriptinfo$plotrow)+2),xaxs="i")
    for(i in 1:nrow(mylist[["chromstatflt"]])){
      
      color_chromhmm = rgb(as.numeric(strsplit(mylist[["chromstatflt"]]$itemRgb[i],",")[[1]])[1],as.numeric(strsplit(mylist[["chromstatflt"]]$itemRgb[i],",")[[1]])[2],as.numeric(strsplit(mylist[["chromstatflt"]]$itemRgb[i],",")[[1]])[3],max=255)

      rect(mylist[["chromstatflt"]]$start[i], max(transcriptinfo$plotrow)+1, mylist[["chromstatflt"]]$end[i], max(transcriptinfo$plotrow)+2,col=color_chromhmm,border=NA)
    }
  }
  
  geneflt = mylist[["geneflt"]]
  if (nrow(transcriptinfo) != 0){
    for (row in (1:nrow(transcriptinfo))) {
      ID = transcriptinfo$names[row]
      subgeneinfo  = geneflt[which(geneflt$ID == ID),]
      if (transcriptinfo$strand[row]=="+"){
        #color = "green"
        color = rgb(23/255,128/255,15/255)
      } else {
        color = rgb(21/255,29/255,252/255)
        #color = "blue"
      }
      rect(subgeneinfo$start, transcriptinfo$plotrow[row] + 0, subgeneinfo$end, transcriptinfo$plotrow[row] + 0 + 0.5 ,col=color,border=NA)
      segments(transcriptinfo$starts[row], transcriptinfo$plotrow[row] + 0.25 , transcriptinfo$stops[row], transcriptinfo$plotrow[row] + 0.25 ,col=color)
    }
  }  
  par(cex=0.5, mai=c(0,0.2,0.1,0.1))
  par(fig=c(x_pic1,x_pic2,parlist[5],parlist[6]), new=TRUE)
  gexbedgraphflt = mylist[["scATACGEXCOVflt"]]
  plot(c(0,0),type='n',xlim=c(start,end), ylim=c(0,yaxis[1]),xaxt="n",xaxs="i",yaxs="i")
  #color = "purple"
  color = "brown"
  if (nrow(gexbedgraphflt) != 0 && showGex){
    for (row in (1:nrow(gexbedgraphflt))) {
      xleft = gexbedgraphflt$start[row]
      ybottom = 0
      xright = gexbedgraphflt$end[row]
      ytop = gexbedgraphflt$value[row]
      #color = "blue"
      color = rgb(21/255,29/255,252/255)
      if(xleft < start) xleft = start
      if(xright > end) xright = end
      rect(xleft, ybottom, xright, ytop,col=color,border=color)
    }
  }

  par(fig=c(x_pic1,x_pic2,parlist[6],parlist[7]), new=TRUE)
  arcbedgraphflt = mylist[["scATACARCCOVflt"]]
  plot(c(0,0),type='n',xlim=c(start,end), ylim=c(yaxis[2],0),xaxt="n",xaxs="i",yaxs="i")
  #color = "purple"

  if (nrow(arcbedgraphflt) != 0 && showArc){
    for (row in (1:nrow(arcbedgraphflt))) {
      xleft = arcbedgraphflt$start[row]
      ybottom = 0
      xright = arcbedgraphflt$end[row]
      ytop = arcbedgraphflt$value[row] 
      #color = "blue"
      #color = rgb(21/255,29/255,252/255)
      color = color_DNA
      if(xleft < start) xleft = start
      if(xright > end) xright = end
      rect(xleft, ybottom, xright, ytop,col=color,border=color)
    }
  }
  # color = "black"
  # peakflt = mylist[["scATACpeakflt"]]
  # rect(peakflt$start, 0.75, peakflt$end, 1.5,col=color,border=color)

  if (twosup) {
    if(givenAnchor){
      loopflt = mylist[["scATACloopgivenanchorflt"]]
    }else{
      loopflt = mylist[["scATAClooptwopeakflt"]]
    }
  } else {
    loopflt = mylist[["scATACloopflt"]]
  }

  
  par(cex=0.5, mai=c(0,0.2,0.1,0.1))

  par(fig=c(x_pic1,x_pic2,parlist[7],parlist[8]), new=TRUE)
  plot(c(0,0),type='n',xlim=c(start,end), ylim=c(0,yaxis[3]),xaxt="n",xaxs="i",yaxs="i")

    #if (nrow(loopflt) > 0){
    if ((!(is.null(loopflt) || nrow(loopflt) == 0)) && showLoop){
      loopflt$count = log10(loopflt$count + 1)
    for (row in (1:nrow(loopflt))) {
      x1 = (loopflt$start1[row] + loopflt$end1[row])/2
      x2 = (loopflt$start2[row] + loopflt$end2[row])/2
      if(twosup&&givenAnchor){
        if(showPETnum){
          height = loopflt$count[row]
        }else{
          height = loopflt$score[row]
        }
      }else{
        height = loopflt$count[row]
      }
      # vert_y  = height
      vert_y  = height
      vert_x= (x1 + x2) / 2
      a = vert_y / ((vert_x - x1) * (vert_x - x2))
      #color = "magenta"
      #color = rgb(126/255,3/255,127/255)
      color = color_DNA
      if(xleft < start) xleft = start
      if(xright > end) xright = end
      curve(a * (x - x1) * (x - x2), from=x1, to=x2, add=TRUE, col=color,lwd = 0.5)
    }

  }

}

adjustBin = function(scATACARCCOVflt,chrom,start,end){
  
  step = as.integer((end-start)/2000)
  bin <- data.frame(chr=chrom,start=seq(start,end-step,step))
  bin$end <- bin$start+step
  bin$queryHits <- seq(1,nrow(bin),by=1)
  bin.gr <- makeGRangesFromDataFrame(bin)

  names(scATACARCCOVflt) <- c("chr","start","end","value")
  scATACARCCOVflt$subjectHits <- seq(1,nrow(scATACARCCOVflt),by=1)

  scATACARCCOVflt.gr <- makeGRangesFromDataFrame(scATACARCCOVflt)
  x <- findOverlaps(bin.gr,scATACARCCOVflt.gr)
  x.data <- as.data.frame(x)
  x.merge.data <- merge(x.data,scATACARCCOVflt,by="subjectHits")
  x.final.data <- x.merge.data %>% group_by(queryHits) %>%  summarise(average = sum(value * (end-start))/step)
  
  x.bin.final.data <- merge(bin,x.final.data,by="queryHits")
  
  scATACARCCOVflt_adjustBin <- x.bin.final.data[,c("chr","start","end","average")]
  names(scATACARCCOVflt_adjustBin) <- c("chrom","start","end","value")
  return(scATACARCCOVflt_adjustBin)
}

### for boxplot
summaryandplot <- function(chrom,start,mid,end){
  library(ggpubr)
  library(ggplot2)
  library(patchwork)
  dis = 8000
  PETcount = 2 
  condition = "PETDNARNA"
  sample_num = -1
  clean = T
  VCscore=0.001
  
  mylist = filterData(data_list[[1]][["genelist"]],data_list[[1]][["scATACfragmentlist"]],data_list[[1]][["scATACpeaklist"]],data_list[[1]][["scATACARCCOVrpkmlist"]],data_list[[1]][["scATACPETlist"]],data_list[[1]][["scATAClooplist"]],data_list[[1]][["scATAClooptwopeaklist"]],data_list[[1]][["scATACloopgivenanchorlist"]],data_list[[1]][["scATACGEXdupCOVrpkmlist"]],data_list[[1]][["scRNAfragmentlist"]],data_list[[1]][["chromstatlist"]],chrom,start,end,dis=dis,PETcount = PETcount,VCscore=VCscore,condition=condition,sample_num = sample_num,cleanPET=clean)

  Total_FRAGMENT = mylist[["scATACfragmentflt"]][mylist[["scATACfragmentflt"]]$end<mid,] # fragment
  Total_PET = mylist[["scATACPETflt"]][mylist[["scATACPETflt"]]$end2<mid,] # PET
  Total_GEX = mylist[["scRNAfragmentflt"]][mylist[["scRNAfragmentflt"]]$end<mid,] # GEX
  Total_FRAGMENT2 = mylist[["scATACfragmentflt"]][mylist[["scATACfragmentflt"]]$start>mid,] # fragment
  Total_PET2 = mylist[["scATACPETflt"]][mylist[["scATACPETflt"]]$start1>mid,] # PET
  Total_GEX2 = mylist[["scRNAfragmentflt"]][mylist[["scRNAfragmentflt"]]$start>mid,] # GEX

  FRAGMENT = data.frame(barcode=names(table(Total_FRAGMENT$barcode)),barcodecount=as.vector(table(Total_FRAGMENT$barcode)))
  FRAGMENT2 = data.frame(barcode=names(table(Total_FRAGMENT2$barcode)),barcodecount=as.vector(table(Total_FRAGMENT2$barcode)))
  PET = data.frame(barcode=names(table(Total_PET$barcode)),barcodecount=as.vector(table(Total_PET$barcode)))
  PET2 = data.frame(barcode=names(table(Total_PET2$barcode)),barcodecount=as.vector(table(Total_PET2$barcode)))
  GEX = data.frame(barcode=names(table(Total_GEX$barcode)),barcodecount=as.vector(table(Total_GEX$barcode)))
  GEX2 = data.frame(barcode=names(table(Total_GEX2$barcode)),barcodecount=as.vector(table(Total_GEX2$barcode))) 
  print(sprintf("There are %d dots in 1st list", sum(FRAGMENT$barcodecount)))
  print(sprintf("There are %d loops in 1st list", sum(PET$barcodecount)))
  print(sprintf("There are %d GEX in 1st list", sum(GEX$barcodecount)))
  print(sprintf("There are %d cells in 1st list", length(GEX$barcode)))
  print(sprintf("There are %d dots in 2st list", sum(FRAGMENT2$barcodecount)))
  print(sprintf("There are %d loops in 2st list", sum(PET2$barcodecount)))
  print(sprintf("There are %d GEX in 2st list", sum(GEX2$barcodecount)))
  print(sprintf("There are %d cells in 2st list", length(GEX2$barcode)))

  m=22
  theme = theme(title=element_text(size=12,hjust=0.2,lineheight=0.2),
              axis.title.x=element_text(size=m,hjust=0.5),
              axis.title.y=element_text(size=m,hjust=0.5),
              axis.text.x=element_text(size=m),
              axis.text.y=element_text(size=m),
              axis.ticks = element_line(size = 1,),
              axis.line = element_line(colour = "black"),
              panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(),
              panel.border = element_blank(), #legend.position="none",
              legend.text = element_text(size=m),
              legend.title = element_text(size=m),
              plot.title = element_text(colour = "black", face = "bold", size = m, vjust = 1,hjust = 0.5))

  tmp <- FRAGMENT
  tmp$group <- "high"
  tmp2 <- FRAGMENT2
  tmp2$group <- "low"
  tmp3 <- rbind(tmp,tmp2)
  #p1 <- ggboxplot(tmp3,y="barcodecount",x="group",add="mean",fill = "group",outlier.shape = NA) + stat_compare_means(comparisons=list(c("high","low")),method="wilcox.test",label = "p.signif",size=8,label.y=11.5) + xlab("")+ylab("Fragment")+labs(fill="") + theme + theme(strip.text.x = element_text(size = 12))+ coord_cartesian (ylim = c (0, 14))

  p1 <- ggplot(tmp3, aes(x = group, y = barcodecount, fill = group)) +
      geom_violin(width=1.2) +
      geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
      #stat_summary(fun = "mean", geom = "point", shape = 8, size = 2.5, color = "black", alpha = 0.7) +
      #stat_summary(fun.data = function(x) data.frame(y = y.umi, label = paste("Median=", round(median(x), 3))), geom = "text", size = 5) +
      stat_compare_means(comparisons=list(c("high","low")),method="wilcox.test",label = "p.format",size=4) + 
      stat_compare_means(comparisons=list(c("high","low")),method="wilcox.test",label = "p.signif",size=8) + 
      scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF")) +
      labs(title = "", x= "", y = "Fragment") +
      theme +
      theme(legend.position = "none") + 
      ylim(0,max(tmp3$barcodecount) * 1.1)

  tmp <- PET
  tmp$group <- "high"
  tmp2 <- PET2
  tmp2$group <- "low"
  tmp3 <- rbind(tmp,tmp2)
  tmp3 <- tmp3[tmp3$barcodecount!=1,]

  p4 <- ggplot(data=tmp3,aes(x=barcodecount,fill=group)) + 
  geom_histogram(binwidth=1,position = position_dodge2(padding = 0.2)) + 
  xlab("PET number") + 
  ylab("Count") + 
  labs(fill="") +
  scale_x_continuous(breaks = seq(1,10,by=1)) + 
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF")) + 
  theme + 
  theme(legend.position = "none")

  tmp <- GEX
  tmp$group <- "high"
  tmp2 <- GEX2
  tmp2$group <- "low"
  tmp3 <- rbind(tmp,tmp2)
  tmp3 <- tmp3[tmp3$barcodecount!=0,]
  #p3 <- ggboxplot(tmp3,y="barcodecount",x="group",add="mean",fill = "group",outlier.shape = NA) + stat_compare_means(comparisons=list(c("high","low")),method="wilcox.test",label = "p.signif",size=8,label.y=11.5) + xlab("")+ylab("GEX")+labs(fill="") + theme + theme(strip.text.x = element_text(size = 20)) + coord_cartesian (ylim = c (0, 14))

  p3 <- ggplot(tmp3, aes(x = group, y = log10(barcodecount), fill = group)) +
      geom_violin(width=1.2) +
      geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
      #stat_summary(fun = "mean", geom = "point", shape = 8, size = 2.5, color = "black", alpha = 0.7) +
      #stat_summary(fun.data = function(x) data.frame(y = y.umi, label = paste("Median=", round(median(x), 3))), geom = "text", size = 5) +
      #scale_fill_manual(values = c("#FF0000","#1E1EFF")) +
      stat_compare_means(comparisons=list(c("high","low")),method="wilcox.test",label = "p.format",size=4) + 
      stat_compare_means(comparisons=list(c("high","low")),method="wilcox.test",label = "p.signif",size=8) +
      scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF")) +
      labs(title = "", x= "", y = expression(paste("Gene expression (",log[10],"count)"))) +
      theme +
      theme(legend.position = "right") + 
      ylim(0,max(log10(tmp3$barcodecount)) * 1.1)

  p <- p1 | p4 | p3
  ggsave("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/high_low_boxplot.pdf",p,width = 15,height = 6)

 #"#FF0000","#1E1EFF"
}
### for counting cell number for each type
countTotalNumber <- function(data_list=data_list){
  name_tran <- c()
  for(i in 1:length(data_list)){
    name_tran[i] <- sub("all","",strsplit(names(data_list),"\\.")[[i]][2])
  }
  count_data <- data.frame(type=name_tran,count=0)
  for(i in 1:length(data_list)){
      count_data$count[i] <- length(unique(data.frame(do.call(rbind,data_list[[i]][["scATACfragmentlist"]]))$barcode))
  }
  write.csv(count_data,file="/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/celltyperesult/count_data.csv",quote=FALSE,row.names=FALSE)
}

downsample <- function(file=c("mouseBrain.Astroependymal.all.data.RData","mouseBrain.Excitatory_Neurons.all.data.RData","mouseBrain.Inhibitory_Neurons.all.data.RData","mouseBrain.Microglia.all.data.RData","mouseBrain.Oligodendrocytes.all.data.RData","mouseBrain.OPC.all.data.RData","mouseBrain.Vascular.all.data.RData"),path="/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitBigClass1/",sample_num=4500,out_path="/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/downsample_data/"){
    
  for(i in file){
      data_list <- list()
      chromstatlist <- list()
      load(paste0(path,i))
      set.seed(1234)

      barcode <- sample(unique(data.frame(do.call(rbind,scATACfragmentlist))$barcode),size=sample_num)

      scRNAfragmentlist_downsample <- list()
      for(j in names(scRNAfragmentlist)){
        scRNAfragmentlist_downsample[[j]] <- scRNAfragmentlist[[j]][scRNAfragmentlist[[j]]$barcode %in% barcode,]
      }

      scATACfragmentlist_downsample <- list()
      for(j in names(scATACfragmentlist)){
        scATACfragmentlist_downsample[[j]] <- scATACfragmentlist[[j]][scATACfragmentlist[[j]]$barcode %in% barcode,]
      }

      scATACPETlist_downsample <- list()
      for(j in names(scATACPETlist)){
        scATACPETlist_downsample[[j]] <- scATACPETlist[[j]][scATACPETlist[[j]]$barcode %in% barcode,]
      }

      scATACfragmentlist <- scATACfragmentlist_downsample
      scATACPETlist <- scATACPETlist_downsample 
      scRNAfragmentlist <- scRNAfragmentlist_downsample

      save(genelist,scATACfragmentlist,scATACpeaklist,scATACARCCOVlist,scATACARCCOVrpkmlist,scATACPETlist,scATAClooplist,scATAClooptwopeaklist,scATACloopgivenanchorlist,scATACGEXCOVlist,scATACGEXCOVrpkmlist,scATACGEXdupCOVlist,scATACGEXdupCOVrpkmlist,scRNAfragmentlist,scRNAdupfragmentlist,chromstatlist, file=paste0(out_path,i))

  }

}


loadData_allele <- function(file=c("patski.allele.all.G1.data.RData"),path="/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineCellCycleFinal/"){
  data_list <- list()
  for(i in file){
      chromstatlist <- list()
      load(paste0(path,i))
      
      tmp <- gsub(".allele.all",".allele.P",i)
      data_list[[tmp]] <- list(genelist,scATACfragmentPlist,scATACpeakPlist,scATACARCCOVPlist,scATACARCCOVPrpkmlist,scATACPETPlist,scATACloopPlist,scATAClooptwopeakPlist,scATACloopgivenanchorPlist,scATACGEXCOVPlist,scATACGEXCOVPrpkmlist,scATACGEXdupCOVPlist,scATACGEXdupCOVPrpkmlist,scRNAfragmentPlist,scRNAdupfragmentPlist,chromstatlist)
      names(data_list[[tmp]]) <- c("genelist","scATACfragmentlist","scATACpeaklist","scATACARCCOVlist","scATACARCCOVrpkmlist","scATACPETlist","scATAClooplist","scATAClooptwopeaklist","scATACloopgivenanchorlist","scATACGEXCOVlist","scATACGEXCOVrpkmlist","scATACGEXdupCOVlist","scATACGEXdupCOVrpkmlist","scRNAfragmentlist","scRNAdupfragmentlist","chromstatlist")

      tmp <- gsub(".allele.all",".allele.M",i)
      data_list[[tmp]] <- list(genelist,scATACfragmentMlist,scATACpeakMlist,scATACARCCOVMlist,scATACARCCOVMrpkmlist,scATACPETMlist,scATACloopMlist,scATAClooptwopeakMlist,scATACloopgivenanchorMlist,scATACGEXCOVMlist,scATACGEXCOVMrpkmlist,scATACGEXdupCOVMlist,scATACGEXdupCOVMrpkmlist,scRNAfragmentMlist,scRNAdupfragmentMlist,chromstatlist)
      names(data_list[[tmp]]) <- c("genelist","scATACfragmentlist","scATACpeaklist","scATACARCCOVlist","scATACARCCOVrpkmlist","scATACPETlist","scATAClooplist","scATAClooptwopeaklist","scATACloopgivenanchorlist","scATACGEXCOVlist","scATACGEXCOVrpkmlist","scATACGEXdupCOVlist","scATACGEXdupCOVrpkmlist","scRNAfragmentlist","scRNAdupfragmentlist","chromstatlist")
  }
  return(data_list)
}




