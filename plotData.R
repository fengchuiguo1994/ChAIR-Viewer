#### for mm10
####  linux
tmux
srun -N 1 -p ruan_cpu -n 5 --mem=200G --pty /bin/bash
/data/home/ruanlab/huangxingyu/miniconda3/envs/R4/bin/R
#### R

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("k562.all.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/"
data_list <- loadData(file=file,path=path)

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("patski.all.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/"
data_list <- loadData(file=file,path=path,readbulkChIATAC=T)

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("patski.all.allele.P.data.RData","patski.all.allele.M.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/"
data_list <- loadData(file=file,path=path)

### downsampled to 5m data
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("patski.G1.all.data.RData","patski.S.all.data.RData","patski.G2M.all.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellCycle/"
data_list <- loadData(file=file,path=path)

### downsampled to 5m data
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("k562.G1.all.data.RData","k562.S.all.data.RData","k562.G2M.all.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellCycle/"
data_list <- loadData(file=file,path=path)

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitBigClass1/alldata/"
file = c("mouseBrain.Astroependymal.chr16.data.RData","mouseBrain.Excitatory_Neurons.chr16.data.RData","mouseBrain.Inhibitory_Neurons.chr16.data.RData","mouseBrain.Microglia.chr16.data.RData","mouseBrain.Oligodendrocytes.chr16.data.RData","mouseBrain.OPC.chr16.data.RData","mouseBrain.Vascular.chr16.data.RData")
data_list <- loadData(file=file,path=path)
## downsample(file=file,path=path)

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitBigClass1/"
file = c("mouseBrain.P2_Excitatory_Neurons.all.data.RData","mouseBrain.P11_Excitatory_Neurons.all.data.RData","mouseBrain.P95_Excitatory_Neurons.all.data.RData","mouseBrain.P365_Excitatory_Neurons.all.data.RData","mouseBrain.P730_Excitatory_Neurons.all.data.RData")
data_list <- loadData(file=file,path=path)


source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByBrainCellTypeFinal/"
file = list.files(path=path, pattern="*.data.RData" ,all.files=FALSE,full.names=FALSE)
data_list <- loadData(file=file,path=path)

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitCellType2/"
file = list.files(path=path, pattern="*.data.RData" ,all.files=FALSE,full.names=FALSE)
data_list <- loadData(file=file,path=path)

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitCellType2/"
file = c("mouseBrain.CBGRC.all.data.RData","mouseBrain.TEGLU.all.data.RData","mouseBrain.MOL.all.data.RData","mouseBrain.TEINH.all.data.RData","mouseBrain.MSN.all.data.RData","mouseBrain.ACTE.all.data.RData","mouseBrain.OPC.all.data.RData","mouseBrain.MGL.all.data.RData","mouseBrain.VLMC.all.data.RData")
data_list <- loadData(file=file,path=path)


## CBGRC
# source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
# path = "/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/CBGRCauto/"
# file = c("mouseBrain.Young_CBGRC.all.data.RData","mouseBrain.P95_CBGRC.all.data.RData","mouseBrain.P365_CBGRC.all.data.RData","mouseBrain.P730_CBGRC.all.data.RData")
# # file = c("mouseBrain.P11_CBGRC.all.data.RData","mouseBrain.P95_CBGRC.all.data.RData","mouseBrain.P365_CBGRC.all.data.RData","mouseBrain.P730_CBGRC.all.data.RData")
# data_list <- loadData(file=file,path=path)
# ### downsampled CBGRC
# downsample(file=file,path=path,sample_num=7500)
# path = "/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/downsample_data/"
# data_list <- loadData(file=file,path=path)

## CBGRC2
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/CBGRCauto2_P95rep2/"
file = c("mouseBrain.Young_CBGRC.all.data.RData","mouseBrain.P95_Rep2_CBGRC.all.data.RData","mouseBrain.P365_CBGRC.all.data.RData","mouseBrain.P730_CBGRC.all.data.RData")
downsample(file=file,path=path,sample_num=7500,out_path="/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/CBGRCauto2_P95rep2/downsampled/")
path = "/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/CBGRCauto2_P95rep2/downsampled/"
data_list <- loadData(file=file,path=path)

## CBGRC3
# source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
# path = "/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/CBGRCauto3_Old/"
# file = c("mouseBrain.Young_CBGRC.all.data.RData","mouseBrain.Old_CBGRC.all.data.RData")
# downsample(file=file,path=path,sample_num=7500,out_path="/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/CBGRCauto3_Old/downsampled/")
# path = "/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/CBGRCauto3_Old/downsampled/"
# data_list <- loadData(file=file,path=path)



region = " chr1:92455000-92972000 " # chr16 92455098 92972441/chr16 92,455,098 92,972,441/chr16:92455098-92972441/chr16:92,455,098-92,972,441
region = "chr19:20996056-21660000"

region = "chr13:97100000-97200000"

region = "chr19:8400000-10300000" ##2M

region = "chr19:8670000-9130000" ## 0.5M

region = "chr19:8700000-8800000" ## 0.1M

region = "chr15:61950000-62700000" ## Myc

region = "chr15:61950000-62420000" ## micro Myc

region = "chr19:5783000-5855000" ## Neat1 Malat1
gene = "Neat1"
gene = "Malat1"

region = "chrX:103450000-103635000" ## Xist
gene = "Xist"
#cut_region = "chrX:103479533-103487250"

region = "chr18:36868000-38033000" ## Pcdh

region = "chr7:6650000-6950000" ## Peg3

region = "chrX:10466000-11032000" ## Mid1

region = "chrX:152120000-152420000" ## Kdm5c

region = "chrX:13104000-13770000" ## Ddx3x

region = "chrX:50516000-50676000" ## Firre

region = "chrX:105030000-105150000" ## Magee1

region = "chr7:142200000-142700000" ## H19

region = "chr17:47550000-48000000" ## Foxp4

region = "chr18:36868000-38033000" ## Pcdh

region = "chr17:35456000-35956000" ## high
region = "chr17:35857000-35956000" ## micro high

region = "chr11:114622000-114866000" ## Rpl38

region = "chr19:8397000-10274000" ## high and low GEX
region = "chr19:9846000-11144000" ## high and low GEX
region = "chr19:9846000-10513532" ## high and low GEX ### the newest
region = "chr19:9960000-10460000"  ### the newest

region = "chr19:10809887-11098897" ## 20240110
region = "chr19:10832909-10989547"
region = "chr19:10888505-10966564"
# mouse brain
region = "chr10:55584009-57195848"
region = "chr10:55750573-58085052"
gene = "Gja1" #Astroependymal

region = "chr18:15226341-15572148"
region = "chr18:14707632-16090863"
gene = "Aqp4" #Astroependymal

region = "chr16:91201643-91340074" 
region = "chr16:91065636-91519137"
gene = "Olig1"	 #Oligodendrocytes

region = "chr9:119911370-120419833" 
gene = "Mobp"	 #Oligodendrocytes

region = "chr2:33868930-34675677"
gene = "Pbx3" ##OBINH2 OBNBL3

region = "chr2:115761679-116167230"
gene = "Meis2" ##OBINH2 OBNBL3

region = "chr2:127598022-127691905"
gene = "Mal" ## OPC MOL1 MOL2

region = "chr15:8316803-9217422"
gene = "Slc1a3" ## ACBG ACNT ACTE

region = "chr1:56637665-57590088"
gene = "Satb2" ## DGNBL TEGLU(marker) TEINH 

region = "chr9:111944885-112610504"
gene = "Arpp21" ## DGNBL TEGLU(marker) TEINH 

region = "chr7:114620000-116890000"
gene = "Sox6" ## OPC(marker) MOL

region = "chr9:119843000-120352000"
gene = "Mobp" ## OPC MOL(marker)

region = "chr1:159163000-159979000"
gene = "Tnr" ## OPC(marker) MOL

region = "chr7:114340023-115406260"
gene = "Insc" ##

region = "chr12:24420000-25010000"
gene = "Rrm2" ## S

region = "chr10:90960000-91310000"
gene = "Tmpo" # G2M

region = "chr18:71000000-72600000"
region = "chr18:70000000-75000000"
region = "chr18:60000000-80000000"
gene = "Dcc" # excit compartmentRNA gene
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

region = "chr18:66048000-73286000"
cut_region1="chr18:66413541-66984481"
run(region,file=file,cut_region1=cut_region1,cut_region2=cut_region2,cut=T,gene=gene,direction =F,showAllPET=F,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

#### Patj
region = "chr4:98187000-98987000"
gene = "Patj"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

region = "chr4:96000000-106000000"
gene = "Patj"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

region = "chr4:96000000-106000000"
cut_region1="chr4:97543237-99166296"
cut_region2="chr4:104248337-105215077"
run(region,file=file,cut_region1=cut_region1,cut_region2=cut_region2,cut=T,gene=gene,direction =F,showAllPET=F,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

#### for hg38
#### linux
tmux
srun -N 1 -p ruan_cpu -n 5 --mem=200G --pty /bin/bash
/data/home/ruanlab/huangxingyu/miniconda3/envs/R4/bin/R
#### R
load("k562.all.data.rds")

region = "chr9:136475000-136650000" # NOTCH1

region = "chrX:48772000-48950000" # GATA1

region = "chr8:127649000-128249000" # MYC

region = "chr11:65402000-65522000" #NEAT1 & MALAT1

region = "chr5:163370000-163520000"
gene = "HMMR" ## g2m gene

region = "chr12:56610000-56900000"
gene = "PRIM1" ## s gene

region = "chr11:61716000-62149000"
gene = "FEN1" ## s gene

##########################################
#### read P-M data
# aa1 = filterData(genelist,scATACfragmentPlist,scATACpeaklist,scATACARCCOVPrpkmlist,scATACPETPlist,scATACloopPlist,scATAClooptwopeakPlist,scATACGEXCOVPrpkmlist,scRNAfragmentPlist,chrom,start,end,dis=dis,PETcount = PETcount,PETDNARNA=PETDNARNA)
# aa2 = filterData(genelist,scATACfragmentMlist,scATACpeaklist,scATACARCCOVMrpkmlist,scATACPETMlist,scATACloopMlist,scATAClooptwopeakMlist,scATACGEXCOVMrpkmlist,scRNAfragmentMlist,chrom,start,end,dis=dis,PETcount = PETcount,PETDNARNA=PETDNARNA)
# aa_list <- list(aa1,aa2)
###################################

### for high and low 

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("k562.all.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/"
for(read in c("bulkChIATAC","bulkATAC","bulkATAC","bulkRNA","")){
data_list <- loadData(file=file,path=path,addLoad="/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/k562.bulk.all.data.RData",readOtherData=read)
region = "chr11:65367000-65517000"
if(read==""){
    run(region,file=file,direction =F,cell_number=800,PETcount=3,VCscore=0.001,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=T,readOtherData=read,outfmt="pdf",output=paste0("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/k562.","scChaiT",".pdf"))
}else{
    run(region,file=file,direction =F,cell_number=800,PETcount=3,VCscore=0.001,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=T,readOtherData=read,outfmt="pdf",output=paste0("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/k562.",read,".pdf"))
}
}

read=""
region = "chr11:65367000-65517000"
data_list <- loadData(file=file,path=path,addLoad="/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/k562.bulk.all.data.RData",readOtherData=read)
region = "chr11:65367000-65517000"
run(region,file=file,direction =F,cell_number=800,PETcount=3,VCscore=0.001,type="PET",condition="PETDNARNA",sample_num = -1,clean=T,showPETnum=T,readOtherData=read,outfmt="pdf",output=paste0("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/k562.tss","scChaiT",".pdf"))

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("patski.all.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/"
for(read in c("bulkChIATAC","bulkATAC","bulkATAC","bulkRNA","")){
data_list <- loadData(file=file,path=path,addLoad="/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/patski.bulk.all.data.RData",readOtherData=read)
region = "chr19:9960000-10460000"  ### the newest
if(read==""){
    run(region,file=file,direction =F,cell_number=800,PETcount=3,VCscore=0.001,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=T,readOtherData=read,outfmt="pdf",output=paste0("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/patski.","scChaiT",".pdf"))
}else{
   run(region,file=file,direction =F,cell_number=800,PETcount=3,VCscore=0.001,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=T,readOtherData=read,outfmt="pdf",output=paste0("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/patski.",read,".pdf")) 
}
}

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("patski.chr19.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/alldata/"
read=""
data_list <- loadData(file=file,path=path,addLoad="/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/patski.bulk.all.data.RData",readOtherData=read)
region = "chr19:9960000-10460000"  ### the newest
run(region,file=file,direction =F,cell_number=800,PETcount=3,VCscore=0.001,type="PET",condition="PETDNARNA",sample_num = -1,clean=T,showPETnum=T,readOtherData=read,outfmt="pdf",output=paste0("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/patski.tss","scChaiT",".pdf"))

### for patski cell cycle
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("patski.all.G1.data.RData","patski.all.S.data.RData","patski.all.G2M.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineCellCycleFinal/"
data_list <- loadData(file=file,path=path)

region = "chr10:90960000-91310000"
gene = "Tmpo" # G2M
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 7500,clean=F,showPETnum=F,outfmt="pdf")
region = "chr12:24420000-25010000"
gene = "Rrm2" ## S
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 7500,clean=F,showPETnum=F,outfmt="pdf")

region = "chrX:103450000-103635000" ## Xist
gene = "Xist" ## P most
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 7500,clean=F,showPETnum=F,outfmt="pdf")

### for k562 cell cycle 
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("k562.all.G1.data.RData","k562.all.S.data.RData","k562.all.G2M.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineCellCycleFinal/"
data_list <- loadData(file=file,path=path)

region = "chr5:163370000-163520000"
gene = "HMMR" ## g2m gene
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 10000,clean=F,showPETnum=F,outfmt="pdf")
region = "chr12:56610000-56900000"
gene = "PRIM1" ## s gene
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 10000,clean=F,showPETnum=F,outfmt="pdf")
region = "chr11:61716000-62149000"
gene = "FEN1" ## s gene
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 10000,clean=F,showPETnum=F,outfmt="pdf")

### for P and M allele
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("patski.all.allele.P.data.RData","patski.all.allele.M.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/"
data_list <- loadData(file=file,path=path)

region = "chrX:103450000-103635000" ## Xist
gene = "Xist" ## P most
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",,sameCell=T,sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

region = "chrX:73850000-74511000"
gene = "Flna" ## M most and around chrX
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",,sameCell=T,sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

#### for P and M allele also for G1 S G2M

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
region = "chrX:103450000-103635000" ## Xist
gene = "Xist" ## P most

for(k in c("G1","S","G2M")){
file = c(paste0("patski.allele.all.",k,".data.RData"))
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineCellCycleFinal/"
data_list <- loadData_allele(file=file,path=path)
file = c(paste0("patski.allele.P.",k,".data.RData"),paste0("patski.allele.M.",k,".data.RData"))
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=3,VCscore=0.001,type="PET",condition="PETDNARNA",,sameCell=T,sample_num = 7980,clean=F,showPETnum=F,outfmt="png")
}



# region = "chr4:140625122-141340116"
# gene = "Sdhb" ## M most
# region = "chrX:95293747-96319027"
# gene = "Zc3h12b" ## M most and around chrX
# region = "chrX:105888000-106334000"
# gene = "Cox7b" ## M most
# region = "chrX:98735604-99787538"
# gene = "Efnb1" ## M most and around Xist
# region = "chrX:6746192-8073965"
# gene = "Clcn5" ## M most and around chrX

## 7 big class list
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitBigClass1/alldata/"


# region_gene_list <- list(c("chr18:15334000-15501000","Aqp4"),c("chr7:45121000-45218000","Slc17a7"),c("chr2:70464000-70671000","Gad1"),c("chr4:136812000-136975000","C1qa"),c("chr16:91201000-91340000","Olig1"),c("chr5:74302000-75441000","Pdgfra"),c("chr16:18721000-18872000","Cldn5"),c("chr7:51465000-52094000","Slc17a6"))

region_gene_list <- list(c("chr7:45121000-45218000","Slc17a7"),c("chr2:70464000-70671000","Gad1"),c("chr18:15334000-15501000","Aqp4"),c("chr16:91201000-91340000","Olig1"),c("chr5:74302000-75441000","Pdgfra"),c("chr16:18721000-18872000","Cldn5"),c("chr4:136812000-136975000","C1qa"))

for(i in 1:length(region_gene_list)){
    region = region_gene_list[[i]][[1]]
    gene = region_gene_list[[i]][[2]]
    region = gsub("-|:"," ",region)
    regionarr = strsplit(trimws(region)," |\t")
    chrom = regionarr[[1]][1]
    start = as.numeric(regionarr[[1]][2])
    end = as.numeric(regionarr[[1]][3])

    file = c(paste0("mouseBrain.Excitatory_Neurons.",chrom,".data.RData"),paste0("mouseBrain.Inhibitory_Neurons.",chrom,".data.RData"),paste0("mouseBrain.Astroependymal.",chrom,".data.RData"),paste0("mouseBrain.Oligodendrocytes.",chrom,".data.RData"),paste0("mouseBrain.OPC.",chrom,".data.RData"),paste0("mouseBrain.Vascular.",chrom,".data.RData"),paste0("mouseBrain.Microglia.",chrom,".data.RData"))  
    file = file[i]
    data_list <- loadData(file=file,path=path)
    run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")
    #run(region,file=file,gene=gene,direction =F,cell_number=800,PETcount=2,VCscore=0.001,type="PET",condition="PETDNARNA",sample_num = -1,clean=T,showPETnum=F,outfmt="pdf")
}


##### Dcc and Patj
## one gene
region = "chr18:71000000-73100000"
gene = "Dcc"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

region = "chr8:22846000-23285000"
gene = "Ank1"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

region = "chr8:23500000-24497000"
gene = "Zmat4"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")


## compartment
# region = "chr18:60000000-80000000"
# gene = "Dcc"
# run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

# region = "chr8:20000000-40000000"
# gene = "Ank1"
# run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")


##cut region
## DCC
region = "chr18:60000000-80000000"
cut_region1= c("chr18:70705100-73090908")
cut_region2= c("chr18:66119734-66953436","chr18:63005530-63554160","chr18:74314856-74705099","chr18:75803069-75997855","chr18:76402444-76871617")
run(region,file=file,cut_region1=cut_region1,cut_region2=cut_region2,cut=T,gene=gene,direction =F,showAllPET=F,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf",showBarcode=T,showGex=F,showArc=F,showLoop=F)
## Ank1
region = "chr8:20000000-40000000"
cut_region1= c("chr8:22403547-23246118")
cut_region2= c("chr8:25028825-26997782","chr8:31005464-31295160","chr8:33623060-36567627")
run(region,file=file,cut_region1=cut_region1,cut_region2=cut_region2,cut=T,gene=gene,direction =F,showAllPET=F,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf",showBarcode=T,showGex=F,showArc=F,showLoop=F)
## Zmat4
region = "chr8:20000000-40000000"
cut_region1= c("chr8:23352550-24886917")
cut_region2= c("chr8:27308204-30882482","chr8:31325942-33179600","chr8:36824834-40000000")
run(region,file=file,cut_region1=cut_region1,cut_region2=cut_region2,cut=T,gene=gene,direction =F,showAllPET=F,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf",showBarcode=T,showGex=F,showArc=F,showLoop=F)
### linage 
 
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
# OPC AND OLIGO
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitBigClass1/"
file = c("mouseBrain.OPC.all.data.RData","mouseBrain.Oligodendrocytes.all.data.RData")
data_list <- loadData(file=file,path=path)
region="chr5:74877549-75319997"
gene="Pdgfra"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")
region="chr2:127457197-127832732"
gene="Mal"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

# region="chr5:74877549-75319997"
# gene="Pdgfra"
# run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 5650,clean=F,showPETnum=F,outfmt="png")
# region="chr2:127457197-127832732"
# gene="Mal"
# run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 5650,clean=F,showPETnum=F,outfmt="png")

# CBNBL2 AND CBGRC
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitCellType/"
file = c("mouseBrain.CBNBL2.all.data.RData","mouseBrain.CBGRC.all.data.RData")
data_list <- loadData(file=file,path=path)
# region="chr2:79284180-79634387"
# gene="Neurod1"
region = "chr19:46154336-46743519"
gene = "Sufu"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")
region="chr11:42079568-42547951"
gene="Gabra6"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

region = "chr19:46154336-46743519"
gene = "Sufu"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 5680,clean=F,showPETnum=F,outfmt="png")
region="chr11:42079568-42547951"
gene="Gabra6"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 5680,clean=F,showPETnum=F,outfmt="png")


# OBINH2 AND OBNBL3
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitCellType/"
file = c("mouseBrain.OBINH2.all.data.RData","mouseBrain.OBNBL3.all.data.RData")
data_list <- loadData(file=file,path=path)
region="chr15:81996730-82422265"
gene="Shisa8"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")
region="chr17:3019417-3838508"
gene="Tiam2"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

region="chr15:81996730-82422265"
gene="Shisa8"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 2275,clean=F,showPETnum=F,outfmt="png")
region="chr17:3019417-3838508"
gene="Tiam2"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 2275,clean=F,showPETnum=F,outfmt="png")

region="chr18:45060280-45894611"
gene="Kcnn2"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")
region="chr7:67826891-68388576"
gene="Igf1r"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

region="chr18:45060280-45894611"
gene="Kcnn2"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 2275,clean=F,showPETnum=F,outfmt="png")
region="chr7:67826891-68388576"
gene="Igf1r"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 2275,clean=F,showPETnum=F,outfmt="png")


# DGNBL AND DGGRC
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitCellType2/"
file = c("mouseBrain.DGNBL.all.data.RData","mouseBrain.DGGRC.all.data.RData")
data_list <- loadData(file=file,path=path)
region="chr13:83533470-84589541"
gene="Gm17750"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")
region="chr11:90146275-90581162"
gene="Hlf"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

# region="chr13:83533470-84589541"
# gene="Gm17750"
# run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 2150,clean=F,showPETnum=F,outfmt="png")
# region="chr11:90146275-90581162"
# gene="Hlf"
# run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = 2150,clean=F,showPETnum=F,outfmt="png")

# DGNBL1 DGNBL1 AND DGGRC2
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitCellType/"
file = c("mouseBrain.DGNBL1.all.data.RData","mouseBrain.DGNBL2.all.data.RData","mouseBrain.DGGRC1.all.data.RData")
data_list <- loadData(file=file,path=path)
region="chr11:60024850-63015319"
gene="Mfap4"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")
region="chr12:26274492-28402363"
gene="Sox11"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")
region="chr11:3902659-6215874"
gene="Rasl10a"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=2,VCscore=0.001,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")


### TEGLU layer

# L23	Nell2	    chr15:94364098-95695661
# L4	Egr1	    chr18:34761289-35036782
# L56	Neto2       chr8:85532065-86659090
### TEGLU layer
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitAuto/"
file = c("mouseBrain.Layer23.all.data.RData","mouseBrain.Layer4.all.data.RData","mouseBrain.Layer56.all.data.RData")
data_list <- loadData(file=file,path=path)

region = "chr15:94364000-95695000"
gene = "Nell2"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = 5720,clean=F,showPETnum=F,outfmt="pdf")

region = "chr18:34761000-35036000"
gene = "Egr1"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = 5720,clean=F,showPETnum=F,outfmt="pdf")

region = "chr8:85532000-86659000"
gene = "Neto2"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = 5720,clean=F,showPETnum=F,outfmt="pdf")


### CBGRC 
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitCellType2/"
file = c("mouseBrain.CBGRC.all.data.RData")
data_list <- loadData(file=file,path=path)

region = "chr6:19403360-23403360"
gene = "Kcnd2"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")


region = "chr18:32770000-33870000"
gene = "Nrep"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

region = "chr8:86902000-87702000"
gene = "Cbln1"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")

region = "chr12:37610000-40280000"
gene = "Etv1"
run(region,file=file,gene=gene,direction =T,showAllPET=T,,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")


## TEGLU
source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.mouseBrainSplitCellType2/"
file = c("mouseBrain.P11_TEGLU.all.data.RData","mouseBrain.P95_TEGLU.all.data.RData","mouseBrain.P365_TEGLU.all.data.RData","mouseBrain.P730_TEGLU.all.data.RData")
data_list <- loadData(file=file,path=path)
downsample(file=file,path=path,sample_num=2230,out_path="/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/TEGLUauto/downsampled/")
path="/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/TEGLUauto/downsampled/"
data_list <- loadData(file=file,path=path)

region = "chr18:71000000-73100000"
gene = "Dcc"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

region = "chr8:22846000-23285000"
gene = "Ank1"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

region = "chr8:23500000-24497000"
gene = "Zmat4"
run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PETDNARNA",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png")

##cut region
## DCC
region = "chr18:60000000-80000000"
cut_region1= c("chr18:70705100-73090908")
cut_region2= c("chr18:66119734-66953436","chr18:63005530-63554160","chr18:74314856-74705099","chr18:75803069-75997855","chr18:76402444-76871617")
run(region,file=file,cut_region1=cut_region1,cut_region2=cut_region2,cut=T,gene=gene,direction =F,showAllPET=F,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png",showBarcode=T,showGex=F,showArc=F,showLoop=F)
## Ank1
region = "chr8:20000000-40000000"
cut_region1= c("chr8:22403547-23246118")
cut_region2= c("chr8:25028825-26997782","chr8:31005464-31295160","chr8:33623060-36567627")
run(region,file=file,cut_region1=cut_region1,cut_region2=cut_region2,cut=T,gene=gene,direction =F,showAllPET=F,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png",showBarcode=T,showGex=F,showArc=F,showLoop=F)
## Zmat4
region = "chr8:20000000-40000000"
cut_region1= c("chr8:23352550-24886917")
cut_region2= c("chr8:27308204-30882482","chr8:31325942-33179600","chr8:36824834-40000000")
run(region,file=file,cut_region1=cut_region1,cut_region2=cut_region2,cut=T,gene=gene,direction =F,showAllPET=F,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="png",showBarcode=T,showGex=F,showArc=F,showLoop=F)



### condition determine this region's barcode character
### type means choose barcode with this type and sort the 500(cell_number) highest type 
### clean means only show gene promoter region's PET
### ORDER: sample - condition - type - clean
countTotalNumber(data_list)

summaryandplot("chr19",9960000,10210000,10460000)

# region1 = "chr19:9960000-10210000" 
# region1 = gsub("-|:"," ",region)
# regionarr = strsplit(trimws(region1)," |\t")
# chrom1 = regionarr[[1]][1]
# start1 = as.numeric(regionarr[[1]][2])
# end1 = as.numeric(regionarr[[1]][3])
# region2 = "chr19:10210000-10460000" 
# region2 = gsub("-|:"," ",region2)
# regionarr = strsplit(trimws(region2)," |\t")
# chrom2 = regionarr[[1]][1]
# start2 = as.numeric(regionarr[[1]][2])
# end2 = as.numeric(regionarr[[1]][3])
# summaryandplot("chr19",9960000,10210000,"chr19",10210000,10460000)


###### plot scatter

### gene count
# source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
# path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByBrainBigClass1Final/"
# file = list.files(path=path, pattern="*.data.RData" ,all.files=FALSE,full.names=FALSE)
# data_list <- loadData(file=file,path=path)
# scATACPETlist_list <- list()
# scATACfragmentlist_list <- list()
# scRNAfragmentlist_list <- list()
# scATACPETlist <- list()
# scATACfragmentlist <- list()
# scRNAfragmentlist <- list()
# for(j in names(data_list[["mouseBrain.allAstroependymal.data.RData"]][["scATACPETlist"]])){
#     for(i in file){
#     scATACPETlist_list[[j]][[i]] <- data_list[[i]][["scATACPETlist"]][[j]]
#     scATACfragmentlist_list[[j]][[i]] <- data_list[[i]][["scATACfragmentlist"]][[j]]
#     scRNAfragmentlist_list[[j]][[i]] <- data_list[[i]][["scRNAfragmentlist"]][[j]]
#     }
#     scATACPETlist[[j]] <- as.data.frame(do.call(rbind,scATACPETlist_list[[j]]))
#     scATACfragmentlist[[j]] <- as.data.frame(do.call(rbind,scATACfragmentlist_list[[j]]))
#     scRNAfragmentlist[[j]] <- as.data.frame(do.call(rbind,scRNAfragmentlist_list[[j]]))
# }

# genelist <- data_list[["mouseBrain.allAstroependymal.data.RData"]][["genelist"]]
# test <- read.table("/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/8.mergeSampleBrain/mouseBrain.RNA.final.anno.counts.txt")
# rowsum <- rowSums(test)

#############

source("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/source.r")
file = c("patski.all.data.RData")
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/9.splitByCellLineFinal/"
data_list <- loadData(file=file,path=path)

library(Signac)
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(clustree)
library(patchwork)
path = "/data/home/ruanlab/huangxingyu/Haoxi20230215/finalresult/finalresult/8.mergeSample"
setwd(path)
set.seed(1234)
SCG = readRDS(file="patski.final.cg.rds")
# SCGS = subset(SCG, subset=(Phase=="S"))
# SCGG1 = subset(SCG, subset=(Phase=="G1"))
# SCGG2M = subset(SCG, subset=(Phase=="G2M"))
test <- as.data.frame(SCG@assays$RNA@counts)
rowsum <- rowSums(test)

#########################
genelist=data_list[[1]][["genelist"]]

scATACPETlist=data_list[[1]][["scATACPETlist"]]
scATACfragmentlist=data_list[[1]][["scATACfragmentlist"]]
scRNAfragmentlist=data_list[[1]][["scRNAfragmentlist"]]

gene = as.data.frame(do.call(rbind,genelist))
geneflt = gene[gene$type == "exon", ]
### gene length
numberofgeneinfo = length(names(table(geneflt$ID)))
namesgeneinfo = names(table(geneflt$ID))
starts = c()
stops  = c()
sizes  = c()
strands = c()
chrom = c()
for (i in (1:numberofgeneinfo)) {
    subgeneinfo  = geneflt[which(geneflt$ID == namesgeneinfo[i]),]
    starts = c(starts,min(subgeneinfo$start,subgeneinfo$end))
    stops  = c(stops, max(subgeneinfo$start,subgeneinfo$end))
    sizes  = c(sizes, stops[i] - starts[i])
    strands = c(strands,subgeneinfo$strand[1])
    chrom = c(chrom,subgeneinfo$chrom[1])
}
transcriptinfo = data.frame(names=namesgeneinfo,chrom=chrom,starts=starts,stops=stops,sizes=sizes,strand=strands)
transcriptinfo = transcriptinfo[order(transcriptinfo$names),]
for(i in 1:nrow(transcriptinfo)){
    transcriptinfo$names_new[i] <- strsplit(transcriptinfo[i,1],"-")[[1]][1]
}
transcriptinfo <- transcriptinfo[order(transcriptinfo$names_new,-transcriptinfo$sizes),]
for(i in nrow(transcriptinfo):2){
    if(transcriptinfo$names_new[i]==transcriptinfo$names_new[i-1] && transcriptinfo$sizes[i]<=transcriptinfo$sizes[i-1]){
        transcriptinfo <- transcriptinfo[-i,]
    }
}
transcriptinfo$count <- 0
count=0
for(i in 1:length(rowsum)){
    if(sum(transcriptinfo$names_new %in% names(rowsum)[i])==1){
        count=count+1
        transcriptinfo[transcriptinfo$names_new %in% names(rowsum)[i],]$count <- rowsum[i]
    }
}
### cal FPKM
totalcount <- sum(transcriptinfo$count)
transcriptinfo$fpkm <- 10^9 * transcriptinfo$count / (transcriptinfo$sizes * totalcount)

### count PETDNARNA number of gene and PET of gene's number
transcriptinfo$tsspet <- 0
transcriptinfo$PETDNARNA <- 0
for(i in 1:nrow(transcriptinfo)){
    if(transcriptinfo$strand[i]=="+"){
        transcriptinfo$pos[i]=transcriptinfo$starts[i]
    }else{
        transcriptinfo$pos[i]=transcriptinfo$stops[i]
    }
    transcriptinfo$tsspet[i] <- sum((scATACPETlist[[transcriptinfo$chrom[i]]]$start1>transcriptinfo$pos[i]-5000 & scATACPETlist[[transcriptinfo$chrom[i]]]$end1<transcriptinfo$pos[i]+5000)|(scATACPETlist[[transcriptinfo$chrom[i]]]$start2>transcriptinfo$pos[i]-5000 & scATACPETlist[[transcriptinfo$chrom[i]]]$end2<transcriptinfo$pos[i]+5000))

    scATACPETlist_tmp <- scATACPETlist[[transcriptinfo$chrom[i]]][(scATACPETlist[[transcriptinfo$chrom[i]]]$start1>transcriptinfo$pos[i]-5000 & scATACPETlist[[transcriptinfo$chrom[i]]]$end1<transcriptinfo$pos[i]+5000)|(scATACPETlist[[transcriptinfo$chrom[i]]]$start2>transcriptinfo$pos[i]-5000 & scATACPETlist[[transcriptinfo$chrom[i]]]$end2<transcriptinfo$pos[i]+5000),]

    scATACfragmentlist_tmp <- scATACfragmentlist[[transcriptinfo$chrom[i]]][(scATACfragmentlist[[transcriptinfo$chrom[i]]]$start>transcriptinfo$pos[i]-5000 & scATACfragmentlist[[transcriptinfo$chrom[i]]]$end<transcriptinfo$pos[i]+5000),]

    scRNAfragmentlist_tmp <- scRNAfragmentlist[[transcriptinfo$chrom[i]]][(scRNAfragmentlist[[transcriptinfo$chrom[i]]]$start>transcriptinfo$starts[i] & scRNAfragmentlist[[transcriptinfo$chrom[i]]]$end<transcriptinfo$stops[i]),]    

    transcriptinfo$PETDNARNA[i] <- length(intersect(scATACPETlist_tmp$barcode,intersect(scATACfragmentlist_tmp$barcode,scRNAfragmentlist_tmp$barcode)))
}

transcriptinfo_fpkm1 <- transcriptinfo[transcriptinfo$fpkm>1,]
transcriptinfo_fpkm1$log10fpkm <- log10(transcriptinfo_fpkm1$fpkm)
library(ggplot2)
library(ggpubr)
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
p <- ggplot(data=transcriptinfo_fpkm1,aes(x=log10(fpkm),y=tsspet)) + geom_hex(bins=300,color = NA) + theme + stat_cor(method="spearman",size=10) + stat_smooth(method="lm",se=TRUE) + ylab("PET number of promoter regions") + xlab(expression(log[10](fpkm))) + theme(legend.position=c(0.9,0.9)) + scale_fill_gradient(low = "#4DBBD5FF", high = "#3C5488FF")
ggsave(filename="/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/tsspet_fpkm_spearman.pdf",p,width=10,height=10)

p <- ggplot(data=transcriptinfo_fpkm1,aes(x=log10(fpkm),y=PETDNARNA)) + geom_hex(bins=300,color = NA) + theme + stat_cor(method="spearman",size=10) + stat_smooth(method="lm",se=TRUE) + ylab("Detected cells") + xlab(expression(log[10](fpkm))) + theme(legend.position=c(0.9,0.9))
ggsave("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/PETDNARNA_fpkm_spearman.pdf",p,width=12,height=10) + scale_fill_gradient(low = "#4DBBD5FF", high = "#3C5488FF")

p <- ggplot(data=transcriptinfo_fpkm1,aes(x=log10(fpkm),y=tsspet)) + geom_hex(bins=300,color = NA) + theme + stat_cor(method="pearson",size=10) + stat_smooth(method="lm",se=TRUE) + ylab("PET number of promoter regions") + xlab(expression(log[10](fpkm))) + theme(legend.position=c(0.9,0.9)) + scale_fill_gradient(low = "#4DBBD5FF", high = "#3C5488FF")             
ggsave("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/tsspet_fpkm_pearson.pdf",p,width=12,height=10)

p <- ggplot(data=transcriptinfo_fpkm1,aes(x=log10(fpkm),y=PETDNARNA)) + geom_hex(bins=300,color = NA) + theme + stat_cor(method="pearson",size=10) + stat_smooth(method="lm",se=TRUE) + ylab("Detected cells") + xlab(expression(log[10](fpkm))) + theme(legend.position=c(0.9,0.9)) + scale_fill_gradient(low = "#4DBBD5FF", high = "#3C5488FF")
ggsave("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/PETDNARNA_fpkm_pearson.pdf",p,width=12,height=10) 

transcriptinfo_fpkm1$group <- 4
transcriptinfo_fpkm1[transcriptinfo_fpkm1$fpkm<=quantile(transcriptinfo_fpkm1$fpkm)[4],]$group <- 3
transcriptinfo_fpkm1[transcriptinfo_fpkm1$fpkm<=quantile(transcriptinfo_fpkm1$fpkm)[3],]$group <- 2
transcriptinfo_fpkm1[transcriptinfo_fpkm1$fpkm<=quantile(transcriptinfo_fpkm1$fpkm)[2],]$group <- 1
transcriptinfo_fpkm1$group <- factor(transcriptinfo_fpkm1$group)

p <- ggplot(transcriptinfo_fpkm1, aes(x = group, y = log10fpkm, fill = group)) +
      geom_violin() +
      geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
      stat_compare_means(comparisons=list(c("1","2"),c("2","3"),c("3","4")),method="wilcox.test",label = "p.format",size=4,step.increase = c(0.05,0.05,0.05)) +
      stat_compare_means(comparisons=list(c("1","2"),c("2","3"),c("3","4")),method="wilcox.test",label = "p.signif",size=8,step.increase = c(0.05,0.05,0.05)) +
      scale_fill_manual(values = c("#91D1C2FF","#4DBBD5FF","#FF6F00FF","#E64B35FF")) +
      labs(title = "", x= "Gene expression (FPKM) group", y = "Gene expression (FPKM)") +
      theme +
      labs(fill="") + 
      theme(legend.position = "top") + 
      ylim(0,max(transcriptinfo_fpkm1$log10fpkm) * 1.2)
ggsave("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/fpkm_boxplot.pdf",p,width=6,height=8)

p <- ggplot(transcriptinfo_fpkm1, aes(x = group, y = tsspet, fill = group)) +
      geom_violin() +
      geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
      stat_compare_means(comparisons=list(c("1","2"),c("2","3"),c("3","4")),method="wilcox.test",label = "p.format",size=4,step.increase = c(0.05,0.05,0.05)) +
      stat_compare_means(comparisons=list(c("1","2"),c("2","3"),c("3","4")),method="wilcox.test",label = "p.signif",size=8,step.increase = c(0.05,0.05,0.05)) +
      scale_fill_manual(values = c("#91D1C2FF","#4DBBD5FF","#FF6F00FF","#E64B35FF")) +
      labs(title = "", x= "Gene expression (FPKM) group", y = "PET number of promoter regions") +
      theme +
      labs(fill="") + 
      theme(legend.position = "top") + 
      ylim(0,max(transcriptinfo_fpkm1$tsspet) * 1.2)
ggsave("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/tsspet_boxplot.pdf",p,width=6,height=8)

p <- ggplot(transcriptinfo_fpkm1, aes(x = group, y = PETDNARNA, fill = group)) +
      geom_violin() +
      geom_boxplot(fill = "white", width = 0.1, outlier.shape = NA) +
      stat_compare_means(comparisons=list(c("1","2"),c("2","3"),c("3","4")),method="wilcox.test",label = "p.format",size=4,step.increase = c(0.05,0.05,0.05)) +
      stat_compare_means(comparisons=list(c("1","2"),c("2","3"),c("3","4")),method="wilcox.test",label = "p.signif",size=8,step.increase = c(0.05,0.05,0.05)) +
      scale_fill_manual(values = c("#91D1C2FF","#4DBBD5FF","#FF6F00FF","#E64B35FF")) +
      labs(title = "", x= "Gene expression (FPKM) group", y = "Detected cells") +
      theme +
      labs(fill="") + 
      theme(legend.position = "top") + 
      ylim(0,max(transcriptinfo_fpkm1$PETDNARNA) * 1.2)
ggsave("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/cells_boxplot.pdf",p,width=6,height=8)

region_gene <- c("Best1","Dagla","Fads1","Fads2","Fads3","Fen1","Fth1","Gm10143","Gm30042","Lrrc10b","Myrf","Rab3il1","Syt7","Tmem258")
p <- ggplot(data=transcriptinfo[transcriptinfo$names_new%in%region_gene,],aes(x=log10(fpkm+1e-5),y=tsspet)) + geom_point(size=8) + theme + stat_cor(method="spearman",size=10) + stat_smooth(method="lm",se=TRUE) + ylab("PET number of promoter regions") + xlab(expression(log[10](fpkm)))          
ggsave(filename="/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/tsspet_fpkm_spearman_regiongene.pdf",p,width=10,height=10)
# p <- ggplot(data=transcriptinfo[transcriptinfo$names_new%in%region_gene,],aes(x=log10(fpkm+1e-5),y=tsspet)) + geom_point() + theme + stat_cor(method="pearson",size=10) + stat_smooth(method="lm",se=TRUE) + ylab("PET number of promoter regions")              
# ggsave("/data/home/ruanlab/huangjiaxiang/scChIATAC_new/001.plot/tsspet_fpkm_pearson_regiongene.png",p,width=12,height=10)