source("source.r")
file = c("k562.all.data.RData")
path = "./"
data_list <- loadData(file=file,path=path)

region = "chr19:8400000-10300000" ##2M
region = "chr19:8670000-9130000" ## 0.5M
region = "chr19:8700000-8800000" ## 0.1M
region = "chr15:61950000-62700000" ## Myc
region = "chr15:61950000-62420000" ## micro Myc
region = "chr19:5783000-5855000" ## Neat1 Malat1
gene = "Neat1"
gene = "Malat1"

run(region,file=file,gene=gene,direction =T,showAllPET=T,cell_number=800,PETcount=0,VCscore=0,type="PET",condition="PETDNARNA",sample_num = -1,clean=F,showPETnum=F,outfmt="pdf")
