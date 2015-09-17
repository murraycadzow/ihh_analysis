get_rsb = function(pop1, pop2, chr ){
  if(pop1 != pop2){
    combined= merge(x=read.table(paste0(pop1,"chr",chr,".ihh")), y=read.table(paste0(pop2,"chr",chr,".ihh")), by="POSITION")

    a= ies2rsb(combined[,c("CHR.x","POSITION","FREQ_a.x","IHHa.x","IHHd.x","IES.x")], combined[,c("CHR.y","POSITION","FREQ_a.y","IHHa.y","IHHd.y","IES.y")],
               popname1=pop1, popname2=pop2, method="bilateral")
    b=a$res.rsb[!(is.na(a$res.rsb[3])),]
    b$POP1= as.factor(pop1)
    b$POP2= as.factor(pop2)
    b$Test = as.factor("Rsb")
    names(b)=c("CHR","POSITION","Std_iHH","Pvalue", "POP1","POP2","Test")
  } else{
    b= get_ihs(pop1, chr)
  }
  return(b)

}

get_ihs = function(pop1,chr){
  a = ihh2ihs(read.table(paste0(pop1,"chr",chr,".ihh")))
  b=a$res.ihs[!(is.na(a$res.ihs[3])),]
  b$POP1= as.factor(pop1)
  b$POP2= as.factor(pop1)
  b$Test = as.factor("iHS")
  names(b)= c("CHR","POSITION","Std_iHH","Pvalue","POP1","POP2","Test")
  return(b)
}


df=data.table::data.table()
for(i in 1:22){
  df=rbindlist(list(df,read_delim(file=paste0("/run/user//1001/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/extract/fixed_vcf/snpEff_classic/rs_genes",i,".txt"), delim="\t", col_names=c("CHR","start","POSITION","rs_ids","gene"))) )
}





for( i in 1:22){
  pop1 <- "AXIOM"
  rsb=data.frame()
  for( pop2 in c("AXIOM","OMNI","CEU","CHB","CHS","GBR","YRI")){
    print(pop2)
    if(file.exists(file = paste0(pop1,"chr",i,".ihh")) & file.exists(file = paste0(pop2,"chr",i,".ihh"))){
      rsb=rbind(rsb, get_rsb(pop1=pop1, pop2=pop2, chr=i))

    }
  }

  rsb <-rsb %>% group_by(POP2) %>% mutate(Test_rank= rank(abs(Std_iHH), ties.method="random"))
  rsb <- rsb %>% group_by(POP2) %>% arrange(Test_rank)  %>% mutate(new_rank = order(Test_rank, decreasing = TRUE))
  rsb= data.frame( chrom=rsb$CHR, chrom_start=rsb$POSITION-1, chrom_end=rsb$POSITION, Std_iHH=rsb$Std_iHH, neglogPvalue=rsb$Pvalue, Std_iHH_rank=rsb$new_rank, POP1 = rsb$POP1, POP2=rsb$POP2, Test=rsb$Test)

  write.table(rsb, file=paste0("~/",pop1,i,"_rsb.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}


for( i in 1:22){
  pop1 <- "OMNI"
  rsb=data.frame()
  for( pop2 in c("AXIOM","OMNI","CEU","CHB","CHS","GBR","YRI")){
    print(pop2)
    if(file.exists(file = paste0(pop1,"chr",i,".ihh")) & file.exists(file = paste0(pop2,"chr",i,".ihh"))){
      rsb=rbind(rsb, get_rsb(pop1=pop1, pop2=pop2, chr=i))

    }
  }
  rsb<-rsb %>% group_by(POP2) %>% mutate(Test_rank= rank(abs(Std_iHH), ties.method="random"))
  rsb <- rsb %>% group_by(POP2) %>% arrange(Test_rank)  %>% mutate(new_rank = order(Test_rank, decreasing = TRUE))
  rsb= data.frame( chrom=rsb$CHR, chrom_start=rsb$POSITION-1, chrom_end=rsb$POSITION,Std_iHH=rsb$Std_iHH, neglogPvalue=rsb$Pvalue, Std_iHH_rank=rsb$new_rank, POP1 = rsb$POP1, POP2=rsb$POP2, Test= rsb$Test)

  write.table(rsb, file=paste0("~/",pop1,i,"_rsb.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}


library(RMySQL)
drv = dbDriver("MySQL")
db = dbConnect(drv, user="murray", host="127.0.0.1", dbname="selection")

dbGetQuery(db, "CREATE TABLE `axiom_rsb`(`chrom` int(31),`chrom_start` int(10),`chrom_end` int(10), `Std_iHH` float, `neglogPvalue` float, `Std_iHH_rank` int(10),`Population1` varchar(20), `Population2` varchar(20), `Test` varchar(20) );")
for(i in 1:22){
dbGetQuery(db,paste0("load data infile '/home/murraycadzow/AXIOM",i,"_rsb.txt' into table axiom_rsb fields terminated by '\t' lines terminated by '\n' ignore 1 rows;"))
}

dbGetQuery(db, "CREATE TABLE `omni_rsb`(`chrom` int(31),`chrom_start` int(10),`chrom_end` int(10), `Std_iHH` float, `neglogPvalue` float, `Std_iHH_rank` int(10),`Population1` varchar(20), `Population2` varchar(20), `Test` varchar(20) );")
for(i in 1:22){
  dbGetQuery(db,paste0("load data infile '/home/murraycadzow/OMNI",i,"_rsb.txt' into table omni_rsb fields terminated by '\t' lines terminated by '\n' ignore 1 rows;"))
}



####
#### iHS
####

library(dplyr)
library(rehh)
library(qqman)
library(readr)
library(data.table)
read_ihs=function(POP){
  if(!file.exists(paste0(POP,"_ihs.RData"))){
    ihh=data.frame()
    for( i in 1:22){
      if(file.exists(file = paste0(POP,"chr",i,".ihh"))){
        ihh=rbind(ihh,read.table(file = paste0(POP,"chr",i,".ihh"), header=TRUE))
      } else {
        print(paste0(POP,"chr",i,".ihh is missing"))
      }
    }
    ihs = ihh2ihs(ihh)
    ihs = ihs$res.ihs
    ihs = ihs[!is.na(ihs$iHS),]
    write.table(ihs, file=paste0(POP,"_ihs.txt"), row.names=FALSE, quote=FALSE,col.names = TRUE, sep="\t")
  }else{
    ihs = read_delim(file=paste0(POP,"_ihs.txt"), col_names=TRUE, delim="\t")
  }
  return(ihs)
}

#pops AXIOM OMNI CEU GBR CHB CHS YRI
a<-read_ihs(POP="YRI")
a= data.frame(marker=rownames(a), chrom=a$CHR, chrom_start=a$POSITION-1, chrom_end=a$POSITION, iHS=a$iHS, Population = "YRI_1KGP", neglogPvalue=a$Pvalue)
a <-a %>% group_by(chrom) %>% mutate(rank= rank(abs(iHS), ties.method="random"))
a <- a %>% group_by(chrom) %>% arrange(rank)  %>% mutate(new_rank = order(rank, decreasing = TRUE))
write.table(a[,c("chrom", "chrom_start", "chrom_end", "marker", "iHS", "new_rank", "neglogPvalue","Population" )], file="~/yri_ihs.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
rm(a)

##mysql
#CREATE TABLE `omni_ihs`(`chrom` int(31),`chrom_start` int(10),`chrom_end` int(10),`marker` varchar(40), `iHS` float, `iHS_rank` int(10),`neglogPvalue` float,`Population` varchar(20) );

#load data infile '/home/murraycadzow/axiom_ihs.txt' into table axiom_ihs fields terminated by '\t' lines terminated by "\n" ignore 1 rows;
