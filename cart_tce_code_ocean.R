
library(bannerCommenter)
library(stringr)
library(dplyr)
library(plyr)
library(deconstructSigs)
library(dndscv)
library(mmsig)
library(plotrix)
library(igraph)
library(quantsmooth)
library(dplyr)
library(stringr)
require(plyr)
library('GenomicFeatures')
library(IRanges)
library("GenomicRanges")
library(GenomicRanges)
library(IRanges)
library(mmsig)
library(pheatmap)
library(deconstructSigs)
library(survminer)
library(survival)
library(reshape2)

## The analysis was done on hg38 genome
library(BSgenome.Hsapiens.UCSC.hg38)

## This represents the theme we use for the CCF plot (panel 3)
gg_th <- (panel_theme = theme_bw() + theme(
  panel.border = element_blank(),
  text = element_text(color = "black"),
  axis.text = element_text(color = "black"),
  legend.position = c(0.3,1.15), #'top', #need to increase top margin 
  legend.key.size = unit(5, 'mm'),
  legend.title = element_blank(),
  legend.direction = 'horizontal',
  legend.text = element_text(size = 8),
  #legend.justification = c(0, -5),
  legend.spacing.x = unit(0.01, "mm"),
  plot.subtitle = element_text(hjust = 0.5, size = 8),
  plot.title = element_text(face = 'bold', hjust = 0, vjust = -2, size = 12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  axis.text.x = element_text(size = 10),
  axis.title = element_text(size = 10),
  axis.line = element_line(),
  plot.margin = unit(c(1,0,0,0.1), 'cm'))) 

##################################################################
##                                                              ##
##            Title: Timing antigenic escape analysis           ##
##                                                              ##
##################################################################

## we set the working directory to the folder were all the test files are located
wd = "/Users/mxp2050/Desktop/cart_tce_code_ocean/"
setwd(wd)

## This test code is prepared for sample MM-15 (P3298, P2510)
samples = c("P3298", "P2510")

##############
############## adding files into the environment
##############

##### snv
snv = read.delim("data/all_snv_test.txt")
unique(snv$SAMPLE)
# [1] "P3298" "P2510"

## keeping high confidence calls
snv = snv %>%
  dplyr::filter(HighConfidence == "True")
## adding a code column for downstream analysis
snv$code = paste(snv$CHROM, snv$POS)

snv1 = snv[snv$SAMPLE %in% samples,]

## we keep only single nucleotide variants 
snv1 = snv1 %>%
  dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1)

table(table(snv1$code))
# 
# 1     2 
# 1236 12673 

## this shows that 12673 SNVs are shared by both samples

##### cnv
cnv = read.delim("data/cnv_test.txt")
unique(cnv$sample)
# [1] "P3298" "P2510"

cnv1 = cnv[cnv$sample %in% samples,]

######################################################################
###
### #15 - VN - P2510-P3298
###
#####################################################################

all_snv_branch = snv1

## adding the dp output into the environment
MUT_PTS = read.delim("data/VN_2000iters_1000burnin_bestConsensusAssignments_test.bed")
MUT_PTS$code = paste(MUT_PTS$chr, MUT_PTS$end)
unique(MUT_PTS$cluster)
# [1]  2  1  4  6  3 NA

## merge the snv file with the dp output file
all_snv_branch_dp = join(all_snv_branch,MUT_PTS, by = "code" )
all_snv_branch_dp2<- unique(all_snv_branch_dp)

## some mutations have NA values for likelihood belonging into any cluster
dim(all_snv_branch_dp[!is.na(all_snv_branch_dp$likelihood),])
table(all_snv_branch_dp2$cluster, all_snv_branch_dp2$SAMPLE)
# 
# P2510 P3298
# 1   129   144
# 2 11731 12173
# 3   495   490
# 4    15   204
# 6   127    87

clusters = read.delim("data/VN_2000iters_1000burnin_bestClusterInfo_test.txt")
colnames(clusters)
# [1] "cluster.no"      "VNP2510T"        "VNP3298T"        "no.of.mutations"

## changing the names of the 2 columns to match the samples we use
colnames(clusters)[2:3] = c("P2510", "P3298")

## this is to check the frequency of each mutation per each cluster
kk<- as.data.frame(table(all_snv_branch_dp2$code, all_snv_branch_dp2$cluster))
colnames(kk)[1]<-"code"
all_snv_branch_dp22<- unique(merge(all_snv_branch_dp2, kk, by ="code"))
table(all_snv_branch_dp22$Freq, all_snv_branch_dp22$cluster)

missing<- all_snv_branch_dp[is.na(all_snv_branch_dp$likelihood),]

all_snv_branch_dp2<- unique(all_snv_branch_dp[!is.na(all_snv_branch_dp$cluster),])

## vaf boxplot per each cluster
boxplot(all_snv_branch_dp2$t_VAF~all_snv_branch_dp2$cluster)
# cluster per each chromosome
unique(all_snv_branch_dp2$chr)
# [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18"
# [19] "19" "20" "21" "22" "X" 
all_snv_branch_dp2$chr = factor(all_snv_branch_dp2$chr,
                                levels = c(1:22, "X"))
barplot(table(all_snv_branch_dp2$cluster, all_snv_branch_dp2$chr), las=2)

## we collapse clusters 3 and 2 to 1 (this is done manually based on prior biological knowledge) --> should be individualized per patient
variants<- all_snv_branch_dp2 %>% dplyr::mutate(cluster = ifelse(cluster== 2 | cluster == 3,1, cluster))

## this represents our newly formed clusters based on the prior collapse
tab=as.data.frame(table(variants$cluster))
tab
# Var1  Freq
# 1    1 25162
# 2    4   219
# 3    6   214

##############
############## FIRST PANEL: TREE
##############

##TREE
clu<- read.delim("data/VN_2000iters_1000burnin_bestClusterInfo_test.txt")
cluname<- "VN"
rownames(clu) = clu$cluster.no

clu<- clu[,c(2,3)] 

## we want to collapse clusters 3 and 2 to 1
clu<- clu[-c(2,3),]
clu[1,1]<- 1
clu[1,2]<- 1
clu
# VNP2510T    VNP3298T
# 1  1.00000 1.000000000
# 4  0.00000 0.985693898
# 6  0.15094 0.001589591

MEDIANS.1<- sapply(clu, function(x) x * 100)
rownames(MEDIANS.1) = row.names(clu)

## add cluster 6
tab=merge(tab,MEDIANS.1,by.x="Var1",by.y="row.names")

## we ascribe if its part of the trunk, part of the PD or the BAS sample based on the CCF for each cluster
tab$Cluster=ifelse(tab$VNP2510T==100 & tab$VNP3298T==100,"Trunk",
                   ifelse(tab$VNP2510T>tab$VNP3298T,"BAS",
                          ifelse(tab$VNP2510T<tab$VNP3298T,"PD",NA)))
# tab$Color=COLORS[1:nrow(tab)] # come back to this for color options

## we assign the colors we want to use for our 3 clusters
tab$Color = ifelse(tab$Cluster == "Trunk",
                   yes = "grey",
                   ifelse(tab$Cluster == "PD",
                          yes = "blue",
                          ifelse(tab$Cluster == "BAS",
                                 yes = "yellow",
                                 no = NA)))

CLUSTERS=unique(tab$Cluster)
CLUSTERS=CLUSTERS[CLUSTERS!="Trunk"]
SUM=NULL
for(j in 1:length(CLUSTERS)){
  S=sum(tab[tab$Cluster==CLUSTERS[j],"Freq"])
  SUM=rbind(SUM,c(CLUSTERS[j],S))
}

SUM=as.data.frame(SUM)
colnames(SUM)=c("Cluster","Sum")

## we need to manually add the cluster 6 from the split
SUM$Sum=as.numeric(as.character(SUM$Sum))
if(SUM[SUM$Cluster=="BAS","Sum"]>SUM[SUM$Cluster=="PD","Sum"]){
  VAL=as.numeric(as.character(SUM[SUM$Cluster=="BAS","Sum"]))
} else if(SUM[SUM$Cluster=="BAS","Sum"]<SUM[SUM$Cluster=="PD","Sum"]){
  VAL=as.numeric(as.character(SUM[SUM$Cluster=="PD","Sum"]))
} else if(SUM[SUM$Cluster=="BAS","Sum"]==SUM[SUM$Cluster=="PD","Sum"]){
  VAL=as.numeric(as.character(unique(SUM$Sum)))
} else{
  print("Boo")
}
SUM
# Cluster Sum
# 1      PD 219
# 2     BAS 214

yMAX=round(tab[tab$Cluster=="Trunk","Freq"]+VAL+((tab[tab$Cluster=="Trunk","Freq"]+VAL)/10))

pdf(paste(wd, "output/", cluname, ".tree.pdf",
          sep=""),
    width=4,height=5) #(4x5)
{
  {
    plot(1.5,
         type = "n", 
         frame = FALSE, 
         xlab = "", 
         ylab = "Number of Variants",
         main= cluname, #unique(all_snv_branch_dp2$sampleID)[1],
         las=1,
         xaxt="n",
         xlim=c(0.9,2.1),
         ylim=c(yMAX,0), #c(yMAX,0),
         lwd=2, cex.axis = 0.8, cex.lab = 0.8)
    lines(c(1.5,1.5),c(0,tab[tab$Cluster=="Trunk","Freq"]), 
          pch = 19, 
          col = c("grey70", tab[tab$Cluster=="Trunk","Color"]),  
          type = "o",
          lwd=10)
    points(c(1.5),0,col="grey70",lwd=10,pch=19,cex =1) #bg = "grey70",cex =2) 
  }
  {
    ADD.p=purrr::accumulate(rep(0.1, 20), ~.x/2)
    ADD.p=ADD.p-0.1
    ADD.p=ADD.p*(-1)
    ADD.p=ADD.p+1.5
    #ADD.p[1]=ADD.p[1]-0.1
    ADD.p=c(0,ADD.p)
    ADD.n=purrr::accumulate(rep(0.1, 20), ~.x/2)
    ADD.n=ADD.n-0.1
    ADD.n=ADD.n*(-1)
    ADD.n=1.5-ADD.n
    ADD.n=c(0,ADD.n)
    GROUP=unique(tab$Cluster)
    GROUP=GROUP[GROUP!="Trunk"]
    GROUP=GROUP[order(GROUP)]
  }
  for(zz in 1:length(GROUP)){
    tab.sub=tab[tab$Cluster==GROUP[zz],]
    rownames(tab.sub)=1:length(tab.sub$Var1) 
    ROWS=(1:nrow(tab.sub))
    for(zzz in (1:length(ROWS))){
      tab.sub.1=tab.sub[rownames(tab.sub)==c(1:zzz),]
      if(length(tab.sub.1$Var1)==1){
        if(tab.sub.1$Cluster=="BAS"){
          lines(c((ADD.n[zzz+1]),ADD.n[zzz+2]),
                c(tab[tab$Cluster=="Trunk","Freq"],
                  (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))), 
                pch = 19, 
                col = tab.sub$Color[zzz], 
                type = "o",
                lwd=10)
          points(c(ADD.n[zzz+2]),
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq)),
                 col=tab[tab$Cluster=="BAS","Color"], #"grey70",
                 lwd=10,
                 pch=19) 
        } 
        else if(tab.sub.1$Cluster=="PD"){
          lines(c((ADD.p[zzz+1]),ADD.p[zzz+2]),
                c(tab[tab$Cluster=="Trunk","Freq"],
                  (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))), 
                pch = 19, 
                col = tab.sub$Color[zzz], 
                type = "o",
                lwd=10)
          points(c(ADD.p[zzz+2]),
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq)),
                 col=tab[tab$Cluster=="PD","Color"],#"grey70",
                 lwd=10,
                 pch=19)
          
        } #added PD2 seg
        else if(tab.sub.1$Cluster=="PD2"){
          lines(c((ADD.p[zzz+1]),ADD.p[zzz+6]),
                c(tab[tab$Cluster=="Trunk","Freq"],
                  (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))), 
                pch = 19, 
                col = tab.sub$Color[zzz], 
                type = "o",
                lwd=10)
          points(c(ADD.p[zzz+6]),
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq)),
                 col=tab[tab$Cluster=="PD2","Color"],#"grey70",
                 lwd=10,
                 pch=19)
        }
        else{
          print("Bo!!!")
        }
      } 
      else if(length(tab.sub.1$Var1)==2){
        if(tab.sub.1$Cluster=="BAS"){
          lines(c((ADD.n[zzz+1]),ADD.n[zzz+2]),
                c((tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))), 
                pch = 19, 
                col = tab.sub$Color[zzz], #"grey70", 
                type = "o",
                lwd=10)
          #  points(c(ADD.p[zzz+1]),
          #        (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-1],
          #       col= tab[tab$Cluster=="BAS","Color"][1],#"grey70",
          #      lwd=10,
          #     pch=19)
          #  points(c(ADD.p[zzz+2]),
          #        (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz],
          #       col=tab[tab$Cluster=="BAS","Color"][2],#"grey70",
          #      lwd=10,
          #     pch=19) 
          points(c(ADD.n[zzz+1]), ###added this
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-1], 
                 col=tab[tab$Cluster=="BAS","Color"][1],#"grey70",
                 lwd=10,
                 pch=19)
        } 
        else if(tab.sub.1$Cluster=="PD"){
          lines(c((ADD.p[zzz+1]),ADD.p[zzz+2]),
                c((tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))), 
                pch = 19, 
                col = tab.sub$Color[zzz], #"grey70",
                type = "o",
                lwd=10)
          points(c(ADD.p[zzz+1]),
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-1],
                 col=tab[tab$Cluster=="PD","Color"],#"grey70",
                 lwd=10,
                 pch=19)
          points(c(ADD.p[zzz+2]),
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz],
                 col=tab[tab$Cluster=="PD","Color"][2], #"grey70",
                 lwd=10,
                 pch=19)
        } #added PD2
        else if(tab.sub.1$Cluster=="PD2"){
          lines(c((ADD.p[zzz+6]),ADD.p[zzz+10]),
                c((tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))), 
                pch = 19, 
                col = tab.sub$Color[zzz], #"grey70",
                type = "o",
                lwd=10)
          points(c(ADD.p[zzz+6]),
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-1],
                 col=tab[tab$Cluster=="PD2","Color"],#"grey70",
                 lwd=10,
                 pch=19)
          points(c(ADD.p[zzz+10]),
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz],
                 col=tab[tab$Cluster=="PD2","Color"][2], #"grey70",
                 lwd=10,
                 pch=19)
        }
        else{
          print("Bo!!!")
        }
      } else if(length(tab.sub.1$Var1)>2){
        if(tab.sub.1$Cluster=="BAS"){
          lines(c((ADD.n[zzz]),ADD.n[zzz+1], ADD.n[zzz+2]), 
                c((tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))), 
                pch = 19, 
                col = tab.sub$Color[zzz], #"grey70",
                type = "o",
                lwd=10)
          #  points(c(ADD.p[zzz+1]),
          #        (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-1],
          #       col= tab[tab$Cluster=="BAS","Color"][1],#"grey70",
          #      lwd=10,
          #     pch=19)
          #  points(c(ADD.p[zzz+2]),
          #        (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz],
          #       col=tab[tab$Cluster=="BAS","Color"][2],#"grey70",
          #      lwd=10,
          #     pch=19) 
          points(c(ADD.n[zzz]), ###added this
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-2], 
                 col=tab[tab$Cluster=="BAS","Color"][1],#"grey70",
                 lwd=10,
                 pch=19)
          points(c(ADD.n[zzz+1]), ###added this
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-1], ####here need to fix the position
                 col=tab[tab$Cluster=="BAS","Color"][2],#"grey70",
                 lwd=10,
                 pch=19)
          points(c(ADD.n[zzz+2]), ###added this
                 (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz], ####here need to fix the position
                 col=tab[tab$Cluster=="BAS","Color"][3],#"grey70",
                 lwd=10,
                 pch=19)
        } 
        else if(length(tab.sub.1$Var1)>2){
          if(tab.sub.1$Cluster=="PD"){
            lines(c((ADD.p[zzz]),ADD.p[zzz+1], ADD.p[zzz+2]), 
                  c((tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))), 
                  pch = 19, 
                  col = tab.sub$Color[zzz], #"grey70",
                  type = "o",
                  lwd=10)
            #  points(c(ADD.p[zzz+1]),
            #        (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-1],
            #       col= tab[tab$Cluster=="BAS","Color"][1],#"grey70",
            #      lwd=10,
            #     pch=19)
            #  points(c(ADD.p[zzz+2]),
            #        (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz],
            #       col=tab[tab$Cluster=="BAS","Color"][2],#"grey70",
            #      lwd=10,
            #     pch=19) 
            points(c(ADD.p[zzz]), ###added this
                   (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-2], ####here need to fix the position
                   col=tab[tab$Cluster=="PD","Color"][1],#"grey70",
                   lwd=10,
                   pch=19)
            points(c(ADD.p[zzz+1]), ###added this
                   (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz-1], ####here need to fix the position
                   col=tab[tab$Cluster=="PD","Color"][2],#"grey70",
                   lwd=10,
                   pch=19)
            points(c(ADD.p[zzz+2]), ###added this
                   (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab.sub.1$Freq))[zzz], ####here need to fix the position
                   col=tab[tab$Cluster=="PD","Color"][3],#"grey70",
                   lwd=10,
                   pch=19)
          } }
        
        
        else{
          print("Bo!!!")
        }
      }
      
      else{
        print("Boo!!")
      }
      
    }
  }
  
  {
    points(c(1.5),tab[tab$Cluster=="Trunk","Freq"],col=tab[tab$Cluster=="Trunk","Color"],lwd=10,pch=19)
    #  text((ADD.n[nrow(tab[tab$Cluster=="BAS",])+1]),
    #      (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab[tab$Cluster=="BAS","Freq"])+
    #        (tab[tab$Cluster=="Trunk","Freq"]+cumsum(tab[tab$Cluster=="BAS","Freq"]))/10),
    #    label="BAS",cex=1)
    # text((ADD.p[nrow(tab[tab$Cluster=="PD",])+1]),
    #     (tab[tab$Cluster=="Trunk","Freq"]+sum(tab[tab$Cluster=="PD","Freq"])+
    #       (tab[tab$Cluster=="Trunk","Freq"]+sum(tab[tab$Cluster=="PD","Freq"]))/10),
    #     label="PD",cex=1)
  }
}
dev.off()

## We can manually decrease afterwards the y-axis for the mutations in the trunk.

##############
############## SECOND PANEL: MUTATIONAL SIGNATURES
##############

## Mutational signatures were generated per each cluster using mmsig fitting SBS1, SBS2, SBS5, SBS8, SBS9, SBS13, SBS18 and SBS99(SBS-MM1) in samples exposed to melphalan
## according to the work of EH Rustad et al. 2021

mat_sig2<-read.delim("data/sbs_test.txt")
colnames(mat_sig2)[9] = "SBS99"

df = mat_sig2

# Colors
colors=c("SBS1"="#4D9B75",
         "SBS2" = "#D95F02",
         "SBS5"="#7570B3",
         "SBS8"="#E7298A",
         "SBS9"="#66A61E",
         "SBS13" = "#E6AB02",
         "SBS18"="#A6761D",
         "SBS99" = "grey")

# Combine SBS columns into a matrix for barplot
mat <- as.matrix(df[, c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS18", "SBS99")])
rownames(mat) = df$sample

# Set up a 1x1 plotting area
par(mfrow=c(1,1))

pdf(paste(wd, "output/", cluname, ".sigs.pdf",
          sep=""),
    width=4,height=5) #(4x5)

barplot(t(mat), beside = FALSE, col = colors, ylim = c(0, 1), las = 1,
        ylab = "Proportion", width = rep(1, nrow(mat)), border = "black")
dev.off()


##############
############## THIRD PANEL: CCF PLOT
##############

meltab<- tab
meltab$Var1<- paste("Cluster", tab$Var1)
meltab<- melt(meltab[,c(1,3,4,6 )])
meltab$variable<- as.character(meltab$variable)

meltab<- meltab %>%
  mutate(variable = fct_reorder(variable, Var1))

ggplot(data = meltab, aes(x = variable, y = value , fill = Var1, group = Var1)) + geom_line(aes(color = Var1), size = 2) + geom_point(size = 4, shape = 21) +
  scale_color_manual(guide = "none", values =  meltab$Color[1:(length(meltab$Color)/2)], labels = meltab$Var1[1:(length(meltab$Var1)/2)]) + 
  scale_fill_manual(values = meltab$Color[1:(length(meltab$Color)/2)], labels = meltab$Var1[1:(length(meltab$Var1)/2)]) + 
  labs(fill = "Cluster", y = "Cancer Cell Fraction", x = NULL) +gg_th   + theme(axis.text.x=element_text(angle=90,hjust=1))

ggsave(paste(wd , "output/",cluname, ".CCF.pdf",
             sep=""),
       plot = last_plot(), height = 2.7, width = 2.2, units = "in")  

library(readr)
snv = read_tsv("/Users/mxp2050/Library/CloudStorage/Box-Box/CoMMpass IA21/Final Files/6262ea84d1211200016ce44c_Somatic Mutation Files - SNV and INDEL_MMRF_CoMMpass_IA21_combined_vcfmerger2_All_Canonical_Variants.tsv")

