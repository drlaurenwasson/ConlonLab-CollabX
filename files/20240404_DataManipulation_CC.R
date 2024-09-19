library(shiny)
library(readxl)
library(ggplot2)
library(plyr)
library(tidyverse)
library(ggpattern)

#Load in the table for collaborative cross
s1c <- read_excel("NIHMS1806954-supplement-Table_S1__Differentially_expressed_genes_or_proteins_in_specific_experiments_or_analyses___S1A-P_.xlsx", sheet = "S1C")
View(s1c)
genes<- make.names(s1c$Gene, unique = TRUE)
s1c<- s1c[,24:59]
rownames(s1c)<- genes

f6<- s1c[,1:9]
f6norm<- (f6)/f6$`F6_ REF`

f7<- s1c[,10:18]
f7norm<- (f7)/f7$`F7_ REF`

f9<- s1c[,19:27]
f9norm<- (f9)/f9$`F9_ REF`

f10<- s1c[,28:36]
f10norm<- (f10)/f10$`F10_ REF`

#PWK
pwk<- f6norm[,c(2,8)]
rownames(pwk)<- genes
a<- f7norm[,2]
pwk<- cbind(pwk,a)
colnames(pwk)<- c("F6_M_PWK", "F6_F_PWK", "F7_M_PWK")
a<- f10norm[,3]
pwk<- cbind(pwk,a)
colnames(pwk)<- c("F6_M_PWK", "F6_F_PWK", "F7_M_PWK", "F10_F_PWK")
pwk2<-pwk
pwk2$M_PWK_avg<- rowMeans(subset(pwk, select = c(F6_M_PWK, F7_M_PWK)), na.rm = TRUE)
pwk2$F_PWK_avg<- rowMeans(subset(pwk, select = c(F6_F_PWK, F10_F_PWK)), na.rm = TRUE)

pwk2$FCMvsF_PWK<- round((pwk2$M_PWK_avg/pwk2$F_PWK_avg), digits = 2)
pwk2$Log2MvsF_PWK<- round(log2(pwk2$M_PWK_avg/pwk2$F_PWK_avg), digits = 2)

#C57
c57<- as.data.frame(f6norm[,3])
rownames(c57)<- genes
a<- f9norm[,c(3,9)]
c57<- cbind(c57,a)
colnames(c57)<- c("F6_F_C57", "F9_F_C57", "F9_M_C57")
a<- f10norm[,9]
c57<- cbind(c57,a)
colnames(c57)<- c("F6_F_C57", "F9_F_C57", "F9_M_C57", "F10_M_C57")
c572<- c57
c572$M_C57_avg<- rowMeans(subset(c57, select = c(F9_M_C57, F10_M_C57)), na.rm = TRUE)
c572$F_C57_avg<- rowMeans(subset(c57, select = c(F6_F_C57, F9_F_C57)), na.rm = TRUE)
c572$FCMvsF_c57<- round((c572$M_C57_avg/c572$F_C57_avg), digits = 2)
c572$Log2MvsF_c57<- round(log2(c572$M_C57_avg/c572$F_C57_avg), digits = 2)

#A/J
aj<- as.data.frame(f6norm[,4])
rownames(aj)<- genes
a<- f7norm[,c(9)]
aj<- cbind(aj,a)
a<- f9norm[,c(2,7)]
aj<- cbind(aj,a)
colnames(aj)<- c("F6_F_AJ", "F7_M_AJ", "F9_F_AJ", "F9_M_AJ")
aj2<- aj
aj2$M_aj_avg<- rowMeans(subset(aj, select = c(F7_M_AJ, F9_M_AJ)), na.rm = TRUE)
aj2$F_aj_avg<- rowMeans(subset(aj, select = c(F6_F_AJ, F9_F_AJ)), na.rm = TRUE)
aj2$FCMvsF_aj<- round((aj2$M_aj_avg/aj2$F_aj_avg), digits = 2)
aj2$Log2MvsF_aj<- round(log2(aj2$M_aj_avg/aj2$F_aj_avg), digits = 2)

#129
s129<- as.data.frame(f6norm[,5])
a<- f7norm[,c(3,7)]
s129<- cbind(s129,a)
rownames(s129)<- genes
a<- f10norm[,6]
s129<- cbind(s129,a)
colnames(s129)<- c("F6_M_129", "F7_M_129", "F7_F_129", "F10_F_129")
s1292<- s129
s1292$M_s129_avg<- rowMeans(subset(s129, select = c(F6_M_129, F7_M_129)), na.rm = TRUE)
s1292$F_s129_avg<- rowMeans(subset(s129, select = c(F7_F_129, F10_F_129)), na.rm = TRUE)
s1292$FCMvsF_s129<- round((s1292$M_s129_avg/s1292$F_s129_avg), digits = 2)
s1292$Log2MvsF_s129<- round(log2(s1292$M_s129_avg/s1292$F_s129_avg), digits = 2)


#WSB
wsb<- as.data.frame(f6norm[,6])
rownames(wsb)<- genes
a<- f9norm[,c(8)]
wsb<- cbind(wsb,a)
a<- f10norm[,c(4,8)]
wsb<- cbind(wsb,a)
colnames(wsb)<- c("F6_F_WSB", "F9_M_WSB", "F10_F_WSB", "F10_M_WSB")
wsb2<- wsb
wsb2$M_wsb_avg<- rowMeans(subset(wsb, select = c(F9_M_WSB, F10_M_WSB)), na.rm = TRUE)
wsb2$F_wsb_avg<- rowMeans(subset(wsb, select = c(F6_F_WSB, F10_F_WSB)), na.rm = TRUE)
wsb2$FCMvsF_wsb<- round((wsb2$M_wsb_avg/wsb2$F_wsb_avg), digits = 2)
wsb2$Log2MvsF_wsb<- round(log2(wsb2$M_wsb_avg/wsb2$F_wsb_avg), digits = 2)

#CAS
cas<- as.data.frame(f6norm[,7])
rownames(cas)<- genes
a<- f7norm[,c(6)]
cas<- cbind(cas,a)
a<- f10norm[,c(2,5)]
cas<- cbind(cas,a)
colnames(cas)<- c("F6_F_CAS", "F7_M_CAS", "F10_F_CAS", "F10_M_CAS")
cas2<- cas
cas2$M_cas_avg<- rowMeans(subset(cas, select = c(F7_M_CAS, F10_M_CAS)), na.rm = TRUE)
cas2$F_cas_avg<- rowMeans(subset(cas, select = c(F6_F_CAS, F10_F_CAS)), na.rm = TRUE)
cas2$FCMvsF_cas<- round((cas2$M_cas_avg/cas2$F_cas_avg), digits = 2)
cas2$Log2MvsF_cas<- round(log2(cas2$M_cas_avg/cas2$F_cas_avg), digits = 2)

#NZO
nzo<- as.data.frame(f6norm[,9])
rownames(nzo)<- genes
a<- f7norm[,c(8)]
nzo<- cbind(nzo,a)
a<- f9norm[,c(4)]
nzo<- cbind(nzo,a)
a<- f10norm[,7]
nzo<- cbind(nzo,a)
colnames(nzo)<- c("F6_M_NZO", "F7_F_NZO", "F9_M_NZO", "F10_F_NZO")
nzo2<- nzo
nzo2$M_nzo_avg<- rowMeans(subset(nzo, select = c(F6_M_NZO, F9_M_NZO)), na.rm = TRUE)
nzo2$F_nzo_avg<- rowMeans(subset(nzo, select = c(F7_F_NZO, F10_F_NZO)), na.rm = TRUE)
nzo2$FCMvsF_nzo<- round((nzo2$M_nzo_avg/nzo2$F_nzo_avg), digits = 2)
nzo2$Log2MvsF_nzo<- round(log2(nzo2$M_nzo_avg/nzo2$F_nzo_avg), digits = 2)

#NOD
nod<- as.data.frame(f7norm[,c(4,5)])
rownames(nod)<- genes
a<- f9norm[,c(5,6)]
nod<- cbind(nod,a)
colnames(nod)<- c("F7_M_NOD", "F7_F_NOD", "F9_M_NOD", "F9_F_NOD")
nod2<- nod
nod2$M_nod_avg<- rowMeans(subset(nod, select = c(F7_M_NOD, F9_M_NOD)), na.rm = TRUE)
nod2$F_nod_avg<- rowMeans(subset(nod, select = c(F7_F_NOD, F9_F_NOD)), na.rm = TRUE)
nod2$FCMvsF_nod<- round((nod2$M_nod_avg/nod2$F_nod_avg), digits = 2)
nod2$Log2MvsF_nod<- round(log2(nod2$M_nod_avg/nod2$F_nod_avg), digits = 2)

#
#Merge them
cc_df<- cbind(pwk,c57,aj,s129,wsb,cas,nzo,nod )
cc_df2<- as.data.frame(cbind(pwk2[,8],c572[,8],aj2[,8],s1292[,8],wsb2[,8],cas2[,8],nzo2[,8],nod2[,8]))
cc_df3<- as.data.frame(cbind(pwk2[,7],c572[,7],aj2[,7],s1292[,7],wsb2[,7],cas2[,7],nzo2[,7],nod2[,7]))

rownames(cc_df2)<- genes
colnames(cc_df2)<- c("PWK", "C57", "AJ", "129", "WSB", "CAS", "NZO", "NOD")
write.table(cc_df, file = "files/cc_df.txt", sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
write.table(cc_df2, file = "cc_df2.txt", sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)

rownames(cc_df3)<- genes
colnames(cc_df3)<- c("PWK", "C57", "AJ", "129", "WSB", "CAS", "NZO", "NOD")
write.table(cc_df3, file = "cc_df3.txt", sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)


#Normalize the regular fc m/f so that C57 is 1
df3test<- as.data.frame(t(cc_df3[rownames(cc_df3)==gene,]))
df3$strain<- c(rep("PWK",4), rep("C57", 4), rep("AJ", 4), rep("129", 4), rep("WSB", 4), rep("CAS", 4), rep("NZO", 4), rep("NOD", 4))







#Beginnings of app plots
#Make the data look like we want it to
gene<- "Myh4"
df<- as.data.frame(t(cc_df[rownames(cc_df)==gene,]))
df$strain<- c(rep("PWK",4), rep("C57", 4), rep("AJ", 4), rep("129", 4), rep("WSB", 4), rep("CAS", 4), rep("NZO", 4), rep("NOD", 4))
colnames(df)[1]<-"Expression"
df$gene<- gene
df$sex<- c("M", "F", "M", "F", "F", "F", "M",
           "M", "F", "M", "F", "M", "M", "M",
           "F", "F", "F", "M", "F", "M", "F", 
           "M", "F", "M", "M", "F", "M", "F",
           "M", "F", "M", "F")
df$ss<- paste(df$strain,df$sex)

#Plot just the abundance values
p2<- df %>%
  ggplot(aes(x=strain, y=Expression, group = ss, fill = ss)) + 
  geom_boxplot() 
p2<- p2+labs(title=paste0("Normalized Protein Expression- ", gene), x="Strain", y = "Avg. Abundance (Normalized)")+
  theme_classic() 

#Plot the log2 fold change
df3<- as.data.frame(t(cc_df2[rownames(cc_df2)==gene,]))
df3$strain<- colnames(cc_df2)
colnames(df3)[1]<-"Log2FC"
p3<- ggplot(df3,aes(x=strain, y=Log2FC))+
  geom_bar(stat="identity", fill =c("#ff0000") ) 
p3<- p3+labs(title=paste0("Normalized Protein Expression- ", gene), x="Strain", y = "Log2FC M/F")+
  theme_classic() 







