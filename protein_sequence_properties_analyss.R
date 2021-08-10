# Arabidopsis Proteome sequence/structure properties analysis

# 1. convert TAIR to Uniprot accession
# 2. protein sequence properties analysis : 
# 2.1 bacis analysis: length	Mass	pI	charge	alphaHelixIndex	BetaSheet	hdrophobicity	instabilityIndex
# 2.2 transmembrane property analysis: TM	TMregion	expAA	PredHel	aveTM
# 3. protein 3D structure property analysis:
# 3.1 PDB id download
# 3.2 PBD structure download
# 3.3 homology structure download from SWISS-MODEL
# 3.4 protein solvent accessible area calculation 
# 4. statistical analysis and results visualization



install.packages("Peptides")
install.packages("seqinr")
library(Peptides)
library(seqinr)

library(tidyverse)
# GET all proteins fasta sequence
#method 1. from tair bulk sequence downloading function
#method 2. 

TAIR10<-read.fasta(file = "TAIR10.fasta", seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
pro.raw<-read.csv("Protein_raw.csv",stringsAsFactors = FALSE)
pro.tair.name<-pro.raw$Protein.Acc
pro.seq<-TAIR10[names(TAIR10) %in% pro.tair.name]
length(pro.seq)
length(TAIR10)
write.fasta(sequences = pro.seq, names = names(pro.seq),file.out = "pro.seq.fasta")


##-----------------GET protein properties: length, mw, pi, hydrophobicity, etc--------------

my_seq<-pro.seq
#length
z<-data.frame()
for (i in 1:2753){
  z<-append(z,lengthpep(my_seq[[i]]))
}

dat <- do.call("rbind.data.frame", z)
dat$Protein.Acc <- names(my_seq)
properties<-merge(pro.raw,dat,by='Protein.Acc')
names(properties)[length(names(properties))]<-"length" 

#mass
z<-data.frame()
for (i in 1:2753){
  z<-append(z,mw(my_seq[[i]],monoisotopic = FALSE))
}

dat <- do.call("rbind.data.frame", z)
dat$Protein.Acc <- names(my_seq)
properties<-merge(properties,dat,by='Protein.Acc')
names(properties)[length(names(properties))]<-"Mass" 


#isoelectic point (pI) 
z<-data.frame()
for (i in 1:2753){
  z<-append(z,pI(my_seq[[i]],pKscale = "EMBOSS" ))
}

dat = do.call("rbind.data.frame", z)
dat$Protein.Acc <- names(my_seq)
properties<-merge(properties,dat,by='Protein.Acc')
names(properties)[length(names(properties))]<-"pI" 

#charge
z<-data.frame()
for (i in 1:2753){
  z<-append(z,charge(my_seq[[i]],pH = 7.8, pKscale = "Lehninger"))
}

dat = do.call("rbind.data.frame", z)
dat$Protein.Acc <- names(my_seq)
properties<-merge(properties,dat,by='Protein.Acc')
names(properties)[length(names(properties))]<-"charge" 

# ALPHA HELIX VALUE
z<-data.frame()
for (i in 1:2753){
  z<-append(z,hmoment(my_seq[[i]],angle = 100, window = 11))
}
dat = do.call("rbind.data.frame", z)
dat$Protein.Acc <- names(my_seq)
properties<-merge(properties,dat,by='Protein.Acc')
names(properties)[length(names(properties))]<-"alphaHelixIndex" 

# BETA SHEET VALUE
z<-data.frame()
for (i in 1:2753){
  z<-append(z,hmoment(my_seq[[i]],angle = 160, window = 11))
}
dat = do.call("rbind.data.frame", z)
dat$Protein.Acc <- names(my_seq)
properties<-merge(properties,dat,by='Protein.Acc')
names(properties)[length(names(properties))]<-"BetaSheet" 

# Hydrophobicity
z<-data.frame()
for (i in 1:2753){
  z<-append(z,hydrophobicity(my_seq[[i]], scale = "KyteDoolittle"))
}
dat = do.call("rbind.data.frame", z)
dat$Protein.Acc <- names(my_seq)
properties<-merge(properties,dat,by='Protein.Acc')
names(properties)[length(names(properties))]<-"hdrophobicity" 

# instability index
z<-data.frame()
for (i in 1:2753){
  z<-append(z,instaIndex(my_seq[[i]]))
}
dat = do.call("rbind.data.frame", z)
dat$Protein.Acc <- names(my_seq)
properties<-merge(properties,dat,by='Protein.Acc')
names(properties)[length(names(properties))]<-"instabilityIndex" 

write_csv(properties,"pro.basicSeq.properties.csv")

#Transmembran analysis-------------------------------
#download membrane_protein.tair10.txt file from tair, TM total length prediction 
#was analyzed using TMHMM server (v.2.0) (http://www.cbs.dtu.dk/services/TMHMM/)

tair.tm<-read_delim("Membrane_Proteins.tai10.txt",delim="\t",col_names = FALSE)
names(tair.tm) <-c("Protein.Acc","TM","TMregion")

properties<- properties%>% left_join(tair.tm,by="Protein.Acc")

tm.raw<-properties %>% filter(TM>=1) %>% select(c(1,14,15))

TMHMM<-read_csv("TMHMM.csv")
TMHMM<-TMHMM%>% separate(expAA,into = c("temp","expAA"),sep = "=") %>%
                separate(PredHel,into = c("temp1","PredHel"),sep = "=")

tm.raw <-tm.raw %>% left_join(TMHMM[,c(1,4,7)],by=c("Protein.Acc"="accession")) 

tm.raw[,c(4:5)]<-sapply(tm.raw[,4:5], as.numeric)
tm.raw$aveTM<- tm.raw$expAA/tm.raw$PredHel

properties <- properties %>% left_join(tm.raw,by="Protein.ACC")
                    
write_csv(properties,"pro.basicSeq.properties.TM.csv")          

# PDB polar, apolar,surface, buried analysis------------------------------
# get protein PDB id from uniprot
#download PDB file from PDB website

BiocManager::install("bio3d")
library(bio3d)

#get.pdb("3gcb")

pdb_id<-read.csv("pdb-id.csv")
pdb_id1<-pdb_id[,c(1,7)]

for(i in 1:135){
  pdb_id2 <- as.character(pdb_id1[i,2])
  get.pdb(pdb_id2)
}

#for proteins without available 3D structure, the protein homology model -----------
#with a normalized QMEAN at least 0.6 was downloaded from SWISS-MODEL 
#(https://swissmodel.expasy.org/).

library("jsonlite")

#swiss.index<-read_delim("INDEX",delim = '\t')
swiss.indexjson<-read_json("INDEX.json")
swiss.indexjson.df<-as.data.frame(do.call(rbind, lapply(swiss.indexjson$index, as.vector)))

pro.swiss<-merge(pro.raw,swiss.indexjson.df,by.x = "uniprot",by.y = "uniprot_ac",all.x = TRUE,all.y = FALSE)

pro.swiss <- apply(pro.swiss,2,as.character)
write.csv(pro.swiss,"pro.swiss.csv")

#--download all models data from website, then try to find the target ones and copy them-----
# to new directory

fname<-read.csv("swiss-id.csv")

f = list.files("./")
f<- as.data.frame(f)

fappend<-data.frame()

for (i in 1:10) {
  f1<- as.character(fname[i,])
  fsave <- as.data.frame(f[grep(f1,f$f),])
  fappend <- rbind(fappend,fsave)
  
}

for (i in 1:10){
  file.copy(as.character(fappend[i,]), "C:/Users/liup1/Documents/Programming-projects/up-low-phase-properties/3")
}

# The protein solvent accessible surface area was calculated using GETAREA (1.0 beta)---------- 
#(http://curie.utmb.edu/area_man.html).using multiple pdb files, using "Getarea_Remote.pl"  
# after running getarea_remote.pl, get pdb.txt files.now organize these files.
f = list.files(path = "./areatxt")
f

dat = lapply(f, function(i){
  x = read.csv(i, header = TRUE)
  # Get the columns you want, e.g. 1, 3, 5
  # x = x[, c(1, 3, 5)]
  # You may want to add a column to say which file they're from
  x$file = i
  # Return your data
  x })

dat = do.call("rbind.data.frame", dat)
write.csv(dat,file="areaall2.csv")

dat2<-read.csv("areaall_anlaysis.csv")
library(tidyr)

data_wide <-  pivot_wider(dat2, names_from=variable, values_from=value)

write.csv(data_wide,"areaall_analysis2.csv")

##statistical analysis,  Kruskal-Wallis Test was used------------ 
#for comparison the difference among groups. and generate ggplot figures.

dat2<-read.csv("charge.csv")

library(ggpubr)

myGRAVY<-data.frame(dat2$code,dat2$GRAVY)

pGRAVY<-myGRAVY %>%
  mutate(dat2.code = factor(dat2.code, levels=c("<0.67", "0.67-1.5", ">1.5"))) %>%
  ggplot( aes(x=dat2.code, y=dat2.GRAVY,fill=dat2.code)) +
  geom_boxplot()+
  scale_fill_brewer(palette="RdBu")+
  labs(title="Hydrophobicity",x="Up/Low Ratio", y = "GRAVY")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank())

pGRAVY

KRGRAVY <- kruskal.test(dat2.GRAVY ~ dat2.code, data = myGRAVY)

pvalue2 <- data.frame("sample"=c("GRAVY","aIndex","charge","aHelix","bsheet","instability"),
                      "pvalue"=c(KRGRAVY$p.value,KRaIndex$p.value,KRcharge$p.value,KRaHelix$p.value,
                                 KRbsheet$p.value,KRinsta$p.value))
pvalue<-rbind(pvalue,pvalue2)
pvalue<-arrange(pvalue,desc(pvalue))

PKRlength<-pairwise.wilcox.test(my_mass$dat.Length, my_mass$dat.code,
                                p.adjust.method = "BH")
library(ggpubr)

ggarrange(pTM, plength, pmass, paHelix,pGRAVY,
          pbsheet, pIEP, pcharge, pinsta, paIndex,
          ncol=5, nrow=2, common.legend = TRUE, legend="bottom")

write.csv(pvalue,file="pvalue.csv")

##subcellular localization
# Protein subcellular localization prediction was based on the TAIR10 subcellular 
# prediction database from TAIR and Subcellular Localization Database for Arabidopsis Proteins (SUBA). 
# stacked plot for subcellular
mydata<-read.csv("stacked-barchart-sub.csv")
mydata<-mydata[,c(1,2,3,4)]

longdata <- melt(mydata, id.vars = c("location"))
write.csv(longdata,"longdata.csv")

mydata<-read.csv("longdata.csv")

mydata<-mydata%>%rename(Ratio=rating)

tiff("total_subcellular_percentage_stacked.tiff",width = 7, height = 4, units = 'in', res = 300, compression = 'lzw')

pdf("total_subcellular_percentage_stacked2.pdf",width = 7,height = 4)
mydata%>%
  mutate(Ratio = factor(Ratio, levels=c("Ratio(>1.5)", "Ratio(0.67-1.5)", "Ratio(<0.67)"))) %>%
  mutate(location=factor(location,levels=c("PM","Extracellular","Nucleus","Golgi",
                                           "Cytosol","Unknown","Plastid","Others","Vacuole","ER","Peroxisome","Mitochondrion"))) %>%
  ggplot(aes(y=number, x=location, fill=Ratio)) + 
  scale_fill_manual(values=c("#FAA582", "#BABABA", "#92C5DE"))+
  geom_bar(stat="identity") + 
  xlab('') + ylab('Percentage')+
  theme(axis.text.x = element_text(angle = 45,hjust=0, vjust=0,size = 14,color = "black"),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12))+
  scale_x_discrete(position = "top") 


dev.off()

