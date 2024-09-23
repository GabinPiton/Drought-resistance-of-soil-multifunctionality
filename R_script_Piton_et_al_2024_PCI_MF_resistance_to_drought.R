
library(ggplot2)
library(lsmeans)
library(multcomp)
library(vegan)

# upload data from github
D=read.csv("https://raw.githubusercontent.com/GabinPiton/Drought-resistance-of-soil-multifunctionality/refs/heads/main/Piton_et_al_2024_PCI_MF_response_to_drought.csv",sep=";")


########### ANOVA3

#select variables tested with ANOVA 3
varlist<- c("MF_Mic_mean","MBC","MBN","MBP","SIR","EEA","EEC","EEN","EEP","NMP","NEA","DEA","H","DOC","logNO3","logNH4","logDON","PO4","logCN.Mic","CP.Mic","NP.Mic","DCN","DCP","DNP","DCN.MicCN","DCP.MicCP","DNP.MicNP","GPGN","Bact_Richness","Bact_shannon","Fungi_Richness","Fungi_shannon")

#Loop to report ANOVA3 statistics in a table for all variable in varlist

TAB=t(rep(NA,times=35))


for (i in varlist) {
  MDL=eval(substitute(aov(j~ Clim*Ferti*Plant, data=D,na.action=na.omit), list(j = as.name(i))))
  M=anova(MDL)
  m=as.data.frame(M)
  tab=data.frame(NA)
  
  for (y in 1:7){
    m2=m[y,]
    rownames(m2)= as.name(i)
    tab= cbind(tab,m2)
  }
  
  tab=tab[,-1]
  colnames(TAB)= colnames(tab)
  TAB= rbind(TAB,tab)
 
 
}

colnames(TAB)=c(c(paste(c(rep(row.names(m)[1],4)),(colnames(TAB)[1:5]))),
  c(paste(c(rep(row.names(m)[2],4)),(colnames(TAB)[1:5]))),
  c(paste(rep(row.names(m)[3],4),(colnames(TAB)[1:5]))),
  c(paste(rep(row.names(m)[4],4),(colnames(TAB)[1:5]))),
  c(paste(rep(row.names(m)[5],4),(colnames(TAB)[1:5]))),
  c(paste(rep(row.names(m)[6],4),(colnames(TAB)[1:5]))),
  c(paste(rep(row.names(m)[7],4),(colnames(TAB)[1:5]))))

TAB=TAB[-1,]

#Export table
write.table(TAB, "TABLE_S2_all_stat.txt")



### Loop to make BOXPLOTS on all variables tested with the ANOVA, plots are export on working repository
D$Grp<-paste(D$Plant,D$Ferti,D$Clim)
DATA<-D

dev.off()

for (i in varlist) {
  MDL=eval(substitute(lm(log(j+1)~ Clim*Plant*Ferti,data=DATA,na.action=na.omit), list(j = as.name(i))))
  M=anova(MDL)
  M
  
  cld_group<-cld(lsmeans(MDL,~Clim*Plant*Ferti),Letters=c("abcdef"))
  cld_group<-cld_group[order(cld_group$Plant,cld_group$Ferti,cld_group$Clim),]
  
  Y<-DATA[,as.character(list(j = as.name(i)))]
  png(paste(as.character(list(j = as.name(i))),"D_barplot.png",sep="_"), res=600, width=600*4 , height=450*4)
  
  p <- ggplot(DATA, aes(x=Grp, y=Y,fill=Grp)) + 
    theme_classic(base_size = 12)+
    scale_fill_manual(values=c("gray35","gray35","snow3","snow3","olivedrab4","olivedrab4","olivedrab2","olivedrab2"))+
    theme(axis.title.y  =element_text(size=9),axis.title.x  =element_text(size=8),
          axis.text=element_text(colour="black",angle=90,hjust=1),
          panel.border = element_rect(fill = NA, colour = "black", size=0.8),legend.position='none') +
    labs(x ="Plant*Ferti*Clim", y =as.character(list(j = as.name(i))))+
    geom_boxplot()+
    # stat_summary(geom = 'text', label = cld_group$.group, fun.y = max, vjust = -0.5,size=4)+
    ylim(c(min(Y),max(Y)+0.02*max(Y)))+
    scale_x_discrete(limits=c("NP NF Normal","NP NF Drought",
                              "NP F Normal","NP F Drought",
                              "P NF Normal","P NF Drought",
                              "P F Normal","P F Drought"))
  
  
  print(p)
  dev.off()
  
}

################
#PERMANOVA on Fungal and Bacteria OTU reported in Table 1

#Fungi
Fungal_OTU=read.csv("https://raw.githubusercontent.com/GabinPiton/Drought-resistance-of-soil-multifunctionality/refs/heads/main/Fungi_OTU.csv",sep=",")
#OTU=read.csv("Fungi_OTU.csv",sep=",")


treatments <- Fungal_OTU[,853:857]#Fungi

rownames(Fungal_OTU)<-Fungal_OTU$Num
Fungal_OTU<- Fungal_OTU[,3:851]#Fungi

Fungal_OTU_hell<-decostand(Fungal_OTU,"hellinger")
Fungal_OTU_dist<-vegdist(Fungal_OTU_hell,method="bray")

Fungi_perMANOVA<-adonis2(Fungal_OTU_dist ~ Plant*Ferti*Clim , data = treatments, permutations = 9999)
Fungi_perMANOVA


#Bacteria
Bact_OTU=read.csv("https://raw.githubusercontent.com/GabinPiton/Drought-resistance-of-soil-multifunctionality/refs/heads/main/Bacteria_OTU.csv",sep=",")

treatments <- Bact_OTU[,3574:3576]#bacteria

rownames(Bact_OTU)<-Bact_OTU$Num
Bact_OTU<-   Bact_OTU[,3:3572]#bacteria
Bact_OTU_hell<-decostand(Bact_OTU,"hellinger")
Bact_OTU_dist<-vegdist(Bact_OTU_hell,method="bray")

Bact_perMANOVA<-adonis2(Bact_OTU_dist ~ Plant*Ferti*Clim , data = treatments, permutations = 9999)
Bact_perMANOVA

### Correlation plots of Figure 2

#create a function for IC plot
IC <- function(x) 1.96*sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

#create group for IC and color on plot

D$CPF<-paste(D$Clim,D$Plant,D$Ferti)
D$PF<-paste(D$Plant,D$Ferti)

# Model fit
MD<-lm(MF_Mic_mean~ H, data =  D,na.action=na.omit)
summary(MD)
anova(MD)

# Model prediction line

predslm = predict(MD, interval = "confidence")
predslm<-cbind(predslm,D)

#IC calculation for the plot
STD_MF_Mic_mean = aggregate(MF_Mic_mean ~ CPF, data = D, IC)
colnames(STD_MF_Mic_mean)= c("CPF", "IC_MF_Mic_mean")
STD_H = aggregate(H ~ CPF, data = D, IC)
colnames(STD_H)= c("CPF", "IC_H")
STD<-merge(STD_MF_Mic_mean,STD_H,"CPF")
MOY_MF_Mic_mean = aggregate(MF_Mic_mean ~ CPF, data = D, mean)
colnames(MOY_MF_Mic_mean)= c("CPF", "Moy_MF_Mic_mean")
MOY_H = aggregate(H ~ CPF, data = D, mean)
colnames(MOY_H)= c("CPF", "Moy_H")
MOY<-merge(MOY_MF_Mic_mean,MOY_H,"CPF")
MOY_STD<-merge(MOY,STD,"CPF")


#Split CFP for color on plot
split<-c(strsplit(as.character(MOY_STD$CPF), " "))
N<-nrow(MOY_STD)

for (i in 1:N) {MOY_STD$Clim[i]<-split[[i]][1]
MOY_STD$Plant[i]<-split[[i]][2]
MOY_STD$Ferti[i]<-split[[i]][3]
}

MOY_STD$PF<-paste(MOY_STD$Plant,MOY_STD$Ferti)


jpeg(filename = paste("MF_H_IC","2.jpg"),
     width = 800*2.8, height = 600*2.8, units = "px", pointsize = 6,res=250)

p <- ggplot()+
  geom_ribbon(data=predslm, aes(x=H,ymin = lwr, ymax = upr, fill = "black", color = NULL), alpha = .2, show.legend = F) +
  geom_line(data=D,aes(y=predict(MD),x=H, colour="black"),linewidth=1, show.legend = F)+
  geom_errorbar(data=MOY_STD,aes(x=Moy_H,ymin=Moy_MF_Mic_mean-IC_MF_Mic_mean, ymax=Moy_MF_Mic_mean+IC_MF_Mic_mean), width=0.05*mean(MOY_STD$Moy_H))+
  geom_errorbarh(data=MOY_STD,aes(y=Moy_MF_Mic_mean,xmin=Moy_H-IC_H, xmax=Moy_H+IC_H), height=0.5*mean(abs(MOY_STD$Moy_MF_Mic_mean)))+
  geom_point(data=MOY_STD,aes(y=Moy_MF_Mic_mean,x=Moy_H,fill=PF, shape=Clim, colour="black"), size = 6, alpha = 1, show.legend = F, stroke = 1)+
  scale_shape_manual(values=c(22,21,24,23)) +
  scale_fill_manual(values=c("black","gray30","gray80","olivedrab4","olivedrab2")) +
  scale_color_manual(values=c("black","gray30","gray80","olivedrab4","olivedrab2")) +
  scale_linetype_manual(values=c("solid","solid","solid","solid"))+
  xlab(expression(atop(paste("Moisture"),paste("%"))))+
  ylab(expression("Microbial multifunctionality"))+
  theme_classic(base_size=22) 
p

dev.off()
print(p)






#Same with stoichiometric imbalance


MD<-lm(MF_Mic_mean~ DCN.MicCN, data =  D,na.action=na.omit)
summary(MD)
anova(MD)


predslm = predict(MD, interval = "confidence")
predslm<-cbind(predslm,D)


STD_MF_Mic_mean = aggregate(MF_Mic_mean ~ CPF, data = D, IC)
colnames(STD_MF_Mic_mean)= c("CPF", "IC_MF_Mic_mean")
STD_H = aggregate(DCN.MicCN ~ CPF, data = D, IC)
colnames(STD_H)= c("CPF", "IC_H")
STD<-merge(STD_MF_Mic_mean,STD_H,"CPF")
MOY_MF_Mic_mean = aggregate(MF_Mic_mean ~ CPF, data = D, mean)
colnames(MOY_MF_Mic_mean)= c("CPF", "Moy_MF_Mic_mean")
MOY_H = aggregate(DCN.MicCN ~ CPF, data = D, mean)
colnames(MOY_H)= c("CPF", "Moy_H")
MOY<-merge(MOY_MF_Mic_mean,MOY_H,"CPF")
MOY_STD<-merge(MOY,STD,"CPF")


split<-c(strsplit(as.character(MOY_STD$CPF), " "))
N<-nrow(MOY_STD)

for (i in 1:N) {MOY_STD$Clim[i]<-split[[i]][1]
MOY_STD$Plant[i]<-split[[i]][2]
MOY_STD$Ferti[i]<-split[[i]][3]
}

MOY_STD$PF<-paste(MOY_STD$Plant,MOY_STD$Ferti)


jpeg(filename = paste("MF_DCN.MicCN_IC","2.jpg"),
     width = 800*2.8, height = 600*2.8, units = "px", pointsize = 6,res=250)

p <- ggplot()+
  geom_ribbon(data=predslm, aes(x=DCN.MicCN,ymin = lwr, ymax = upr, fill = "black", color = NULL), alpha = .2, show.legend = F) +
  geom_line(data=D,aes(y=predict(MD),x=DCN.MicCN, colour="black"),linewidth=1, show.legend = F)+
  geom_errorbar(data=MOY_STD,aes(x=Moy_H,ymin=Moy_MF_Mic_mean-IC_MF_Mic_mean, ymax=Moy_MF_Mic_mean+IC_MF_Mic_mean), width=0.05*mean(MOY_STD$Moy_H))+
  geom_errorbarh(data=MOY_STD,aes(y=Moy_MF_Mic_mean,xmin=Moy_H-IC_H, xmax=Moy_H+IC_H), height=0.5*mean(abs(MOY_STD$Moy_MF_Mic_mean)))+
  geom_point(data=MOY_STD,aes(y=Moy_MF_Mic_mean,x=Moy_H,fill=PF, shape=Clim, colour="black"), size = 6, alpha = 1, show.legend = F, stroke = 1)+
  scale_shape_manual(values=c(22,21,24,23)) +
  scale_fill_manual(values=c("black","gray30","gray80","olivedrab4","olivedrab2")) +
  scale_color_manual(values=c("black","gray30","gray80","olivedrab4","olivedrab2")) +
  scale_linetype_manual(values=c("solid","solid","solid","solid"))+
  xlab(expression(atop(paste("C:N imbalance"),paste("(DOC:TDN)/(Microbial-C:N)"))))+
  ylab(expression("Microbial multifunctionality"))+
  theme_classic(base_size=22) 
p

dev.off()
print(p)






#Same with GPGN


MD<-lm(MF_Mic_mean~ GPGN, data =  D,na.action=na.omit)
summary(MD)
anova(MD)
predslm = predict(MD, interval = "confidence")
predslm<-cbind(predslm,D)


STD_MF_Mic_mean = aggregate(MF_Mic_mean ~ CPF, data = D, IC)
colnames(STD_MF_Mic_mean)= c("CPF", "IC_MF_Mic_mean")
STD_H = aggregate(GPGN ~ CPF, data = D, IC)
colnames(STD_H)= c("CPF", "IC_H")
STD<-merge(STD_MF_Mic_mean,STD_H,"CPF")
MOY_MF_Mic_mean = aggregate(MF_Mic_mean ~ CPF, data = D, mean)
colnames(MOY_MF_Mic_mean)= c("CPF", "Moy_MF_Mic_mean")
MOY_H = aggregate(GPGN ~ CPF, data = D, mean)
colnames(MOY_H)= c("CPF", "Moy_H")
MOY<-merge(MOY_MF_Mic_mean,MOY_H,"CPF")
MOY_STD<-merge(MOY,STD,"CPF")


split<-c(strsplit(as.character(MOY_STD$CPF), " "))
N<-nrow(MOY_STD)

for (i in 1:N) {MOY_STD$Clim[i]<-split[[i]][1]
MOY_STD$Plant[i]<-split[[i]][2]
MOY_STD$Ferti[i]<-split[[i]][3]
}

MOY_STD$PF<-paste(MOY_STD$Plant,MOY_STD$Ferti)


jpeg(filename = paste("MF_GPGN_IC","2.jpg"),
     width = 800*2.8, height = 600*2.8, units = "px", pointsize = 6,res=250)

p <- ggplot()+
  geom_ribbon(data=predslm, aes(x=GPGN,ymin = lwr, ymax = upr, fill = "black", color = NULL), alpha = .2, show.legend = F) +
  geom_line(data=D,aes(y=predict(MD),x=GPGN, colour="black"),linewidth=1, show.legend = F)+
  geom_errorbar(data=MOY_STD,aes(x=Moy_H,ymin=Moy_MF_Mic_mean-IC_MF_Mic_mean, ymax=Moy_MF_Mic_mean+IC_MF_Mic_mean), width=0.05*mean(MOY_STD$Moy_H))+
  geom_errorbarh(data=MOY_STD,aes(y=Moy_MF_Mic_mean,xmin=Moy_H-IC_H, xmax=Moy_H+IC_H), height=0.5*mean(abs(MOY_STD$Moy_MF_Mic_mean)))+
  geom_point(data=MOY_STD,aes(y=Moy_MF_Mic_mean,x=Moy_H,fill=PF, shape=Clim, colour="black"), size = 6, alpha = 1, show.legend = F, stroke = 1)+
  scale_shape_manual(values=c(22,21,24,23)) +
  scale_fill_manual(values=c("black","gray30","gray80","olivedrab4","olivedrab2")) +
  scale_color_manual(values=c("black","gray30","gray80","olivedrab4","olivedrab2")) +
  scale_linetype_manual(values=c("solid","solid","solid","solid"))+
  xlab(expression("Gram positive : Gram negative"))+
  ylab(expression("Microbial multifunctionality"))+
  theme_classic(base_size=22) 
p

dev.off()
print(p)


#Same with bact shannon



MD<-lm(MF_Mic_mean~ Bact_shannon, data =  D,na.action=na.omit)
summary(MD)
anova(MD)
predslm = predict(MD, interval = "confidence")
predslm<-cbind(predslm,D)


STD_MF_Mic_mean = aggregate(MF_Mic_mean ~ CPF, data = D, IC)
colnames(STD_MF_Mic_mean)= c("CPF", "IC_MF_Mic_mean")
STD_H = aggregate(Bact_shannon ~ CPF, data = D, IC)
colnames(STD_H)= c("CPF", "IC_H")
STD<-merge(STD_MF_Mic_mean,STD_H,"CPF")
MOY_MF_Mic_mean = aggregate(MF_Mic_mean ~ CPF, data = D, mean)
colnames(MOY_MF_Mic_mean)= c("CPF", "Moy_MF_Mic_mean")
MOY_H = aggregate(Bact_shannon ~ CPF, data = D, mean)
colnames(MOY_H)= c("CPF", "Moy_H")
MOY<-merge(MOY_MF_Mic_mean,MOY_H,"CPF")
MOY_STD<-merge(MOY,STD,"CPF")


split<-c(strsplit(as.character(MOY_STD$CPF), " "))
N<-nrow(MOY_STD)

for (i in 1:N) {MOY_STD$Clim[i]<-split[[i]][1]
MOY_STD$Plant[i]<-split[[i]][2]
MOY_STD$Ferti[i]<-split[[i]][3]
}

MOY_STD$PF<-paste(MOY_STD$Plant,MOY_STD$Ferti)


jpeg(filename = paste("MF_Bact_shannon_IC","2.jpg"),
     width = 800*2.8, height = 600*2.8, units = "px", pointsize = 6,res=250)

p <- ggplot()+
  geom_ribbon(data=predslm, aes(x=Bact_shannon,ymin = lwr, ymax = upr, fill = "black", color = NULL), alpha = .2, show.legend = F) +
  geom_line(data=D,aes(y=predict(MD),x=Bact_shannon, colour="black"),linewidth=1, show.legend = F)+
  geom_errorbar(data=MOY_STD,aes(x=Moy_H,ymin=Moy_MF_Mic_mean-IC_MF_Mic_mean, ymax=Moy_MF_Mic_mean+IC_MF_Mic_mean), width=0.003*mean(MOY_STD$Moy_H))+
  geom_errorbarh(data=MOY_STD,aes(y=Moy_MF_Mic_mean,xmin=Moy_H-IC_H, xmax=Moy_H+IC_H), height=0.5*mean(abs(MOY_STD$Moy_MF_Mic_mean)))+
  geom_point(data=MOY_STD,aes(y=Moy_MF_Mic_mean,x=Moy_H,fill=PF, shape=Clim, colour="black"), size = 6, alpha = 1, show.legend = F, stroke = 1)+
  scale_shape_manual(values=c(22,21,24,23)) +
  scale_fill_manual(values=c("black","gray30","gray80","olivedrab4","olivedrab2")) +
  scale_color_manual(values=c("black","gray30","gray80","olivedrab4","olivedrab2")) +
  scale_linetype_manual(values=c("solid","solid","solid","solid"))+
  xlab(expression(atop(paste("Bacterial diversity"),paste("(Shannon index)"))))+
  ylab(expression("Microbial multifunctionality"))+
  theme_classic(base_size=22) 
p

dev.off()
print(p)




