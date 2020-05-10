##############################################################
##Limma stats
##############################################################
library(limma)

TS<-factor(colnames(Norm), levels = c("M0",
                                      "M1",
                                      "M2",
                                      "H0",
                                      "H1",
                                      "H2"))
design<-model.matrix(~0+TS) #based on page 45 of the limma user guide
colnames(design)<-levels(TS)
design
fit<-lmFit(Norm, design)

#lets make a contrast matrix for the MCF7 data
#we don't need to actually run all three
#same genes list is if you add M1-M2
cont.M<-makeContrasts(
  M0vM1=M0-M1,
  M0vM2=M0-M2,
  H0vH1=H0-H1,
  H0vH2=H0-H2,
  DiffG0vG1=(M0-M1)-(H0-H1),
  DiffG0vG2=(M0-M2)-(H0-H2),
  levels=design)

fit2<-contrasts.fit(fit, cont.M)
fit2<-eBayes(fit2)
topTable(fit2)

results<-decideTests(fit2)
#get a tabular table of results
summary(results)

write.fit(fit2, results, "Fucci.write.fit.out.txt", adjust = "BH")
###lets add the limma stats to the stats table just
#for ftest and then run a scatter plot to see how they look 
AnnotStats$Limma.Ftest.pvalue<-fit2$F.p.value
#lets plot the two ftest statistics because the limma seems way 
plot(AnnotStats$F.mt.rawp, AnnotStats$Limma.Ftest.pvalue)
#lets adjust the AnnotStats ftest pvalue not that BH is the same as fdr
AnnotStats$Limma.F.p.BH<-p.adjust(AnnotStats$Limma.Ftest.pvalue,
                                  method = "BH")

#lets add the other important stats 
#apparently all the pvalues in fit2 are 
#accessible by in the $p.value table make
#sure to use fit2 not fit1
colnames(fit2$p.value)
AnnotStats$Lm.p.M0vM1<-fit2$p.value[,1]
AnnotStats$Lm.p.H0vH1<-fit2$p.value[,3]
AnnotStats$Lm.p.M0vM2<-fit2$p.value[,2]
AnnotStats$Lm.p.H0vH2<-fit2$p.value[,4]
AnnotStats$Lm.p.Int.G0vG1<-fit2$p.value[,5]
AnnotStats$Lm.p.Int.G0vG2<-fit2$p.value[,6]

#interesting plot looking at Limma pvalues vs rowftest pvalues
plot(AnnotStats$G0vG1m.pvalue,AnnotStats$Lm.p.M0vM1)


write.table(AnnotStats,file = "AnnotStats.txt", sep = "\t")