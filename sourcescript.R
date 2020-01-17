library(BIEN)
library(rotl)
library(geiger)
library(ape)
library(readr)
library(caper)
library(phytools)

TV <- read_csv("~/R/TVEL.csv")
head(TVEL)
#LA <- na.omit(LA$trait_value)
TV <- TV[rev(order(TV$velocity)),]
TV <- TV[!duplicated(TV$species),]
TV <- TV[!duplicated(TV),]
#TV <- TV[sample(1:nrow(TV), 250, replace=TRUE),]
TV <- TV[complete.cases(TV[ , 2]),]
TV <- TV[complete.cases(TV[ , 1]),]
write.csv(TV, file = "TV.csv",row.names=TRUE)
as.matrix(TV)
rownames(TV)<-TV[,2]
as.data.frame(TV)
row.names(TV)<-TV$species
tibble::column_to_rownames(as.data.frame(TV$species))
View(TV)
TV<-TV[,-2]

library("brranching")
tr <- phylomatic(TV$species, taxnames = TRUE, get = "POST", informat = "newick",
                 method = "phylomatic", storedtree = "zanne2014", treeuri = NULL,
                 taxaformat = "slashpath", outformat = "newick", clean = TRUE,
                 db = "apg", verbose = TRUE)
TVtree <- tr
tr <- TVtree
tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores=TRUE)
plot(TVtree)
plot(tr)
is.binary.tree(tr)#use next line if FALSE
tr <- multi2di(tr)
row.names(tr) <- TV$species
matches <- match(TV$species, tr$tip.label, nomatch = 0)
TV <- subset(TV, matches != 0)
TV
tr <- drop.tip(tr, setdiff(tr$tip.label,TV$species))
name.check(tr, TV)
tr$tip.label <- Hmisc::capitalize(tr$tip.label)
tr<-pruned.tree
TV <- TV[tr$tip.label, ]
TV <- TV[match(tr$tip.label,rownames(TV)),]
setdiff(tr$tip.label, TV$species)
setdiff(TV$species, tr$tip.label)
TV$species
tr$tip.label




library("caper")
comparative.data(tr, TV, species, vcv=T, na.omit=TRUE,
                 force.root=FALSE, warn.dropped=FALSE, scope=NULL)
library(phytools)
contMap(pruned.tree, TV$velocity)
fastAnc(pruned.tree, TV$velocity)


#adding a geological time scale !!!tree$root.time is missing, check tree is time scaled!!!
library(strap) #loading the library
tr$root.time <- max(nodeHeights(tr))
geoscalePhylo(tr, cex.ts=0.6, cex.tip=0.6)

PICreconstruction <- ace(TV$velocity, tr, type="continuous", method="pic")
PICreconstruction
traitvaluefinal <- c(TV$trait_value, PICreconstruction$ace)
traitvaluefinal
plot(tr, show.tip.label=FALSE)
tiplabels(pch = 21, cex=(TV$velocity))
nodelabels(pch = 21, cex=(PICreconstruction$ace)*1000)

BMfit <- fitContinuous(tr, TV$velocity, model="BM")
FBM <-fastBM(tr)
aa <- fastAnc(tr, TVELO)

TVELO <- as.vector(as.numeric(TV$log.velocity))
names(TVELO) <- as.vector(TV$species)
TVELO
Anc <- fastAnc(tr, TVELO)
Anc
tr$tip.label

TVELO
#phenogram
phenogram(tr, TVELO, spread.labels = TRUE)
par(mai=c(2.02,1.82,0.82,0.42))

#density traitgram (95% CI)
fancyTree(tr, type = "phenogram95", x = TVELO, spread.labels = TRUE)
par(mai=c(2.02,1.82,0.82,0.42))

TVELO <- as.vector(as.numeric(TV$velocity))
names(TVELO) <- as.vector(TV$species)
TVELO
Anc <- fastAnc(tr, TVELO)
Anc
tr$tip.label
TV$species
##create contMap
XX <- contMap(tr, TVELO, fsize = 0.1, legend=FALSE)
plot(XX, fsize = c(0.4, 0.5))
plot(XX, type = "fan", fsize = c(0.1, 0.5))
add.color.bar(40,obj$cols,title="Seed Terminal Velocity (m/s)",
              lims=obj$lims,digits=3,prompt=FALSE,x=0,
              y=1-0.08*(Ntip(obj$tree)-1),lwd=4,fsize=1,subtitle="")

## simulate a third character
z <- fastBM(tree, sig2 = 2)
XYZ <- cbind(XY, z)
fancyTree(tree, X = XYZ, type = "scattergram")

branchtimes <- branching.times(tr)
branchtimes
datum <- cbind(branchtimes, PICreconstruction$ace)
datum
write.csv(datum, file = "TVdates.csv",row.names=TRUE)
TVdates <- read_csv("TVdates.csv")
View(TVdates)

plot(TVdates$X3 ~ TVdates$branchtimes, ylab="Reconstucted Trait Value at Node", xlab= "Node Date")
fit<-lm(TVdates$X3 ~ TVdates$branchtimes)
abline(fit)

library(mgcv)
z <- gam(TVdates$X3 ~ s(TVdates$branchtimes), data = TVdates, family = gaussian,
         method = "GCV.Cp")
residuals(z, type = "response")
summary(z)
plot(z)
init_fin <- par("fin") #gives 7 7
init_mai <- par("mai") #gives 1.02 0.82 0.82 0.42

par(mai=c(2.02,1.82,0.82,0.42))

plot(z, se = 1, seWithMean=TRUE)
plot(z, se = 1, seWithMean = TRUE, rug = TRUE, shift = mean(predict(z)),
     trans = exp)  

VELO<-as.data.frame(TVELO)
VELO<-rownames_to_column(VELO)
colnames(VELO) <- c("species","velocity")
VELO$species->rownames(VELO)
VELO<-as.data.frame(VELO)
View(VELO)
tr$node.label
Tree <- tr
tr
tr<-makeLabel(tr)
tr$node.label
tr$tip.label

remove.packages("caper")

name.check(tr, VELO)
cdat <- comparative.data(data=TVV, phy=tr, names.col="species",vcv=T)

pglsModel <- gls(TVELO ~ 1, correlation = corBrownian(phy = tr),
                 data = TVELO, method = "ML")
summary(pglsModel)
est.lambda <- pgls(VELO$velocity ~ 1, data = cdat, lambda = "ML")
summary(est.lambda)
plot(est.lambda)

####
library(picante)
Kcalc(TVV[tr$tip.label], tr) #0.007676882
phylosignal(TVV[tr$tip.label], tr, reps = 1000) #K = 0.007676882, PIC.variance.obs = 0.05601669, PIC.variance.rnd.mean = 0.2358978, PIC.variance.P = 0.000999001 PIC.variance.Z = -0.8077313
phylosig(tr, TVV, method="K",test=TRUE) #$lambda[1] 0.939423; $logL[1] 104.8503; $logL0[1] -315.0623 ;$P[1] 1.184583e-184

Kcalc(TVVA[trA$tip.label], trA) #0.0734311
phylosignal(TVVA[trA$tip.label], trA, reps = 1000) #K = 0.0734311, PIC.variance.obs = 0.05370388, PIC.variance.rnd.mean = 0.236527, PIC.variance.P = 0.000999001 PIC.variance.Z = -2.318462
phylosig(trA, TVVA, method="lambda",test=TRUE) #$lambda[1] 0.8681716; $logL[1] -31.28509; $logL0[1] -51.4735 ;$P[1] 2.094121e-10
phylosig(trA, TVVA, method="lambda",test=TRUE) #$lambda[1] 0.8681716; $logL[1] -31.28509; $logL0[1] -51.4735 ;$P[1] 2.094121e-10

###
library(ape)
richness.yule.test(TVVA, trA)
mcconwaysims.test(trA)
yule.cov(tr, formula, data = NULL)
yule.time(trA, TVVA, BIRTH = NULL, root.time = 0, opti = "nlm", start = 0.01)

diversi.gof((branching.times(tr)), null = "exponential", z = NULL) 
tra<-extract.clade(tr, 1763)
diversi.gof((branching.times(tra)), null = "exponential", z = NULL) 

#OUTPUT (FULL)
##Tests of Constant Diversification Rates:
##Data: (branching.times(tr)) 
##Number of branching times: 1019 
##Null model: exponential
##Anderson-Darling test: A2 = 156.333   P < 0.01
##Cramer-von Mises test: W2 = 109.493   P < 0.01

#OUTPUT (Asteraceae)
#Tests of Constant Diversification Rates
#Data: (branching.times(tra)) 
#Number of branching times: 134 
#Null model: exponential
#Cramer-von Mises test: W2 = 13.134   P < 0.01
#Anderson-Darling test: A2 = 14.328   P < 0.01

diversi.time((branching.times(tr)), census = NULL, censoring.codes = c(1, 0), Tc = NULL)
diversi.time((branching.times(tra)), census = NULL, censoring.codes = c(1, 0), Tc = NULL)

?diversi.time

library(phylobase)
install.packages(dispeRsal)
library("taxize")
library("ncbit", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("wikitaxa", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("taxa", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("taxonomizr")
library("taxizedb")
library("Taxonstand")
library("readr")

View(TVEL)
taxa<-TPL(TVEL$species, genus = NULL, species = NULL, infrasp = NULL,
          infra = TRUE, corr = TRUE, diffchar = 2, max.distance = 1,
          version = "1.1", encoding = "UTF-8", author = TRUE,
          drop.lower.level = FALSE, file = "", silent = TRUE, repeats = 6)
tx <- taxa
tx <- cbind(tx, TVEL)
View(tx)

boxplot(tx$velocity ~ tx$Family)
boxplot(tx$velocity ~ tx$Genus)

tx <- tx[rev(order(tx$velocity)),]
tx <- tx[!duplicated(tx$species),]
tx <- tx[!duplicated(tx),]
#TV <- TV[sample(1:nrow(TV), 250, replace=TRUE),]
tx <- TV[complete.cases(tx[ , 2]),]
tx <- TV[complete.cases(tx[ , 1]),]
write.csv(tx, file = "tx.csv",row.names=TRUE)

devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup)
txt<-lookup_table(tx$species, by_species=TRUE,include_counts = TRUE)
knitr::kable(head(add_higher_order(txt)))
txt<-tibble::rownames_to_column(txt)
colnames(txt)[1] <- "species"
#txt<-txt[,-2]
View(txt)

txtx <- plyr::join(tx, txt, by = c("species"), type = "right")
View(txtx)

plot(txtxt$number.of.accepted.species ~ txtxt$log.velocity)
abline(lm(txtxt$number.of.accepted.species ~ txtxt$log.velocity))
plot(lm(txtxt$number.of.accepted.species ~ txtxt$log.velocity))

plot(txtxt$number.of.accepted.species ~ txtxt$velocity)
abline(lm(txtxt$number.of.accepted.species ~ txtxt$velocity))
plot(lm(txtxt$number.of.accepted.species ~ txtxt$velocity))

colnames(dispersal_type)[1] <- "species"

txtxt <- plyr::join(dispersal_type, txtx, by = c("species"), type = "right")
write.csv(txtxt, file = "txtxt.csv",row.names=TRUE)

plot(Asteraceae$number.of.accepted.species ~ Asteraceae$log.velocity)
abline(lm(Asteraceae$number.of.accepted.species ~ Asteraceae$log.velocity))
plot(lm(Asteraceae$number.of.accepted.species ~ Asteraceae$log.velocity))

boxplot(Asteraceae$log.velocity ~ Asteraceae$Genus)

A <- mean(Asteraceae$log.velocity[Asteraceae$Genus])

r1<-tapply(Asteraceae$log.velocity, Asteraceae$Genus, mean)
r2<-tapply(Asteraceae$number.of.accepted.species, Asteraceae$Genus, mean)

u1<-tapply(Asteraceae$log.velocity, Asteraceae$Genus, mean)
u2<-tapply(Asteraceae$number.of.accepted.and.unresolved.species, Asteraceae$Genus, mean)


plot(r2 ~ r1)
abline(lm(r2 ~ r1))
plot(lm(r2 ~ r1))
mdl<-lm(r2 ~ r1) ### taraxacum, senecio, hieracium
summary(mdl) ###

plot(u2 ~ u1)
abline(lm(u2 ~ u1))
plot(lm(u2 ~ u1))
mudl<-lm(u2 ~ u1)
summary(mudl)

library(taxonlookup)
DDDE<-lookup_table(DispersalDistanceDataEDIT$Species, by_species=TRUE,include_counts = TRUE)
View(DDDE)
DDDE<-tibble::rownames_to_column(DDDE)
colnames(DDDE)[1] <- "species"
colnames(DispersalDistanceDataEDIT)[1] <- "species"
#txt<-txt[,-2]
View(DDDE)
DDD <- plyr::join(DispersalDistanceDataEDIT, DDDE, by = c("species"), type = "right")

View(DDD)

ecdf.x <- ecdf(DDD$`Seed_terminal_velocity_(m/s)`)
ecdf<-ecdf.x(DDD$`Seed_terminal_velocity_(m/s)`)


dr1<-tapply(DDD$`Maximum_dispersal_distance_analysis_(m)`, DDD$Genus, mean)
dr2<-tapply(DDD$number.of.accepted.species, DDD$Genus, mean)
dr3<-tapply(ecdf, DDD$Genus, mean)

du1<-tapply(DDD$`Maximum_dispersal_distance_analysis_(m)`, DDD$Genus, mean)
du2<-tapply(DDD$number.of.accepted.and.unresolved.species, DDD$Genus, mean)

plot(dr2 ~ dr3)
abline(lm(dr2 ~ dr1))
plot(lm(dr2 ~ dr1))
mdl<-lm(dr2 ~ dr1) ### taraxacum, euphorbia, hieracium +4 more to identify()
summary(mdl) ####

plot(du2 ~ du1)
abline(lm(du2 ~ du1))
plot(lm(du2 ~ du1))
mudl<-lm(du2 ~ du1) ##Hieracium ##Rubus #Taraxacum ##Senecio
summary(mudl)
identify(du1)

###########################################################
TV<-DDD
as.matrix(TV)
rownames(TV)<-TV[,1]
as.data.frame(TV)
rownames(TV)<-TV$species
column_to_rownames(as.data.frame(TV$species))
View(TV)

library(readr)
TV <- read_csv("TV.csv")
View(TV)

library("brranching")
tr <- phylomatic(TV$species, taxnames = TRUE, get = "POST", informat = "newick",
                 method = "phylomatic", storedtree = "zanne2014", treeuri = NULL,
                 taxaformat = "slashpath", outformat = "newick", clean = TRUE,
                 db = "apg", verbose = TRUE)
TVtree <- tr
tr <- TVtree
tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores=TRUE)
plot(TVtree)
plot(tr)
is.binary.tree(tr)#use next line if FALSE
#binarytree <- multi2di(tr)
row.names(tr) <- TV$species
#matches <- match(TV$species, tr$tip.label, nomatch = 0)
#TV <- subset(TV, matches != 0)
#TV
#pruned.tree<-drop.tip(tr,tr$tip.label[-match(TV$species, tr$tip.label)])
#pruned.tree<-drop.tip(tr, tr$tip.label[-na.omit(match(TV$species, tr$tip.label))])
as.matrix(TV)
rownames(TV)<-TV[,2]
as.data.frame(TV)
rownames(TV)<-TV$species
column_to_rownames(as.data.frame(TV$species))
View(TV)

library(readr)
TV <- read_csv("TV.csv")
View(TV)

library("brranching")
tr <- phylomatic(TV$species, taxnames = TRUE, get = "POST", informat = "newick",
                 method = "phylomatic", storedtree = "zanne2014", treeuri = NULL,
                 taxaformat = "slashpath", outformat = "newick", clean = TRUE,
                 db = "apg", verbose = TRUE)
TVtree <- tr
tr <- TVtree
tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores=TRUE)
plot(TVtree)
plot(tr)
is.binary.tree(tr)#use next line if FALSE
#binarytree <- multi2di(tr)
row.names(tr) <- TV$species
#matches <- match(TV$species, tr$tip.label, nomatch = 0)
#TV <- subset(TV, matches != 0)
#TV
#pruned.tree<-drop.tip(tr,tr$tip.label[-match(TV$species, tr$tip.label)])
#pruned.tree<-drop.tip(tr, tr$tip.label[-na.omit(match(TV$species, tr$tip.label))])
name.check(tr, TV)
tr<-pruned.tree
TV <- TV[tr$tip.label, ]
TV <- TV[match(tr$tip.label,rownames(TV)),]
setdiff(tr$tip.label, TV$species)
setdiff(TV$species, tr$tip.label)
TV$species
tr$tip.label<-gsub("\\_", " ", tr$tip.label)

tr2 <- drop.tip(tr, setdiff(tr$tip.label,TV$species))
tr<-tr2
tr$tip.label<- capitalize(tr$tip.label)
tr$tip.label<-sub("^Achillea erba-rottas.moschata$", "Achillea erba-rotta s. moschata", tr$tip.label)
tr$tip.label<-sub("^Alnus viridissubsp.crispa$", "Alnus viridis subsp. crispa", tr$tip.label)
tr$tip.label<-sub("^Arabis hirsutasubsp.hirsuta$", "Achillea erba-rotta subsp. moschata", tr$tip.label)
tr$tip.label<-sub("^Cnidoscolus urenssubsp.stimulosus$", "Cnidoscolus urens subsp. stimulosus", tr$tip.label)
tr$tip.label<-sub("^Viola pubescenssubsp.scabriuscula$", "Viola pubescens subsp. scabriuscula", tr$tip.label)
tr$tip.label<-sub("^Viola sororiasubsp.sororia$", "Viola sororia subsp. sororia", tr$tip.label)
tr$tip.label<-sub("^Ficus crassirameasubsp.stupenda$", "Ficus crassiramea subsp. stupenda", tr$tip.label)
tr$tip.label<-sub("^Penstemon laevigatussubsp.digitalis$", "Penstemon laevigatus subsp. digitalis", tr$tip.label)
tr$tip.label<-sub("^Picris angustifoliasubsp.angustiolia$", "Picris angustifolia subsp. angustiolia", tr$tip.label)
tr$tip.label<-sub("^Erigeron acrissubsp.angulosus$", "Erigeron acris subsp. angulosus", tr$tip.label)
tr$tip.label<-sub("^Tripolium pannonicumsubsp.tripolium$", "Tripolium pannonicum subsp. tripolium", tr$tip.label)
tr$tip.label<-sub("^Hieracium murorumaggr.$", "Hieracium murorum aggr.", tr$tip.label)
tr$tip.label<-sub("^Halesia tetrapterasubsp.monticola$", "Halesia tetraptera subsp. monticola", tr$tip.label)
tr$tip.label<-sub("^Atriplex prostratasubsp.calotheca$", "Atriplex prostrata subsp. calotheca", tr$tip.label)
tr$tip.label<-sub("^Physalis longifoliasubsp.subglabrata$", "Physalis longifolia subsp. subglabrata", tr$tip.label)
tr$tip.label<-sub("^Plantago patagonicasubsp.aristata$", "Plantago patagonica subsp. aristata", tr$tip.label)
tr$tip.label<-sub("^Solidago virgaureasubsp.alpestris$", "Solidago virgaurea subsp. alpestris", tr$tip.label)
tr$tip.label<-sub("^Tilia platyphyllossubsp.cordifolia$", "Tilia platyphyllos subsp. cordifolia", tr$tip.label)
tr$tip.label<-sub("^Arceuthobium vaginatumsubsp.cryptopodum$", "Arceuthobium vaginatum subsp. cryptopodum", tr$tip.label)
tr$tip.label<-sub("^Heterotheca subaxillarissubsp.latifolia$", "Heterotheca subaxillaris subsp. latifolia", tr$tip.label)
tr$tip.label<-sub("^Silene latifoliasubsp.alba$", "Silene latifolia subsp. alba", tr$tip.label)

tr$tip.label<-sub("Achillea erba-rottas.moschata", "Achillea erba-rotta s. moschata", tr$tip.label)
tr$tip.label<-sub("Armeria maritimas.alpina", "Armeria maritima s. alpina", tr$tip.label)
tr$tip.label<-sub("^Senecio aquaticuss.barbareifolius$", "Senecio aquaticus s. barbareifolius", tr$tip.label)
tr$tip.label<-sub("^Pulsatilla vulgariss.grandis$", "Pulsatilla vulgaris s. grandis", tr$tip.label)
tr$tip.label<-sub("^Taraxacum laevigatumag.$", "Taraxacum laevigatum ag.", tr$tip.label)
tr$tip.label<-sub("^Taraxacum sec.erythrosperma$", "Taraxacum Sec. Erythrosperma", tr$tip.label)
tr$tip.label<-sub("^Taraxacum sec.ruderalia$", "Taraxacum Sec. Ruderalia", tr$tip.label)
tr$tip.label<-sub("^Spartina xtownsendii$", "Spartina x townsendii", tr$tip.label)
tr$tip.label<-sub("^Tragopogon pratensiss.orientalis$", "Tragopogon pratensis s. orientalis", tr$tip.label)
tr$tip.label<-sub("^Tragopogon pratensiss.pratensis$", "Tragopogon pratensis s. pratensis", tr$tip.label)
tr$tip.label<-sub("^Vicia sativas.nigra$", "Vicia sativa s. nigra", tr$tip.label)
tr$tip.label<-sub("^Thymus glabrescenss.decipiens$", "Thymus glabrescens s. decipiens", tr$tip.label)
tr$tip.label<-sub("^Odontites vernas.serotina$", "Odontites verna s. serotina", tr$tip.label)
tr$tip.label<-sub("^Leontodon taraxacoidessubsp.taraxacoides$", "Leontodon taraxacoides subsp. taraxacoides", tr$tip.label)
tr$tip.label<-sub("^Centaurea phrygiasubsp.pseudophrygia$", "Centaurea phrygia subsp. pseudophrygia", tr$tip.label)
tr$tip.label<-sub("^Medicago sativas.falcata$", "Medicago sativa s. falcata", tr$tip.label)
tr$tip.label<-sub("^Molinia caeruleas.arundinacea$", "Molinia caerulea s. arundinacea", tr$tip.label)
TV

tr$tip.label<-sub("Achillea erba-rottas.moschata", "Achillea erba-rotta s. moschata", tr$tip.label)
matches <- match(TV$species, tr$tip.label, nomatch = 0)
TVX <- subset(TV, matches != 0)
TV<-TVX

trx<-TV[!duplicated(TV$species), ]
TV<-trx
duplicated(TV$species)

read

source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
biocLite("treeio")

pruned.tree<-drop.tip(tr,tr$tip.label[-match(TV$species, tr$tip.label)])
tr<-pruned.tree

install.packages('caper', repos='http://r-forge.r-project.org')
library("caper")
tr<-makeLabel(tr)
tr$node.label
tr$node.label <- with(tr, ((max(edge)-Nnode) +1):max(edge))
TRV <- comparative.data(tr, TV, species, vcv=T, na.omit=F,
                        force.root=FALSE, warn.dropped=T, scope=NULL)


TV$velocity<-log10(TV$velocity)
as.data.frame(TV)

TRV$dropped$tips
TRV$dropped$unmatched.rows

model.pgls <- pgls(log(TV$`Maximum_dispersal_distance_analysis_(m)`) ~ log(TV$`Seed_release_height_(m)`),
                   data = TRV, lambda = "ML")
plot(log(TV$`Maximum_dispersal_distance_analysis_(m)`) ~ log(TV$`Seed_release_height_(m)`), col=(as.numeric(as.factor(TV$Growth_form))+1), pch=c(as.factor(TV$Dispersal_syndrome))+14, data=TRV$data)
abline(model.pgls)

model.pglsB <- pgls((TV$`Seed_terminal_velocity_(m/s)`) ~ log(TV$`Seed_release_height_(m)`),
                    data = TRV, lambda = "ML")
plot(log(TV$`Seed_terminal_velocity_(m/s)`) ~ log(TV$`Seed_release_height_(m)`), col=(as.numeric(as.factor(TV$Growth_form))+1), pch=c(as.factor(TV$Dispersal_syndrome))+14, data=TRV$data)
abline(model.pglsB)
lambda.profile <- pgls.profile(model.pgls, "lambda")
plot(lambda.profile)
par(mfrow = c(2, 2))
plot(model.pgls)
par(mfrow = c(1, 1))
plot(log(TV$`Maximum_dispersal_distance_analysis_(m)`) ~ log(TV$`Seed_terminal_velocity_(m/s)`), data=TRV$data)

model.pgls2 <- pgls(log(TV$`Maximum_dispersal_distance_analysis_(m)`) ~ log(TV$`Seed_terminal_velocity_(m/s)`),
                    data = TRV, lambda = "ML")
plot(log(TV$`Maximum_dispersal_distance_analysis_(m)`) ~ log(TV$`Seed_terminal_velocity_(m/s)`), col=(as.numeric(as.factor(TV$Growth_form))+1), pch=c(as.factor(TV$Dispersal_syndrome))+14, data=TRV$data)
abline(model.pgls2)
lambda.profile <- pgls.profile(model.pgls2, "lambda")
plot(lambda.profile)
par(mfrow = c(2, 2))
plot(model.pgls2)
summary(model.pgls2)

lambda.profile <- pgls.profile(model.pgls2, "lambda")
plot(lambda.profile)
par(mfrow = c(2, 2))
plot(model.pgls2)

plot(log(TV$`Maximum_dispersal_distance_analysis_(m)`) ~ log(TV$number.of.accepted.species), pch=c(as.factor(TV$Growth_form)), col=c(as.factor(TV$Dispersal_syndrome)), data=TRV$data)
plot(log(TV$number.of.accepted.species)~log(TV$`Seed_terminal_velocity_(m/s)`), pch=c(as.factor(TV$Growth_form)), col=c(as.factor(TV$Dispersal_syndrome)), data=TRV$data)


TV2 <- subset(TV, Dispersal_syndrome == "wind (special)"|"wind (none)")
View(TV2)
tr2 <- drop.tip(tr, setdiff(tr$tip.label,TV2$species))
TRV2 <- comparative.data(tr2, TV2, species, vcv=T, na.omit=F,
                         force.root=FALSE, warn.dropped=T, scope=NULL)

plot(log(TV2$`Maximum_dispersal_distance_analysis_(m)`) ~ log(TV2$number.of.accepted.and.unresolved.species), pch=c(as.factor(TV2$Growth_form)), col=c(as.factor(TV2$Dispersal_syndrome)), data=TRV2$data)
plot(log(TV2$number.of.accepted.and.unresolved.species)~log(TV2$`Seed_terminal_velocity_(m/s)`), pch=c(as.factor(TV2$Growth_form)), col=c(as.factor(TV2$Dispersal_syndrome)), data=TRV2$data)

library(phytools)
contMap(tr, TV$velocity)
fastAnc(tr, TV$velocity)

model.pgls <- pgls(log(TV$`Maximum_dispersal_distance_analysis_(m)`) ~ log(TV$`Seed_terminal_velocity_(m/s)`),
                   data = primate, lambda = "ML")

PICreconstruction <- ace(TV$`Maximum_dispersal_distance_analysis_(m)`, tr, type="continuous", method="pic")
PICreconstruction
traitvaluefinal <- c(TV$trait_value, PICreconstruction$ace)
traitvaluefinal
plot(tr, show.tip.label=FALSE)
tiplabels(pch = 21, cex=(TV$velocity))
nodelabels(pch = 21, cex=(PICreconstruction$ace)*1000)

BMfit <- fitContinuous(tr, TVELO, model="BM")
FBM <-fastBM(tr)
aa <- fastAnc(tr, TVELO)

TVELO <- as.vector(as.numeric(TV$`Seed_terminal_velocity_(m/s)`))
names(TVELO) <- as.vector(TV$species)
TVELO
Anc <- fastAnc(tr, TVELO)
Anc
tr$tip.label

##create contMap
XX <- contMap(tr, TVELO, fsize = 0.1, legend=FALSE)
plot(XX, fsize = c(0.4, 0.5))
plot(XX, type = "fan", fsize = c(0.1, 0.5))
add.color.bar(40,obj$cols,title="Seed Terminal Velocity (m/s)",
              lims=obj$lims,digits=3,prompt=FALSE,x=0,
              y=1-0.08*(Ntip(obj$tree)-1),lwd=4,fsize=1,subtitle="")

phylosig(tr,TVV,test=TRUE) #$K[1] 0.0119373, $P[1] 0.035
lambda <- phylosig(tr,TVV,test=TRUE, method = "lambda") #$lambda[1] 0.7723352, $logL[1] -884.8534, $logL0 [1] -958.095, $P[1] 1.017874e-33
lambda
lambda.profile <- pgls.profile(model.pgls, "lambda")
plot(lambda.profile)

library(diversitree)

tr2 <- force.ultrametric(tr2)
tr <-multi2di(tr)
tr <- force.ultrametric(tr)
TVV <- as.vector(TV$velocity)
names(TVV) <- as.vector(TV$species)
make.quasse(tr, TVV, 0.93, 0.27,sampling.f=(1020/250000))
tree.quasse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE, x0=NA,
            single.lineage=TRUE, verbose=FALSE)

p<-starting.point.quasse(tr,TVV)
set.seed(1)
phy <- tree.quasse(p, max.taxa=2000, x0=0,
                   single.lineage=FALSE)
TVV
tr$tip.label

citation("brranching")

TVV<-log10(TVV)
mass.sd <- 0.05
xr <- range(TVV) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
make.primates <- function(lambda, mu)
  make.quasse(tr, TVV, mass.sd, lambda, mu)
f.c <- make.primates(constant.x, constant.x)
f.l <- make.primates(linear.x, constant.x)
f.s <- make.primates(sigmoid.x, constant.x)
f.h <- make.primates(noroptimal.x, constant.x)
f.h

control <- list(parscale=.1, reltol=0.001)
nodrift <- function(f)
  constrain(f, drift ~ 0)
mle.c <- find.mle(nodrift(f.c), p, lower=0, control=control, verbose=0)
p.c <- mle.c$par
p.l <- c(p.c[1], l.m=0, p.c[2:3])
p.s <- p.h <- c(p.c[1], p.c[1], mean(xr), 1, p.c[2:3])
names(p.s) <- argnames(nodrift(f.s))
names(p.h) <- argnames(nodrift(f.h))
mle.l <- find.mle(nodrift(f.l), p.l, control=control, verbose=0)
mle.s <- find.mle(nodrift(f.s), p.s, control=control, verbose=0)
mle.h <- find.mle(nodrift(f.h), p.h, control=control, verbose=0)
c(linear=coef(mle.l),
  sigmoidal=coef(mle.s),
  hump=coef(mle.h))
anova(mle.c, linear=mle.l, sigmoidal=mle.s, hump=mle.h)
mle.d.l <- find.mle(f.l, coef(mle.l, TRUE), control=control, verbose=0)
mle.d.s <- find.mle(f.s, coef(mle.s, TRUE), control=control, verbose=0)
mle.d.h <- find.mle(f.h, coef(mle.h, TRUE), control=control, verbose=0)
c(linear=coef(mle.d.l)[["drift"]],
  sigmoidal=coef(mle.d.s)[["drift"]],
  hump=coef(mle.d.h)[["drift"]])
anova(mle.c, linear=mle.l, sigmoidal=mle.s, hump=mle.h,
      drift.linear=mle.d.l, drift.sigmoidal=mle.d.s,
      drift.hump=mle.d.h)

coef(mle.d.h)
#       l.y0        l.y1      l.xmid        l.s2         m.c       drift   diffusion 
#0.001804145 0.196609048 0.298228185 0.438370075 0.161860121 0.000964086 0.007222570 

nodelabels(tr)

plot(tr, cex=0.049, edge.width = 0.5, label.offset = 1)
plot(mytree, cex=0.049, edge.width = 0.5, label.offset = 1)

identify.phylo(tr) #1763
getDescendants(tr, 1763)
identify(mytree, tips= T, labels = T)

sample.f
f.cc <- make.quasse.split(tr, TVV, mass.sd, constant.x, constant.x,
                          1763, Inf, sampling.f=(135/24000))
argnames(f.cc)
g.cc <- constrain(f.cc, drift.1 ~ 0, drift.2 ~ 0,
                  diffusion.2 ~ diffusion.1)
argnames(g.cc)
p.cc <- c(p.c, p.c[1:2])
names(p.cc) <- argnames(g.cc)
mle.c$lnLik - g.cc(p.cc)

mle.cc <- find.mle(g.cc, p.cc, control=control, lower=0, verbose=0)

f.ll <- make.quasse.split(tr, TVV, mass.sd, linear.x, constant.x,
                          1763, Inf,sampling.f=(135/24000))
g.ll <- constrain(f.ll, drift.1 ~ 0, drift.2 ~ 0,
                  diffusion.2 ~ diffusion.1)
g.lc <- constrain(g.ll, l.m.2 ~ 0)
g.cl <- constrain(g.ll, l.m.1 ~ 0)
#Generate a starting points: start with the function where both speciation rates are linear functions.
p.cc <- coef(mle.cc)
p.ll <- c(p.cc[1], 0, p.cc[2:4], 0, p.cc[5])
names(p.ll) <- argnames(g.ll)
#Run the ML searches for this model:
mle.ll <- find.mle(g.ll, p.ll, control=control, verbose=0)
#Then generate starting points for models with just one of the sections of the tree having a linear speciation function:
p.lc <- c(coef(mle.ll)[1:3], p.ll[c(4, 5, 7)])
p.cl <- c(p.ll[c(1, 3, 4)], coef(mle.ll)[5:7])
#and run the ML search:
mle.lc <- find.mle(g.lc, p.lc, control=control, verbose=0)
mle.cl <- find.mle(g.cl, p.cl, control=control, verbose=0)
#We can then compare the models again:
anova(mle.c, linear=mle.l, sigmoidal=mle.s, hump=mle.h,
      part.constant=mle.cc,
      part.linear.bg=mle.lc,
      part.linear.fg=mle.cl,
      part.linear=mle.ll)

plot.hisse.states(trA, mle.ll)

mle.nd.lc <- find.mle(nodrift(g.lc), coef(mle.lc, TRUE), control=control, verbose=0)
mle.nd.cl <- find.mle(nodrift(g.cl), coef(mle.cl, TRUE), control=control, verbose=0) 
mle.nd.ll <- find.mle(nodrift(g.ll), coef(mle.ll, TRUE), control=control, verbose=0)  

coef(mle.ll)
#   l.c.1        l.m.1        m.c.1  diffusion.1        l.c.2        l.m.2        m.c.2 
#25.869539624  0.004035947 25.859026135  0.005784569 60.264371137 -0.029294614 60.213903458 
########################################ESsim###########################################
essim <- function(phy, trait, nsim = 1000, is) {
  
  require(ape)
  require(mvtnorm)
  
  if(missing(is)) { # If inverse equal splits statistics not provided, calculate it
    rootnode <- length(phy$tip.label) + 1
    is <- numeric(length(phy$tip.label))
    for (i in 1:length(is)){
      node <- i
      index <- 1
      qx <- 0
      while (node != rootnode){
        el <- phy$edge.length[phy$edge[,2] == node]
        node <- phy$edge[,1][phy$edge[,2] == node]			
        qx <- qx + el* (1 / 2^(index-1))			
        index <- index + 1
      }
      is[i] <- 1/qx
    }		
    names(is) <- phy$tip.label
  }
  
  is <- log(is[phy$tip.label]) # log transform
  trait <- trait[phy$tip.label]
  
  # Pearson's correlation between splits statistic and trait
  res <- cor.test(is, trait, method="pearson")
  
  # Fit Brownian motion model to get diffusion rate and root state estimates
  vv <- vcv.phylo(as.phylo(phy))
  onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
  root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
  rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))
  
  # Brownian simulations 
  sims <- t(rmvnorm(nsim, sigma=rate*vv))
  rownames(sims) <- rownames(vv)
  
  # Pearson's correlations of simulated datasets
  sim.r <- sapply(1:nsim, function(x) cor.test(is[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
  
  # Calculate the two-tailed p value
  corr <- res$estimate
  upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
  lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
  pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed
  
  result <- as.vector(c(corr, pval))
  names(result) <- c("rho", "P Value")
  return(result)
  
}


ES-sim/R/essim.R 
##es <- log(es[phy$tip.label]) # log transform
trait <- trait[phy$tip.label]
write.tree(tra, file = "tra.tre")
trA<-read.tree("tra.tre")
write.tree(tr, file = "trp.tre")
trp<-read.tree("trp.tre")
library(rncl)
plot(trp)
trp<-as.phylo(trp)
trp$tip.label
TV$species
esim<-essim(trp, TVV, nsim= 10000)  #Named num [1:2] -0.0164 0.9535 - attr(*, "names")= chr [1:2] "rho" "P Value"
esima<-essim(trA, TVVA, nsim= 10000)
tra$tip.label
TVVa$rowname
TVVa<-as.data.frame(TVV)
TVVa<-tibble::rownames_to_column(TVVa)
row.names(TVVa) <- TVVa$rowname
matchesa <- match(TVVa$rowname, tra$tip.label, nomatch = 0)
TVVa <- subset(TVVa, matchesa != 0)
setdiff(trA$tip.label, TVVa$rowname)
setdiff(TVVa$rowname, trA$tip.label)
View(TVVA)
esim # rho = -0.01644765; P Value =  0.95350465
esima # rho = 0.0394; P Value = 0.8921
TVVA <- as.vector(as.numeric(TVVa$TVV))
names(TVVA) <- as.vector(TVVa$rowname)
tra$tip.label

