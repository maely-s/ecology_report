#####################################################
# HAV701B Ecology: concepts, tools and applications 
#####################################################

############################################################################################
# Field workshops - Plant community diversity in the north of Montpellier
# Lecturers : Guillaume Papuga & Bastien Mérigot  
# TP November 2022
############################################################################################

# Load the community data with the read.table() function

x<-read.table(file.choose(), sep=";", dec=",", header=T)

# 1.	Checking the loaded data

edit(x) # if manual changes are needed: x<-edit(x)

### close the window before using the console again

dim(x)

names(x)
attach(x)
head(x)
structure(x)

summary(x)

plot(longitude,latitude)
identify(longitude,latitude)
# in Rstudio press esc to display the data after selecting it on the graph

# Table with species only 
z<-x[1:65,11:139]
dim(z)

# Create two separate tables for each slope with the subset() function

x<-x[1:65,1:139]

slope<-versant
south<-subset(x, slope=="S")
north<-subset(x, slope=="N")

dim(south)
dim(north)

# Table with species only 
zsud<-south[,11:ncol(south)]
znord<-north[,11:ncol(north)]

dim(zsud)
dim(znord)


require(ade4)

# Graphical representation of the data table 
table.value(z, clegend=0) 
# Ditto ordered with respect to a factor variable
slope<-x$versant
table.value(z[order( slope),], row.labels =slope[order( slope)], clegend=0 ) 

#ditto by sub-array of slope data
table.value(zsud, clegend=0 ) 
table.value(znord, clegend=0 ) 


# 2.	Exhaustivity of sampling



# Make 1 curve for the total area + 1 curve per slope; method without replacement (number of permutations must be equal to the number of samples)
require(vegan)
sp2 <- specaccum(z, "random", permutations=65)
plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

# Total richness (gamma)
round(sp2$richness[65],2)
round(sp2$sd[65],2)

rsud<-sud[,11:ncol(sud)]
sp3 <- specaccum(zsud, "random", permutations=42)
plot(sp3, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

# Rarified richness at n=17 quadra
round(sp3$richness[17],2)
round(sp3$sd[17],2)

rnord<-nord[,11:ncol(nord)]

sp4 <- specaccum(znord, "random", permutations=30)
plot(sp4, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

round(sp4$richness[17],2)
round(sp4$sd[17],2)

par(mfrow=c(1,3))
plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
plot(sp3, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
plot(sp4, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")


# Was the sampling sufficiently comprehensive for the method used?


# Frequency rank diagram

require(BiodiversityR)
x$versant<-as.factor(x$versant)
slope<-x$versant

# or versant

rankabuncomp(z, y=x, factor="slope", scale="proportion", legend=T)
rankabuncomp(z, y=x, factor="slope", scale="abundance", legend=T)
rankabuncomp(z, y=x, factor="slope", scale="logabun", legend=T)

# Models fitting data with radfit()

nzsud<-apply(zsud,2,sum)
nznord<-apply(znord,2,sum)

radfit(nzsouth)
plot(radfit(nzsud))

radfit(nznord)
plot(radfit(nznordd))


#3. Diversity indices 

# 3.1) Calculation of indices (alpha)

#Specific richness S
sr<-specnumber(z,MARGIN=1)

# Total number of individuals N, and abundance of dominant species Nmax, per quadra
N<-apply(z,1,sum)
Nmax<-apply(z,1, max)

# Berger-Parker index
bg<-Nmax/N

# proportion p of each species per sample
p<-z/N
p<-as.matrix(p)

# Simpson's heterogeneous (RS+evenness) index (unbiased) PIE: (n/n-1) x 1-D
simpc <- (N/(N-1))*diversity(z, "simpson")

# Simpson 1-D / (1- 1/S) Evenness
esimp<-diversity(z, "simpson")/(1-(1/sr))

# Functional diversity

# Load functional trait data

trait<-read.table(file.choose(), sep= , dec=, header=, row.names= )

# average the quantitative trait value for the 3 individuals measured from one species (use correct variable names)
trait1=aggregate(trait$chlo ~ species, data = trait, mean)

# CWM (Community Weighted Mean) index
require(weimea)
icwm<-cwm(z,trait1)

# Functional distances dij between species measured from the two traits
require(FD)
d<-gowdis(trait1)

#Functional classification of species (complementarity-redundancy) with CAH (mean linkage) on Gower
plot(hclust(d, "average"), hang=-1) 

# Quadratic entropy (extension of 1-D with dij)
require(adiv)
qe<-QE(z,d)

# Put the 6 indices in the same matrix 
indices_matrix <- cbind(sr, bg, simpc, esimp, qe, icwm)


# Make a Draftsman plot to study their relationship (ie empirical redundancy or not)

your_plotting_function(indices_matrix)  #replace 'your_plotting_function' with the actual plotting function

#3.2) Effect of explanatory variables

# Create 4 sub figures (ie 4 indices, or 6 with functional indices) for the 3 following effects. For tests recover DESINF and EVA EU scripts):

# 3.2.1) Effect of slope
# boxplot() slope + t.test() (check hypotheses and orient the test more precisely)

par(mfrow=c(2,2))

boxplot(sr~slope)
your_statistical_test_function(sr, slope) #(replace 'your_plotting_function' and 'your_statistical_test_function' with actual functions)

boxplot(bg~slope)
your_statistical_test_function(bg, slope)

boxplot(simpc~slope)
your_statistical_test_function(simpc, slope)

boxplot(esimp~slope)
your_statistical_test_function(esimp, slope)

# Mean and sd data for each modality of a qualitative variable
round(tapply(sr,slope, mean),2)
round(tapply(sr,slope, sd), 2)

round(tapply(bg,slope, mean),2)
round(tapply(bg,slope, sd), 2)

round(tapply(simpc,slope, mean),2)
round(tapply(simpc,slope, sd), 2)

round(tapply(esimp,slope, mean),2)
round(tapply(esimp,slope, sd), 2)

# Statistical test of comparison of 2 means (cf script UEs DESINF and/or EVA)


# 3.2.2) Effect of tree cover
#plot() vs tree cover + lm and correlation if linear link

plot(sr ~ rec_espece_dominant_arboree)
scatter.smooth(sr ~ rec_espece_dominant_arboree, span = 2/3, degree = 2)
correlation_result <- cor.test(sr, rec_espece_dominant_arboree)
correlation_result

# Also to be done for green oak and Aleppo pine overlays


# 3.2.3) Effect of cover compared on the south and north slopes

#Indices vs. cover according to slope -> ANCOVA
#etc.

# 3.2.4) Effect of cover on dominant trees

#boxplot S vs. dominant trees + anova() (check hypotheses and orient the test more precisely)
#etc.



#How is alpha diversity structured according to the explanatory variables?




#3.2.1.1) Overall
#Gamma decomposition into beta and average alpha to determine beta for the area, by slope, between each of the 3 replicates for a given group.
#Does the overall species composition within the study area, by slope

gamma_result <- gamma.diversity(z, method = "multiplicative")



#3.2) Beta diversity

#3.2.1) Taxonomic beta diversity

#3.2.1.1) Global
#Gamma decomposition into beta and average alpha to determine beta, for the area, by slope, between each of the 3 replicates for a given group.
#Is the species composition overall within the study area, by slope

beta_global_result <- beta.diversity(z, method = "multiplicative")

# Perform the calculation with the multiplicative model


#3.2.1.2) Inter-quadra
#Presence-absence data: Jaccard
jac <- vegdist(z, method = "jaccard")

# Jaccard-PCoA (MDS)
require(vegan)
jac<-vegdist(z, "jaccard")

library(plotly)
gom<-as.matrix(jac)
gom
plot_ly(x=colnames(gom), y=rownames(gom), z = gom, type = "heatmap")
plot_ly(x=colnames(gom), y=rownames(gom), z = gom, type = "heatmap", colorscale= "Earth")

#PCoA
pcbc<-dudi.pco(jac)

# percentages axis
perc<-round((pcbc$eig/sum(pcbc$eig))*100,2)
pourc
cumsum(pourc)

# Projections of samples
s.label(pcbc$li, sub="Jaccard")

# Projections of samples according to factorial variables
par(mfrow = c(1,2))
s.class(pcbc$li, versant, col=c(1:4))
s.class(pcbc$li, zone, col=c(1:3))

# A posteriori projection of variables contribution (species correlations to axes)
require(ape)
pcobc=pcoa(jac)
biplot(pcobc, z)

# Overlay data: same as above with Bray-Curtis method="bray".

# 3.2.1.2) Inter-quadra (continued)
# Perform the calculation with the multiplicative model
beta_inter_quadra_result <- beta.diversity(z, method = "multiplicative")


#3.2.2) Functional beta division

#3.2.2.1) Global

# Quadratic entropy 

require(adiv)
discomQE(z, d, structrures=versant) 

# see also
require(picante)
raoD(z, d) 

# 4) Sampling bias
# boxplot RS on 3 surveys for different areas depending on the  morning /afternoon group
# Assuming 'your_data' is the relevant dataset and 'morning_afternoon' is a categorical variable
boxplot(your_data$RS ~ your_data$area * your_data$morning_afternoon,
        xlab = "Area and Time of Day",
        ylab = "Richness (RS)",
        main = "Sampling Bias Analysis",
        col = c("lightblue", "lightgreen"),
        border = "black",
        notch = TRUE)