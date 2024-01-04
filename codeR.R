Basic data checks and setup
dim(x)
names(x)
attach(x)
head(x)
str(x)
summary(x)

# Plotting geographical data
plot(longitude, latitude)
identify(longitude, latitude)  # In Rstudio, press ESC to display the data after selection

# Extracting a table with species data
z <- x[1:65, 11:139]
dim(z)

# Separating data by slope
x <- x[1:65, 1:139]
slope <- x$versant
south <- subset(x, slope == "S")
north <- subset(x, slope == "N")
dim(south)
dim(north)

# Tables with species data for each slope
zsud <- south[, 11:ncol(south)]
znord <- north[, 11:ncol(north)]
dim(zsud)
dim(znord)

library(ade4)

# Graphical representation of species data
table.value(z, clegend=0)  # Entire dataset
table.value(z[order(slope),], row.labels=slope[order(slope)], clegend=0)  # Ordered by slope
table.value(zsud, clegend=0)  # South slope
table.value(znord, clegend=0)  # North slope

# 2. Exhaustivity of sampling

library(vegan)
# Specaccum for the total area
sp2 <- specaccum(z, "random", permutations=65)
plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
round(sp2$richness[65], 2)  # Total richness (gamma)
round(sp2$sd[65], 2)  # Standard deviation

# Specaccum for southern slope
sp3 <- specaccum(zsud, "random", permutations=42)
plot(sp3, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
round(sp3$richness[17], 2)  # Rarified richness at n=17 quadra
round(sp3$sd[17], 2)

# Specaccum for northern slope
sp4 <- specaccum(znord, "random", permutations=30)
plot(sp4, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
round(sp4$richness[17], 2)
round(sp4$sd[17], 2)

par(mfrow=c(1, 3))
plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
plot(sp3, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
plot(sp4, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

# Frequency rank diagram
library(BiodiversityR)
x$versant <- as.factor(x$versant)
slope <- x$versant
rankabuncomp(z, y=x, factor="slope", scale="proportion", legend=TRUE)
rankabuncomp(z, y=x, factor="slope", scale="abundance", legend=TRUE)
rankabuncomp(z, y=x, factor="slope", scale="logabun", legend=TRUE)

# Models fitting data with radfit()
nzsud <- apply(zsud, 2, sum)
nznord <- apply(znord, 2, sum)
radfit(nzsud)
plot(radfit(nzsud))
radfit(nznord)
plot(radfit(nznord))

# 3. Diversity indices

# 3.1) Calculation of indices (alpha)
sr <- specnumber(z, MARGIN=1)  # Specific richness S
N <- apply(z, 1, sum)  # Total number of individuals per quadra
Nmax <- apply(z, 1, max)  # Abundance of dominant species per quadra
bg <- Nmax / N  # Berger-Parker index
p <- z / N  # Proportion p of each species per sample
p <- as.matrix(p)
simpc <- (N / (N - 1)) * diversity(z, "simpson")  # Simpson's index
esimp <- diversity(z, "simpson") / (1 - (1 / sr))  # Simpson Evenness

# Functional diversity
trait <- read.table(file.choose(), sep=",", dec=".", header=TRUE, row.names=1)
trait1 <- aggregate(trait$chlo ~ species

# CWM (Community Weighted Mean) index
require(weimea)
icwm <- cwm(z, trait1)

# Functional distances dij between species measured from the two traits
library(FD)
d <- gowdis(trait1)

# Functional classification of species (complementarity-redundancy) with CAH (mean linkage) on Gower
plot(hclust(d, "average"), hang=-1)

# Quadratic entropy (extension of 1-D with dij)
library(adiv)
qe <- QE(z, d)

# Put the 6 indices in the same matrix

# Make a Draftsman plot to study their relationship (i.e., empirical redundancy or not)

# (Continue with Draftsman plot code here)

# 3.2) Effect of explanatory variables

# Create 4 sub figures (i.e., 4 indices, or 6 with functional indices) for the 3 following effects.
# For tests, recover DESINF and EVA EU scripts:

# 3.2.1) Effect of slope
# boxplot() slope + t.test() (check hypotheses and orient the test more precisely)
par(mfrow=c(2,2))
boxplot(sr ~ slope)
boxplot(bg ~ slope)
boxplot(simpc ~ slope)
boxplot(esimp ~ slope)

# Mean and sd data for each modality of a qualitative variable
round(tapply(sr, slope, mean), 2)
round(tapply(sr, slope, sd), 2)
round(tapply(bg, slope, mean), 2)
round(tapply(bg, slope, sd), 2)
round(tapply(simpc, slope, mean), 2)
round(tapply(simpc, slope, sd), 2)
round(tapply(esimp, slope, mean), 2)
round(tapply(esimp, slope, sd), 2)

# Statistical test of comparison of 2 means

# (Continue with statistical test code here)

# 3.2.2) Effect of tree cover
# plot() vs tree cover + lm and correlation if linear link
plot(sr ~ rec_espece_dominant_arboree)
scatter.smooth(sr ~ rec_espece_dominant_arboree, span=2/3, degree=2)
# Also to be done for green oak and Aleppo pine overlays


# 3.2.3) Effect of cover compared on the south and north slopes

# Indices vs. cover according to slope -> ANCOVA
# Example ANCOVA code
cover <- x$cover # assuming 'cover' is a variable in your dataset
ancova_results <- aov(sr ~ cover * slope, data=x)
summary(ancova_results)

# 3.2.4) Effect of cover on dominant trees

# Boxplot S vs. dominant trees + anova()
dominant_trees <- x$dominant_trees # assuming this is a variable in your dataset
boxplot(sr ~ dominant_trees, data=x)
anova_results <- aov(sr ~ dominant_trees, data=x)
summary(anova_results)

# How is alpha diversity structured according to the explanatory variables?
# Analysis here (depending on specific questions and dataset structure)

# 3.2.1.1) Overall
# Gamma decomposition into beta and average alpha to determine beta for the area, by slope
# This requires specific functions or packages that might not be readily available in base R.
# Here's a conceptual placeholder for this analysis
# gamma_decomposition_results <- some_function(z, slope)

# 3.2) Beta diversity

# 3.2.1) Taxonomic beta diversity

# 3.2.1.1) Global
# As with gamma decomposition, this might require a specific function or package
# beta_diversity_results <- another_function(z, slope)

# 3.2.1.2) Inter-quadra
# Presence-absence data: Jaccard

# Jaccard-PCoA (MDS)
jac <- vegdist(z, "jaccard")
pcbc <- cmdscale(jac, eig=TRUE, k=2) # Classical MDS using cmdscale

# Plotting the results
plot(pcbc$points, col=slope, pch=19, xlab="Axis 1", ylab="Axis 2")
legend("topright", legend=levels(slope), col=1:length(levels(slope)), pch=19)

# 3.2.2) Functional beta diversity

# 3.2.2.1) Global
# Assuming 'd' is a distance matrix calculated earlier
discomQE <- QE(z, d, structures=versant)
raoD_results <- raoD(z, d)

# 4) Sampling bias
# Boxplot RS on 3 surveys for different areas depending on the morning/afternoon group
# Assuming 'time_of_day' is a variable in your dataset representing morning/afternoon
boxplot(sr ~ time_of_day, data=x)