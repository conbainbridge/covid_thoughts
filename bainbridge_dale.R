########################################
####### Streams of consciousness #######
######## students during covid #########
### Constance Bainbridge & Rick Dale ###
########################################

##########################
### Initialize & clean ###
##########################
library(jsonlite)
library(stringr)
library(glmnet)
library(ggplot2)
library(ggpubr)
library(igraph)
library(scico)

# Import study1 LIWC csv, study2 LIWC csv, combine into one dataset with columns for subject id and condition added in, plus which dataset (fall 2021 or winter 2022)
setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/data", sep=""))
dir = getwd( )
setwd(dir)

fall21_LIWC = read.csv("study1_LIWC_all.csv")
winter22_LIWC = read.csv("study2_LIWC_all.csv")

### Add subject ID, condition, prompt type as variables (fall quarter 2021)
for (i in 1:length(fall21_LIWC$Filename)) {
  # Subject ID
  id = fall21_LIWC$Filename[i]
  id = str_replace(id, "(_.*).txt$", "")
  fall21_LIWC$id[i] = id
  
  # Condition
  if (str_detect(fall21_LIWC$Filename[i], "PrFP", negate = FALSE)) {
    fall21_LIWC$condition[i] = "PrFP"
  } else if (str_detect(fall21_LIWC$Filename[i], "PFPr", negate = FALSE)) {
    fall21_LIWC$condition[i] = "PFPr"
  }
  
  # Prompt type
  if (str_detect(fall21_LIWC$Filename[i], "present", negate = FALSE)) {
    fall21_LIWC$prompt[i] = "present"
  } else if (str_detect(fall21_LIWC$Filename[i], "past", negate = FALSE)) {
    fall21_LIWC$prompt[i] = "past"
  } else if (str_detect(fall21_LIWC$Filename[i], "future", negate = FALSE)) {
    fall21_LIWC$prompt[i] = "future"
  }
}

# Add subject ID, condition, prompt type as variables (winter quarter 2022)
for (i in 1:length(winter22_LIWC$Filename)) {
  # Subject ID
  id = winter22_LIWC$Filename[i]
  id = str_replace(id, "(_.*).txt$", "")
  winter22_LIWC$id[i] = id
  
  # Condition
  if (str_detect(winter22_LIWC$Filename[i], "PrFP", negate = FALSE)) {
    winter22_LIWC$condition[i] = "PrFP"
  } else if (str_detect(winter22_LIWC$Filename[i], "PFPr", negate = FALSE)) {
    winter22_LIWC$condition[i] = "PFPr"
  }
  
  # Prompt type
  if (str_detect(winter22_LIWC$Filename[i], "present", negate = FALSE)) {
    winter22_LIWC$prompt[i] = "present"
  } else if (str_detect(winter22_LIWC$Filename[i], "past", negate = FALSE)) {
    winter22_LIWC$prompt[i] = "past"
  } else if (str_detect(winter22_LIWC$Filename[i], "future", negate = FALSE)) {
    winter22_LIWC$prompt[i] = "future"
  }
}
#####



################
### DF setup ###
################
as.matrix(winter22_LIWC)
as.matrix(fall21_LIWC)
winter22_f_data = winter22_LIWC[winter22_LIWC$prompt == "future",]
fall21_f_data = fall21_LIWC[fall21_LIWC$prompt == "future",]
#####



#################
### WINTER-22 ###
#################
#################
### PCA setup ###
#################
winter_f_pcaData = winter22_f_data[,3:119] # Keeping just LIWC scores
winter_f_pcaData = winter_f_pcaData[,colSums(winter_f_pcaData)>0] # Remove any LIWC categories where none had a value
winter_f_pcaModel = prcomp(scale(winter_f_pcaData))
summary(winter_f_pcaModel)

### Fig. 1, cumulative variance for PCs
png("fig1.png", units="px", width=7000, height=6125, res=1000) # Then converted to 300 dpi in Adobe Photoshop
fig1_title = strwrap("Cumulative variance explained for 91 total Principal Component", width = 40)
plot(cumsum(winter_f_pcaModel$sdev^2)/sum(winter_f_pcaModel$sdev^2),ylim=c(0,1), ylab="Cumulative proportion", xlab="Principal Component", main="Cumulative variance explained for 91 total Principal Components")
abline(v=20,col="blue") 
abline(h=0.69368,col="darkgreen") # Cumulative prop. variance including up to PC20
dev.off()

dim(winter_f_pcaModel$x)

initialComponents = c(1:20)
winter_initial_model = glm(winter22_f_data$condition=='PrFP'~.,
                     data=data.frame(winter_f_pcaModel$x[,initialComponents]),family='binomial')
summary(winter_initial_model)

selectedComponents = c(1,4,5,10,11,13,18) # These are the ones that are significant/approach significance out of the first 20

winter_f_model = glm(winter22_f_data$condition=='PrFP'~.,
               data=data.frame(winter_f_pcaModel$x[,selectedComponents]),family='binomial')
summary(winter_f_model)

### Fix two PCs that have negative slope to simplify interpretation of loading scores
winter_f_pcaModel$rotation[,5]=-1*winter_f_pcaModel$rotation[,5]
winter_f_pcaModel$rotation[,10]=-1*winter_f_pcaModel$rotation[,10]
#####



#################
### PCA stats ###
#################
### Check top variables in each PC
top_indices = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116)

loading_scores_PC1 = winter_f_pcaModel$rotation[,1]
sorted_PC1 = sort(loading_scores_PC1) # Higher loading score, for example health is typical of during pandemic
PCA_top1 = sorted_PC1[top_indices]
var_scores_PC1 = abs(loading_scores_PC1)
var_scores_ranked_PC1 = sort(var_scores_PC1, decreasing=TRUE)
top_10_vars_PC1 = names(var_scores_ranked_PC1[1:10])
top_10_vars_PC1

loading_scores_PC4 = winter_f_pcaModel$rotation[,4]
sorted_PC4 = sort(loading_scores_PC4)
PCA_top4 = sorted_PC4[top_indices]
var_scores_PC4 = abs(loading_scores_PC4)
var_scores_ranked_PC4 = sort(var_scores_PC4, decreasing=TRUE)
top_10_vars_PC4 = names(var_scores_ranked_PC4[1:10])
top_10_vars_PC4

loading_scores_PC5 = winter_f_pcaModel$rotation[,5]
sorted_PC5 = sort(winter_f_pcaModel$rotation[,5])
PCA_top5 = sorted_PC5[top_indices]
var_scores_PC5 = abs(loading_scores_PC5)
var_scores_ranked_PC5 = sort(var_scores_PC5, decreasing=TRUE)
top_10_vars_PC5 = names(var_scores_ranked_PC5[1:10])
top_10_vars_PC5

loading_scores_PC10 = winter_f_pcaModel$rotation[,10]
sorted_PC10 = sort(loading_scores_PC10)
PCA_top10 = sorted_PC10[top_indices]
var_scores_PC10 = abs(loading_scores_PC10)
var_scores_ranked_PC10 = sort(var_scores_PC10, decreasing=TRUE)
top_10_vars_PC10 = names(var_scores_ranked_PC10[1:10])
top_10_vars_PC10

loading_scores_PC11 = winter_f_pcaModel$rotation[,11]
sorted_PC11 = sort(loading_scores_PC11)
PCA_top11 = sorted_PC11[top_indices]
var_scores_PC11 = abs(loading_scores_PC11)
var_scores_ranked_PC11 = sort(var_scores_PC11, decreasing=TRUE)
top_10_vars_PC11 = names(var_scores_ranked_PC11[1:10])
top_10_vars_PC11

loading_scores_PC13 = winter_f_pcaModel$rotation[,13]
sorted_PC13 = sort(loading_scores_PC13)
PCA_top13 = sorted_PC13[top_indices]
var_scores_PC13 = abs(loading_scores_PC13)
var_scores_ranked_PC13 = sort(var_scores_PC13, decreasing=TRUE)
top_10_vars_PC13 = names(var_scores_ranked_PC13[1:10])
top_10_vars_PC13

loading_scores_PC18 = winter_f_pcaModel$rotation[,18]
sorted_PC18 = sort(loading_scores_PC18)
PCA_top18 = sorted_PC18[top_indices]
var_scores_PC18 = abs(loading_scores_PC18)
var_scores_ranked_PC18 = sort(var_scores_PC18, decreasing=TRUE)
top_10_vars_PC18 = names(var_scores_ranked_PC18[1:10])
top_10_vars_PC18

# Top 5 positive and negative loading scores and relevant LIWC category per PC, all together
winter_f_PCA_top1 = data.frame(loading = PCA_top1, PC = rep(c("PC1"),each=20))
winter_f_PCA_top4 = data.frame(loading = PCA_top4, PC = rep(c("PC4"),each=20))
winter_f_PCA_top5 = data.frame(loading = PCA_top5, PC = rep(c("PC5"),each=20))
winter_f_PCA_top10 = data.frame(loading = PCA_top10, PC = rep(c("PC10"),each=20))
winter_f_PCA_top11 = data.frame(loading = PCA_top11, PC = rep(c("PC11"),each=20))
winter_f_PCA_top13 = data.frame(loading = PCA_top13, PC = rep(c("PC13"),each=20))
winter_f_PCA_top18 = data.frame(loading = PCA_top18, PC = rep(c("PC18"),each=20))
winter_f_PCA_top = rbind(winter_f_PCA_top1, winter_f_PCA_top4, winter_f_PCA_top5, winter_f_PCA_top10, winter_f_PCA_top11, winter_f_PCA_top13, winter_f_PCA_top18)
winter_f_PCA_top

# Fig. 2 - PC5 violin plot
png("fig2.png", units="px", width=6125, height=7000, res=1000) # Then converted to 300 dpi in Adobe Photoshop
violin_PC5 <- ggplot(winter22_f_data, aes(x=condition=='PrFP', y=winter_f_pcaModel$x[,5], fill=condition)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  labs(title="PC5 scores per condition", x="Condition", y="PC score") +
  scale_x_discrete(labels = c("Pre-pandemic", "During pandemic")) +
  ylim(NA, 8) +
  stat_compare_means(aes(label=..p.signif..), method = "t.test", label.y = 7, label.x = 1.5, size=10) +
  theme(legend.position = "none", text = element_text(size = 16), axis.text = element_text(size = 16)) +
  scale_fill_manual(values=c("#bf615e", "#5ca84f"))
violin_PC5 # p = 0.011
dev.off()

# Look at other PCs in violin plots - slot in selected component values to build each. 4 is the only other one with a significant difference between violins. For reference, full list of selected components: 1, 4, 5, 10, 11, 13, 18
violin_PC4 <- ggplot(winter22_f_data, aes(x=condition=='PrFP', y=winter_f_pcaModel$x[,4], fill=condition)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  labs(title="PC4 scores per condition", x="Condition", y="PC score") +
  scale_x_discrete(labels = c("Pre-pandemic", "During pandemic")) +
  ylim(NA, 8) +
  stat_compare_means(aes(label=..p.signif..), method = "t.test", label.y = 7, label.x = 1.5, size=10) +
  theme(legend.position = "none", text = element_text(size = 14)) +
  scale_fill_manual(values=c("#5ca84f", "#bf615e"))
violin_PC4 # p = 0.025

summary(lm(winter_f_pcaModel$x[,5]~winter22_f_data$condition=='PrFP'))
summary(lm(winter22_f_data$focusfuture~winter22_f_data$condition=='PrFP'))
#####



####################
### Network code ###
####################
topScoring = order(rowSums(abs(winter_f_pcaModel$rotation[,selectedComponents])),decreasing=T)[1:50]
presentScore = rowSums(abs(winter_f_pcaModel$rotation[topScoring,selectedComponents]))
theDists = as.matrix(dist(winter_f_pcaModel$rotation[topScoring,selectedComponents],method='maximum'))
hist(theDists,100)
g = graph_from_adjacency_matrix(theDists<.15,mode="undirected")
g = simplify(g)
length(V(g))
length(E(g))
set.seed(42)
l = layout_with_fr(g)
V(g)$size = .1
mxScale = (presentScore-min(presentScore))/(max(presentScore)-min(presentScore))
V(g)$label.color = rgb(mxScale,0,mxScale)

### Fig. 3 - LIWC network
png("fig3.png", units="px", width=5299, height=4323, res=700) # Then converted to 300 dpi and white space trimmed using Adobe Photoshop
plot(g, layout=l, vertex.label.cex=1, edge.width=3) # Sizing is for ~1211 x 988px for exporting figure
dev.off()

hist(degree(g))
V(g)$name[order(degree(g),decreasing=T)][1:10]

### Network stats
winter_f_LIWCdegree = degree(g, v = V(g))
winter_f_LIWCbetweenness <- betweenness(g)
winter_f_LIWCtransitivity = transitivity(g)
winter_f_meanDegree = mean(winter_f_LIWCdegree)
winter_f_meanBetween = mean(winter_f_LIWCbetweenness)
#####



########################
### Permutation test ###
########################
winter_f_subData = winter22_f_data[,3:119]
# Scale data
for (i in 1:117) {
  winter_f_subData[,i] = scale(winter_f_subData[,i])
}

is.nan.data.frame <- function(x)
  do.call( cbind, lapply(x, is.nan))
winter_f_subData[is.nan(winter_f_subData)] <- 0

permDegreeDist = NULL
permBetweenDist = NULL
permCoefDist = NULL
shufRow = NULL
shufCol = NULL

set.seed(42)

for (i in 1:10000) { # Loop through determined number of permutations - 10,000 optimal
  shufData = winter_f_subData
  for (j in 1:117) {
    # Shuffle columns
    shufCol = sample(shufData[,j])
    shufData[,j] = shufCol
    shufCol = NULL
  }
  for (j in 1:91) {
    # Shuffle rows
    shufRow = sample(unname(shufData[j,]))
    shufData[j,] = shufRow
    shufRow = NULL
  }
  shufData = shufData[,colSums(abs(shufData))>0]
  shufModel = prcomp(shufData)

  # Run network stats
  shufTopScoring = order(rowSums(abs(shufModel$rotation[,selectedComponents])),decreasing=T)[1:50]
  shufDists = as.matrix(dist(shufModel$rotation[shufTopScoring,selectedComponents],method='maximum'))
  shufg = graph_from_adjacency_matrix(shufDists<.15,mode="undirected")
  shufg = simplify(shufg)
  length(V(shufg))
  length(E(shufg))
  l = layout_with_kk(shufg)
  V(shufg)$size = .1
  shufDegree = degree(shufg, v = V(shufg))
  shufBetween <- betweenness(shufg)
  shufTrans = transitivity(shufg)
  permDegreeDist[i] = mean(shufDegree)
  permBetweenDist[i] = mean(shufBetween)
  permCoefDist[i] = shufTrans
}

### Compare real data to distributions from permutations
hist(permDegreeDist)
abline(v=winter_f_meanDegree,col="blue")
1-pnorm((winter_f_meanDegree-mean(permDegreeDist))/sd(permDegreeDist))
# Last run w/ 10,000 permutations: 0.038/0.04

hist(permBetweenDist)
abline(v=winter_f_meanBetween,col="blue")
1-pnorm((winter_f_meanBetween-mean(permBetweenDist))/sd(permBetweenDist))
# Last run w/ 10,000 permutations: 0.185



hist(permCoefDist)
abline(v=winter_f_LIWCtransitivity,col="blue")
# p value
1-pnorm((winter_f_LIWCtransitivity-mean(permCoefDist))/sd(permCoefDist))
# Last run w/ 10,000 permutations: 0.0039/0.004
#####



###############
### FALL-21 ###
###############
#################
### PCA setup ###
#################
fall_f_pcaData = fall21_f_data[,3:119] # Keeping just LIWC scores
fall_f_pcaData = fall_f_pcaData[,colSums(fall_f_pcaData)>0] # Remove any LIWC categories where none had a value
fall_f_pcaModel = prcomp(scale(fall_f_pcaData))
summary(fall_f_pcaModel)

### Fig. 1, cumulative variance for PCs
plot(cumsum(fall_f_pcaModel$sdev^2)/sum(fall_f_pcaModel$sdev^2),ylim=c(0,1), ylab="Cumulative proportion", xlab="Principal Component", main="Cumulative variance explained for 91 total Principal Components")

dim(fall_f_pcaModel$x)

initialComponents = c(1:20)
fall_initial_model = glm(fall21_f_data$condition=='PrFP'~.,
                           data=data.frame(fall_f_pcaModel$x[,initialComponents]),family='binomial')
summary(fall_initial_model)

## Determine top ordering for fall21 single significant PC4
topScoring = order(-abs(fall_f_pcaModel$rotation[,4]))[1:10]
presentScore = fall_f_pcaModel$rotation[topScoring,4]


### n = 91 random subset test
fall_f_subData = fall21_f_data[sample(nrow(fall21_f_data), 91), ]
fall_f_subpcaData = fall_f_subData[,3:119] # Keeping just LIWC scores
fall_f_subpcaData = fall_f_subpcaData[,colSums(fall_f_pcaData)>0] # Remove any LIWC categories where none had a value
fall_f_subpcaModel = prcomp(scale(fall_f_subpcaData))
summary(fall_f_subpcaModel)

plot(cumsum(fall_f_subpcaModel$sdev^2)/sum(fall_f_subpcaModel$sdev^2),ylim=c(0,1), ylab="Cumulative proportion", xlab="Principal Component", main="Cumulative variance explained for 91 total Principal Components")

dim(fall_f_subpcaModel$x)

initialComponents = c(1:20)
fall_initial_model = glm(fall_f_subData$condition=='PrFP'~.,
                         data=data.frame(fall_f_subpcaModel$x[,initialComponents]),family='binomial') #
summary(fall_initial_model)

## Determine top ordering for fall21 single significant PC4
topScoring = order(-abs(fall_f_subpcaModel$rotation[,4]))[1:10]
presentScore = fall_f_subpcaModel$rotation[topScoring,4]
#####

