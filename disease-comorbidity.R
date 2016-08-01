#### PROLOGUE ####
## libraries 
library(igraph)
library(ggplot2)
library(leaflet)

## setwd, infile, outfile
networkfile <- "./data/allnet3.net"

#### DATA MANAGEMENT ####

## read data
data <- read.delim(networkfile, header = F)
names(data) <- c("d1", "d2", "p1", "p2", "co", "rr", "rrcil", "rrcir", "phi", "t")

## add disease grouping
    unique_d1 <- unique(data$d1)
    unique_d2 <- unique(data$d2)
    diseases <- unique(c(unique_d1, unique_d2))


diseaseGroup <- data.frame(icd9code = diseases, group = character(numDisease), stringsAsFactors = F)


# assign the icd9 groupings (Lists compiled from Wikipedia: https://en.wikipedia.org/wiki/List_of_ICD-9_codes)
categories <- "List of ICD-9 codes 001-139: infectious and parasitic diseases
List of ICD-9 codes 140-239: neoplasms
List of ICD-9 codes 240-279: endocrine, nutritional and metabolic diseases, and immunity disorders
List of ICD-9 codes 280-289: diseases of the blood and blood-forming organs
List of ICD-9 codes 290-319: mental disorders
List of ICD-9 codes 320-359: diseases of the nervous system
List of ICD-9 codes 360-389: diseases of the sense organs
List of ICD-9 codes 390-459: diseases of the circulatory system
List of ICD-9 codes 460-519: diseases of the respiratory system
List of ICD-9 codes 520-579: diseases of the digestive system
List of ICD-9 codes 580-629: diseases of the genitourinary system
List of ICD-9 codes 630-679: complications of pregnancy, childbirth, and the puerperium
List of ICD-9 codes 680-709: diseases of the skin and subcutaneous tissue
List of ICD-9 codes 710-739: diseases of the musculoskeletal system and connective tissue
List of ICD-9 codes 740-759: congenital anomalies
List of ICD-9 codes 760-779: certain conditions originating in the perinatal period
List of ICD-9 codes 780-799: symptoms, signs, and ill-defined conditions
List of ICD-9 codes 800-999: injury and poisoning"

write(categories, file = "categories-temp.txt") # convert to textfile for easy conversion to vector
categories <- readLines("categories-temp.txt")

category <- gsub("List of ICD-9 codes [0-9]{3}-[0-9]{3}: ", "", categories)
uppercode <- as.numeric(gsub("^.*([0-9]{3}).*", "\\1", categories))

## assign the categories
for (i in seq(length(uppercode), 1, -1)){
  upperBound <- uppercode[i]
  diseaseGroup$group[diseaseGroup$icd9code <= upperBound] <- category[i]
}
#diseaseGroup$group <- diseaseGroup$group)

diseaseGroup$group <- as.numeric(as.factor(diseaseGroup$group))

diseaseGroup <- diseaseGroup[order(diseaseGroup$group),]

## add prevalence
uniqueD1 <- data[!duplicated(data$d1), c("d1", "p1")]
uniqueD2 <- data[!duplicated(data$d2), c("d2", "p2")]
names(uniqueD1) <- c("icd9code", "prevalence")
names(uniqueD2) <- c("icd9code", "prevalence")
prevalences <- unique(rbind(uniqueD1, uniqueD2))

diseaseGroup <- merge(diseaseGroup, prevalences, by = "icd9code")

## color assignment for the network graph
colorFun <- colorNumeric("Set", 1:18)
colors <- colorFun(1:18)

numReps <- split(diseaseGroup$icd9code, diseaseGroup$group)
numReps <- sapply(numReps, length)

colorList <- rep(colors, numReps)
diseaseGroup$color <- colorList

## descriptives
  
# diseases
    # number of diseases
    numDisease <- length(diseases)
    numDisease     # 994 unique diseases


# co-occurence
hist(log(data$co))
boxplot(data$co)

summary(data$co) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1       4      21    1473     164 1403000 

head(data[order(data$co, decreasing = T), c("d1", "d2")], 20)     # shows top 20 disease co-occurence

## replication of graphs

# distribution of RR between all disease pairs
p <- ggplot(data = data) + geom_density(aes(x = log(rr,2), ..count..)) + geom_vline(xintercept = 0) + xlim(-4,6)

d <- ggplot_build(p)$data[[1]]  # gets the x,y points of the density curve calculated by geom_density

p <- p + geom_area(data = subset(d, x > 0), aes(x = x, y = y), fill = "green") +
  geom_area(data = subset(d, x < 0), aes(x = x, y = y), fill = "red")
  
p

# scatterplot between RR and Phi
p <- ggplot(data = data) + geom_point(aes(x = phi, y = log(rr)))
p       # this takes a while since there are a LOT of points

# prevalence distribution for all diseases

#############
#       NOT SURE ABOUT THIS. is this cumulative or no? shows long tail distribution in the paper, but not in the replicated version
# prevalence <- hist(c(log(data$p1), log(data$p2)))
#############


# number of links and nodes at different values of phi

pairs <- data[order(data$phi, decreasing = F),c("d1", "d2", "phi")]
pairs <- pairs[0.0001 <= pairs$phi & pairs$phi < 1, ]
pairs$num <- seq(nrow(pairs), 1, -1)


phis <- c(0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

numNodes <- sapply(phis, function(x){
  remaining <- pairs[pairs$phi >= x, ]
  length(unique(c(remaining$d1, remaining$d2)))
})

numEdges <- sapply(phis, function(x){
  remaining <- pairs[pairs$phi >= x, ]
  nrow(remaining)
})


plot(numEdges/nrow(data) ~ phis, type = "l", col = "blue", ylim = c(0,1), log = "x", ylab = "Fraction of Nodes \n Fraction of links", xlab = "phi*")
lines(numNodes/numDisease ~ phis, type = "l", col = "red", log = "x")

# replicate Fig S 2
cutoffs <- c(0.04, 0.06, 0.1)

par(mfrow = c(1,3))
for (cutoff in cutoffs) {
  data_sub <- data[data$phi > cutoff & data$t >= 2.58,]   # significance at 1% level, t >= 2.58, at 5% level, t >= 1.96
  numDiseasesAtCutOff <- length(unique(c(data_sub$d1, data_sub$d2)))
  numEdgesAtCutOff <- nrow(data_sub)
  graph <- graph.data.frame(data_sub, T, vertices = diseaseGroup)
  
  degree <- as.vector(degree(graph, mode = "total"))
  graph <- delete.vertices(graph, which(degree == 0))    # remove isolates
  layout <- layout_with_fr(graph)
  
  plot(graph,
       vertex.label = NA, 
       edge.arrow.size = 0.001, 
       vertex.size = 6, 
       layout = layout)
  mtext(paste("phi* = ", cutoff, "N =",numDiseasesAtCutOff, "L =", numEdgesAtCutOff), side = 1)
}



#### graph of diseases ####

## phi
data_sub <- data[data$phi > 0.06 & data$t >= 2.58, ]     # significance at 1% level, t >= 2.58, at 5% level, t >= 1.96

prevalencePerc <- diseaseGroup$prevalence*100 / sum(diseaseGroup$prevalence)
vertexSize <- rep(3, length(prevalencePerc))
vertexSize[prevalencePerc > 1] <- 5
vertexSize[prevalencePerc > 30] <- 6

diseaseGroup$vertexSize <- vertexSize


graph <- graph.data.frame(data_sub, T, vertices = diseaseGroup)

edgeWidth <- rep(0.001, nrow(data_sub))
edgeWidth[data_sub$phi > 0.1] <- 0.05
edgeWidth[data_sub$phi > 0.2] <- 2

edgeColor <- rep("grey85", nrow(data_sub))
edgeColor[data_sub$phi > 0.1] <- "steelblue"
edgeColor[data_sub$phi > 0.2] <- "coral"

E(graph)$width <- edgeWidth
E(graph)$color <- edgeColor


edgeList <- get.edgelist(graph, names = T)
fromGroup <- diseaseGroup$group[as.numeric(edgeList[,1])]
toGroup <- diseaseGroup$group[as.numeric(edgeList[,2])]
is.SameGroup <- fromGroup == toGroup

# adjust edge weight according to groups; nodes of disease within groups will be closer than those that are not

E(graph)$weight <- rep(1, length(E(graph)))
E(graph)$weight[is.SameGroup] <- 10000000


degree <- as.vector(degree(graph, mode = "total"))
graph <- delete.vertices(graph, which(degree == 0))   # delete isolates
layout <- layout_with_fr(graph)

pdf("graph1-phi.pdf", height = 20, width = 20)
plot(graph, 
     vertex.label.cex = 1,
     vertex.label.color = "black",
     vertex.label.family = "Helvetica",
     vertex.label.font = 2,
     edge.arrow.size = 0.001, 
     vertex.size = V(graph)$vertexSize, 
     layout = layout,
     edge.width = E(graph)$width,
     edge.color = E(graph)$color
     )
dev.off()

## biggest subgraph
comps <- components(graph)   # returns components statistics
whichGroup <- which.max(comps$csize)   # returns the index of the biggest group

subgraph <- induced_subgraph(graph, which(comps$membership == whichGroup))

degree <- as.vector(degree(subgraph, mode = "total"))
subgraph <- delete.vertices(subgraph, which(degree == 0))   # delete isolates
layout <- layout_with_fr(subgraph, grid = "grid")

pdf("subgraph1-phi.pdf", height = 20, width = 20)
plot(subgraph, 
     vertex.label.cex = 1,
     vertex.label.color = "black",
     vertex.label.family = "Helvetica",
     vertex.label.font = 2,
     edge.arrow.size = 0.001, 
     vertex.size = V(subgraph)$vertexSize, 
     layout = layout,
     edge.width = E(subgraph)$width,
     edge.color = E(subgraph)$color
     )
dev.off()


## RR
data_sub <- data[data$rr > 20 & data$t >= 2.58, ]     # significance at 1% level, t >= 2.58, at 5% level, t >= 1.96

prevalencePerc <- diseaseGroup$prevalence*100 / sum(diseaseGroup$prevalence)
vertexSize <- rep(3, length(prevalencePerc))
vertexSize[prevalencePerc > 1] <- 5
vertexSize[prevalencePerc > 30] <- 6

diseaseGroup$vertexSize <- vertexSize


graph <- graph.data.frame(data_sub, T, vertices = diseaseGroup)

edgeWidth <- rep(0.001, nrow(data_sub))
edgeWidth[data_sub$phi > 0.1] <- 0.05
edgeWidth[data_sub$phi > 0.2] <- 2

edgeColor <- rep("grey85", nrow(data_sub))
edgeColor[data_sub$phi > 0.1] <- "steelblue"
edgeColor[data_sub$phi > 0.2] <- "coral"

E(graph)$width <- edgeWidth
E(graph)$color <- edgeColor


edgeList <- get.edgelist(graph, names = T)
fromGroup <- diseaseGroup$group[as.numeric(edgeList[,1])]
toGroup <- diseaseGroup$group[as.numeric(edgeList[,2])]
is.SameGroup <- fromGroup == toGroup

# adjust edge weight according to groups; nodes of disease within groups will be closer than those that are not

E(graph)$weight <- rep(1, length(E(graph)))
E(graph)$weight[is.SameGroup] <- 10000000


degree <- as.vector(degree(graph, mode = "total"))
graph <- delete.vertices(graph, which(degree == 0))   # delete isolates
layout <- layout_with_fr(graph)

pdf("graph1-rr.pdf", height = 20, width = 20)
plot(graph, 
     vertex.label.cex = 1,
     vertex.label.color = "black",
     vertex.label.family = "Helvetica",
     vertex.label.font = 2,
     edge.arrow.size = 0.001, 
     vertex.size = V(graph)$vertexSize, 
     layout = layout,
     edge.width = E(graph)$width,
     edge.color = E(graph)$color
     )
dev.off()

## biggest subgraph
comps <- components(graph)   # returns components statistics
whichGroup <- which.max(comps$csize)   # returns the index of the biggest group

subgraph <- induced_subgraph(graph, which(comps$membership == whichGroup))

degree <- as.vector(degree(subgraph, mode = "total"))
subgraph <- delete.vertices(subgraph, which(degree == 0))   # delete isolates
layout <- layout_with_fr(subgraph, grid = "grid")

pdf("subgraph1-rr.pdf", height = 20, width = 20)
plot(subgraph, 
     vertex.label.cex = 0.3*V(subgraph)$vertexSize,
     vertex.label.color = "black",
     vertex.label.family = "Helvetica",
     vertex.label.font = 2,
     edge.arrow.size = 0.001, 
     vertex.size = V(subgraph)$vertexSize, 
     layout = layout,
     edge.width = E(subgraph)$width,
     edge.color = E(subgraph)$color
     )
dev.off()



### SANDBOX
    # some functions from stackoverflow hunts
# delete.isolates <- function(graph, mode = 'all') {
#   isolates <- which(degree(graph, mode = mode) == 0) - 1
#   delete.vertices(graph, isolates)
# }
# 
# weight.community=function(row,membership,weigth.within,weight.between){
#   if(as.numeric(membership[which(names(membership)==row[1])])==as.numeric(membership[which(names(membership)==row[2])])){
#     weight=weigth.within
#   } else {
#     weight=weight.between
#   }
# return(weight)
# }
# 
# networkgraph <- delete.isolates(networkgraph)
# E(networkgraph)$weight=apply(get.edgelist(networkgraph),1,weight.community,membership,10,1)
#layout <- layout.fruchterman.reingold(graph)

# create a logical vector: true if fromVertex and toVertex are from the same category, false otherwise

# edgeList <- get.edgelist(biggestSubgraph, names = T)
# fromGroup <- diseaseGroup$group[as.numeric(edgeList[,1])]
# toGroup <- diseaseGroup$group[as.numeric(edgeList[,2])]
# is.SameGroup <- fromGroup == toGroup
# 
# # adjust edge weight according to groups; nodes of disease within groups will be closer than those that are not
# 
# E(biggestSubgraph)$weight <- 1000
# E(biggestSubgraph)$weight[is.SameGroup] <- 10

# # remove isolated nodes:
# degree <- as.vector(degree(graph, mode = "total"))
# graphNoIsolates <- delete.vertices(graph, which(degree == 0))
# 
# # partition graph into components to get the biggest connected subgraph
# comps <- components(graphNoIsolates)   # returns components statistics
# whichGroup <- which.max(comps$csize)   # returns the index of the biggest group
# 
# biggestSubgraph <- induced_subgraph(graphNoIsolates, which(comps$membership == whichGroup))
# 
# 
# 
# ## create grid layout
# # this layout draws the nodes in space such that icd9 codes of same groups get 'clustered'
# baseX <- rep(c(5,10,15,20,25), 4)[1:18]
# baseY <- rep(c(5,10,15,20), c(5,5,5,5))[1:18]
# 
# x <- rep(baseX, numReps)[which(degree != 0)]
# y <- rep(baseY, numReps)[which(degree != 0)]
# 
# edgeList <- get.edgelist(graphNoIsolates, names = T)
# fromGroup <- diseaseGroup$group[as.numeric(edgeList[,1])]
# toGroup <- diseaseGroup$group[as.numeric(edgeList[,2])]
# is.SameGroup <- fromGroup == toGroup
# 
# # adjust edge weight according to groups; nodes of disease within groups will be closer than those that are not
# 
# weight <- rep(1, length(E(graphNoIsolates)))
# weight[is.SameGroup] <- 1000
# 
# 
# ## graph: no isolates
# #layout
# delta <- 2
# layout <- layout_with_fr(graphNoIsolates, minx = x - delta, maxx = x + delta, 
#                          miny = y - delta, maxy = y + delta, weights = weight)
# 
# layout <- layout_with_fr(graphNoIsolates, weights = weight, coords = matrix(c(x,y), ncol = 2))
# #envelopes 
# envelopes <- split(V(graphNoIsolates)$name, V(graphNoIsolates)$group)
# 
# plot(graphNoIsolates, 
#      vertex.label = NA, 
#      edge.arrow.size = 0.001, 
#      vertex.size = 4,
#      layout = layout, 
#      rescale = T
#      ,
#      mark.groups = envelopes,
#      mark.col = rainbow(length(envelopes), alpha = 0.05),
#      mark.shape = -1
#      )
# 
# # graph: biggest subgraph
# 
# 
# 
# 
# edgeList <- get.edgelist(biggestSubgraph, names = T)
# fromGroup <- diseaseGroup$group[as.numeric(edgeList[,1])]
# toGroup <- diseaseGroup$group[as.numeric(edgeList[,2])]
# is.SameGroup <- fromGroup == toGroup
# weight <- rep(0.0001, length(E(biggestSubgraph)))
# weight[is.SameGroup] <- 1
# layout <- layout_with_fr(biggestSubgraph, weights = weight)
# 
# x <- rep(baseX, numReps)
# y <- rep(baseY, numReps)
# 
# layout <- matrix(data = c(jitter(x, 2.5)*100, jitter(y, 2.5)*100), ncol = 2)[which(degree != 0),][which(comps$membership == whichGroup),]
# 
# envelopes <- split(V(biggestSubgraph)$name, V(biggestSubgraph)$group)
# 
# plot(biggestSubgraph, 
#      vertex.label = NA, 
#      edge.arrow.size = 0.001, 
#      vertex.size = 4,
#      edge.width = weight,
#      layout = layout,
#      mark.groups = envelopes,
#      mark.col = rainbow(length(envelopes), alpha = 0.05),
#      mark.shape = -1)
# 
# #
# ## create grid layout
# # this layout draws the nodes in space such that icd9 codes of same groups get 'clustered'
# baseX <- rep(c(1,2,3,4,5), 4)[1:18]
# baseY <- rep(c(1,2,3,4), c(5,5,5,5))[1:18]
# numReps <- aggregate(icd9code ~ group, FUN = length, data = diseaseGroup)[,2][order(categories)]
# x <- rep(baseX, numReps)
# y <- rep(baseY, numReps)
# layout <- matrix(data = c(jitter(x, 2)*100, jitter(y, 2)*100), ncol = 2)
# 
# # export graph as .gml / can be read by cytoscape for manual visualization
# write.graph(graph, file = "pdn.gml", format = "gml")
