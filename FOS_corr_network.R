# Introduction ####
# This script is used to generate FOS correlation networks from fold change in FOS expression data (steps 1 and 3),
# visualize these networks (step 4), and compute graph theory metrics (step 5). In addition, positive Pearson's r values
# are Fisher z-transformed for comparison across treatment conditions by one-way ANOVA (step 2).

# Load packages ####
library(corrplot) #to generate correlation matrices
library(igraph) #to visualize networks and compute graph theory metrics
library(DescTools) #to perform z transformation

# 1. Generate FOS correlation matrices ####
naive <- read.csv("~/Desktop/FCbyRegion_Naive.csv")
h24 <- read.csv("~/Desktop/FCbyRegion_24h.csv")
wk4 <- read.csv("~/Desktop/FCbyRegion_4wk.csv")

cor_naive<-cor(naive,method="pearson")
cor_h24<-cor(h24,method="pearson")
cor_wk4<-cor(wk4,method="pearson")

# visualize correlation matrix (repeat this code for each treatment condition)
col <- colorRampPalette(c("#224E7A","#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444","#670000"))
corrplot(cor_naive, method="color", type="full", col=col(200),  
         #addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.cex=0.8, tl.srt = 45, #Text label color, size and angle
         diag=TRUE 
)

# 2. Fisher z transform positive Pearson's r values (repeat this code for each treatment condition) ####
# first transform the data from a matrix into a dataframe
ds_naive = data.frame(cor = as.vector(cor_naive), 
                id1=rownames(cor_naive), 
                id2=rep(rownames(cor_naive), each=nrow(cor_naive)))
ds_naive = ds_naive[as.vector(upper.tri(cor_naive, TRUE)),] #remove duplicates
ds_naive = ds_naive[!ds_naive$cor == 1.00000000, ] #exclude self-correlations
ds_naive = ds_naive[!ds_naive$cor < 0, ] #exclude negative correlations
# then z transform the correlation values in the dataframe
zs_naive <- FisherZ(ds_naive$cor)
# the z-transformed Pearson's r values for each treatment condition were compared by one-way ANOVA using Prism software

# 3. Make weighted FOS correlation network graphs ####
#set up the graph from the correlation matrix
g1<-graph_from_adjacency_matrix(cor_naive,mode="undirected",weighted=TRUE,diag=FALSE)
g2<-graph_from_adjacency_matrix(cor_h24,mode="undirected",weighted=TRUE,diag=FALSE)
g3<-graph_from_adjacency_matrix(cor_wk4,mode="undirected",weighted=TRUE,diag=FALSE)
#remove multiples and loops
g1<- simplify(g1, remove.multiple=TRUE, remove.loops=TRUE)
g2<- simplify(g2, remove.multiple=TRUE, remove.loops=TRUE)
g3<- simplify(g3, remove.multiple=TRUE, remove.loops=TRUE)
#remove negative edge weights
g1 <- delete_edges(g1, E(g1)[which(E(g1)$weight<0)])
g2 <- delete_edges(g2, E(g2)[which(E(g2)$weight<0)])
g3 <- delete_edges(g3, E(g3)[which(E(g3)$weight<0)])
#create vectors from the edge weights for each graph
edgeweightsg1 <- E(g1)$weight
edgeweightsg2 <- E(g2)$weight
edgeweightsg3 <- E(g3)$weight

##examine features of the network graph (repeat this code for each treatment condition)
vcount(g1) #number of nodes
ecount(g1) #numer of edges

# 4. Plot each weighted network graph (repeat this code for each treatment condition) ####
#assign node shape
V(g1)$shape <- "circle"
#color the nodes according to their anatomical group
V(g1)$color[V(g1)$name %in% c("dACC", "vACC", "AId","AIv","Cla")] <- "palegreen3"
V(g1)$color[V(g1)$name %in% c("CPu", "NAc")] <- "lightblue3"
V(g1)$color[V(g1)$name %in% c("VP","BNST")] <- "rosybrown1"
V(g1)$color[V(g1)$name %in% c("BLA","CeA")] <- "yellow3"
V(g1)$color[V(g1)$name %in% c("DG")] <- "lightcyan2"
V(g1)$color[V(g1)$name %in% c("MHb","LHb","PVT")] <- "darksalmon"
V(g1)$color[V(g1)$name %in% c("PAG","SNc","SNr","VTA")] <- "plum2"
#plot the network
par(family="Arial")
plot(g1,
        layout=layout.fruchterman.reingold,
        edge.curved=FALSE,
        vertex.size=18,
        #vertex.label.dist=-0.5,
        vertex.label.color="black",
        vertex.label.cex=0.5,
        asp=FALSE,
        vertex.label.cex=0.5,
        edge.width=edgeweightsg1,
        edge.arrow.mode=0
)

# 5. Compute graph theory metrics ####
#weighted degree
g1_strength<-graph.strength(g1,weights=edgeweightsg1) #weighted degree
g2_strength<-graph.strength(g2,weights=edgeweightsg2)
g3_strength<-graph.strength(g3,weights=edgeweightsg3)
ks.test(g1_strength, g2_strength) #repeat this for each pair of conditions
#bonferroni correct across p-values from the three pairwise comparisons
p=c(0.0007141,0.009224,0.3057)
p.adjust(p, method = "bonferroni", n = length(p))

#betweenness centrality
g1_bet<-betweenness(g1,weights=edgeweightsg1)
g2_bet<-betweenness(g2,weights=edgeweightsg2)
g3_bet<-betweenness(g3,weights=edgeweightsg3)
ks.test(g1_bet, g2_bet) #repeat this for each pair of conditions
#bonferroni correct across p-values from the three pairwise comparisons
p=c(0.2997,0.5262,0.5262)
p.adjust(p, method = "bonferroni", n = length(p))