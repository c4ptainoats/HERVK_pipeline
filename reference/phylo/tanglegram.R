library(ape)
library(phytools)
library(dendextend)
library(viridis)
library(dplyr)
library(phylogram)
# Loads trees
tree1 <- read.nexus(file = "1.tree")
# Assumes root at the smallest distance from any of the tips 
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "2.tree")
tree2 <- midpoint.root(tree2)
# Gets branch lengths
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
# Draws the tree
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
# Runs an algorithm to align the trees together optimally - so as many tips are faceing each other as possible
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)
# Computes tanglegram lines with formatting
dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
# Produces pdf output
dev.copy(pdf, '12.pdf', width = 10, height = 11)
dev.off()

###### Lines below correspond to the methods above, only follow comparison between the consecutive sequence fragments.
tree1 <- read.nexus(file = "2.tree")
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "3.tree")
tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
#dend12_corrected <- dendextend::dendlist(tree1, tree2)
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)

dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
dev.copy(pdf, '23.pdf', width = 10, height = 11)
dev.off()

######
tree1 <- read.nexus(file = "3.tree")
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "4.tree")
tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
#dend12_corrected <- dendextend::dendlist(tree1, tree2)
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)

dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
dev.copy(pdf, '34.pdf', width = 10, height = 11)
dev.off()
######
tree1 <- read.nexus(file = "4.tree")
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "5.tree")
tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
#dend12_corrected <- dendextend::dendlist(tree1, tree2)
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)

dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
dev.copy(pdf, '45.pdf', width = 10, height = 11)
dev.off()
######
tree1 <- read.nexus(file = "5.tree")
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "6.tree")
tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
#dend12_corrected <- dendextend::dendlist(tree1, tree2)
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)

dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
dev.copy(pdf, '56.pdf', width = 10, height = 11)
dev.off()
######
tree1 <- read.nexus(file = "6.tree")
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "7.tree")
tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
#dend12_corrected <- dendextend::dendlist(tree1, tree2)
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)

dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
dev.copy(pdf, '67.pdf', width = 10, height = 11)
dev.off()
######
tree1 <- read.nexus(file = "7.tree")
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "8.tree")
tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
#dend12_corrected <- dendextend::dendlist(tree1, tree2)
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)

dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
dev.copy(pdf, '78.pdf', width = 10, height = 11)
dev.off()
######
tree1 <- read.nexus(file = "8.tree")
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "9.tree")
tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
#dend12_corrected <- dendextend::dendlist(tree1, tree2)
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)

dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
dev.copy(pdf, '89.pdf', width = 10, height = 11)
dev.off()
######
tree1 <- read.nexus(file = "9.tree")
tree1 <- midpoint.root(tree1)
tree2 <- read.nexus(file = "10.tree")
tree2 <- midpoint.root(tree2)
tree1 <- compute.brlen(tree1)
tree2 <- compute.brlen(tree2)
tree1<- as.dendrogram(tree1)
tree1
tree2<- as.dendrogram(tree2)
#dend12_corrected <- dendextend::dendlist(tree1, tree2)
dend12_corrected<-untangle_step_rotate_2side(tree1,tree2)

dendextend::tanglegram(dend12_corrected,common_subtrees_color_lines = FALSE, margin_inner = 5.5, lab.cex = 0.3, lwd = 
                         0.2, edge.lwd = 0.2, type = "r")
dev.copy(pdf, '910.pdf', width = 10, height = 11)
dev.off()