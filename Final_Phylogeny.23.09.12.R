
library(phangorn)
library(adegenet)
options(scipen=999)
par(xpd=TRUE)


#################
dna <- adegenet::fasta2DNAbin(file="MLST_muscle_alignment_renamed_pick2.fasta")
D <- dist.dna(dna, model = "TN93")
Dm <- as.matrix(D)

annot <- NULL
annot$id <- row.names(Dm)
annot <- as.data.frame(annot)

annot$new <- "old"
annot$new[substr(annot$id, 1, 2) %in% c("S0", "S1")] <- "new"

annot$col <- "white"
annot$col[substr(annot$id, 1, 2) %in% c("S0", "S1")] <- "gold"

tre <- nj(D)
tre <- ladderize(tre)

set.seed(42)
myBoots <- boot.phylo(tre, dna, function(e) nj(dist.dna(e, model = "TN93")), B = 100)

par(mfrow = c(1, 2))

{plot(tre, cex = 0.5)
  title(paste("NJ tree\n", "MLST Genes"),
        #       "\n PERMANOVA p value =", res$`Pr(>F)`[1]),
        sub = "TN93")
  nodelabels(myBoots,  frame = "none", font = 4,
             cex=.5, adj = c(1.1, -.5), col = "darkred", bg = "lightblue")
  tiplabels(tre$tip.label, adj = c(0, 0.5), bg=annot$col,cex=.5, fg="transparent")
  axisPhylo()}

temp <- tre
N <- length(tre$tip.label)
toCollapse <- match(which(myBoots<25)+N, temp$edge[,2])
temp$edge.length[toCollapse] <- 0
tre <- di2multi(temp, tol=0.00001)

set.seed(42)
myBoots <- boot.phylo(tre, dna, function(e) nj(dist.dna(e, model = "TN93")), B = 100)

#{setEPS(bg = "white", family = "Times", width = 6, pointsize = 12)
# postscript("MLST_tree.eps", onefile = TRUE)

{tiff(filename = "MLST_tree.tiff", res = 500, width = 2500, height = 2500 )
  {plot(tre, cex = 0.5)
    title("NJ tree MLST Genes",
          sub = "TN93 \n Nodes with bootstrab < 25 collapsed")
    nodelabels(myBoots,  frame = "none", font = 4,
               cex=.5, adj = c(1.4, -.5), col = "darkred", bg = "lightblue")
    
    tiplabels(tre$tip.label, adj = c(0, 0.5), bg=annot$col,cex=.5, fg="transparent")
    axisPhylo()}
  dev.off()}
###################


dna.core <- adegenet::fasta2DNAbin(file="CORE_muscle_alignment_renamed_pick2.fasta")
D <- dist.dna(dna.core, model = "TN93")
Dm <- as.matrix(D)

annot <- NULL
annot$id <- row.names(Dm)
annot <- as.data.frame(annot)

annot$new <- "old"
annot$new[substr(annot$id, 1, 2) %in% c("S0", "S1")] <- "new"

annot$col <- "white"
annot$col[substr(annot$id, 1, 2) %in% c("S0", "S1")] <- "gold"

tre <- nj(D)
tre <- ladderize(tre)

set.seed(42)
myBoots <- boot.phylo(tre, dna.core, function(e) nj(dist.dna(e, model = "TN93")))

{tiff(filename = "Core_tree.tiff", res = 500, width = 2500, height = 2500 )
  
  {plot(tre, cex = 0.5)
    title(paste("NJ tree\n", "Core Genome"),
          #       "\n PERMANOVA p value =", res$`Pr(>F)`[1]),
          sub = "TN93")
    nodelabels(myBoots,  frame = "none", font = 4,
               cex=.5, adj = c(1.4, -.5), col = "darkred", bg = "lightblue")
    tiplabels(tre$tip.label, adj = c(0, 0.5), bg=annot$col,cex=.5, fg="transparent")
    axisPhylo()}
  dev.off()}

