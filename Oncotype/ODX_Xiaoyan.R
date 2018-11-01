# OncotypeDX scoring
# modified from 4) Run RS November 9 2015.R by Nick
# 2016-9-21


# Make a recurrence score function
oneToFifteen <- function(x) {
  
  x <- x - min(x)
  x <- x * (15/max(x))
  x
}

GRB7groupscore <- function(arr) {
  
  arr["GRB7", ] <- arr["GRB7", ] * 0.9
  arr["HER2", ] <- arr["HER2", ] * 0.1
  
  gs <- apply(arr[c("GRB7", "HER2"), ], 2, sum)
  gs[gs < 8] <- 8
  gs
}

ERgroupscore <- function(arr) {
  
  arr["ER", ] <- arr["ER", ] * 0.2
  arr["PGR", ] <- arr["PGR", ] * 0.3
  arr["BCL2", ] <- arr["BCL2", ] * 0.25
  arr["SCUBE2", ] <- arr["SCUBE2", ] * 0.25
  
  gs <- apply(arr[c("ER", "PGR", "BCL2", "SCUBE2"), ], 2, sum)
  gs
}

proliferationgroupscore <- function(arr) {
  
  arr["BIRC5", ] <- arr["BIRC5", ] * 0.2
  arr["Ki67", ] <- arr["Ki67", ] * 0.2
  arr["MYBL2", ] <- arr["MYBL2", ] * 0.2
  arr["CCNB1", ] <- arr["CCNB1", ] * 0.2
  arr["STK15", ] <- arr["STK15", ] * 0.2
  
  gs <- apply(arr[c("BIRC5", "Ki67", "MYBL2", "CCNB1", "STK15"), ], 2, 
              sum)
  gs[gs < 6.5] <- 6.5
  gs
}

invasiongroupscore <- function(arr) {
  
  arr["CTSL2", ] <- arr["CTSL2", ] * 0.5
  arr["MMP11", ] <- arr["MMP11", ] * 0.5
  
  gs <- apply(arr[c("CTSL2", "MMP11"), ], 2, sum)
  gs
}

RS <- function(arr) {
  
  arr["GRB7gs", ] <- arr["GRB7gs", ] * 0.47
  arr["ERgs", ] <- arr["ERgs", ] * -0.34
  arr["proliferationgs", ] <- arr["proliferationgs", ] * 1.04
  arr["invasiongs", ] <- arr["invasiongs", ] * 0.1
  arr["CD68", ] <- arr["CD68", ] * 0.05
  arr["GSTM1", ] <- arr["GSTM1", ] * -0.08
  arr["BAG1", ] <- arr["BAG1", ] * -0.07
  
  RSu <- apply(arr, 2, sum)
  RS <- 20 * (RSu - 6.7)
  RS[RS < 0] <- 0
  RS[RS > 100] <- 100
  
  class <- rep("intermediate", length(RS))
  class[RS >= 31] <- "high"
  class[RS < 18] <- "low"
  rbind(RSu, RS, class)
}

paikfunction <- function(mat, group) {
  
  reference <- apply(mat[group, ], 2, mean)
  matr <- t(apply(mat, 1, function(x) x - reference))
  # print(matr[,1:5])
  
  matr <- apply(matr, 1, function(x) oneToFifteen(x))
  matr <- t(matr)
  # print(matr[,1:25])
  
  GRB7gs <- GRB7groupscore(matr)
  ERgs <- ERgroupscore(matr)
  proliferationgs <- proliferationgroupscore(matr)
  invasiongs <- invasiongroupscore(matr)
  
  mat2 <- rbind(GRB7gs, ERgs, proliferationgs, invasiongs, CD68 = matr["CD68", ],
                GSTM1 = matr["GSTM1", ],
                BAG1 = matr["BAG1", ])
  
  # print(mat2)[1:5,1:5]
  o <- RS(mat2)
  
  rbind(mat2[c("GRB7gs", "ERgs", "proliferationgs", "invasiongs"), ], o)
}

# Load in paik annotation object and clean
load("E:/GitLocal/iss-analysis/Subtyping/paik110504.RData")
tf1 = duplicated(paik$symbol) 
paik.clean = paik[!tf1,] # remove duplicates
paik.clean = paik.clean[,1:4] # remove superfluous information

colnames(paik.clean)[4] = "EntrezGene.ID"
levels(paik.clean$symbol)[levels(paik.clean$symbol)=="Survivin"] <- "BIRC5"


####### DESseq ########
library(DESeq2)

# get command line input
args = commandArgs(trailingOnly = FALSE)

# read in and set genes as rownames
# library(gdata)
# orig.RNA.seq = read.csv("C:/Users/qxyyx/OneDrive/work/windowsize_5000.csv")
orig.RNA.seq = read.csv(args[7])
rownames(orig.RNA.seq) = orig.RNA.seq[,"GeneName"]
orig.RNA.seq = orig.RNA.seq[,-1]

table(rownames(orig.RNA.seq) %in% paik.clean[,"symbol"]) # TRUE 21
countdata = orig.RNA.seq

colnames(countdata) = sub("X","",colnames(countdata)) 


# clean incorrect names
# already corrected in the input file /Xiaoyan
# rownames(orig.RNA.seq)[3] = "RPLPO"

# count data slot
countdata <- orig.RNA.seq + 1

coldata = data.frame(1:dim(countdata)[2])
rownames(coldata) <- colnames(countdata)



ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ 1) # this desgin is fine for explotatory analysis but not for DEG

rld <- rlog(ddsMat, blind=TRUE)

# # some plots to see what the transformation did
# par( mfrow = c( 1, 2 ) )
# dds <- estimateSizeFactors(ddsMat)
# plot(log2(counts(ddsMat, normalized=FALSE)[,1:2] + 1),
#      pch=16, cex=0.3)
# plot(assay(rld)[,1:2],
#      pch=16, cex=0.3)

# # any samples similar?
# sampleDists <- dist( t( assay(rld) ) )
# sampleDists
# 
# # heatmap
# library("RColorBrewer")
# library("pheatmap")
# sampleDistMatrix <- as.matrix( sampleDists )
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix, col=colors)

# run oncotypeDX on this object
gr = paik.clean[match(rownames(rld),paik.clean$symbol),"group"]=="ref"
rld_rs = paikfunction(assay(rld), gr)
rld_rs.df = t(data.frame(rld_rs))

# same again but mean center first (as range of control genes differ so much)
library(Biobase)
mc.rld = rld

# rowmedians and center
rms = rowMedians(assay(mc.rld))
mc.data =  assay(mc.rld) -  rms

# add back in to assay data
assay(mc.rld) = mc.data

# # re-run OncotypeDx
# mc.gr = paik.clean[match(rownames(mc.rld),paik.clean$symbol),"group"]=="ref"
# mc.rld_rs = paikfunction(assay(mc.rld), gr)
# mc.rld_rs.df = t(data.frame(mc.rld_rs))
# 
# table(mc.rld_rs.df[,"class"], rld_rs.df[,"class"]) # identical
# plot(mc.rld_rs.df[,"RS"], rld_rs.df[,"RS"]) # identical

# write out results
library(xlsx)

xlsx.writeMultipleData <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i])
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE)
  }
}


xlsx.writeMultipleData(args[8], rld_rs.df)

