#install and call all necessary libraries for analysis

library(limma)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(org.Dm.eg.db)

counts <- read.delim(file = "counts_Drosophila.txt")
targets <- read.delim(file = "SampleInfo_Drosophila.txt")

head(counts)
table(targets$Group)

mycpm <- cpm(counts)
plot(counts[,1], mycpm[,1], xlim = c(0,20), ylim = c(0,5))
abline(v = 10, col = 2)
abline(h = 2, col = 4)

threshold <- mycpm > 2
keep <- rowSums(threshold) >= 3
table(keep)

counts.keep <- counts[keep,]
dim(counts.keep)

y <- DGEList(counts.keep)
#barplot(y$samples$lib.size)

par(mfrow = c(1,1))
logcpm <- cpm(y$counts, log = TRUE)
boxplot(logcpm, xlab ="", ylab = "Log2 counts per million", las =2, outline = FALSE)
abline(h = median(logcpm), col = "blue")
title("Boxplots of logCPMs (unnormalized)")

par(mfrow = c(1,2), oma = c(2,0,0,0))
group.col <- c("red", "blue")[targets$Group]
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by groups)",cex.main=0.8)

lib.col <- c("light pink","light green")[targets$Library]
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2, col=lib.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by library prep)",cex.main=0.8)

#MDS plot
par(mfrow = c(1,1))
plotMDS(y)

par(mfrow=c(1,2))
plotMDS(y,col=group.col)
legend("topright",legend=levels(targets$Group),fill=c("red","blue"))
plotMDS(y,col=lib.col)
legend("topleft",legend=levels(targets$Library),fill=c("light pink","light green"))

#Hierarchical clustering 

logcounts <- cpm(y, log = TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing = TRUE))[1:500]
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)

#plot heatmap

heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=group.col,scale="row",margins=c(10,5))

heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=lib.col,scale="row",margins=c(10,5))

#Normalization
norm.y <- calcNormFactors(y)
norm.y$samples 

par(mfrow=c(1,2))
plotMD(logcounts,column=2)
abline(h=0,col="grey")
plotMD(norm.y,column = 2)
abline(h=0,col="grey")



abline(h=median(norm.y),col=4)
title("Boxplots of log counts \n(after normalization)",cex.main=0.8)
       
#Differential expression analysis

design <- model.matrix(~targets$Library + targets$Group)
design

colnames(design) <- c("Int", "SEvsPE", "UVsT")

#voom transform the data

par(mfrow = c(1,1))
v <- voom(norm.y, design, plot =TRUE)
v
par(mfrow=c(1,2))
boxplot(logcounts)
abline(h=median(logcounts),col=4)

boxplot(v$E)
abline(h=median(v$E),col=4)
title("Boxplots of log counts \n(after normalization)",cex.main=0.8)
#test for differential expression

fit <- lmFit(v, design)
fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)

gen <- topTable(fit, coef = 3, sort.by = "p", n = 5)
gen <- write.csv(gen,"gen.csv" )
rownames(fit)
#add annotation 
columns(org.Dm.eg.db)

ann <- select(org.Dm.eg.db, keys =rownames(fit), columns = c("FLYBASE", "ENTREZID", "SYMBOL","GENENAME"),keytype = "FLYBASE" )
head(ann)

table(ann$FLYBASE==rownames(fit))
fit$genes <- ann

topTable(fit, coef = 3, sort.by = "p")
ls("package:org.Dm.eg.db")
fly <- toTable(org.Dm.egFLYBASE)
head(fly)

symbol <- toTable(org.Dm.egSYMBOL)
genename <- toTable(org.Dm.egGENENAME)

ann1 <- merge(fly, symbol, by = "gene_id")
head(ann1)

ann2 <- merge(ann1, genename, by = "gene_id")
head(ann2)

m <- match(rownames(fit),ann2$flybase_id)
table(is.na(m)) # check for unmatched rows

ann3 <- ann2[m[!is.na(m)],]
head(ann3)

head(fit$genes)
topTable(fit, coef = 3, sort.by = "p")

#check expression of pasilla
ps <- grep("pasilla", fit$genes$GENENAME)
topTable(fit[ps,], coef = 3)

#plots after testing for DE
par(mfrow = c(1,2))
plotMD(fit, coef = 3, status = results[,"UVsT"])
volcanoplot(fit,coef=3,highlight=100,names=fit$genes$SYMBOL)
stripchart(v$E["FBgn0025111",]~targets$Group)
stripchart(v$E["FBgn0261552",]~targets$Group)

#testing relative to a threshold

fit.treat <- treat(fit, lfc = 1)
res.treat <- decideTests(fit.treat)
summary(res.treat)

topTreat(fit.treat, coef = 3)

plotMD(fit.treat,coef=3,status=res.treat[,"UVsT"])
abline(h=0,col="grey")

#Gene set testing with goana
go <- goana(fit, coef="UVsT",species = "Dm", geneid="ENTREZID")
topGO(go, n=10)
