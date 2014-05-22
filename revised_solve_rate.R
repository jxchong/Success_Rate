library(RColorBrewer)
colors = brewer.pal(12, "Set3")
# blue, orange, purple, gray, green
bestcolors=colors[c(9,7,5,6,10)]
par.default <- par()

color.mapping_success <- hsv(h=rgb2hsv(col2rgb(bestcolors[2]))['h',], s=rgb2hsv(col2rgb(bestcolors[2]))['s',], v=.6)
color.mapping_unsuccess <- hsv(h=rgb2hsv(col2rgb(bestcolors[1]))['h',], s=rgb2hsv(col2rgb(bestcolors[1]))['s',], v=.7)
colors.AD_comb <- c(bestcolors[2], color.mapping_success, bestcolors[1], color.mapping_unsuccess)
colors.novelty <- c(hsv(h=rgb2hsv(col2rgb(bestcolors[3]))['h',], s=rgb2hsv(col2rgb(bestcolors[3]))['s',], v=.6), bestcolors[3])



#########
# read in data
#########

setwd("/Users/jxchong/Dropbox/Postdoc/CMG/UWCMG/Solve rate")
# data <- read.table("CMG_Progress_Report_2013_12_solve_rate.txt", head=TRUE, sep="\t")
# data <- read.table("CMG_Progress_Report_2014_02_success_rate.txt", head=TRUE, sep="\t")
data <- read.table("CMG_Progress_Report_2014_05_success_rate.txt", head=TRUE, sep="\t")
filedate <- "2014-05"


modelnames <- c("AD", "AR", "AR\nconsang.", "de novo", "unknown", "X-linked", "Total")

include <- subset(data, subset=(overall=="first_pass_done"|overall=="complete"))
nocohorts <- subset(include, subset=cohort=="N")



######################################################################################## 
# Success rate (by-phenotype/gene without inheritance model)
# 1) for each phenotype, separate success and unsuccessful kindreds
# 2) group successful kindreds by causal gene and count genes
# 3) group all unsolved kindreds as a single unsuccessful phenotype

# NOTE: need to manually add in the double hits
######################################################################################## 


success_v_model <- table(subset(nocohorts, select=c("Solved", "Model")))

unsuccess_mapping <- table(subset(nocohorts, Solved=="N", select=c("Mapping_Data", "Model")))
success_mapping <- table(subset(nocohorts, Solved=="Y", select=c("Mapping_Data", "Model")))
total.unsuccess_mapping <- table(subset(nocohorts, Solved=="N", select=c("Mapping_Data")))
total.success_mapping <- table(subset(nocohorts, Solved=="Y", select=c("Mapping_Data")))

success_v_model.mapping <- data.frame( matrix(data=0, nrow=4, ncol=length(modelnames)), row.names=c("Solved, -mapping data", "Solved, +mapping data", "Unsolved, -mapping data", "Unsolved, +mapping data"))
names(success_v_model.mapping) <- modelnames
success_v_model.mapping[c("Unsolved, -mapping data","Unsolved, +mapping data"),] <- cbind(unsuccess_mapping, total.unsuccess_mapping)
success_v_model.mapping[c("Solved, -mapping data","Solved, +mapping data"),] <- cbind(success_mapping, total.success_mapping)

success_v_model.mapping.ordered <- success_v_model.mapping[,rev(c(7,5,1,4,2,3,6))]
labels <- c("Solved, -mapping data", "Solved, +mapping data", "Unsolved, -mapping data", "Unsolved, +mapping data")

######## adjust for multiple hits
success_v_model.mapping.ordered["Solved, -mapping data","AR"] = success_v_model.mapping.ordered["Solved, -mapping data","AR"]+1
########

pdf(paste0("success_v_model.",filedate,".pdf"), width=10, height=6)
par(xpd=T, mar=c(5,4,4,10))
bars.x <- barplot(as.matrix(success_v_model.mapping.ordered), beside=FALSE, col=colors.AD_comb, ylab="# of Phenotypes (no cohorts)", xlab="Inheritance Model", border=NA, cex.main=0.95, cex.names=0.9)
legend(8.5, 20, col=colors.AD_comb, legend=labels, border=NA, fill=colors.AD_comb, bty="n", cex=0.8, xpd=TRUE)
text(x=bars.x, y=colSums(success_v_model.mapping.ordered), labels=paste0(colSums(success_v_model.mapping.ordered[c(1:2),]), ":", as.vector(colSums(success_v_model.mapping.ordered)), " (", round(colSums(success_v_model.mapping.ordered[c(1:2),])/colSums(success_v_model.mapping.ordered), digits=2), ")"), xpd=TRUE, pos=3, cex=0.8)
dev.off()




## with NHLBI cohorts
success_v_model.withNHBLI <- table(subset(include, select=c("Solved", "Model")))

unsuccess_mapping.withNHBLI <- table(subset(include, Solved=="N", select=c("Mapping_Data", "Model")))
success_mapping.withNHBLI <- table(subset(include, Solved=="Y", select=c("Mapping_Data", "Model")))
total.unsuccess_mapping.withNHBLI <- table(subset(include, Solved=="N", select=c("Mapping_Data")))
total.success_mapping.withNHBLI <- table(subset(include, Solved=="Y", select=c("Mapping_Data")))

success_v_model.mapping.withNHBLI <- data.frame( matrix(data=0, nrow=4, ncol=length(modelnames)), row.names=c("Solved, -mapping data", "Solved, +mapping data", "Unsolved, -mapping data", "Unsolved, +mapping data"))
names(success_v_model.mapping.withNHBLI) <- modelnames
success_v_model.mapping.withNHBLI[c("Unsolved, -mapping data","Unsolved, +mapping data"),] <- cbind(unsuccess_mapping.withNHBLI, total.unsuccess_mapping.withNHBLI)
success_v_model.mapping.withNHBLI[c("Solved, -mapping data","Solved, +mapping data"),] <- cbind(success_mapping.withNHBLI, total.success_mapping.withNHBLI)

success_v_model.mapping.ordered.withNHBLI <- success_v_model.mapping.withNHBLI[,rev(c(7,5,1,4,2,3,6))]
labels.withNHBLI <- c("Solved, -mapping data", "Solved, +mapping data", "Unsolved, -mapping data", "Unsolved, +mapping data")

# adjust for multiple hits
success_v_model.mapping.ordered.withNHBLI["Solved, -mapping data","AR"] = success_v_model.mapping.ordered.withNHBLI["Solved, -mapping data","AR"]+1
#

pdf(paste0("success_v_model.withNHLBI.",filedate,".pdf"), width=10, height=6)
par(xpd=T, mar=c(5,4,4,10))
bars.x <- barplot(as.matrix(success_v_model.mapping.ordered.withNHBLI), beside=FALSE, col=colors.AD_comb, ylab="# of Phenotypes", xlab="Inheritance Model", border=NA, cex.main=0.95, cex.names=0.9)
legend(8.5, 20, col=colors.AD_comb, legend=labels.withNHBLI, border=NA, fill=colors.AD_comb, bty="n", cex=0.8, xpd=TRUE)
text(x=bars.x, y=colSums(success_v_model.mapping.ordered.withNHBLI), labels=paste0(colSums(success_v_model.mapping.ordered.withNHBLI[c(1:2),]), ":", as.vector(colSums(success_v_model.mapping.ordered.withNHBLI)), " (", round(colSums(success_v_model.mapping.ordered.withNHBLI[c(1:2),])/colSums(success_v_model.mapping.ordered.withNHBLI), digits=2), ")"), xpd=TRUE, pos=3, cex=0.8)
dev.off()



######################################################################################## 
# NOVEL gene success rate (by-phenotype/gene with inheritance model)
# 1) for each phenotype, separate success and unsuccess kindreds
# 2) group successful kindreds by causal gene and count genes
# 3) group all unsovled kindreds as a single unsuccess phenotype

# NOTE: need to manually add in the double hits
######################################################################################## 

novelty_v_model <- table(subset(nocohorts, select=c("isNovel", "Model")))
novel.totals <- rowSums(novelty_v_model)

summary.novelty_v_model <- data.frame( matrix(data=cbind(novelty_v_model, novel.totals), nrow=2, ncol=length(modelnames)), row.names=c("Not Novel", "Novel"))
names(summary.novelty_v_model) <- modelnames

summary.novelty_v_model.ordered <- summary.novelty_v_model[order(nrow(summary.novelty_v_model):1),rev(c(7,5,1,4,2,3,6))]

# adjust for multiple hits
summary.novelty_v_model.ordered["Novel", "AR"] = summary.novelty_v_model.ordered["Novel", "AR"]+1
#


pdf(paste0("novelty.successes_v_model.",filedate,".pdf"), width=10, height=6)
par(xpd=T, mar=c(5,4,4,6), oma=c(0,0,0,4))
bars.x <- barplot(as.matrix(summary.novelty_v_model.ordered), beside=FALSE, col=colors.novelty, ylab="# of Solved Phenotypes (no cohorts)", xlab="Inheritance Model", border=NA, cex.main=0.95, cex.names=0.9)
legend(8.5, 20, col=colors.novelty, legend=rownames(summary.novelty_v_model.ordered), border=NA, fill=colors.novelty, bty="n", cex=0.8, xpd=TRUE)
text(x=bars.x, y=colSums(summary.novelty_v_model.ordered), labels=paste0(summary.novelty_v_model.ordered["Novel",], " of ", colSums(summary.novelty_v_model.ordered), " (", round(summary.novelty_v_model.ordered["Novel",]/colSums(summary.novelty_v_model.ordered), digits=2), ")"), xpd=TRUE, pos=3)
dev.off()

## with NHLBI cohorts
novelty_v_model.withNHBLI <- table(subset(include, select=c("isNovel", "Model")))
novel.totals.withNHBLI <- rowSums(novelty_v_model.withNHBLI)

summary.novelty_v_model.withNHBLI <- data.frame( matrix(data=cbind(novelty_v_model.withNHBLI, novel.totals.withNHBLI), nrow=2, ncol=length(modelnames)), row.names=c("Not Novel", "Novel"))
names(summary.novelty_v_model.withNHBLI) <- modelnames

summary.novelty_v_model.ordered.withNHBLI <- summary.novelty_v_model.withNHBLI[order(nrow(summary.novelty_v_model.withNHBLI):1),rev(c(7,5,1,4,2,3,6))]

# adjust for multiple hits
summary.novelty_v_model.ordered.withNHBLI["Novel", "AR"] = summary.novelty_v_model.ordered.withNHBLI["Novel", "AR"]+1

pdf(paste0("novelty.successes_v_model.withNHLBI.",filedate,".pdf"), width=10, height=6)
par(xpd=T, mar=c(5,4,4,6), oma=c(0,0,0,4))
bars.x <- barplot(as.matrix(summary.novelty_v_model.ordered.withNHBLI), beside=FALSE, col=colors.novelty, ylab="# of Solved Phenotypes", xlab="Inheritance Model", border=NA, cex.main=0.95, cex.names=0.9)
legend(8.5, 20, col=colors.novelty, legend=rownames(summary.novelty_v_model.ordered.withNHBLI), border=NA, fill=colors.novelty, bty="n", cex=0.8, xpd=TRUE)
text(x=bars.x, y=colSums(summary.novelty_v_model.ordered.withNHBLI), labels=paste0(summary.novelty_v_model.ordered.withNHBLI["Novel",], " of ", colSums(summary.novelty_v_model.ordered.withNHBLI), " (", round(summary.novelty_v_model.ordered.withNHBLI["Novel",]/colSums(summary.novelty_v_model.ordered.withNHBLI), digits=2), ")"), xpd=TRUE, pos=3)
dev.off()



######################################################################################## 
# Successful diagnostic rate (by kindred)
########################################################################################

# with cohorts
# nkindreds.withcohorts.success <- sum(subset(include, select=numerator.kindreds, subset=Solved=="Y"))
nkindreds.withcohorts.unsuccess <- sum(subset(include, select=numerator.kindreds, subset=Solved=="N"))
nkindreds.withcohorts.success.novel <- sum(subset(include, select=numerator.kindreds, subset=(Solved=="Y"&isNovel=="Y")))
nkindreds.withcohorts.success.notnovel <- sum(subset(include, select=numerator.kindreds, subset=(Solved=="Y"&isNovel=="N")))
nkindreds.withcohorts <- as.matrix(rbind(nkindreds.withcohorts.success.novel, nkindreds.withcohorts.success.notnovel, nkindreds.withcohorts.unsuccess))

# without cohorts
# nkindreds.nocohorts.success <- sum(subset(nocohorts, select=numerator.kindreds, subset=Solved=="Y"))
nkindreds.nocohorts.unsuccess <- sum(subset(nocohorts, select=numerator.kindreds, subset=Solved=="N"))
nkindreds.nocohorts.success.novel <- sum(subset(nocohorts, select=numerator.kindreds, subset=(Solved=="Y"&isNovel=="Y")))
nkindreds.nocohorts.success.notnovel <- sum(subset(nocohorts, select=numerator.kindreds, subset=(Solved=="Y"&isNovel=="N")))
nkindreds.nocohorts <- as.matrix(rbind(nkindreds.nocohorts.success.novel, nkindreds.nocohorts.success.notnovel, nkindreds.nocohorts.unsuccess))

successrate.bykindred <- as.matrix(cbind(nkindreds.withcohorts, nkindreds.nocohorts))
rownames(successrate.bykindred) <- c("Novel", "Not Novel", "Unsolved")
successrate.bykindred.ratios <- round(as.matrix(rbind(successrate.bykindred[1,]/colSums(successrate.bykindred),	successrate.bykindred[2,]/colSums(successrate.bykindred))), 2)
successrate.bykindred.ratios.overall <- round(colSums(successrate.bykindred[1:2,])/colSums(successrate.bykindred), 2)

successrate.bykindred.barplot.blocks.ypos <- apply(successrate.bykindred, 2, cumsum)
successrate.bykindred.barplot.blocks.ypos <- successrate.bykindred.barplot.blocks.ypos - successrate.bykindred/2
# successrate.bykindred.barplot.blocks.ypos <- t(successrate.bykindred.barplot.blocks.ypos)
	
	
# make graph
pdf(paste0("success_diagnostic_rate.bykindred.",filedate,".pdf"), width=4, height=6)
par(xpd=T, mar=c(5,4,4,6))
bars.x <- barplot(successrate.bykindred, col=c(colors.novelty, colors.AD_comb[3]), names.arg=c("With Cohorts", "No Cohorts"), ylab="# of Kindreds", xlab="Diagnostic Rate (by kindred)", border=NA, cex.main=0.85, cex.names=0.7)
legend(2.7, max(colSums(successrate.bykindred))/2, col=c(colors.novelty, colors.AD_comb[3]), legend=rownames(successrate.bykindred), border=NA, fill=c(colors.novelty, colors.AD_comb[3]), bty="n", cex=0.8, xpd=TRUE)
text(x=bars.x, y=colSums(successrate.bykindred), labels=paste0(colSums(successrate.bykindred[1:2,]), " of ", colSums(successrate.bykindred), " (", successrate.bykindred.ratios.overall, ")"), xpd=TRUE, pos=3, cex=0.8)
text(x=bars.x, y=t(successrate.bykindred.barplot.blocks.ypos[1:2,]), labels=paste0("N=", t(successrate.bykindred[1:2,]), " (", t(successrate.bykindred.ratios), ")"), cex=0.7)
dev.off()



