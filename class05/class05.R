#' ---
#' title: "Class 5: Introfuction to R graphics"
#' author: "Karen Li"
#' date: "January 22, 2019"
#' output: github_document
#' ---

#Section 2A
read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
plot(weight, typ="o", pch=15, cex=1.5, lwd=2, 
     col="red", ylim=c(2,10), xlab="Age (months)",
     ylab="Weight (kg)", main="Weight vs. Age")

# Section 2B
feat <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, 
                   sep="\t")
par(mar=c(5,11,4,2))
barplot(feat$Count, horiz=TRUE, names.arg=feat$Feature,
        main="Number of Features in the Mouse GRCm38 Genome", las=1, 
        col=heat.colors(12))

#Section 2C


# Section 3A
gender <- read.table("bimm143_05_rstats/male_female_counts.txt", 
                     header=TRUE, sep="\t")
par(mar=c(7,4,4,4))
barplot(gender$Count, names.arg=gender$Sample, las=2, 
        col=rainbow(nrow(gender)))

# Section 3B
genes <- read.table("bimm143_05_rstats/up_down_expression.txt", header=TRUE)
table(genes$State)
palette(c("red", "gray", "blue"))
plot(genes$Condition1, genes$Condition2, col=genes$State)



