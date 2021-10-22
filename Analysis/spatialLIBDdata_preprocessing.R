#-------------------------------------------------
# LIBD Human Dorsolateral Pre-frontal Cortex Data
#-------------------------------------------------

library(spatialLIBD)
library(scry)

ehub <- ExperimentHub::ExperimentHub()
if (!exists("sce")) sce <- fetch_data(type = "sce", eh = ehub)

# ---select the area of interest and the sample
minX <- 250
maxX <- max(sce$imagerow)
minY <- min(sce$imagecol)
maxY <- 370
selectedCells <- sce$imagerow >= minX &
    sce$imagerow <= maxX &
    sce$imagecol >= minY &
    sce$imagecol <= maxY
sce <- sce[,selectedCells]

# ---select the subject
selectedSample <- "151673"
sce <- sce[,sce$sample_name == selectedSample]

# ---plot the spots in the selected area
spe <- sce_to_spe(sce)
vis_clus(spe = spe,
         sampleid = "151673",
         clustervar = "layer_guess_reordered",#"cur.Ds",
         colors = libd_layer_colors,#mypalette,
         spatial = T)+ggplot2::ggtitle("")

# ---selection of the genes using the scry package
# 1. sort the genes according to the deviance
sce <- devianceFeatureSelection(sce, assay="counts", sorted=TRUE)
plot(rowData(sce)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="deviance", main="Gene Selection with Deviance")
abline(v=1000, lty=2, col="red")
abline(v=200, lty=2, col=4)
legend(2000,30000, legend = c("Ideal cut-off","Actual cut-off"), lty = c(2,2), col = c(4,2))

# 2. take the deviance residuals and selects the first 1000 informative genes
sce <- nullResiduals(sce, assay="counts", fam = "binomial", type = "deviance")
sce <- sce[1:1000,]

# 3. save the results
Data <- list()
Data$x <- as.matrix(sce@assays@data$binomial_deviance_residuals)
Data$coordinates <- cbind(sce$imagecol, -sce$imagerow)
Data$gene.names <- rownames(sce@assays@data$binomial_deviance_residuals)
Data$Layers <- sce@colData$layer_guess_reordered
save(Data, file = ".../LIBD_subj151673_1585spots.Rdata")
