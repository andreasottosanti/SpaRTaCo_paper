#------------------------------------------------
# LIBD Human Dorsolateral Prefrontal Cortex Data
#------------------------------------------------

library(spatialLIBD)
library(scry)

ehub <- ExperimentHub::ExperimentHub()
if (!exists("sce")) sce <- fetch_data(type = "sce", eh = ehub)

# ---select the subject
selectedSample <- "151673"
sce <- sce[,sce$sample_name == selectedSample]

# ---pre-selection based on the number of counts in genes and cells
minPerGene <- 300
sce <- sce[which(rowSums(counts(sce)) >= minPerGene),]

# ---plot the spots in the selected area
spe <- sce_to_spe(sce)
vis_clus(spe = spe,
         sampleid = "151673",
         clustervar = "layer_guess_reordered",
         colors = libd_layer_colors,
         spatial = T)+ggplot2::ggtitle("")

# ---selection of the genes using the scry package
# 1. sort the genes according to the deviance
sce <- devianceFeatureSelection(sce, assay="counts", sorted=TRUE)
plot(rowData(sce)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="deviance", main="Gene Selection with Deviance")
abline(v=500, lty=2, col="red")
abline(v=200, lty=2, col=4)
legend(4000,70000, legend = c("Ideal cut-off","Actual cut-off"), lty = c(2,2), col = c(4,2))

# 2. take the deviance residuals and selects the first 500 informative genes
sce <- nullResiduals(sce, assay="counts", fam = "binomial", type = "deviance")
sce <- sce[1:500,]

# 3. save the results
Data <- list()
Data$x <- as.matrix(sce@assays@data$binomial_deviance_residuals)
Data$counts <- as.matrix(sce@assays@data$counts)
Data$coordinates <- cbind(sce$imagerow, sce$imagecol)
Data$gene.id <- rowData(sce)$gene_id
Data$gene.names <- rowData(sce)$gene_name
Data$layers <- sce@colData$layer_guess
save(Data, file = ".../LIBD_subj151673.Rdata")
