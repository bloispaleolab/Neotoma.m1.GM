###README (last update 3-31-2025)
  #This is the code to reproduce data shown in Tables 2-4, Figures 3-4, and
  #Supplementary Figures 1-2 of Fox and Blois (2025) "A geometric 
  #morphometric approach to identifying recent and fossil woodrat molars with 
  #remarks on late Pleistocene Neotoma macrotis from Rancho La Brea". 
  #Please refer to the manuscript for additional context.

###load required packages

##Un-comment to install packages
  #install.packages('geomorph')
  #install.packages('RRPP')
  #install.packages('adegenet')
library(geomorph)
library(RRPP)
library(adegenet)

###Data Preparation 
  #read landmark data files, perform GPA-transformation, add group 
  #classifiers, and reformat data for Discriminant Analysis of Principal
  #Components (DAPC) 

##First landmark repetition with recent specimens, fossils, and N. cinerea 
  #specimens east of California 
data1 <- readland.tps("InputData/Neotoma_combined_R1.tps", specID = "imageID")
  #Examine data before GPA transformation
plotAllSpecimens(data1, mean=TRUE, label=TRUE)
plotOutliers(data1)
  #Perform GPA
Neotoma.gpa_1 <- gpagen(data1)
  #Reexamine data after GPA transformation
plotAllSpecimens(Neotoma.gpa_1$coords, label= TRUE)
plotOutliers(Neotoma.gpa_1$coords)

##Second landmark repetition with recent, fossil, and east N. cinerea specimens
data2 <- readland.tps("InputData/Neotoma_combined_R2.tps", specID = "imageID")
plotAllSpecimens(data2, mean=TRUE, label=TRUE)
plotOutliers(data2)
Neotoma.gpa_2 <- gpagen(data2)
plotAllSpecimens(Neotoma.gpa_2$coords, label= TRUE)
plotOutliers(Neotoma.gpa_2$coords)

##First landmark repetition with unworn and extremely worn teeth removed 
data3 <- readland.tps("InputData/WearVetted_R1.tps", specID = "imageID")
plotAllSpecimens(data3, mean=TRUE, label=TRUE)
plotOutliers(data3)
Neotoma.gpa_3 <- gpagen(data3)
plotAllSpecimens(Neotoma.gpa_3$coords, label= TRUE)
plotOutliers(Neotoma.gpa_3$coords)

##Second landmark repetition with unworn and extremely worn teeth removed
data4 <- readland.tps("InputData/WearVetted_R2.tps", specID = "imageID")
plotAllSpecimens(data4, mean=TRUE, label=TRUE)
plotOutliers(data4)
Neotoma.gpa_4 <- gpagen(data4)
plotAllSpecimens(Neotoma.gpa_4$coords, label= TRUE)
plotOutliers(Neotoma.gpa_4$coords)


##Add classifier for all recent, fossil, and east N. cinerea specimens
classifier <- read.csv("InputData/DAPC_classifier.csv", header=T)

##Add classifier for wear-vetted specimens 
Wear_classifier <- read.csv("InputData/Wear_classifier.csv", header=T)

##Reformat first landmark data repetition for DAPC
Neotoma.coords1<-as.data.frame(Neotoma.gpa_1$coords)
Neotoma.coords1<-t(Neotoma.coords1)
Neotoma.coords1<-cbind(Neotoma.coords1[c(TRUE, FALSE), ],
                       Neotoma.coords1[c(FALSE, TRUE), ])
MergedData1 <- cbind(classifier, Neotoma.coords1)

##Reformat second landmark data repetition for DAPC
Neotoma.coords2<-as.data.frame(Neotoma.gpa_2$coords)
Neotoma.coords2<-t(Neotoma.coords2)
Neotoma.coords2<-cbind(Neotoma.coords2[c(TRUE, FALSE), ],
                       Neotoma.coords2[c(FALSE, TRUE), ])
MergedData2 <- cbind(classifier, Neotoma.coords2)

##Reformat first wear-vetted landmark data repetition for DAPC
Neotoma.coords3<-as.data.frame(Neotoma.gpa_3$coords)
Neotoma.coords3<-t(Neotoma.coords3)
Neotoma.coords3<-cbind(Neotoma.coords3[c(TRUE, FALSE), ],
                       Neotoma.coords3[c(FALSE, TRUE), ])
MergedData3 <- cbind(Wear_classifier, Neotoma.coords3)

##Reformat second wear-vetted landmark data repetition for DAPC
Neotoma.coords4<-as.data.frame(Neotoma.gpa_4$coords)
Neotoma.coords4<-t(Neotoma.coords4)
Neotoma.coords4<-cbind(Neotoma.coords4[c(TRUE, FALSE), ],
                       Neotoma.coords4[c(FALSE, TRUE), ])
MergedData4 <- cbind(Wear_classifier, Neotoma.coords4)

##Save output files for DAPC
  #First landmark repetition with all specimens 
write.csv(MergedData1[1:199,], 'OutputData/Extant_R1.csv') 
write.csv(MergedData1[200:228,], 'OutputData/Unknowns_R1.csv')

  #Second landmark repetition with all specimens 
write.csv(MergedData2[1:199,], 'OutputData/Extant_R2.csv') 
write.csv(MergedData2[200:228,], 'OutputData/Unknowns_R2.csv')

  #First landmark repetition with wear-vetted specimens 
write.csv(MergedData3[1:171,], 'OutputData/Extant_wear_R1.csv') 

  #Second landmark repetition with wear-vetted specimens 
write.csv(MergedData4[1:171,], 'OutputData/Extant_wear_R2.csv')


###Measurement error (ME) analysis

##read-in LM data, excluding fossil and east N. cinerea 'unknowns'
data5 <- readland.tps("InputData/Neotoma_extant_R1.tps", specID = "imageID")
extant.gpa_1 <- gpagen(data5)
data6 <- readland.tps("InputData/Neotoma_extant_R2.tps", specID = "imageID")
extant.gpa_2 <- gpagen(data6)
  #format shape data
shape1 <- extant.gpa_1$coords
shape1<-two.d.array(shape1)
shape2 <- extant.gpa_2$coords
shape2<-two.d.array(shape2)
  #Match specimen row names of first and second landmark repetitions
specimen_names_shape1 <- rownames(shape1)
shape2_reordered <- matrix(NA, nrow = nrow(shape2), ncol = ncol(shape2))
matched_row_names <- character(length = length(specimen_names_shape1))
for (i in 1:length(specimen_names_shape1)) {
  matching_index <- match(specimen_names_shape1[i], rownames(shape2))
  if (!is.na(matching_index)) {
    shape2_reordered[i, ] <- shape2[matching_index, ]
    matched_row_names[i] <- rownames(shape2)[matching_index]
  }
}
shape2_reordered <- shape2_reordered[complete.cases(shape2_reordered), ]
  #Combine the reordered data
coords <- rbind(shape1, shape2_reordered)
rownames(coords) <- c(rownames(shape1), matched_row_names)

##add parameters for ME analysis
RRPP_classifer <- read.csv("InputData/RRPP_classifier.csv", header=T)
groups<-RRPP_classifer$species
reps<-RRPP_classifer$reps
subj<-RRPP_classifer$subj
  #Create RRPP list
rdf <- rrpp.data.frame(coords = coords, groups = groups, reps = reps, subj 
  = subj)

##ME analysis without species groups
ME1 <- measurement.error(
  Y = "coords",
  subjects = "subj",
  replicates = "reps",
  data = rdf)
  #Run ANOVA and ICC without species groups (This is the top part of 
  #Fox and Blois Table 2)
ICCstats(ME1, subjects = "Subjects", with_in = "Systematic ME")

##ME analysis with species groups
ME2 <- measurement.error(
  Y = "coords",
  subjects = "subj",
  replicates = "reps",
  groups = "groups",
  data = rdf)
  #Run ANOVA and ICC with species groups (This is the bottom part of 
  #Fox and Blois Table 2)
ICCstats(ME2, subjects = "Subjects",
         with_in = "Systematic ME", groups = "groups")


###DAPC

##read outputs for extant species training data generated above (lines 96-107)
data7<-read.table('OutputData/extant_R1.csv', header=TRUE, row.names = 1, sep=",")
DAPC_extant_T1 <- data7[c(2:29)]
data8<-read.table('OutputData/extant_R2.csv', header=TRUE, row.names = 1, sep=",")
DAPC_extant_T2 <- data8[c(2:29)]

##read outputs for classifying fossil and east N. cinerea "unknowns"
data9<-read.table('OutputData/Unknowns_R1.csv', header=TRUE, row.names = 1, sep=",")
Unknowns_T1 <- data9[c(2:29)]
data10<-read.table('OutputData/Unknowns_R2.csv', header=TRUE, row.names = 1, sep=",")
Unknowns_T2 <- data10[c(2:29)]

##read outputs for wear-vetted data
data11<-read.table('OutputData/extant_wear_R1.csv', header=TRUE, row.names = 1, sep=",")
DAPC_wear_T1 <- data11[c(2:29)]
data12<-read.table('OutputData/extant_wear_R2.csv', header=TRUE, row.names = 1, sep=",")
DAPC_wear_T2 <- data12[c(2:29)]

##Run cross-validation to identify the number of PCs to retain in DAPC
##(change to dataset of interest directly above)
xvalDapc(DAPC_extant_T1, grp = data7$species, n.pca.max = 24, training.set = 0.9,
         result = "groupMean", center = TRUE, scale = FALSE,
         n.pca = NULL, n.rep = 30, xval.plot = TRUE)
  #Our observed range of lowest MSE output values was 6-22 after dozens of
  #iterations across both data repetitions. We then ran the 'dapc' function 
  #below with a 'n.pca=' value of 6-22. We find that classification stats 
  #("summary(dapc1) $assign.prop") improve with each additional PC retained up 
  #to ~16 and plateau thereafter. "Mean Successful Assignment by Number of PCs 
  #of PCA" in xvalDapc shows similar trends of plateauing successful assignments 
  #around 16 PCs.16 PCs were therefore selected in the final 'dapc' models below.
  #Please see the manuscript and Jombart and Collins (2022) reference within for
  #more information. 

##run DAPC with 16 PCs retained
  #model 1
dapc1 <- dapc(DAPC_extant_T1, grp = data7$species, n.da = 4, n.pca=16)
summary(dapc1)
  #model 2
dapc2 <- dapc(DAPC_extant_T2, grp = data8$species, n.da = 4, n.pca=16)
summary(dapc2)
  #wear-vetted model 1
dapc3 <- dapc(DAPC_wear_T1, grp = data11$species, n.da = 4, n.pca=16)
summary(dapc3)
  #wear-vetted model 2
dapc4 <- dapc(DAPC_wear_T2, grp = data12$species, n.da = 4, n.pca=16)
summary(dapc4)

  #visualize per-specimen misclassification by species (change to the dataset of 
  #interest in directly above)
assignplot(dapc1, subset=1:199)

  #Plot classification table (change to the dataset of interest above)
predicted_groups <- predict(dapc1)
predicted_assignments <- predicted_groups$assign
true_groups <- data7$species #make sure to change this to the appropriate dataset
  #Create a confusion matrix
conf_matrix <- table(true_groups, predicted_assignments)
  #Print the confusion matrix
  #this is part of Fox and Blois Table 3 (rerun with 'dapc2' iteration to 
  #generate data for the bottom half)
print(conf_matrix)

##Predict species memberships of fossil unknowns and generate posterior 
##probabilities of extant species-group membership (Fox and Blois Table 4)
  #Model 1 fossil classification
pred.Unknown_1 <- predict.dapc(dapc1, newdata=Unknowns_T1)
pred.Unknown_1$assign
round(pred.Unknown_1$posterior[1:29, 1:5],2) 

  #Model 2 fossil classification
pred.Unknown_2 <- predict.dapc(dapc2, newdata=Unknowns_T2)
pred.Unknown_2$assign
round(pred.Unknown_2$posterior[1:29, 1:5],2) 

##DAPC scatterplots with fossil and east N. cinerea unknowns
##(This is Fox and Blois Figure 3)
extract_groups <- function(pred_data, ucmvz_id = NULL) {
  # Extract row names
  row_names <- rownames(pred_data$ind.scores)
  # Indices for each group based on prefix
  P23_indices <- grepl("^X\\.P23", row_names)
  SDSM_indices <- grepl("^X\\.SDSM", row_names)
  UCMVZ_indices <- if (!is.null(ucmvz_id)) {
    row_names == ucmvz_id
  } else {
    FALSE
  }
  # Subset the data for each group
  list(
    P23 = list(
      ind.scores = pred_data$ind.scores[P23_indices, , drop = FALSE],
      posterior = pred_data$posterior[P23_indices, , drop = FALSE],
      assign = pred_data$assign[P23_indices]
    ),
    SDSM = list(
      ind.scores = pred_data$ind.scores[SDSM_indices, , drop = FALSE],
      posterior = pred_data$posterior[SDSM_indices, , drop = FALSE],
      assign = pred_data$assign[SDSM_indices]
    ),
    UCMVZ = if (any(UCMVZ_indices)) {
      list(
        ind.scores = pred_data$ind.scores[UCMVZ_indices, , drop = FALSE],
        posterior = pred_data$posterior[UCMVZ_indices, , drop = FALSE],
        assign = pred_data$assign[UCMVZ_indices]
      )
    } else {
      NULL
    }
  )
}
  #For model 1
groups_1 <- extract_groups(pred.Unknown_1, ucmvz_id = "X.UCMVZ-ID-NC-51944L")
pred.Unknown_1_P23 <- groups_1$P23
pred.Unknown_1_SDSM <- groups_1$SDSM
pred.Unknown_1_UCMVZ <- groups_1$UCMVZ
  #For model 2
groups_2 <- extract_groups(pred.Unknown_2, ucmvz_id = "X.UCMVZ-ID-NC-51944L")
pred.Unknown_2_P23 <- groups_2$P23
pred.Unknown_2_SDSM <- groups_2$SDSM
pred.Unknown_2_UCMVZ <- groups_2$UCMVZ
  #Set layout for two plots side-by-side
pdf("Figures/Figure3.pdf", width = 10, height = 5)
  #Plot dapc1
par(mfrow = c(1, 2))  # 1 row, 2 columns
group.col <- c("blue", "plum4", "limegreen", "orange1", "orangered")
scatter(dapc1, scree.da = FALSE, scree.pca = FALSE, col = group.col)
  #Overlay each subset of unknown specimens with different symbols/colors
points(pred.Unknown_1_P23$ind.scores[, 1], pred.Unknown_1_P23$ind.scores[, 2], 
       pch = 18, col = transp("black"), cex = 1.5)
points(pred.Unknown_1_SDSM$ind.scores[, 1], pred.Unknown_1_SDSM$ind.scores[, 2], 
       pch = 15, col = transp("purple"), cex = 1)
points(pred.Unknown_1_UCMVZ$ind.scores[, 1], pred.Unknown_1_UCMVZ$ind.scores[, 2], 
       pch = 15, col = transp("purple"), cex = 1)
myInset <- function() {
  temp <- dapc1$pca.eig
  temp <- 100 * cumsum(temp) / sum(temp)
  plot(temp, col = rep(c("black", "lightgrey"),
                       c(dapc1$n.pca, 24)), ylim = c(0, 100),
       xlab = "PCA axis", ylab = "Variance (%)",
       cex = 1, pch = 20, type = "h", lwd = 2)
}
add.scatter(myInset(), posi = "topright",
            inset = c(-0.068, -0.15), ratio = .2,
            bg = transp("white"))
mtext("A", side = 1, line = 3, cex = 1.5)# Add label 'A' 
  #Plot dapc2
scatter(dapc2, scree.da = FALSE, scree.pca = FALSE, col = group.col)
points(pred.Unknown_2_P23$ind.scores[, 1], pred.Unknown_2_P23$ind.scores[, 2], 
       pch = 18, col = transp("black"), cex = 1.5)
points(pred.Unknown_2_SDSM$ind.scores[, 1], pred.Unknown_2_SDSM$ind.scores[, 2], 
       pch = 15, col = transp("purple"), cex = 1)
points(pred.Unknown_2_UCMVZ$ind.scores[, 1], pred.Unknown_2_UCMVZ$ind.scores[, 2], 
       pch = 15, col = transp("purple"), cex = 1)
myInset <- function() {
  temp <- dapc2$pca.eig
  temp <- 100 * cumsum(temp) / sum(temp)
  plot(temp, col = rep(c("black", "lightgrey"),
                       c(dapc2$n.pca, 24)), ylim = c(0, 100),
       xlab = "PCA axis", ylab = "Variance (%)",
       cex = 1, pch = 20, type = "h", lwd = 2)
}
add.scatter(myInset(), posi = "topright",
            inset = c(-0.068, -0.15), ratio = .2,
            bg = transp("white"))
mtext("B", side = 1, line = 3, cex = 1.5)  # Add label 'B' 

dev.off()

###Ref-to-target vector plots

##Identify rows corresponding to species groups using reformatted data 
##from the first landmark repetition with species classifiers
unknown_group_indices <- which(MergedData1$species == "P23")
eastern_nc_indices <- which(MergedData1$species == "Eastern_Nc")
known_na_indices <- which(MergedData1$species %in% c("Na"))
known_nc_indices <- which(MergedData1$species %in% c("Nc"))
known_nf_indices <- which(MergedData1$species %in% c("Nf"))
known_nl_indices <- which(MergedData1$species %in% c("Nl"))
known_nm_indices <- which(MergedData1$species %in% c("Nm"))
  #Subset merged data for the unknown and known species groups
unknown_group_gpa <- MergedData1[unknown_group_indices, , drop = FALSE]
east_nc_gpa <- MergedData1[eastern_nc_indices, , drop = FALSE]
known_na_gpa <- MergedData1[known_na_indices, , drop = FALSE]
known_nc_gpa <- MergedData1[known_nc_indices, , drop = FALSE]
known_nf_gpa <- MergedData1[known_nf_indices, , drop = FALSE]
known_nl_gpa <- MergedData1[known_nl_indices, , drop = FALSE]
known_nm_gpa <- MergedData1[known_nm_indices, , drop = FALSE]
  #Extract numeric columns for mean computation
numeric_columns <- sapply(unknown_group_gpa, is.numeric)
mean_unknown <- colMeans(unknown_group_gpa[, numeric_columns])
mean_east_nc <- colMeans(east_nc_gpa[, numeric_columns])
mean_na <- colMeans(known_na_gpa[, numeric_columns])
mean_nc <- colMeans(known_nc_gpa[, numeric_columns])
mean_nf <- colMeans(known_nf_gpa[, numeric_columns])
mean_nl <- colMeans(known_nl_gpa[, numeric_columns])
mean_nm <- colMeans(known_nm_gpa[, numeric_columns])
  #Reshape mean data in a 2 column matrix
mean_unknown_matrix <- matrix(mean_unknown, ncol = 2)
mean_east_nc_matrix <- matrix(mean_east_nc, ncol = 2)
mean_na_matrix <- matrix(mean_na, ncol = 2)
mean_nc_matrix <- matrix(mean_nc, ncol = 2)
mean_nf_matrix <- matrix(mean_nf, ncol = 2)
mean_nl_matrix <- matrix(mean_nl, ncol = 2)
mean_nm_matrix <- matrix(mean_nm, ncol = 2)

##Build wireframe links over each species group's mean shape configuration
##Define custom_links matrix connecting landmarks 
custom_links <- matrix(c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,
                         10, 11, 11, 12, 12, 13, 13, 14, 14, 1), ncol = 2, 
                         byrow = TRUE)

  #Define wireframe links in the global environment
  #(make sure to select "finish" and type "n" for each configuration below)
links_data_na <- define.links(mean_na_matrix, ptsize = 2, links = custom_links)
links_data_nc <- define.links(mean_nc_matrix, ptsize = 2, links = custom_links)
links_data_nf <- define.links(mean_nf_matrix, ptsize = 2, links = custom_links)
links_data_nl <- define.links(mean_nl_matrix, ptsize = 2, links = custom_links)
links_data_nm <- define.links(mean_nm_matrix, ptsize = 2, links = custom_links)
links_data_unknown <- define.links(mean_unknown_matrix, ptsize = 2, links = custom_links)
links_data_east_nc <- define.links(mean_east_nc_matrix, ptsize = 2, links = custom_links)

##Plot the wireframe of each extant reference species over the mean shape 
##configuration of fossil unknowns (this is Fox and Blois Figure 4)
pdf("Figures/Figure4.pdf", width=12, height=8)
  #Set up the layout for the panel figure
par(mfrow = c(2, 3), mar = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))  # 2 rows and 3 columns
  #N albigula versus P23
GP1 <- gridPar(pt.size = 1.25)
plotRefToTarget(M1 = as.matrix(mean_unknown_matrix),
                M2 = as.matrix(mean_na_matrix), method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = expression(italic("N. albigula") * " vs P23"), cex.main = 2.5, line = -3)
  #Overlay the wireframe on the scatter plot using wireframe_links
for (i in 1:nrow(links_data_na)) {
  segments(mean_na_matrix[links_data_na[i, 1], 1], mean_na_matrix[links_data_na[i, 1], 2],
           mean_na_matrix[links_data_na[i, 2], 1], mean_na_matrix[links_data_na[i, 2], 2],
           col = "blue", lwd = 1.5)
}

  #N cinerea versus P23
GP1 <- gridPar(pt.size = 1.25)
plotRefToTarget(M1 = as.matrix(mean_unknown_matrix),
                M2 = as.matrix(mean_nc_matrix), method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = expression(italic("N. cinerea") * " vs P23"), cex.main = 2.5, line = -3)
for (i in 1:nrow(links_data_nc)) {
  segments(mean_nc_matrix[links_data_nc[i, 1], 1], mean_nc_matrix[links_data_nc[i, 1], 2],
           mean_nc_matrix[links_data_nc[i, 2], 1], mean_nc_matrix[links_data_nc[i, 2], 2],
           col = "plum4", lwd = 1.5)
}


  #N fuscipes versus P23
GP1 <- gridPar(pt.size = 1.25)
plotRefToTarget(M1 = as.matrix(mean_unknown_matrix),
                M2 = as.matrix(mean_nf_matrix), method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = expression(italic("N. fuscipes") * " vs P23"), cex.main = 2.5, line = -3)
for (i in 1:nrow(links_data_nf)) {
  segments(mean_nm_matrix[links_data_nf[i, 1], 1], mean_nf_matrix[links_data_nf[i, 1], 2],
           mean_nm_matrix[links_data_nf[i, 2], 1], mean_nf_matrix[links_data_nf[i, 2], 2],
           col = "limegreen", lwd = 1.5)
}

  #N lepida versus P23
GP1 <- gridPar(pt.size = 1.25)
plotRefToTarget(M1 = as.matrix(mean_unknown_matrix),
                M2 = as.matrix(mean_nl_matrix), method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = expression(italic("N. lepida") * " vs P23"), cex.main = 2.5, line = -3)
for (i in 1:nrow(links_data_nl)) {
  segments(mean_nl_matrix[links_data_nl[i, 1], 1], mean_nl_matrix[links_data_nl[i, 1], 2],
           mean_nl_matrix[links_data_nl[i, 2], 1], mean_nl_matrix[links_data_nl[i, 2], 2],
           col = "orange1", lwd = 1.5)
}

  #N macrotis versus P23
GP1 <- gridPar(pt.size = 1.25)
plotRefToTarget(M1 = as.matrix(mean_unknown_matrix),
                M2 = as.matrix(mean_nm_matrix), method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = expression(italic("N. macrotis") * " vs P23"), cex.main = 2.5, line = -3)
for (i in 1:nrow(links_data_nm)) {
  segments(mean_nm_matrix[links_data_nm[i, 1], 1], mean_nm_matrix[links_data_nm[i, 1], 2],
           mean_nm_matrix[links_data_nm[i, 2], 1], mean_nm_matrix[links_data_nm[i, 2], 2],
           col = "orangered", lwd = 1.5)
}

  #N macrotis versus N fuscipes
GP1 <- gridPar(pt.size = 1.25, pt.bg = "goldenrod")
plotRefToTarget(M1 = as.matrix(mean_nf_matrix),
                M2 = as.matrix(mean_nm_matrix), method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = expression(italic("N. macrotis") * " vs " * italic("N. fuscipes")), cex.main = 2.5, line = -3)
for (i in 1:nrow(links_data_nm)) {
  segments(mean_nm_matrix[links_data_nm[i, 1], 1], mean_nm_matrix[links_data_nm[i, 1], 2],
           mean_nm_matrix[links_data_nm[i, 2], 1], mean_nm_matrix[links_data_nm[i, 2], 2],
           col = "orangered", lwd = 1.5)
}
  #Reset the plot layout
par(mfrow = c(1, 1))
dev.off()

##Compare mean m1 shapes of California N. cinerea populations with eastern 
##populations from Idaho and South Dakota (Fox and Blois Supplementary Figure 1)

pdf("Figures/Supp_Figure1.pdf", width=6, height=4)
  #Set up the layout for the panel figure
par(mfrow = c(1, 1), mar = c(1, 0, 1, 0), oma = c(0, 0, 0, 0))
GP1 <- gridPar(pt.size = 0.75)
plotRefToTarget(M1 = as.matrix(mean_east_nc_matrix),
                M2 = as.matrix(mean_nc_matrix), method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = expression("California " * italic("N. cinerea") * " vs eastern " * 
  italic("N. cinerea")), cex.main = 1.25)
for (i in 1:nrow(links_data_nc)) {
  segments(mean_nc_matrix[links_data_nc[i, 1], 1], mean_nc_matrix[links_data_nc[i, 1], 2],
           mean_nc_matrix[links_data_nc[i, 2], 1], mean_nc_matrix[links_data_nc[i, 2], 2],
           col = "plum4", lwd = 1)
}
  #Reset the plot layout
par(mfrow = c(1, 1))
dev.off()

##Compare mean m1 shapes of datasets with all recent specimens to datasets with
##wear-vetted specimens (Fox and Blois Supplementary Figure 2)

  #First landmark repetition
consensus_ref1 <- mshape(extant.gpa_1$coords) 
  #calculate mean shape of the reference dataset (all specimens)
consensus_target1 <- mshape(Neotoma.gpa_3$coords)
  #calculate mean shape of the target comparative dataset (wear-vetted specimens)

  #Second landmark repetition
consensus_ref2 <- mshape(extant.gpa_2$coords) 
consensus_target2 <- mshape(Neotoma.gpa_4$coords)

  #(make sure select "finish" and type "n" for each configuration below)
links_target1 <- define.links(consensus_target1, ptsize = 2, links = custom_links)
links_target2 <- define.links(consensus_target2, ptsize = 2, links = custom_links)

pdf("Figures/Supp_Figure2.pdf", width=6, height=8)
  #Set up the layout for the panel figure
par(mfrow = c(2, 1), mar = c(1, 0, 1, 0), oma = c(0, 0, 0, 0))
GP1 <- gridPar(pt.size = 1)
plotRefToTarget(M1 = consensus_ref1,
                M2 = consensus_target1, method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = "All recent specimens vs wear-vetted specimens (R1)", cex.main = 1.25)
for (i in 1:nrow(links_target1)) {
  segments(consensus_target1[links_target1[i, 1], 1], consensus_target1[links_target1[i, 1], 2],
           consensus_target1[links_target1[i, 2], 1], consensus_target1[links_target1[i, 2], 2],
           col = "black", lwd = 1)
}

GP1 <- gridPar(pt.size = 1)
plotRefToTarget(M1 = consensus_ref2,
                M2 = consensus_target2, method = "vector",
                gridPars=GP1, mag = 1, label = FALSE)
title(main = "All recent specimens vs wear-vetted specimens (R2)", cex.main = 1.25)
for (i in 1:nrow(links_target2)) {
  segments(consensus_target2[links_target2[i, 1], 1], consensus_target2[links_target2[i, 1], 2],
           consensus_target2[links_target2[i, 2], 1], consensus_target2[links_target2[i, 2], 2],
           col = "black", lwd = 1)
}
  #Reset the plot layout
par(mfrow = c(1, 1))
dev.off()