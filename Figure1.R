# load packages
library(sf)
library(geodata)
library(RColorBrewer)
library(RgoogleMaps)
library(stringr)
library(tidyverse)
library(grDevices)
library(prettymapr)


# SPECIES DATA ----

# read in the species range maps 
Nalbi_range <- st_read(dsn=path.expand("Map_fig/Input/NatureServ_RangeMaps/neot_albi_pl.shp"))
Nalbi_range <- st_simplify(Nalbi_range, dTolerance=0.01, preserveTopology=FALSE)
Ncine_range <- st_read(dsn=path.expand("Map_fig/Input/NatureServ_RangeMaps/neot_cine_pl.shp"))
Ncine_range <- st_simplify(Ncine_range, dTolerance=0.01, preserveTopology=FALSE)
Nfusc_range <- st_read(dsn=path.expand("Map_fig/Input/NatureServ_RangeMaps/neot_fusc_pl.shp"))
Nfusc_range <- st_simplify(Nfusc_range, dTolerance=0.01, preserveTopology=FALSE)
Nlepi_range <- st_read(dsn=path.expand("Map_fig/Input/NatureServ_RangeMaps/neot_lepi_pl.shp"))
Nlepi_range <- st_simplify(Nlepi_range, dTolerance=0.01, preserveTopology=FALSE)
Nmacr_range <- st_read(dsn=path.expand("Map_fig/Input/NatureServ_RangeMaps/neot_macr_pl.shp"))
Nmacr_range <- st_simplify(Nmacr_range, dTolerance=0.01, preserveTopology=FALSE)

# read in the coordinates for neotoma geomorph occurrences
coords <- read.csv("Map_fig/Input/neotoma_coordinates.csv", stringsAsFactors = F)
coords$Species <- str_trim(coords$Species)
coords_ll <- filter(coords, Latitude != "NA") # note, all occurrences have latitude

# create a coordinates df for each species
Nalbi_coords <- coords_ll %>% filter(Species == "N. albigula")
Ncine_coords <- coords_ll %>% filter(Species == "N. cinerea")
Nfusc_coords <- coords_ll %>% filter(Species == "N. fuscipes")
Nlepi_coords <- coords_ll %>% filter(Species == "N. lepida")
Nmacr_coords <- coords_ll %>% filter(Species == "N. macrotis")

# set colors for the five species ----
color_brewer_palette <- c("blue", "plum3", "limegreen", "orange1", "orangered")
#brewer.pal(6, "BrBG")

Nalbi_col <- color_brewer_palette[1]  # (blue)
Ncine_col <- color_brewer_palette[2]  # (plum3)
Nfusc_col <- color_brewer_palette[3]  # (limegreen)
Nlepi_col <- color_brewer_palette[4]  # (orange1)
Nmacr_col <- color_brewer_palette[5]  # (orangered)

# Set point colors
Nalbi_point_col <- Nalbi_col 
Ncine_point_col <- Ncine_col
Nfusc_point_col <- Nfusc_col  
Nlepi_point_col <- Nlepi_col
Nmacr_point_col <- Nmacr_col

# Set range colors
Nalbi_range_col <- AddAlpha(Nalbi_col,125/255)  
Ncine_range_col <- AddAlpha(Ncine_col,125/255)  
Nfusc_range_col <- AddAlpha(Nfusc_col,150/255)  
Nlepi_range_col <- AddAlpha(Nlepi_col,150/255)  
Nmacr_range_col <- AddAlpha(Nmacr_col,150/255)

# SET LA BREA COORDS ----

labrea <- cbind(-118.3554, 34.0638)
colnames(labrea) <- c('Longitude', 'Latitude') # for some reason, this is not working!
temp <- rbind(Nmacr_coords[, c('Longitude', 'Latitude')], labrea) # this just gets labrea into right df format
labrea <- temp[nrow(temp),]

# BASE MAP DATA ----

# get country data from geodata
NorthAmerica <- gadm(country = country_codes("North America")$ISO3, level = 0, resolution = 2, path = "./Map_fig/Input/geodata_countries")

#check that it plots correctly 
plot(NorthAmerica, xlim = c(-180, -50))

# CREATE THE MAP ELEMENTS----

# Create each map panel separately, then combine them in Illustrator

# central figure
png(file="Map_fig/Output/woodrat-map-all_species_orig.png", height=4, width=4, units="in", res=300)
xlim=c(-142,-97)
ylim=c(23,65)
plot(NorthAmerica, xlim=xlim, ylim=ylim, bty="o")
addscalebar(plotepsg = 4326, unitcategory = "metric", padin=c(0.5,0.15))
addnortharrow(pos="bottomleft", padin = c(0.65, 0.35), scale=0.5)
plot(Nalbi_range, add=T, col=Nalbi_range_col, lty=4, lwd=1.5)
plot(Ncine_range, add=T, col=Ncine_range_col, lty=1)
plot(Nfusc_range, add=T, col=Nfusc_range_col, lty=2, lwd=1.5)
plot(Nlepi_range, add=T, col=Nlepi_range_col, lty=3, lwd=1.5)
plot(Nmacr_range, add=T, col=Nmacr_range_col, lty=5, lwd=1.5)
points(labrea[, c('Longitude', 'Latitude')], cex=1, pch=23, bg="gray", col="black")
dev.off()

#Nalbi
png(file="Map_fig/Output/woodrat-map-n.albigula_orig.png", height=3, width=2, units="in", res=300)
xlim=c(-120,-100) 
ylim=c(24,40)
plot(NorthAmerica, xlim=xlim, ylim=ylim, bty="o")
addscalebar(plotepsg = 4326, unitcategory = "metric", pos="bottomleft", 
            style = "bar", lwd=0.25, padin=c(0.025,0.5), 
            widthhint = 0.1, htin=0.075, labelpadin=0.02, label.cex = 0.4)
plot(Nalbi_range, add=T, col=Nalbi_range_col)
points(labrea[, c('Longitude', 'Latitude')], cex=1.5, pch=23, bg="gray", col="black")
points(Nalbi_coords[, c('Longitude', 'Latitude')], pch=21, cex=0.5, bg=Nalbi_point_col , col="black", lwd=0.5)
mtext(text= ~ italic("Neotoma albigula"), side=3, line=0.5, cex=1)
dev.off()

#Ncine
png(file="Map_fig/Output/woodrat-map-n.cinerea_orig.png", height=3, width=2, units="in", res=300)
xlim=c(-140,-100)
ylim=c(31,66)
plot(NorthAmerica, xlim=xlim, ylim=ylim, bty="o")
addscalebar(plotepsg = 4326, unitcategory = "metric", pos="bottomleft", 
            style = "bar", lwd=0.25, padin=c(0.025,0.5), 
            widthhint = 0.1, htin=0.075, labelpadin=0.02, label.cex = 0.4)
plot(Ncine_range, add=T, col=Ncine_range_col)
points(labrea[, c('Longitude', 'Latitude')], cex=1.5, pch=23, bg="gray", col="black")
points(Ncine_coords[, c('Longitude', 'Latitude')], pch=21, cex=0.5, bg = Ncine_point_col, col="black", lwd=0.5)
mtext(text= ~ italic("Neotoma cinerea"), side=3, line=0.5, cex=1)
dev.off()

#Nfusc
png(file="Map_fig/Output/woodrat-map-n.fuscipes_orig.png", height=3, width=2, units="in", res=300)
xlim=c(-125,-110)
ylim=c(33,46) 
plot(NorthAmerica, xlim=xlim, ylim=ylim, bty="o")
addscalebar(plotepsg = 4326, unitcategory = "metric", pos="bottomright", 
            style = "bar", lwd=0.25, padin=c(0.025,0.5), 
            widthhint = 0.1, htin=0.075, labelpadin=0.02, label.cex = 0.4)
plot(Nfusc_range, add=T, col=Nfusc_range_col)
points(labrea[, c('Longitude', 'Latitude')], cex=1.5, pch=23, bg="gray", col="black")
points(Nfusc_coords[, c('Longitude', 'Latitude')], pch=21, cex=0.5, bg = Nfusc_point_col, col="black", lwd=0.5)
mtext(text= ~ italic("Neotoma fuscipes"), side=3, line=0.5, cex=1)
dev.off()

#Nlepi
png(file="Map_fig/Output/woodrat-map-n.lepida_orig.png", height=3, width=2, units="in", res=300)
xlim=c(-130,-105)
ylim=c(23,45)
plot(NorthAmerica, xlim=xlim, ylim=ylim, bty="o")
addscalebar(plotepsg = 4326, unitcategory = "metric", pos="bottomleft", 
            style = "bar", lwd=0.25, padin=c(0.025,0.5), 
            widthhint = 0.1, htin=0.075, labelpadin=0.02, label.cex = 0.4)
plot(Nlepi_range, add=T, col=Nlepi_range_col)
points(labrea[, c('Longitude', 'Latitude')], cex=1.5, pch=23, bg="gray", col="black")
points(Nlepi_coords[, c('Longitude', 'Latitude')], pch=21, cex=0.5, bg = Nlepi_point_col, col="black", lwd=0.5)
mtext(text= ~ italic("Neotoma lepida"), side=3, line=0.5, cex=1)
dev.off()

#Nmacr
png(file="Map_fig/Output/woodrat-map-n.macrotis_orig.png", height=3, width=2, units="in", res=300)
xlim=c(-125,-112) #16.5
ylim=c(29,40) #11
plot(NorthAmerica, xlim=xlim, ylim=ylim, bty="o")
addscalebar(plotepsg = 4326, unitcategory = "metric", pos="bottomleft", 
            style = "bar", lwd=0.25, padin=c(0.025,0.5), 
            widthhint = 0.1, htin=0.075, labelpadin=0.02, label.cex = 0.4)
plot(Nmacr_range, add=T, col=Nmacr_range_col)
points(labrea[, c('Longitude', 'Latitude')], cex=1.5, pch=23, bg="gray", col="black")
points(Nmacr_coords[, c('Longitude', 'Latitude')], pch=21, cex=0.5, bg = Nmacr_point_col, col="black", lwd=0.5)
mtext(text= ~ italic("Neotoma macrotis"), side=3, line=0.5, cex=1)
dev.off()
