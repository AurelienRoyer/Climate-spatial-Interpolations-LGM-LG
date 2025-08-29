#############################################
# Extract Rodent Species by Locality
# Intersect IUCN range polygons with
# 50 km buffers around localities to generate
# species presence/absence data.
#############################################

# Load necessary libraries
library(sf)
library(terra)
library(lwgeom)

# ------------------------------------------
# 1. Load IUCN Rodent Distribution Shapefile
# A dataset containing European rodents was downloaded from the IUCN Red List of Threatened Species. Version 3, May 2017. https://www.iucnredlist.org. 
# Downloaded on the 23 June 2021
# ------------------------------------------
IUCNfile <- "redlist_species_data/data_0.shp"
rodents <- st_read(IUCNfile)
rodents$BINOMIAL <- as.character(rodents$BINOMIAL)

# Split into list by species
allsp <- unique(rodents$BINOMIAL)
spList <- lapply(allsp, function(sp) rodents[rodents$BINOMIAL == sp,])
names(spList) <- allsp

# Repair invalid geometries
spList <- lapply(spList, function(x) st_make_valid(x))

# ------------------------------------------
# 2. Reproject to Equal-Area (Behrmann)
# ------------------------------------------
Behrmannproj <- "+proj=cea +lon_0=0 +lat_ts=30 +datum=WGS84 +units=m +no_defs"
spListEA <- lapply(spList, function(x) st_transform(x, crs = Behrmannproj))

# ------------------------------------------
# 3. Load Localities and Create 50 km Buffers
# ------------------------------------------
localities <- read.csv("Table C1.csv", sep = ";", dec = ".")
loc_points <- st_as_sf(localities, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
loc_points_proj <- st_transform(loc_points, crs = Behrmannproj)

# Create 50 km radius buffers
buffers <- st_buffer(loc_points_proj, dist = 50000)

# ------------------------------------------
# 4. Spatial Intersection: Species × Locality
# ------------------------------------------
species_matrix <- matrix(0, nrow = length(allsp), ncol = nrow(localities),
                         dimnames = list(allsp, localities$LOCALITY))

for (i in seq_along(spListEA)) {
  dist_poly <- st_geometry(spListEA[[i]])
  for (j in 1:nrow(localities)) {
    buf_poly <- st_geometry(buffers[j, ])
    if (length(st_intersects(buf_poly, dist_poly)[[1]]) > 0) {
      species_matrix[i, j] <- 1
    }
  }
}

# ------------------------------------------
# 5. Save Output
# ------------------------------------------
write.csv2(species_matrix, "species_list_circles50km.csv")
