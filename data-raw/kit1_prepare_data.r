# library(ggplot2)
library(spatstat.geom)
library(sf)
library(geodata)
library(cluster)
library(raster)
library(doParallel)
# library(mgcv)

####################
## input parameters
####################
# distance lags (in km from colony)
lags_sp<- seq(0, 210, by= 10)
# number of weeks to use as predictors
nweeks<- 30
# time lags (in weeks until Jun 30th)
lags_ti<- nweeks:1
# time index (in weeks from Jan 1st)
ID_ti<- 1:nweeks

####################
## Kittiwake data
####################
kitty<- read.csv("../Data_for_BioSS/SMP_kitti_BS.csv", stringsAsFactors= T)

## Supplement for Isle of May
kitty_suppl<- read.csv("../Data_for_BioSS/Bioss_kitti_colonies_subset.csv", stringsAsFactors= T)
kitty_suppl<- kitty_suppl[kitty_suppl$Master.site == "Isle of May SPA", ]

kitty_cols_of_interest<- c("Year", "Site", "Species", "Country", "Master.site", "Subsite", "Unit", "Count", "Fledged.count", "StartLong", "StartLat")
kitty<- rbind(kitty_suppl[, kitty_cols_of_interest], kitty[, kitty_cols_of_interest])

## only keep years for which we have SST data
SST_years<- 1993:2019
kitty<- kitty[kitty$Year %in% SST_years, ]

kitty$SiteYear<- paste(kitty$Site, kitty$Year, sep= ".")

site_by_year<- (table(kitty$Site, kitty$Year) > 0) * 1
hm<- heatmap(site_by_year, Colv= NA, hclustfun = function(d) agnes(d, method = "ward"))
# heatmap(site_by_year, Colv= NA, Rowv= hm$rowInd, col= c(0, 1))
par(las= 1, mar= c(3, 12, 2, 2))
image(t(site_by_year[hm$rowInd, hm$colInd]), col= c(0, 1), axes= F)
axis(1, at= seq(0, 1, l= ncol(site_by_year)), labels= colnames(site_by_year)[hm$colInd])
axis(2, at= seq(0, 1, l= nrow(site_by_year)), labels= rownames(site_by_year)[hm$rowInd], cex.axis= 0.7)

plot(ecdf(rowSums(site_by_year)))
nrow(site_by_year)
# all 136 sites have gaps
# ~ 25% with long-term data (20+)
# ~ 60% with < 10 yr data
# -> loop through available sites * years for extracting environmental data rather than doing all possible combinations

# get unique colony sites
sites<- unique(kitty[, c("Site", "StartLong", "StartLat")])
rownames(sites)<- sites$Site
# create a SpatialPoints object
sites_sp<- SpatialPoints(sites[, c("StartLong", "StartLat")],
    proj4string= CRS("+proj=longlat +datum=WGS84 +no_defs"))

sites_sp_metric<- spTransform(sites_sp, "+init=epsg:3035")

# get unique obs per site*year and store as "kit"
kit<- unique(kitty[, c("SiteYear", "Site", "StartLong", "StartLat", "Year")])
kit<- kit[order(kit$Year, kit$Site), ]

# take sum of nest and fledgling counts per site
kit$AON<- as.vector(tapply(kitty$Count, kitty$SiteYear, sum, na.rm= T)[kit$SiteYear])
kit$Fledg<- as.vector(tapply(kitty$Fledged.count, kitty$SiteYear, sum, na.rm= T)[kit$SiteYear])

plot(table(kit$Year)) # 2013 has the largest number of sites monitored

# check series look sensible
coplot(AON ~ Year | Site, data= kit, show.given = FALSE)
coplot(Fledg ~ Year | Site, data= kit, show.given = FALSE)
coplot(Fledg/AON ~ Year | Site, data= kit, show.given = FALSE)

colony_size<- read.csv("../Data_for_BioSS/Colony_mean_count.csv")
names(colony_size)[2]<- "MeanPop"

# Density of breeders in neighborhood
col_size_full<- read.csv("../Data_for_BioSS/comparative-seabirds-count-dataset-revised-20231213.csv")
col_size_full<- col_size_full[col_size_full$Species == "Black-legged Kittiwake", ]

## add Isle of May data
col_size_full<- rbind(col_size_full, NA)
col_size_full[nrow(col_size_full), c("Master.site", "Site")]<- c("Isle of May SPA", "Isle of May")
col_size_full[nrow(col_size_full), c("Start.Lat", "Start.Long", "S.2000.adjusted.count", "Seabirds.Count.adjusted.count")]<- c(56.1854, -2.5578, 4486, 4486)
# count data above comes from long-term average given by Charlotte

## Distances to other colonies and pop size within distance bands

### function to convert Lon, Lat coordinates from degrees to meters (assumes British National Grid projection by default)
deg_to_meters<- function(coords, origin.CRS= 4326, target.CRS= 27700){
	coords_sf<- st_transform(st_as_sf(coords, coords = 1:2, crs = origin.CRS), crs = target.CRS)
	st_coordinates(coords_sf)
}

## Unique colonies from kit
unikit<- unique(kit[, c("Site", "StartLong", "StartLat")])
## Convert to meters
unikit<- cbind(unikit, deg_to_meters(unikit[, 2:3]))
col_size_full<- cbind(col_size_full, deg_to_meters(col_size_full[, c("Start.Long", "Start.Lat")]))

plot(col_size_full[, c("X", "Y")])
points(unikit[, 4:5], pch= 16, col= rgb(1, 0, 0, 0.5), cex= 3)

plot(unikit[, 4:5], pch= 16, col= 2, cex= 3)
points(col_size_full[, c("X", "Y")])

## Compute pairwise distances
pairwise_distances<- crossdist(unikit$X, unikit$Y, col_size_full$X, col_size_full$Y, squared= F)

## Convert into distance bands (in km)
dist_band<- cut(pairwise_distances, breaks= 1000 * c(-1, lags_sp[-1], 1e7), labels= c(lags_sp[-1], 1e7))
plot(dist_band, ylim= c(0, 1000))
dim(dist_band)<- dim(pairwise_distances)

## Count nests within each (distance band) of kit observations, including own colony

### Using site-mean of Seabird Count census (1998-2002) and Seabirds Count (2015-2021), because many sites have not been visited on both surveys
col_size_mean<- apply(col_size_full[, c("S.2000.adjusted.count", "Seabirds.Count.adjusted.count")], 1, mean, na.rm= T)
col_size_mean[is.na(col_size_mean)]<- 0 # assume that sites never visited had no nests - they should mostly be small colonies anyway

# matplot(log(cbind(col_size_mean, col_size_full[, c("S.2000.adjusted.count", "Seabirds.Count.adjusted.count")])+1))

nest_count<- list()
for(i in 1:nrow(dist_band)){
	nest_count[[i]]<- tapply(col_size_mean, dist_band[i, ], sum, na.rm= T)
}
nest_count<- do.call(rbind, nest_count)[, 1:(length(lags_sp)-1)] # rbind into matrix and remove last distance band
nest_count[is.na(nest_count)]<- 0
rownames(nest_count)<- unikit$Site
head(nest_count)
image(t(nest_count))

unikit$Site[nest_count[, 1] == 0] # some issues, but not over North Sea Coast

## nest counts matrix for models
nest_count_rings<- nest_count[as.character(kit$Site), ]

tail(cbind(kit, nest_count_rings))
# Site subsets (excluding Isle of May, due to lack of open data license)
# Subset of sites including Shetland and Orkney
kit_sub2<- c("Coquet Island RSPB", "Dunbar Coast", "Fair Isle", "Farne Islands", "Foula", "Marwick Head", "Noness", "North Hill RSPB, Papa Westray", "North Sutor Of Cromarty/Castlecraig", "Row Head", "Saltburn Cliffs (Huntcliff)", "Sands of Forvie", "St Abb's Head NNR", "Sumburgh Head")
# Subset of sites excluding Shetland and Orkney
kit_sub1<- c("Coquet Island RSPB", "Dunbar Coast", "Farne Islands", "North Sutor Of Cromarty/Castlecraig", "Saltburn Cliffs (Huntcliff)", "Sands of Forvie", "St Abb's Head NNR")
# Broader subset for North Sea, including Isle of May (does have open license after all, just that data are hosted on different website to the rest)
pol_sub3<- kit$StartLong > -3.0 & kit$StartLat > 54.0 & kit$StartLat < 58.0
with(kit, plot(StartLong, StartLat, col= pol_sub3 + 1))
kit_sub3<- unique(as.character(kit$Site[pol_sub3]))

# subset data to colonies from the North-East
kit2<- kit[kit$Site %in% kit_sub2, ]

kit2$Site<- factor(kit2$Site)

# Add colony size column to the data ("MeanPop")
kit2<- merge(kit2, colony_size)

coplot(AON ~ Year | Site, data= kit2, show.given = FALSE)
coplot(Fledg ~ Year | Site, data= kit2, show.given = FALSE)
coplot(Fledg/AON ~ Year | Site, data= kit2, show.given = FALSE)

## plot kit subset 3
uk_coords<- gadm(country= c("GBR", "IRL"), level= 0, path=tempdir(), resolution= 2)
### Create the contour plot
png("C:/Users/TC44053/Documents/_Research/exposés/2024_ISEC_2024_Swansea/colonies_map.png", width= 480*2, height= 480*3, pointsize= 34)
plot(uk_coords, type = "l", xlab = "Longitude", ylab = "Latitude", col= grey(0.9))
points(col_size_full[, c("Start.Long", "Start.Lat")], pch= 16, col= 4)
points(unikit[unikit$Site %in% kit_sub3, c("StartLong", "StartLat")], pch= 16, col= rgb(0, 1, 0, 0.4), cex= 2)
points(unikit[unikit$Site %in% kit_sub3, c("StartLong", "StartLat")], pch= 1, col= 3, cex= 2)
dev.off()

###############
## SST data
###############
# Geotiffs stored outside of gitlab, as they're a tad large
# import one as template
SST<- stack("../Data_for_BioSS/SST_1993.tif")

# create a metric reprojection template for all future rasters
reproj_template<- projectRaster(SST[[1]], crs= "+init=epsg:3035")
# bug: throws an error on first call -> call again
reproj_template<- projectRaster(SST[[1]], crs= "+init=epsg:3035")

###### Some plots of SST
require(sp)
require(sf)
require(raster)
colony_c<- st_as_sf(data.frame(x= -2.555586, y= 56.185728), coords= c("x", "y"))
st_crs(colony_c)<- 4326

# Define raster cropping extent
e<- extent(-4, -1, 55, 58)
# Crop the raster
SST_cropped<- crop(SST, e)

# Get UK data
uk_data<- raster::getData('GADM', country= 'GBR', level= 0)

for(t in c(1, 6)){
    png(paste0("SSTmap_rings_",t,".png"), 1200, 1200)
    image(SST_cropped[[t]], useRaster= TRUE, col= heat.colors(256), asp= 1)
    for(i in seq(1, 20) / 10){
        symbols(st_coordinates(colony_c), circles= i, inches= F, 
        col= rgb(0.5, 0.5, 0.5, 0.1), cex= 2, add= T, lty= 3)
    plot(uk_data, col= "white", add= T)
    }
    points(st_coordinates(colony_c), pch= 16, col= "blue", cex= 4)
    dev.off()
}

for(t in c(1, 6)){
    png(paste0("SSTbasemap_",t,".png"), 1200, 1200)
    image(SST_cropped[[t]], useRaster= TRUE, col= heat.colors(256), asp= 1)
    plot(uk_data, col= "white", add= T)
    points(st_coordinates(colony_c), pch= 16, col= "blue", cex= 4)
    dev.off()
}

###### Processing SST

# All required computations are within a year, 
# so no need to create a full time x space brick.
# Years can be treated in parallel

# Computing predictors at lags - METHOD #1
# Suitable when predictor needed at sparse set of locations

cl <- makePSOCKcluster(4)
registerDoParallel(cl)

system.time(
SST_extract<- foreach(yr= SST_years, 
                .packages= "raster",
                .errorhandling= c("pass")) %dopar% {

    # identify colonies monitored in the year of interest
    sites_in_year_index<- which(sites$Site %in% kit$Site[kit$Year == yr])
    sites_in_year_names<- sites$Site[sites_in_year_index]
    SiteYear_names<- paste(sites_in_year_names, yr, sep= ".")
	
    sites_in_year_sp_metric<- sites_sp_metric@coords[sites_in_year_index,]

    # construct arrays for storing the outputs
    SST_sum_buffer<- array(NA, dim= c(nweeks,
                                    length(lags_sp[-1]),
                                    sum(kit$Year == yr)), 
                        dimnames= list(week= 1:nweeks,
                                dist= lags_sp[-1],
                                SiteYear= SiteYear_names)
    )

    SST_noNA_buffer<- SST_sum_buffer

    # import SST for year 'yr'
    SSTyr<- stack(paste0("../Data_for_BioSS/SST_", yr, ".tif"))
    SSTyr<- projectRaster(SSTyr, reproj_template)

    # compute SST sum within each distance buffer, for the first 'nweeks' layers
    for(dist in lags_sp[-1]){
        SST_sum_buffer[, as.character(dist), ]<- t(extract(SSTyr[[1:nweeks]], y= sites_in_year_sp_metric, buffer= dist*1000, na.rm= TRUE, layer= 1, n= nweeks, fun= sum))
    }

    # count non-missing SST values within each distance buffer, for the first 'nweeks' layers
    for(dist in lags_sp[-1]){
        SST_noNA_buffer[, as.character(dist), ]<- t(extract(SSTyr[[1:nweeks]], y= sites_in_year_sp_metric, buffer= dist*1000, na.rm= FALSE, layer= 1, n= nweeks, fun= function(z){sum(!is.na(z))}))
    }

    list(year= yr, 
        SiteYear_names= SiteYear_names, 
        SST_sum_buffer= SST_sum_buffer, 
        SST_noNA_buffer= SST_noNA_buffer)
}
)
stopCluster(cl)

# combine yearly outputs into single arrays
SST_SiteYear<- unlist(lapply(SST_extract, FUN= function(x){x$SiteYear_names}))

SST_sum_buffer<- lapply(SST_extract, FUN= function(x){x$SST_sum_buffer}) |>
            unlist() |>
                array(data= _, dim= c(nweeks, length(lags_sp[-1]), nrow(kit)), 
                        dimnames= list(week= 1:nweeks,
                                dist= lags_sp[-1],
                                SiteYear= kit$SiteYear)
                )

SST_noNA_buffer<- lapply(SST_extract, FUN= function(x){x$SST_noNA_buffer}) |>
            unlist() |>
                array(data= _, dim= c(nweeks, length(lags_sp[-1]), nrow(kit)), 
                        dimnames= list(week= 1:nweeks,
                                dist= lags_sp[-1],
                                SiteYear= kit$SiteYear))

# compute rings from difference between successive buffers
SST_sum_ring<- apply(SST_sum_buffer, c(1, 3), FUN= function(x){diff(c(0, x))})
SST_noNA_ring<- apply(SST_noNA_buffer, c(1, 3), FUN= function(x){diff(c(0, x))})

# compute rings mean from SSTsum/Npixels
SST_mean_ring<- SST_sum_ring / SST_noNA_ring


save(SST_sum_buffer, SST_noNA_buffer, SST_mean_ring, kit, file= "kit_SST_v2.RData")
# load("kit_SST_v2.RData")

##############################
## Sandeel density predictions
##############################
# import data
sandeel<- raster("../Data_for_BioSS/nmpwfs_species_distribution_lesser_sandeels_north_sea_predicted_density.tif")
sandeel<- projectRaster(sandeel, reproj_template)
plot(sandeel)

# transform values, to limit the influence of the right skew?
par(mfrow= c(1, 2))
hist(values(sandeel), nclass= 100)
hist(values(sandeel)^(1/3), nclass= 100)

sandeel<- sandeel^(1/3)
par(mfrow= c(1, 1))
plot(sandeel, xlim= c(3400000, 3900000), ylim= c(3250000, 4300000))

# Computing predictors at lags - METHOD #1
# Suitable when predictor needed at sparse set of locations

compute_rings_fixed_time<- function(X, lags= lags_sp, sites.sp){
    # construct arrays for storing the outputs
    site_names<- row.names(sites.sp)
    X_sum_buffer<- array(NA, dim= c(length(lags[-1]),
                                    length(site_names)), 
                        dimnames= list(dist= lags[-1],
                                Site= site_names)
    )

    X_noNA_buffer<- X_sum_buffer

    # compute X sum within each distance buffer
    for(dist in lags[-1]){
        X_sum_buffer[as.character(dist), ]<- t(extract(X, y= sites.sp, buffer= dist*1000, na.rm= TRUE, fun= sum))
    }
    # count non-missing X values within each distance buffer
    for(dist in lags[-1]){
        X_noNA_buffer[as.character(dist), ]<- t(extract(X, y= sites.sp, buffer= dist*1000, na.rm= FALSE, fun= function(z){sum(!is.na(z))}))
    }

    # compute rings from difference between successive buffers
    X_sum_ring<- apply(X_sum_buffer, c(2), FUN= function(x){diff(c(0, x))})
    X_noNA_ring<- apply(X_noNA_buffer, c(2), FUN= function(x){diff(c(0, x))})

    # compute rings mean from Xsum/Npixels
    X_mean_ring<- X_sum_ring / X_noNA_ring

    t(X_mean_ring)
}

sandeel_mean_ring<- compute_rings_fixed_time(sandeel, lags= lags_sp, sites.sp= sites_sp)

save(SST_sum_buffer, SST_noNA_buffer, SST_mean_ring, sandeel_mean_ring, nest_count_rings, kit, kit_sub1, kit_sub2, kit_sub3, kit2, sites_sp, file= "kit_SST_sandeel_v2.RData")

