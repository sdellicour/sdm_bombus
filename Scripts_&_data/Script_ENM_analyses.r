# system("R CMD javareconf"); # rJava_0.9-10.tar.gz didn't work before
# install.packages("rJava_0.9-9.tar.gz", repos=NULL, type="source")
# install.packages("dismo_1.1-4.tar.gz", repos=NULL, type="source")

dyn.load("/Library/Java/JavaVirtualMachines/jdk1.8.0_151.jdk/Contents/Home/jre/lib/server/libjvm.dylib")
library(raster); library(fields); library(RColorBrewer); library(dismo); library(gbm); library(geosphere); library(rgdal); library(pgirmess)
options(java.parameters="-Xmx15000m"); library(rJava)

# 1. Preparing the environmental rasters
	# 1.1. Cropping the WorldClim 2.0 raster files (source: http://worldclim.org/bioclim)
	# 1.2. Preparing the land cover variables (source: https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd12q1)
	# 1.3. Preparing the Natural Earth background raster (visualisation only)
# 2. Variables selection based on pairwise correlations
# 3. Preparing the human population density and null rasters
# 4. Defining the background area		
	# 4.1. With a mask based on D_max
		# 4.1.1. Estimating d_max_list[[h]] with a jackknife
		# 4.1.2. Generating a mask based on d_max_list[[h]]
	# 4.2. Only considering raster cells associated with occurrence data
# 5. Plotting the unique observation points (with or without the mask)
# 6. Maxent analyses with an heterogeneous background
	# 6.1. Sampling pseudo-absences in the background (considering or not human pop. density)
	# 6.2. Removing observations from adjacent cells (optional)
	# 6.3. Keeping only one presence or pseudo-absence point per raster cell (priority = presence points)
	# 6.4. Maxent analysis in itself
	# 6.5. Comparing the results of different Maxent analyses
	# 6.6. Investigating the different response curves
	# 6.7. Extracting and plotting the permutation importances
	# 6.8. Plotting the best SDM result for each species
	# 6.9. Analysing the phylogenetic signal associated with importance
	# 6.10. Analysing the association between species range and bioclimatic importances

datasetsDirectory = "Bombus_observations_2"
resolution = "16"; # resolution = "04"
e_1 = extent(-19, 88, 19, 82); e_2 = extent(-15, 84, 23, 78)
subsampling_adjacent_cells = FALSE; background_cell = TRUE
background_mask = FALSE; human_pop_density_bias = FALSE
cols = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(121)[11:101])
plotWidth = 7.0; plotHeight = 6.2
datasets = list.files(datasetsDirectory); # datasets = datasets[40]
datasets = gsub(".csv","",datasets[which(grepl(".csv",datasets))]) 
species_information = read.csv("Species_information.csv")
subgenera = as.character(unique(species_information[,"subgenera"]))
subgenera_cols = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
"#ddce62","#fdbf6f","#ff7f00","#cab2d6","#b780b6","#6a3d9a","#d3d3d3","#848484")

# 1. Loading the environmental rasters

	# 1.1. Loading the WorldClim 2.0 raster files (source: http://worldclim.org/bioclim)

bioclimatic_variables = as.character(read.csv("Bioclimatic_variables.csv",header=T)[,"short_name"])
files = list.files("WorldClim_2.0_rasters"); rasters_1 = list(); c = 0
for (i in 1:length(files))
	{
		if (grepl(paste0("_",resolution,".asc"),files[i]) == TRUE)
			{
				c = c+1; rast = raster(paste0("WorldClim_2.0_rasters/",files[i])); id = unlist(strsplit(files[i],"_"))[length(unlist(strsplit(files[i],"_")))-1]
				names(rast) = bioclimatic_variables[as.numeric(gsub(".asc","",id))]; crs(rast) = CRS("+proj=longlat +datum=WGS84"); rasters_1[[c]] = rast
			}
	}

	# 1.2. Loading the land cover variables (source: https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd12q1)

rasters_2 = list(); land_cover_variables = c("barren_vegetation","closed_shrubland","croplands","deciduous_broadleaf_forest",
								"deciduous_needleleaf_forest", "evergreen_broadleaf_forest","evergreen_needleleaf_forest","grasslands",
								"mixed_forests","open_shrublands","savannas","snow_ice","urban_areas","wetlands","woody_savannas")
for (i in 1:length(land_cover_variables))
	{
		rast = raster(paste0("Land_cover_rasters/CL_",land_cover_variables[i],"_",resolution,".asc"))
		names(rast) = land_cover_variables[i]; rasters_2[[i]] = rast
	}

	# 1.3. Preparing the Natural Earth background raster (visualisation only)

background = crop(raster("Natural_Earth_GIS_files/Gray_background.tif"), e_1); background[background[]==106] = NA; r = background
cols_background = colorRampPalette(c("grey","white"),bias=1)(max(r[],na.rm=T)-min(r[],na.rm=T))[1:(max(r[],na.rm=T)-min(r[],na.rm=T))]

# 2. Variables selection based on pairwise correlations

observations_list = list(); N_tot_1 = 0; N_tot_2 = 0
species_information = read.csv("Species_information.csv", header=T)
landPolygons_light = rgeos::gSimplify(shapefile("World_land_polygons/Europe_land_polygons.shp"), 0.1)
for (h in 1:length(datasets))
	{
		observations = read.csv(paste0(datasetsDirectory,"/",datasets[h],".csv"), header=T)
		N_tot_1 = N_tot_1 + dim(observations)[1]
		observations = unique(observations[,c("longitude","latitude")])
		N_tot_2 = N_tot_2 + dim(observations)[1]
		observations = observations[!is.na(extract(rasters_1[[1]], observations)),]
		observations = observations[!is.na(extract(rasters_2[[1]], observations)),]
		observations_list[[h]] = observations; cat(datasets[h]); cat("\n")
		if (h == 1)
			{
				xmin = min(observations[,1]); xmax = max(observations[,1])
				ymin = min(observations[,2]); ymax = max(observations[,2])
			}	else		{
				if (xmin > min(observations[,1])) xmin = min(observations[,1])
				if (xmax < max(observations[,1])) xmax = max(observations[,1])
				if (ymin > min(observations[,2])) ymin = min(observations[,2])
				if (ymax < max(observations[,2])) ymax = max(observations[,2])
			}
		hull = chull(observations); hull = c(hull,hull[1])
		p = Polygon(observations[hull,]); ps = Polygons(list(p),1)
		sps = SpatialPolygons(list(ps)); crs(sps) = crs(landPolygons_light)
		range = intersect(sps, landPolygons_light)
		index = which(species_information[,"species"]==datasets[h])
		species_information[index,"area_km2"] = round(sum(areaPolygon(range)/1000000))
	}
write.csv(species_information, "Species_information_NEW.csv", quote=F, row.names=F)

for (h in 1:length(datasets)) # to print the pairs of predictors with a Pearson's r >= 0.7, and the non-variable predictors
	{
		observations = observations_list[[h]]
		for (i in 1:length(rasters_1))
			{
				rast = crop(rasters_1[[i]], e_2); res(rast) = res(rasters_1[[1]])
				crs(rast) = CRS("+proj=longlat +datum=WGS84"); names(rast) = names(rasters_1[[i]]); rasters_1[[i]] = rast
			}
		extractions = matrix(nrow=dim(observations)[1], ncol=length(rasters_1)); colnames(extractions) = bioclimatic_variables
		correlations = matrix(nrow=length(rasters_1), ncol=length(rasters_1)); row.names(correlations) = bioclimatic_variables 
		for (i in 1:length(rasters_1)) extractions[,names(rasters_1[[i]])] = extract(rasters_1[[i]], observations)
		for (i in 1:length(rasters_1))
			{
				for (j in 1:i)
					{
						correlations[i,j] = cor(extractions[,i],extractions[,j]); correlations[j,i] = correlations[i,j]
					}
			}
		for (i in 2:length(rasters_1))
			{
				for (j in 1:(i-1))
					{
						if ((!is.na(correlations[i,j]))&(abs(correlations[i,j]) >= 0.7))
							{
								cat(paste0("	    r >= 0.7 - BIO ",i," - BIO ",j)); cat("\n")
							}
					}
			}
		for (i in 1:length(rasters_2))
			{
				rast = crop(rasters_2[[i]], e_2); res(rast) = res(rasters_1[[1]]); extent(rast) = extent(rasters_1[[1]])
				crs(rast) = CRS("+proj=longlat +datum=WGS84"); names(rast) = names(rasters_2[[i]]); rasters_2[[i]] = rast
			}
		extractions = matrix(nrow=dim(observations)[1], ncol=length(rasters_2)); colnames(extractions) = land_cover_variables
		correlations = matrix(nrow=length(rasters_2), ncol=length(rasters_2)); row.names(correlations) = land_cover_variables 
		for (i in 1:length(rasters_2)) extractions[,names(rasters_2[[i]])] = extract(rasters_2[[i]], observations)
		for (i in 1:length(rasters_2))
			{
				for (j in 1:i)
					{
						correlations[i,j] = cor(extractions[,i],extractions[,j]); correlations[j,i] = correlations[i,j]
					}
			}
		for (i in 2:length(rasters_2))
			{
				for (j in 1:(i-1))
					{
						if ((!is.na(correlations[i,j]))&(abs(correlations[i,j]) >= 0.7))
							{
								cat(paste0("    r >= 0.7 - ",land_cover_variables[i]," - ",land_cover_variables[j])); cat("\n")
							}
					}
			}
		for (i in 1:length(rasters_2))
			{
				if (length(unique(extractions[,i])) == 1)
					{
						cat(paste0("    no variation in ",land_cover_variables[i],", variable number ",i)); cat("\n")
					}
			}
	}
rasters1_to_select = c(1,4,12,15); rasters2_to_select = 1:length(rasters_2); rasters_stacks = list()
for (h in 1:length(datasets))
	{
		if (h == 1)
			{
				rasters = list()
				for (i in rasters1_to_select)
					{
						rast = crop(rasters_1[[i]], e_2); res(rast) = res(rasters_1[[1]])
						crs(rast) = CRS("+proj=longlat +datum=WGS84"); names(rast) = names(rasters_1[[i]])
						rasters[[length(rasters)+1]] = rast
					}
				for (i in rasters2_to_select)
					{
						rast = crop(rasters_2[[i]], e_2); res(rast) = res(rasters_1[[1]]); extent(rast) = extent(rasters_1[[1]])
						crs(rast) = CRS("+proj=longlat +datum=WGS84"); names(rast) = names(rasters_2[[i]])
						rasters[[length(rasters)+1]] = rast
					}
				for (i in 1:length(rasters1_to_select))
					{
						rast1 = rasters[[i]]; rast2 = rasters_2[[1]]; rast1[is.na(rast2[])] = NA; rasters[[i]] = rast1
					}
				variables = rep(NA, length(rasters))
				for (i in 1:length(rasters)) variables[i] = names(rasters[[i]])
				rasters_stacks = stack(rasters)
			}
		observations = observations_list[[h]]
		extractions = matrix(nrow=dim(observations)[1], ncol=length(rasters)); colnames(extractions) = variables
		correlations = matrix(nrow=length(rasters), ncol=length(rasters)); row.names(correlations) = variables 
		for (i in 1:length(rasters)) extractions[,names(rasters[[i]])] = extract(rasters[[i]], observations)
		for (i in 1:length(rasters))
			{
				for (j in 1:i)
					{
						correlations[i,j] = cor(extractions[,i],extractions[,j]); correlations[j,i] = correlations[i,j]
					}
			}
	}

# 3. Preparing the human population density and null rasters

human_pop_density = raster("Human_pop_density.asc"); human_pop_density[human_pop_density[]<0] = 0
human_pop_density_log = human_pop_density; human_pop_density[] = log(human_pop_density[]+1)
rast = crop(human_pop_density, e_2); res(rast) = res(rasters_1[[1]]); extent(rast) = extent(rasters_1[[1]])
crs(rast) = CRS("+proj=longlat +datum=WGS84"); human_pop_density = rast
rast = crop(human_pop_density_log, e_2); res(rast) = res(rasters_1[[1]]); extent(rast) = extent(rasters_1[[1]])
crs(rast) = CRS("+proj=longlat +datum=WGS84"); human_pop_density_log = rast
null_raster = human_pop_density; null_raster[!is.na(null_raster[])] = 1

# 4. Defining the background area
		
	# 4.1. With a mask based on D_max

if (background_mask == TRUE)
	{
		# 4.1.1. Estimating d_max_list[[h]] with a jackknife

		d_max_list = list()
		for (h in 1:length(datasets))
			{	
				nberOfRepetitions = 1000; observations = observations_list[[h]]; d_maxs = c()
				for (i in 1:nberOfRepetitions)
					{
						d_max = 0; nberOfObservationsToDiscard = round(0.2*dim(observations)[1])
						observationsToDiscard = sample(1:dim(observations)[1], nberOfObservationsToDiscard, replace=F)
						observationsToMaintain =  c(1:dim(observations)[1])[!c(1:dim(observations)[1])%in%observationsToDiscard]
						x2 = observations[observationsToMaintain,1:2]
						for (j in 1:length(observationsToDiscard))
							{
								x1 = cbind(observations[observationsToDiscard[j],1], observations[observationsToDiscard[j],2])
								dS = rdist.earth(x1, x2, miles=F)
								if (d_max < min(dS)) d_max = min(dS)
							}
						d_maxs = c(d_maxs, d_max)
					}
				d_max_list[[h]] = mean(d_maxs)
			}
		
		# 4.1.2. Generating a mask based on d_max_list[[h]]
		
			# Note: spatial points have to be transformed with a Mercator* projection in order to be able to use the "gBuffer"
			# function with a radius expressed in meters. (*) The "Web Mercator" or "Pseudo Mercator" projection is a 
			# cylindrical map projection. This is a variant of the regular Mercator projection, except that the computation 
			# is done on a sphere, using the semi-major axis of the ellipsoid (source: proj4.org/operations/projections/webmerc.html)

		background_masks = list()
		for (h in 1:length(datasets))
			{
				df1 = data.frame(observations_list[[h]]); coordinates(df1) = ~ longitude + latitude
				projection(df1) = CRS("+init=epsg:4326"); df2 = spTransform(df1, CRS("+init=epsg:3857")); 
				pol1 = rgeos::gBuffer(df2, width=d_max_list[[h]]*1000)
				pol2 = spTransform(pol1, CRS("+init=epsg:4326")); background_masks[[h]] = pol2
			}
	}

	# 4.2. Only considering raster cells associated with occurrence data

if (background_cell == TRUE)
	{
		cellIDs = c()
		for (h in 1:length(datasets))
			{
				observations = observations_list[[h]]
				cellIDs = c(cellIDs, extract(null_raster, observations, cellnumbers=T))
				cellIDs = unique(cellIDs)
			}
		background_cells = cellIDs
	}

# 5. Plotting the unique observation points (with or without the mask)

for (h in 1:length(datasets))
	{
		if (!file.exists(paste0(datasetsDirectory,"/",datasets[h],".pdf")))
			{
				pdf(paste0(datasetsDirectory,"/",datasets[h],".pdf"), width=plotWidth, height=plotHeight)
				observations = observations_list[[h]] # dev.new(width=7, height=6.2)
				par(mar=c(0,0,0,0), oma=c(0.0,3.6,0.0,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
				plot(crop(background, e_2), col=cols_background, useRaster=T, colNA=NA, box=F, axes=F, legend=F)
				points(observations, col="gray30", pch=3, cex=0.5, lwd=0.5)
				rect(e_2@xmin, e_2@ymin, e_2@xmax, e_2@ymax, lwd=0.2, border="gray30")
				if (background_mask == TRUE) plot(background_masks[[h]], border="red", lwd=0.5, add=T)
				dev.off()
			}
	}

# 6. Maxent analyses with an heterogeneous background

	# 1st series of Maxent analyses based on all the observations ("Maxent_result_files_1") 
		# 1.1: resolution = "04", subsampling_adjacent_cells = F, background_cell = F, background_mask = F, human_pop_density_bias = F
		# 1.2: resolution = "16", subsampling_adjacent_cells = T, background_cell = T, background_mask = F, human_pop_density_bias = F
		# 1.3: resolution = "16", subsampling_adjacent_cells = F, background_cell = T, background_mask = F, human_pop_density_bias = F
	# 2nd series of Maxent analyses based the observations >=1970 ("Maxent_result_files_2")
		# 2.2: resolution = "16", subsampling_adjacent_cells = T, background_cell = T, background_mask = F, human_pop_density_bias = F
		# 2.3: resolution = "16", subsampling_adjacent_cells = F, background_cell = T, background_mask = F, human_pop_density_bias = F

for (h in 1:length(datasets))
	{
		if (!file.exists(paste0("Maxent_result_files_2/",datasets[h],".asc")))
			{
				print(datasets[h])
				
	# 6.1. Sampling pseudo-absences in the background (considering or not human pop. density)
		
				if (background_cell == TRUE)
					{
						null_raster[!(1:length(null_raster[]))%in%background_cells] = NA
						human_pop_density_log[!(1:length(human_pop_density_log[]))%in%background_cells] = NA
					}
				if ((background_mask == TRUE)&(human_pop_density_bias == TRUE))
					{
						background_r = mask(human_pop_density_log, background_masks[[h]]); stack = mask(rasters_stacks, background_masks[[h]])
					}
				if ((background_mask == TRUE)&(human_pop_density_bias == FALSE))
					{
						background_r = mask(null_raster, background_masks[[h]]); stack = mask(rasters_stacks, background_masks[[h]])
					}
				if ((background_mask == FALSE)&(human_pop_density_bias == FALSE))
					{
						background_r = null_raster; stack = rasters_stacks
					}
				observations = observations_list[[h]]
				pseudo_absences = xyFromCell(background_r, sample(which(!is.na(values(background_r))),10000,prob=values(background_r)[!is.na(values(background_r))]))
				presences = cbind(observations, rep(1,dim(observations)[1]))
				absences = cbind(pseudo_absences, rep(0,dim(pseudo_absences)[1]))
				colnames(absences)[1] = "longitude"; colnames(absences)[2] = "latitude"
				colnames(absences)[3] = "response"; colnames(presences)[3] = "response"

	# 6.2. Removing observations from adjacent cells (optional)
					
				if (subsampling_adjacent_cells == T)
					{
						r = rasters_stacks[[1]]; buffer = c()
						presences = presences[which(!is.na(cellFromXY(r, presences[,1:2]))),]
						cellIDs = cellFromXY(r, presences[,1:2])
						for (i in 1:length(unique(cellIDs)))
							{
								if (sum(cellIDs==unique(cellIDs)[i]) > 1)
									{
										tmp = presences[which(cellIDs==unique(cellIDs)[i]),]
										buffer = rbind(buffer, tmp[sample(1:dim(tmp)[1],1),])
									}	else		{
										buffer = rbind(buffer, presences[which(cellIDs==unique(cellIDs)[i]),])
									}
							}
						presences = buffer
						cellIDs1 = cellFromXY(r, presences[,1:2]); cellIDs2 = cellIDs1
						adjacents = adjacent(r, cellIDs1, direction=4)
						adjacents = adjacents[which(adjacents[,1]%in%cellIDs1),]
						adjacents = adjacents[which(adjacents[,2]%in%cellIDs1),]
						while (dim(adjacents)[1] != 0)
							{
								cellIDs_to_remove = unique(c(adjacents[,1],adjacents[,2]))
								maxRepresentations = 0; index = 0
								for (i in 1:length(cellIDs_to_remove))
									{
										representations = sum(c(adjacents[,1],adjacents[,2])==cellIDs_to_remove[i])
										if (maxRepresentations < representations)
											{
												maxRepresentations = representations; index = i
											}
									}
								if (maxRepresentations > 1)
									{
										cellID_to_remove = cellIDs_to_remove[index]			
									}	else		{
										cellID_to_remove = sampl(cellIDs_to_remove,1)
									}
								cellIDs2 = cellIDs2[-which(cellIDs2==cellID_to_remove)]
								adjacents = adjacent(r, cellIDs2, direction=4)
								adjacents = adjacents[which(adjacents[,1]%in%cellIDs2),]
								adjacents = adjacents[which(adjacents[,2]%in%cellIDs2),]
								# print(dim(adjacents))
							}
						presences = presences[which(cellIDs1%in%cellIDs2),]
					}
				data = rbind(presences,absences)
				
	# 6.3. Keeping only one presence or pseudo-absence point per raster cell (priority = presence points)
					
				r = rasters_stacks[[1]]; buffer = c()
				data = data[which(!is.na(cellFromXY(r, data[,1:2]))),]
				cellIDs = cellFromXY(r, data[,1:2])
				for (i in 1:length(unique(cellIDs)))
					{
						if (sum(cellIDs==unique(cellIDs)[i]) > 1)
							{
								tmp = data[which(cellIDs==unique(cellIDs)[i]),]
								if (sum(tmp[,"response"]==1) == 0)
									{
										buffer = rbind(buffer, tmp[sample(1:dim(tmp)[1],1),])
									}
								if (sum(tmp[,"response"]==1) == 1)
									{
										buffer = rbind(buffer, tmp[which(tmp[,"response"]==1),])
									}
								if (sum(tmp[,"response"]==1) >= 2)
									{
										indices = which(tmp[,"response"]==1)
										buffer = rbind(buffer, tmp[sample(indices,1),])
									}
							}	else	{
								buffer = rbind(buffer, data[which(cellIDs==unique(cellIDs)[i]),])
							}
					}
				data = buffer
				
	# 6.4. Maxent analysis in itself
				
				maxent = maxent(x=stack, p=data[which(data[,3]==1),1:2], b=data[which(data[,3]==0),1:2], args=c("randomseed=true","randomtestpoints=20",
					   "replicates=10","jackknife","replicatetype=subsample","noaddsamplestobackground","maximumiterations=10000","writeplotdata=true"))
				write.csv(maxent@results, paste0("Maxent_result_files_2/",datasets[h],".csv"), quote=F)
				saveRDS(maxent, paste0("Maxent_result_files_2/",datasets[h],".RData"))
				map_t0 = predict(maxent, rasters_stacks); writeRaster(mean(map_t0), paste0("Maxent_result_files_2/",datasets[h],".asc"), overwrite=T)
				rast = raster(paste0("Maxent_result_files_2/",datasets[h],".asc"))
				pdf(paste0("Maxent_result_files_2/",datasets[h],".pdf"), width=plotWidth, height=plotHeight)
				par(mar=c(0,0,0,0), oma=c(0.0,3.6,0.0,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
				plot(crop(rast, e_2), col=cols, useRaster=T, colNA="gray90", box=F, axes=F, legend=F)
				plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.815,0.825,0.14,0.34),
					 axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.15, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,0.4,0)), alpha=1)
				dev.off()
			}
	}

	# 6.5. Comparing the results of different Maxent analyses

analyses = c(2:3); AUCs = matrix(nrow=length(datasets), ncol=length(analyses)*2); row.names(AUCs) = datasets
colnames(AUCs) = c("Training_AUC_2","Training_AUC_3","Test_AUC_2","Test_AUC_3")
for (i in 1:length(datasets))
	{
		for (j in 1:length(analyses))
			{
				csv = read.csv(paste0("Maxent_result_files_2/",datasets[i],"_",analyses[j],".csv"), header=T)
				lineIndex = which(csv[,1]=="Training.AUC"); AUCs[i,j] = csv[lineIndex,"species..average."]
				lineIndex = which(csv[,1]=="Test.AUC"); AUCs[i,j+2] = csv[lineIndex,"species..average."]
			}
	}
write.csv(AUCs, "Maxent_AUCs_2.csv", quote=F)

selected_analyses_3 = c("B_armeniacus","B_balteatus","B_fragrans","B_muscorum","B_norvegicus","B_pascuorum","B_rupestris","B_sichelii","B_sulfureus")
analyses = rep("2",length(datasets)); analyses[datasets%in%selected_analyses_3] = "3"

	# 6.6. Investigating the different response curves

predictors = c("grasslands","barren_vegetation","savannas","woody_savannas","open_shrublands","closed_shrubland",
			   "deciduous_broadleaf_forest","deciduous_needleleaf_forest","evergreen_broadleaf_forest","evergreen_needleleaf_forest",
			   "mixed_forests","wetlands","snow_ice","croplands","urban_areas","temperature_seasonality","annual_mean_temperature",
			   "precipitation_seasonality","annual_precipitation"); predictor_values = matrix(nrow=3, ncol=length(predictors))
colnames(predictor_values) = predictors; row.names(predictor_values) = c("median","minV","maxV")
for (i in 1:length(predictors))
	{
		index = which(colnames(head(rasters_stacks))==predictors[i])
		rast = rasters_stacks[[index]]; minV = min(rast[], na.rm=T); maxV = max(rast[], na.rm=T)
		predictor_values[,i] = cbind(median(rast[], na.rm=T), minV, maxV)
	}
dev.new(height=7, width=7.7); par(mfrow=c(5,4), oma=c(1.5,1.5,1,1), mar=c(2,2.3,1,1))
for (i in 1:length(predictors))
	{
		df = data.frame(matrix(nrow=length(seq(predictor_values["minV",i],predictor_values["maxV",i],1)),ncol=length(predictors)))
		colnames(df) = predictors
		for (j in 1:length(predictors))
			{
				if (i == j) df[,predictors[j]] = seq(predictor_values["minV",j],predictor_values["maxV",j],1)
				if (i != j) df[,predictors[j]] = rep(predictor_values["median",j],dim(df)[1])
			}
		for (j in 1:length(datasets))
			{
				maxent = readRDS(paste0("Maxent_result_files_2/",datasets[j],"_",analyses[j],".RData"))
				best_AUC = which(maxent@results["Training.AUC",1:10]==max(maxent@results["Training.AUC",1:10]))[1]
				best_model = maxent@models[[best_AUC]]
				p = predict(best_model, df, type="response")
				if (j == 1)
					{
						minX = min(df[,predictors[i]]); maxX = max(df[,predictors[i]]); minY = min(p); maxY = max(p)
					}	else	{
						if (minX > min(df[,predictors[i]])) minX = min(df[,predictors[i]])
						if (maxX < max(df[,predictors[i]])) maxX = max(df[,predictors[i]])
						if (minY > min(p)) minY = min(p)
						if (maxY < max(p)) maxY = max(p)
					}
			}
		for (j in 1:length(datasets))
			{
				maxent = readRDS(paste0("Maxent_result_files_2/",datasets[j],"_",analyses[j],".RData"))
				best_AUC = which(maxent@results["Training.AUC",1:10]==max(maxent@results["Training.AUC",1:10]))[1]
				best_model = maxent@models[[best_AUC]]
				p = predict(best_model, df, type="response")
				if (j == 1)
					{
						plot(df[,predictors[i]], p, col="gray30", ann=F, axes=F, lwd=0.1, type="l", xlim=c(minX,maxX), ylim=c(minY, maxY))
					}	else	{
						lines(df[,predictors[i]], p, col="gray30", lwd=0.1)
					}
			}
		axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.1,0))
		axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.3,0))
		title(xlab=gsub("_"," ",predictors[i]), cex.lab=0.8, mgp=c(0.9,0,0), col.lab="gray30")
		title(ylab="predicted values", cex.lab=0.8, mgp=c(1.3,0,0), col.lab="gray30")
		box(lwd=0.2, col="gray30")
	}
dev.copy2pdf(file="Maxent_responses_NEW.pdf"); dev.off()

	# 6.7. Extracting and plotting the permutation importances

original_names = c("precipitation_of_wettest_quarter.1","precipitation_of_wettest_quarter.2",
				   "precipitation_of_wettest_quarter.3","precipitation_of_wettest_quarter.4")
corrected_names = c("annual_mean_temperature","temperature_seasonality","annual_precipitation","precipitation_seasonality")
predictors = c("grasslands","barren_vegetation","savannas","woody_savannas","open_shrublands","closed_shrubland",
			   "deciduous_broadleaf_forest","deciduous_needleleaf_forest","evergreen_broadleaf_forest","evergreen_needleleaf_forest",
			   "mixed_forests","wetlands","snow_ice","croplands","urban_areas","temperature_seasonality","annual_mean_temperature",
			   "precipitation_seasonality","annual_precipitation"); predictors = rev(predictors)
colours = c("gold","darkkhaki","goldenrod2","goldenrod3","darkseagreen1","darkseagreen2","darkolivegreen1","darkolivegreen2",
			"darkolivegreen3","darkolivegreen4","forestgreen","darkslategray2","darkslategray3","lemonchiffon2","darkgrey",
			"firebrick1","firebrick3","dodgerblue1","dodgerblue3"); colours = rev(colours)
importances = matrix(nrow=length(datasets), ncol=19); row.names(importances) = datasets; colnames(importances) = predictors
for (i in 1:length(datasets))
	{
		csv = read.csv(paste0("Maxent_result_files_2/",datasets[i],"_",analyses[i],".csv"), header=T)
		csv[,1] = as.character(csv[,1])
		for (j in 1:length(original_names))
			{
				for (k in 1:dim(csv)[1]) csv[k,1] = gsub(original_names[j],corrected_names[j],csv[k,1])
			}
		for (j in 1:length(predictors))
			{
				lineIndex = which(csv[,1]==paste0(predictors[j],".permutation.importance"))
				importances[i,j] = csv[lineIndex,"species..average."]
			}
	}
write.csv(importances, "Maxent_importances.csv", quote=F)
importances = read.csv("Maxent_importances.csv", header=T)
importances = importances[as.character(species_information[,"species"]),]

dev.new(width=12.7, height=3.5); par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(7.0,2.0,1.7,1.5), lwd=0.2)
barplot(t(importances), col=colours, ann=F, axes=F, axisnames=F, lwd=1.3, space=0, width=1, xlim=c(0,77), border="gray30")
axis(side=1, lwd.tick=0, cex.axis=0.75, lwd=0, tck=-0.020, col.axis="gray30", mgp=c(0,0,0), at=seq(0.5,67.5,1), labels=gsub("B_","B. ",row.names(importances)), las=2, line=0.5)
axis(side=2, lwd.tick=0.2, cex.axis=0.70, lwd=0.2, tck=-0.030, col.axis="gray30", mgp=c(0,0.3,0), at=seq(0,100,20), pos=-1)
title(ylab="Importance of the predictor (%)", cex.lab=0.75, mgp=c(-0.2,0,0), col.lab="gray30")
legend(68.5, 105, rev(gsub("_"," ",predictors)), col=rev(colours), text.col="gray30", pch=16, pt.cex=1.6, box.lty=0, cex=0.7, y.intersp=1.3)
dev.copy2pdf(file="Maxent_importances1_NEW.pdf"); dev.off()

library(ade4)
importances = read.csv("Maxent_importances.csv", header=T); colours = c()
pca = dudi.pca(importances, scannf=F, nf=19); lis = pca$li[,1:2]; cos = pca$co
for (i in 1:dim(lis)[1])
	{
		index = which(species_information[,"species"]==row.names(lis)[i])
		subgenus = species_information[index,"subgenera"]
		colours = c(colours, subgenera_cols[which(subgenera==subgenus)])
	}

dev.new(width=7, height=6); par(mar=c(3,3,0,1.5))
plot(lis, col=colours, cex=0.7, pch=16, ann=F, axes=F, xlim=c(-6.3,4.1), ylim=c(-4,9.5))
text(lis[,1], lis[,2], labels=gsub("B_","",row.names(lis)), cex=0.6, col="gray50", pos=4, offset=0.25)
points(lis, col=colours, cex=0.9, pch=16); points(lis, col="gray30", cex=0.9, pch=1, lwd=0.3)
axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.05,0), at=seq(-8,5,1))
axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.30,0), at=seq(-5,9,1))
title(xlab="PCA axis 1", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
title(ylab="PCA axis 2", cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
legend(2.3, 9.3, subgenera, col=subgenera_cols, text.col="gray30", pch=16, pt.cex=1.1, box.lty=0, cex=0.6, y.intersp=1.2)
dev.copy2pdf(file="Maxent_importances2_NEW.pdf"); dev.off()

	# 6.8. Plotting the best SDM result for each species

dev.new(height=7.7, width=5.6); e_3 = extent(-13, 62, 30, 72) # 5 maps on each row
par(mfrow=c(7,5), mar=c(0,0,0,0), oma=c(1.3,1.2,0.30,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in 1:35)
	{
		rast = crop(raster(paste0("Maxent_result_files_2/",datasets[i],"_",analyses[i],".asc")), e_3)
		observations = read.csv(paste0(datasetsDirectory,"/",datasets[i],".csv"), header=T)
		observations = unique(observations[,c("longitude","latitude")])
		observations = observations[!is.na(extract(rast, observations)),]
		plot(rast, col="white", useRaster=T, colNA="gray90", box=F, axes=F, legend=F)
		points(observations, pch=16, cex=0.3, col="gray30")
		rect(e_3@xmin, e_3@ymin, e_3@xmax, e_3@ymax, lwd=0.2, border="gray30")
		mtext(side=1, gsub("B_","B. ",datasets[i]), cex=0.55, col="gray30", line=-0.8)
	}
dev.copy2pdf(file="Observations_map1_NEW.pdf"); dev.off()
dev.new(height=7.7, width=5.6); e_3 = extent(-13, 62, 30, 72)
par(mfrow=c(7,5), mar=c(0,0,0,0), oma=c(1.3,1.2,0.30,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in 36:length(datasets))
	{
		rast = crop(raster(paste0("Maxent_result_files_2/",datasets[i],"_",analyses[i],".asc")), e_3)
		observations = read.csv(paste0(datasetsDirectory,"/",datasets[i],".csv"), header=T)
		observations = unique(observations[,c("longitude","latitude")])
		observations = observations[!is.na(extract(rast, observations)),]
		plot(rast, col="white", useRaster=T, colNA="gray90", box=F, axes=F, legend=F)
		points(observations, pch=16, cex=0.3, col="gray30")
		rect(e_3@xmin, e_3@ymin, e_3@xmax, e_3@ymax, lwd=0.2, border="gray30")
		mtext(side=1, gsub("B_","B. ",datasets[i]), cex=0.55, col="gray30", line=-0.8)
	}
dev.copy2pdf(file="Observations_map2_NEW.pdf"); dev.off()

lakes = shapefile("Natural_Earth_GIS_files/Natural_Earth_lakes.shp")
template = crop(raster(paste0("Maxent_result_files_2/",datasets[1],"_",analyses[1],".asc")), e_3)
template[!is.na(template[])] = 0; template = mask(template, lakes)
cols = colorRampPalette(brewer.pal(9,"YlOrBr"))(121)[1:101]
dev.new(height=7.7, width=5.6); e_3 = extent(-13, 62, 30, 72) # 5 maps on each row
par(mfrow=c(7,5), mar=c(0,0,0,0), oma=c(1.3,1.2,0.30,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in 1:35)
	{
		rast = crop(raster(paste0("Maxent_result_files_2/",datasets[i],"_",analyses[i],".asc")), e_3)
		rast[template[]==0] = NA; plot(rast, col=cols, useRaster=T, colNA="gray90", box=F, axes=F, legend=F)
		rect(e_3@xmin, e_3@ymin, e_3@xmax, e_3@ymax, lwd=0.2, border="gray30")
		mtext(side=1, gsub("B_","B. ",datasets[i]), cex=0.55, col="gray30", line=-0.8)
	}
dev.copy2pdf(file="Maxent_predictions1_NEW.pdf"); dev.off()
dev.new(height=7.7, width=5.6); e_3 = extent(-13, 62, 30, 72) # 5 maps on each row
par(mfrow=c(7,5), mar=c(0,0,0,0), oma=c(1.3,1.2,0.30,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in 36:length(datasets))
	{
		rast = crop(raster(paste0("Maxent_result_files_2/",datasets[i],"_",analyses[i],".asc")), e_3)
		rast[template[]==0] = NA; plot(rast, col=cols, useRaster=T, colNA="gray90", box=F, axes=F, legend=F)
		rect(e_3@xmin, e_3@ymin, e_3@xmax, e_3@ymax, lwd=0.2, border="gray30")
		mtext(side=1, gsub("B_","B. ",datasets[i]), cex=0.55, col="gray30", line=-0.8)
	}
dev.copy2pdf(file="Maxent_predictions2_NEW.pdf"); dev.off()
dev.new(height=7.7, width=5.6); e_3 = extent(-13, 62, 30, 72) # 5 maps on each row
par(mfrow=c(7,5), mar=c(0,0,0,0), oma=c(1.3,1.2,0.30,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in c(29,63,38,45,57,58,9,24))
	{
		rast = crop(raster(paste0("Maxent_result_files_2/",datasets[i],"_",analyses[i],".asc")), e_3)
		rast[template[]==0] = NA; plot(rast, col=cols, useRaster=T, colNA="gray90", box=F, axes=F, legend=F)
		rect(e_3@xmin, e_3@ymin, e_3@xmax, e_3@ymax, lwd=0.2, border="gray30")
		mtext(side=1, gsub("B_","B. ",datasets[i]), cex=0.55, col="gray30", line=-0.8)
	}
dev.copy2pdf(file="Selected_predictions_NEW.pdf"); dev.off()

	# 6.9. Analysing the phylogenetic signal associated with importance

library(ape); library(phytools); library(picante)
tree1 = read.nexus("Bombus_Cameron.tree")[[1]]
tree1$tip.label = gsub("portchins", "portschins", tree1$tip.label)
tree1$tip.label = gsub("sichelli", "sichelii", tree1$tip.label)
tips1 = tree1$tip.label; tips_to_remove = c()
for (i in 1:length(tips1))
	{
		if (sum(grepl(tips1[i],datasets))!=1) tips_to_remove = c(tips_to_remove, tips1[i])
	}
tree2 = drop.tip(tree1, tips_to_remove)
tips2 = tree2$tip.label; tips_missings = c(); tips3 = tips2
for (i in 1:length(datasets))
	{
		sp_found = FALSE
		for (j in 1:length(tips2))
			{
				if (grepl(tips2[j],datasets[i]))
					{
						sp_found = TRUE; tips3[j] = datasets[i]
					}
			}
		if (sp_found == FALSE) tips_missings = c(tips_missings, datasets[i])
	}
tree2$tip.label = tips3

traits = read.csv("Maxent_importances.csv", head=T)
rowNames = traits[1:dim(traits)[1],1]
traits = traits[,2:dim(traits)[2]]
row.names(traits) = rowNames
indices = which(row.names(traits)%in%tree2$tip.label)
traits = traits[indices,]; traits = traits[tree2$tip.label,]

K_pValues = matrix(nrow=dim(traits)[2], ncol=2)
row.names(K_pValues) = colnames(traits)
colnames(K_pValues) = c("K","p-value")
for (i in 1:dim(traits)[2])
	{
		values = matrix(nrow=dim(traits)[1], ncol=1)
		values[,1] = traits[,colnames(traits)[i]]
		row.names(values) = row.names(traits)
		test = phytools::phylosig(tree2, values, test=T, method="K")
		K_pValues[i,] = cbind(test$K, test$P)
	}
K_pValues = round(K_pValues, 3); print(K_pValues)
write.csv(K_pValues, "Phylogenetic_signals.csv", quote=F)

	# to plot the croplands "importance" on tips:

values = matrix(nrow=dim(traits)[1], ncol=1)
values[,1] = traits[tree2$tip.label,"croplands"]
colourScale = colorRampPalette(brewer.pal(9,"YlGn"))(101)
cols = colourScale[(((values-min(values))/(max(values)-min(values)))*100)+1]
tree3 = tree2; tree3$tip.label = gsub("B_", "    B. ", tree3$tip.label)
dev.new(width=7, height=7.5); par(oma=c(0,0,0,0), mar=c(0.8,1.2,0.5,2.2), lwd=0.2)
plot(tree3, show.tip.label=T, show.node.label=T, edge.width=0.5, cex=0.6, align.tip.label=3, col="gray30")
l1 = length(tree3$tip.label)+1; l2 = (2*length(tree3$tip.label))-2
for (i in 1:dim(tree3$edge)[1])
	{
		if (!tree3$edge[i,2]%in%tree3$edge[,1])
			{
				nodelabels(node=tree3$edge[i,2], pch=16, cex=1.2, col=cols[tree3$edge[i,2]])
				nodelabels(node=tree3$edge[i,2], pch=1, cex=1.2, col="gray30", lwd=0.5)
			}
	}
mat = matrix(nrow=1, ncol=2); mat[1,1] = min(values); mat[1,2] = max(values); rast = raster(matrix(mat))
plot(rast, legend.only=T, add=T, col=colorRampPalette(brewer.pal(9,"YlGn"))(101), legend.width=0.5, alpha=1,
	 legend.shrink=0.3, smallplot=c(0.11,0.123,0.15,0.50), legend.args=list(text="", cex=0.8, line=0.5, col="gray30"),
	 axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, tck=-0.5, line=0, mgp=c(0,0.4,0), at=seq(0,35,5),
	 labels=c("0.00","0.05","0.10","0.15","0.20","0.25","0.30","0.35")))
add.scale.bar(length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
title1 = "Importance of"; mtext(title1, side=3, line=-16.1, at=0.10, cex=0.75, font=1, col="gray30")
title2 = "croplands as"; mtext(title2, side=3, line=-16.9, at=0.10, cex=0.75, font=1, col="gray30")
title3 = "a SDM predictor"; mtext(title3, side=3, line=-17.7, at=0.10, cex=0.75, font=1, col="gray30")
dev.copy2pdf(file="MCC_tree_croplands_NEW.pdf"); dev.off()

	# 6.10. Analysing the association between species range and bioclimatic importances

species_information = read.csv("Species_information.csv", header=T)
bioclim_variables = c("annual_precipitation","precipitation_seasonality","annual_mean_temperature","temperature_seasonality")
importances = read.csv("Maxent_importances.csv", header=T); colours = c()
bioclim_importances = rep(NA, dim(species_information)[1])
for (i in 1:dim(species_information)[1])
	{
		index = which(row.names(importances)==species_information[i,"species"])
		bioclim_importances[i] = sum(importances[index,bioclim_variables])
		subgenus = species_information[i,"subgenera"]
		colours = c(colours, subgenera_cols[which(subgenera==subgenus)])
	}
response = bioclim_importances
predictor = species_information[,"area_km2"]
formula = as.formula(paste0("response ~ predictor"))
LR = lm(formula); R2 = summary(LR)$r.squared

dev.new(width=4.5, height=6); par(mar=c(1.7,2.2,0.5,1.5)); options(scipen=9)
plot(species_information[,"area_km2"], bioclim_importances, col=colours, cex=0.7, pch=16, ann=F, axes=F, xlim=c(0,1.25*10^7), ylim=c(0,100))
text(species_information[,"area_km2"], bioclim_importances, labels=gsub("B_","",species_information[,"species"]), cex=0.6, col="gray50", pos=4, offset=0.25)
axis(side=1, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,-0.02,0), at=seq(0,1.5*10^7,0.3*10^7))
axis(side=2, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.010, col.axis="gray30", mgp=c(0,0.18,0), at=seq(0,100,10))
title(xlab="Approximated range size (km2)", cex.lab=0.7, mgp=c(-0.2,0,0), col.lab="gray30")
title(ylab="Total importance of bioclimatic variables (%)", cex.lab=0.7, mgp=c(0.4,0,0), col.lab="gray30")
abline(LR, lwd=1.5, col="gray50", lty=2)
for (i in 1:length(bioclim_importances))
	{
		points(species_information[i,"area_km2"], bioclim_importances[i], col=colours[i], cex=0.9, pch=16)
		points(species_information[i,"area_km2"], bioclim_importances[i], col="gray30", cex=0.9, pch=1, lwd=0.3)
	}
legend(0.93*10^7, 37.0, subgenera, col=subgenera_cols, text.col="gray30", pch=16, pt.cex=1.0, box.lty=0, cex=0.55, y.intersp=1.1)
dev.copy2pdf(file="Bioclim_range_sizes_NEW.pdf"); dev.off()

df = as.data.frame(cbind(bioclim_importances,as.character(species_information[,"pop_trend_IUCN_Europe"])))
colnames(df) = c("bioclim_importances","pop_trend_IUCN_Europe")
df[,"pop_trend_IUCN_Europe"] = as.character(df[,"pop_trend_IUCN_Europe"])
df[df[,"pop_trend_IUCN_Europe"]=="unknown","pop_trend_IUCN_Europe"] = 1
df[df[,"pop_trend_IUCN_Europe"]=="not_assessed","pop_trend_IUCN_Europe"] = 2
df[df[,"pop_trend_IUCN_Europe"]=="decreasing","pop_trend_IUCN_Europe"] = 3
df[df[,"pop_trend_IUCN_Europe"]=="stable","pop_trend_IUCN_Europe"] = 4
df[df[,"pop_trend_IUCN_Europe"]=="increasing","pop_trend_IUCN_Europe"] = 5
df[,"bioclim_importances"] = as.numeric(df[,"bioclim_importances"])
df[,"pop_trend_IUCN_Europe"] = as.numeric(df[,"pop_trend_IUCN_Europe"])
anova = aov(as.formula("bioclim_importances ~ pop_trend_IUCN_Europe"), data=df); summary(anova)

df = as.data.frame(cbind(bioclim_importances,as.character(species_information[,"status_IUCN_Europe"])))
colnames(df) = c("bioclim_importances","status_IUCN_Europe")
df[,"status_IUCN_Europe"] = as.character(df[,"status_IUCN_Europe"])
df[df[,"status_IUCN_Europe"]=="data_deficient","status_IUCN_Europe"] = 1
df[df[,"status_IUCN_Europe"]=="not_assessed","status_IUCN_Europe"] = 2
df[df[,"status_IUCN_Europe"]=="least_concern","status_IUCN_Europe"] = 3
df[df[,"status_IUCN_Europe"]=="vulnerable","status_IUCN_Europe"] = 4
df[df[,"status_IUCN_Europe"]=="near_threatened","status_IUCN_Europe"] = 5
df[df[,"status_IUCN_Europe"]=="endangered","status_IUCN_Europe"] = 6
df[df[,"status_IUCN_Europe"]=="critically_endangered","status_IUCN_Europe"] = 7
df[,"bioclim_importances"] = as.numeric(df[,"bioclim_importances"])
df[,"status_IUCN_Europe"] = as.numeric(df[,"status_IUCN_Europe"])
anova = aov(as.formula("bioclim_importances ~ status_IUCN_Europe"), data=df); summary(anova)

