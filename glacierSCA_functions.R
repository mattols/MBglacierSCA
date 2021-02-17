#
# Glacier Trend Functions
# MODDRFS & MODSCAG data input
# Available:
#     https://snow-data.jpl.nasa.gov/moddrfs
#     https://snow-data.jpl.nasa.gov/modscag
# Date: 11/12/2019
#
# # # # # # # # # # # #

library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)
library(reshape2)
library(dplyr)

##### LOAD PATHS ###############################################################################
load.paths <- function(macpaths=FALSE){
  # Specify whether using mac or not
  if(macpaths){
    rfpath <<- "~/FOLDER_RF/"
    scapath <<- "/FOLDER_SCA/"
  } else{
    rfpath <- 'F:/HiMAT/MATTO/FOLDER_RF/'
    scapath <- 'F:/HiMAT/MATTO/FOLDER_SCA/'
    tfpath <- 'F:/HiMAT/MATTO/TOPO_FORCING/'
  }
}

##### LOAD EXTRA DATA ###############################################################################
load.data <- function(macpaths=FALSE){
  if(macpaths){
    rgi <<- readOGR(dsn="/RGI_GLACIERS", layer="rgi50_josh_mb")
    dpath <- "/DEM_REGION.tif"
    demL <<- raster(dpath)
  } else{
    rgi <<- readOGR(dsn="F:/HiMAT/MATTO/RGI_GLACIERS", layer="rgi50_josh_mb")
    dpath <- "F:/HiMAT/MATTO/DEM_REGION.tif"
    demL <<- raster(dpath)
  }
}

##### LOAD NEW CSV ###############################################################################
mb.csv.new <- function(macpaths=FALSE){
  #
  # LOAD IN MASS BALANCE DATA .CSV
  # adds field for: 
  # rgi <region>, hydro <sub_bas>, and <sub_bas2> (with some basins aggregated)
  #
  if(macpaths){
    mb_path <- "~/MAURER_GLACIERMB.csv"
    bpath = "RGI_SUBREGION.shp"
  } else{
    mb_path <- "F:/HiMAT/MATTO/MAURER_GLACIERMB.csv"
    bpath = "F:/HiMAT/MATTO/RGI_SUBREGION.shp"
  }
  #
  # READ IN MB .CSV
  mb = read.csv(mb_path)
  mb$Object.ID = gsub("'",'', mb$Object.ID) # Fix character string of glacier name
  mb$Object.ID = gsub("[[:punct:]]","_",mb$Object.ID) 
  mb = mb[!is.na(mb$lon),] 
  mbc = mb[mb$category == levels(mb$category)[1],] # Only clean ice glaciers
  mbg = mbc[,c(1,7,10,12,16,18)] # only keep area2016, meanElev, meanAspect, meanELevChj, and geoMassBal
  mbg.coords = mbc[,4:5]
  #
  # EXTRACT REGION | Tibet 8 | West H 12 | Central H 13 | East H 14
  sp.mb = SpatialPointsDataFrame(mbg.coords, mbg, proj4string = CRS("+proj=longlat +datum=WGS84"))
  basin = readOGR(bpath)
  sp.mb2 = spTransform(sp.mb,crs(basin))
  mbg$region = as.factor(extract(basin, sp.mb2)$poly.ID)
  #
  # CHANGE REGION 8 - TIBET 
  mbg$region[mbg$region == '8' & sp.mb$lon > 86] <- 14 # east_H
  mbg$region[mbg$region == '8' & sp.mb$lon < 86] <- 13 # central_H
  mbg$region = as.factor(as.character(mbg$region))
  # | West H 12 | Central H 13 | East H 14 | #
  #
  # get hydrobasin
  hydro <- readOGR("~/HMA_HYDROBASINS.shp")
  mb.bas = over(sp.mb2, hydro)
  mbg$sub_bas = as.factor(mb.bas$SUB_BAS)
  mbg$sub_bas[is.na(mbg$sub_bas)] <- as.factor(61023)
  mbg$sub_bas2 = mbg$sub_bas 
  # aggregate subbasins (higher n)
  # add all 61009 -> 61012 | 45006 -> 45009 | 45052 -> 45042 | 45050 -> 45042 
  # 45036 -> 45037 | 45035 -> 45020 | 45058 -> 45062 | NA -> 61023
  mbg$sub_bas2[mbg$sub_bas2 == "61009"] <- "61012"; mbg$sub_bas2[mbg$sub_bas2 == "45006"] <- "45009"
  mbg$sub_bas2[mbg$sub_bas2 == "45052"] <- "45042"; mbg$sub_bas2[mbg$sub_bas2 == "45050"] <- "45042"
  mbg$sub_bas2[mbg$sub_bas2 == "45036"] <- "45037"; mbg$sub_bas2[mbg$sub_bas2 == "45035"] <- "45020"
  mbg$sub_bas2[mbg$sub_bas2 == "45058"] <- "45062"; mbg$sub_bas2[mbg$sub_bas2 == "45070"] <- "45044"
  mbg$sub_bas2 <- as.factor(mbg$sub_bas2)
  return(mbg)
}


##############################################################################################################################
### TREND FUNCTIONS ################################################################################################################
##############################################################################################################################

## DONE
##############################################################################################################################
glacier.var <- function(glacier_id, varpath, rgi, demL, noZero=TRUE, threshold=NULL){
  #
    # FOR GIVEN GLACIER                                       #
    # create VARIABLES for glacier_i                          #
    # mask with inner polygon                                 #
    #                                                         #
    # Returns:                                                #
    # vstk | velv | vdates 
  #
  # READ IN VARIABLE STACK
  var_match = list.files(varpath, full.names = T, pattern = paste0(glacier_id,".*grd"))
  v_stk <- stack(var_match)
  #
  # FIND RGI SHAPEFILE AND MASK FOR INSIDE POLYGON
  shape = rgi[gsub("[[:punct:]]", "_", rgi@data$RGIId) == glacier_id,]
  cent_distance <- res(v_stk)[1]/2
  suppressWarnings(inner_shape <- buffer(shape,width=-cent_distance))
  #
  # DEM
  dem <- crop(demL,v_stk)
  velv <<- mask(resample(dem, v_stk), inner_shape)
  #
  # MASK VAR FOR ALL YEARS
  vstk <- mask(v_stk, inner_shape)
  #
  # DELETE ALL 0vals
  if(noZero){vstk[vstk==0] <- NA}
  # THRESHOLD
  if (!is.null(threshold)){
    # add threshold
    # cat("! Error threshold for RF not dependent on SCA !\n  SCA threshold valid")
    vstk[vstk < threshold] <- NA
    print(paste("SCA threshold of >", threshold))
  }
  #
  # EXTRACT DATES
  v_n = substring(names(vstk),2,8)
  vdates <- as.Date(v_n,format="%Y%j")
  #
  # EXCLUDE 2017, 2018 VALS
  # didx <- format(vdates, "%Y") != "2018"
  didx <- format(vdates, "%Y") != "2018" & format(vdates, "%Y") != "2017"
  vdates <<- vdates[didx]
  vstk <<- vstk[[which(didx)]]
  #
  # cat("Created variables: \n   vstk | velv | vdates \n <glacier i> \n")  
  # print(glacier_id)
}

## New - 04/22/2020
##############################################################################################################################
glacier.var2 <- function(glacier_id, varpath, rgi, demL, makeZero=TRUE, threshold=NULL){
  #
  # FOR GIVEN GLACIER                                       #
  # create VARIABLES for glacier_i                          #
  # mask with inner polygon                                 #
  #                                                         #
  # Returns:                                                #
  # vstk | velv | vdates | shape
  #
  # READ IN VARIABLE STACK
  var_match = list.files(varpath, full.names = T, pattern = paste0(glacier_id,".*grd"))
  v_stk <- stack(var_match)
  #
  # FIND RGI SHAPEFILE AND MASK FOR INSIDE POLYGON
  shape <<- rgi[gsub("[[:punct:]]", "_", rgi@data$RGIId) == glacier_id,]
  cent_distance <- res(v_stk)[1]/2
  suppressWarnings(inner_shape <- buffer(shape,width=-cent_distance))
  #
  # DEM
  dem <- crop(demL,v_stk)
  velv <<- mask(resample(dem, v_stk), inner_shape)
  #
  # MASK VAR FOR ALL YEARS
  vstk <- mask(v_stk, inner_shape)
  #
  # THRESHOLD
  if (!is.null(threshold)){
    # add threshold
    # cat("! Error threshold for RF not dependent on SCA !\n  SCA threshold valid")
    vstk[vstk < threshold] <- NA
    print(paste("SCA threshold of >", threshold))
  }
  #
  # MAKE NA VALUES WITHIN glacier 0 
  if(makeZero){vstk[is.na(vstk)] <- 0;vstk <- mask(v_stk, inner_shape)}
  #
  # EXTRACT DATES
  v_n = substring(names(vstk),2,8)
  vdates <- as.Date(v_n,format="%Y%j")
  #
  # EXCLUDE 2017, 2018 VALS
  # didx <- format(vdates, "%Y") != "2018"
  didx <- format(vdates, "%Y") != "2018" & format(vdates, "%Y") != "2017"
  vdates <<- vdates[didx]
  vstk <<- vstk[[which(didx)]]
  #
  # cat("Created variables: \n   vstk | velv | vdates \n <glacier i> \n")  
  # print(glacier_id)
}

## DONE
##############################################################################################################################
glacier.var.vals <- function(df_glaciers, varpath, rgi, demL, threshold=NULL){
  strt = Sys.time(); ci=ncol(df_glaciers)
  o_names <- names(df_glaciers)
  df_glaciers[,10:50] <- NA # was 38
  ## LOOP THROUGH EACH GLACIER
  # !!!!
  # quick fix (last name has no SCA)
  for (i in 1:(length(df_glaciers$Object.ID)-1)){
    g_name <- df_glaciers$Object.ID[i]
    glacier.var2(g_name, varpath, rgi, demL, threshold=threshold)
    # have vstk, vdates, velv 
    
    # OUTSIDE FUNCTIONS
    trend_vals = year.season.trends(vstk, vdates, velv)
    
    # add values to data frame
    df_glaciers[i,(ci+1)] <- mean(getValues(vstk),na.rm=T) # mean over all years
    df_glaciers[i,(ci+2)] <- sum(getValues(vstk),na.rm=T) # cumulative sum
    df_glaciers[i,((ci+3):(ci+41))] <- trend_vals # was 29
    
    # if(vvar=="RF"){
    #   # RF specific
    #   
    # } else{
    #   # SCA specific
    # }
    
    print(i)
    # add values to dataframe
    
    # SCA %*n_pix
    # Min_sca*n_pix
    # RF frequency (# maxima)
    #
      
    # cumulative
    
    if (i%%20==0){
      end = Sys.time() - strt 
      print(end)
    }
    
  }
  
  # add names (move to end)
  names(df_glaciers) = c(o_names,c("tot_mean","tot_sum","y_t","y_p","y_r2","mam_t","mam_p","mam_r2","jja_t",
                                   "jja_p","jja_r2","mamj_t","mamj_p","mamj_r2","jaso_t",
                                   "jaso_p","jaso_r2","aso_t","aso_p","aso_r2","amj_t",
                                   "amj_p","amj_r2","na_t","na_p","na_r2","nm_t","nm_p","nm_r2",
                                   "sa_yt","sa_yp","sa_yr2","sa_mamt","sa_mamp","sa_mamr2",
                                   "sa_amjt","sa_amjp","sa_amjr2","sa_jjat","sa_jjap","sa_jjar2"))
  #12 more for sca

  end = Sys.time() - strt 
  print(end)
  return(df_glaciers)
}

## DONE
##############################################################################################################################
year.season.trends <- function(vstk, vdates, velv){
  
  # INDEX YEARS, MONTHS, MAM, JJA
  yid = format(vdates,"%Y"); yidx = as.numeric(yid)-2000
  mid = format(vdates,"%m"); midx = as.numeric(mid)
  mami = midx == 3 | midx == 4 | midx == 5; jjai = midx == 6 | midx == 7 | midx == 8
  mamji = midx == 3 | midx == 4 | midx == 5 | midx == 6; jasoi = midx == 7 | midx == 8 | midx == 9 | midx == 10
  amji = midx == 4 | midx == 5 | midx == 6; asoi = midx == 8 | midx == 9 | midx == 10
  mam = yidx[mami] # March, April, May etc.
  jja = yidx[jjai]
  amj = yidx[amji]
  mamj = yidx[mamji]
  jaso = yidx[jasoi]
  aso = yidx[asoi]
  
  # variables
  # if(ncol(df_glaciers)!=9){cat("! Incorrect column index !\n> ncol != 9 ");stop()}
  
  # MAIN FUNCTIONS
  fun_time = 'mean'     # joining images over time
  fun_spatial = 'mean'  # over glacier area
  #
  # EXTRA FUNCTIONS
  fun_lm_time = function(x) { if (is.na(x[1])){ NA } else {lm(x ~ time)$coefficients[2] }} 
  fun_count_na = function(i, ...) sum(!is.na(i))
  fun_count = function(i, ...) sum((i) < mean(i))
  fun_count_mean = function(x) {
    comp <- x[1]
    main <- x[2:length(x)]
    length(main[main >= comp])
    #sum(main >= comp,na.rm=T)
  }
  
  # YEARLY/SEASONAL STACKS (ALL SAME NLAYERS)
  # 
  yearly_means_stk0 <- stackApply(vstk,indices=yidx, fun_time)
  yearly_means_stk <- cellStats(stackApply(vstk,indices=yidx, fun_time), fun_spatial)
  mam_means_stk <- cellStats(stackApply(vstk[[which(mami)]],indices=mam, fun_time), fun_spatial)
  jja_means_stk <- cellStats(stackApply(vstk[[which(jjai)]],indices=jja, fun_time), fun_spatial)
  mamj_means_stk <- cellStats(stackApply(vstk[[which(mamji)]],indices=mamj, fun_time), fun_spatial)
  jaso_means_stk <- cellStats(stackApply(vstk[[which(jasoi)]],indices=jaso, fun_time), fun_spatial)
  aso_means_stk <- cellStats(stackApply(vstk[[which(asoi)]],indices=aso, fun_time), fun_spatial)
  amj_means_stk <- cellStats(stackApply(vstk[[which(amji)]],indices=amj, fun_time), fun_spatial)
  #
  # odd
  #monthly_means_stk <- stackApply(vstk,indices=midx, mean) # 12 layers (months)
  #year_month_means_stk <- stackApply(vstk,indices=(yidx + midx/100), mean)
  # 
  # NCELL COUNT !NAs
  years_ncell_na <-  stackApply(vstk, indices = yidx, fun_count_na)
  years_ncell_na[years_ncell_na==0] <- NA
  years_ncell_na <- cellStats(years_ncell_na, fun_spatial)
  #
  # NCELLS ABOVE MEAN
  for (i in 1:nlayers(yearly_means_stk0)){
    vmeans = vstk[[which(yidx==i)]]
    #vvm <- calc(vmeans, mean)
    vvm <- yearly_means_stk0[[i]]
    vmeans2 = stack(vvm, vmeans)
    r = calc(vmeans2, fun_count_mean)
    # get rid of NA
    r[is.na(vmeans[[1]])] = NA
    if (i==1){vms = r}else{vms=stack(vms,r)}
  }
  ncells_peryear_over_mean = cellStats(vms,mean)
  
  ## ADD SCA (Should add 'type'=SCA/RF)
  # sca_weight_area <- (cellStats(vstk,mean)/100)*cellStats(vstk,fun_count_na)*(0.5^2)
  sca_weight_area <- cellStats(vstk,mean)
  sca_area_yr <- with(data.frame(value=sca_weight_area,factor=yidx), tapply(value, factor, mean,na.rm=T))
  sca_area_mam <- with(data.frame(value=sca_weight_area[which(mami)],factor=mam), tapply(value, factor, mean,na.rm=T))
  sca_area_amj <- with(data.frame(value=sca_weight_area[which(amji)],factor=amj), tapply(value, factor, mean,na.rm=T))
  sca_area_jja <- with(data.frame(value=sca_weight_area[which(jjai)],factor=jja), tapply(value, factor, mean,na.rm=T))
  
  # TO ADD
  # difference between max and min SCA every year
  # day of min SCA
  
  # TREND VALS
  # create trends and summaries
  ys <- 1:nlayers(yearly_means_stk0)
  vals_list = list()
  #
  fit <- lm(yearly_means_stk~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(mam_means_stk~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(jja_means_stk~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(mamj_means_stk~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(jaso_means_stk~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(aso_means_stk~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(amj_means_stk~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(years_ncell_na~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(ncells_peryear_over_mean~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  # NEW FOR SCA 12
  #
  fit <- lm(sca_area_yr~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(sca_area_mam~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(sca_area_amj~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  #
  fit <- lm(sca_area_jja~ys)
  vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$coefficients[2,4],summary(fit)$r.squared)
  

   return(unlist(vals_list))
  
}



##############################################################################################################################
############# CORRELLATIONS FUNCTIONS ##################

## DONE
##############################################################################################################################
g.match.rf.sca <- function(glacier_id, rgi, demL, noZero=TRUE, threshold = NULL,
                                 rfpath, scapath, tfpath){
  #
    #                                                         #
    # FOR GIVEN GLACIER                                       #
    # create matching SCA and RF data for glacier_i           #
    # mask with inner polygon                                 #
    #                                                         #
    # Returns:                                                #
    # rsub | ssub | elv | match_dates                       #
    #                                                         #
  #
  # READ RF FILE
  rf_match <- list.files(rfpath, full.names = T, pattern = paste0(glacier_id,".*grd"))
  rfstk = stack(rf_match)
  #
  # RGI SHAPEFILE MASK (inner)
  shape <<- rgi[gsub("[[:punct:]]", "_", rgi@data$RGIId) == glacier_id,]
  cent_distance <- res(rfstk)[1]/2
  suppressWarnings(inner_shape <<- buffer(shape,width=-cent_distance))
  #
  # DEM
  dem <- crop(demL,rfstk)
  elv <<- mask(resample(dem, rfstk), inner_shape)
  #
  # find matching SCA for glacier
  sca_match = list.files(scapath, full.names = T, pattern = paste0(glacier_id,".*grd"))
  if (!is.na(sca_match)){scastk = try(stack(sca_match))} else{print("!!! SCA_glacier_match does not exist!")}
  #
  # mask for all years
  rf_shp = mask(rfstk, inner_shape)
  sca_shp = mask(scastk, inner_shape)
  #
  # DELETE ALL 0vals
  if(noZero){rf_shp[rf_shp==0] <- NA;sca_shp[sca_shp==0] <- NA}
  #
  # OVERLAP SCA & RF DATES
  rf_n = substring(names(rf_shp),2,8)
  sca_n = substring(names(sca_shp),2,8)
  #
  # MATCH RF & SCA (Keep all RF)
  # EXCLUDE 2018 DATES (& 2017)
  mixd = intersect(rf_n, sca_n)
  match_dates <- as.Date(mixd,format="%Y%j")
  didx <- format(match_dates, "%Y") != "2018" & format(match_dates, "%Y") != "2017"
  match_dates <<- match_dates[didx]
  mixd <- mixd[didx]
  #
  # SUBSET
  rsub <<- rf_shp[[match(mixd,rf_n)]]
  ssub <<- sca_shp[[match(mixd,sca_n)]]
  #
  if (!is.null(threshold)){
    # add threshold
    if(length(threshold)>1){
      rsub[ssub < threshold[1]] <- NA
      ssub[ssub < threshold[2]] <- NA
      rsub <<- rsub; ssub <<- ssub
      cat("Threshold of >", threshold[1], "SCA to RF \n ","Threshold of >", threshold[2], "SCA to SCA")
    } else{
      rsub[ssub < threshold] <- NA
      ssub[ssub < threshold] <- NA
      rsub <<- rsub; ssub <<- ssub
      print(paste("Threshold of >", threshold, "to all"))
    }
    
  } 
  # 
  # cat("Created variables: \n  rsub | ssub | elv | match_dates ")  
  # print(paste("<glacier_id>",glacier_id))
}

## NEW 06/23/2020
##############################################################################################################################
g.match.rf.sca2 <- function(glacier_id, rgi, demL, noZero=TRUE, threshold = NULL,
                           rfpath, scapath, tfpath){
  #
  #                                                         #
  # FOR GIVEN GLACIER                                       #
  # create matching SCA and RF data for glacier_i           #
  # mask with inner polygon                                 #
  #                                                         #
  # Returns:                                                #
  # rsub | ssub | elv | match_dates                       #
  #                                                         #
  #
  # READ RF FILE
  rf_match <- list.files(rfpath, full.names = T, pattern = paste0(glacier_id,".*grd"))
  rfstk = stack(rf_match)
  #
  # RGI SHAPEFILE MASK (inner)
  shape <<- rgi[gsub("[[:punct:]]", "_", rgi@data$RGIId) == glacier_id,]
  cent_distance <- res(rfstk)[1]/2
  suppressWarnings(inner_shape <<- buffer(shape,width=-cent_distance))
  #
  # DEM
  dem <- crop(demL,rfstk)
  elv <<- mask(resample(dem, rfstk), inner_shape)
  #
  # find matching SCA for glacier
  sca_match = list.files(scapath, full.names = T, pattern = paste0(glacier_id,".*grd"))
  if (!is.na(sca_match)){scastk = try(stack(sca_match))} else{print("!!! SCA_glacier_match does not exist!")}
  #
  # mask for all years
  rf_shp = mask(rfstk, inner_shape)
  sca_shp = mask(scastk, inner_shape)
  #
  # DELETE ALL 0vals
  # if(noZero){rf_shp[rf_shp==0] <- NA;sca_shp[sca_shp==0] <- NA}
  #
  # OVERLAP SCA & RF DATES
  rf_n = substring(names(rf_shp),2,8)
  sca_n = substring(names(sca_shp),2,8)
  #
  # MATCH RF & SCA (Keep all RF)
  # EXCLUDE 2018 DATES (& 2017)
  mixd = intersect(rf_n, sca_n)
  match_dates <- as.Date(mixd,format="%Y%j")
  didx <- format(match_dates, "%Y") != "2018" & format(match_dates, "%Y") != "2017"
  match_dates <<- match_dates[didx]
  mixd <- mixd[didx]
  #
  # SUBSET
  rsub <<- rf_shp[[match(mixd,rf_n)]]
  ssub <<- sca_shp[[match(mixd,sca_n)]]
  #
  if (!is.null(threshold)){
    # add threshold
    if(length(threshold)>1){
      rsub[ssub < threshold[1]] <- NA
      ssub[ssub < threshold[2]] <- NA
      cat("Threshold of >", threshold[1], "SCA to RF \n ","Threshold of >", threshold[2], "SCA to SCA")
    } else{
      rsub[ssub < threshold] <- NA
      ssub[ssub < threshold] <- NA
      print(paste("Threshold of >", threshold, "to all"))
    }
    
  } 
  # MAKE NA VALUES WITHIN glacier 0 (NEW NEW)
  if(noZero){
    rsub[is.na(rsub)] <- 0;rsub <- mask(rsub, inner_shape)
    ssub[is.na(ssub)] <- 0;ssub <- mask(ssub, inner_shape)
    rsub <<- rsub; ssub <<- ssub
    print("...applied 0-mask")
    
  }
  #
  # cat("Created variables: \n  rsub | ssub | elv | match_dates ")  
  # print(paste("<glacier_id>",glacier_id))
}

## DONE
#####################################################################################
glacier.corr.year <- function(df_glaciers, rgi,demL,noZero = TRUE,threshold = NULL,
                              rfpath,scapath,tfpath){
  # 
  # Data frame of cross correlation between RF and SCA 
  # 
  strt = Sys.time(); ci=ncol(df_glaciers)
  o_names <- names(df_glaciers)
  df_glaciers[,10:68] <- NA # make room for new variables
  ## LOOP THROUGH EACH GLACIER
  # !!!!
  # quick fix (last name has no SCA)
  for (i in 1:(length(df_glaciers$Object.ID)-1)){
    #
    g.match.rf.sca(df_glaciers$Object.ID[i],rgi,demL,noZero = TRUE,threshold = threshold,
                   rfpath,scapath,tfpath)
    # ...creates: rsub | ssub | elv | match_dates
    #
    # CROSS-CORRELATION
    rf_vals <- cellStats(rsub,mean)
    sca_vals <- cellStats(ssub,mean)
    continuous_weeks <- as.Date(seq(match_dates[1], tail(match_dates,1), by="7 days"))
    df_corr <- data.frame(time = continuous_weeks)
    df_corr$sca <- df_corr$rf <- NA
    match_weeks_idx <- match(match_dates,continuous_weeks)
    df_corr$rf[match_weeks_idx] <- rf_vals; df_corr$rf <- fill.na.mean(df_corr$rf)
    df_corr$sca[match_weeks_idx] <- sca_vals; df_corr$sca <- fill.na.mean(df_corr$sca)
    #
    ## CORRELATION
    ccvals <- ccf(df_corr$rf,df_corr$sca, plot = FALSE)
    cloess <- loess(ccvals$acf~ccvals$lag, span=0.2)
    lines(cloess$fitted~cloess$x, col='red',lwd=2)
    cf <- min(cloess$fitted)
    lag <- ccvals$lag[which.min(cloess$fitted)]
    #
    # split by years
    yid = format(match_dates,"%Y"); yidx = as.numeric(yid)-2000
    #
    ## AREA
    sca_weight_area <- (cellStats(ssub,mean)/100)*cellStats(ssub,fun_count_na)*(0.5^2)
    rf_weight_area <- (cellStats(rsub,mean)/100)*cellStats(rsub,fun_count_na)*(0.5^2)
    #
    sca_mean_yr <- with(data.frame(value=sca_weight_area,factor=yidx), tapply(value, factor, mean,na.rm=T))
    sca_sd_yr <- with(data.frame(value=sca_weight_area,factor=yidx), tapply(value, factor, sd,na.rm=T))
    sca_min_yr <- with(data.frame(value=sca_weight_area,factor=yidx), tapply(value, factor, min,na.rm=T))
    rf_mean_yr <- with(data.frame(value=rf_weight_area,factor=yidx), tapply(value, factor, mean,na.rm=T))
    rf_max_yr <- with(data.frame(value=rf_weight_area,factor=yidx), tapply(value, factor, max,na.rm=T))
    #
    rsr = cellStats(rsub,mean)*(cellStats(ssub,mean)/100)
    rf_sca_frac <- with(data.frame(value=rsr,factor=yidx), tapply(value, factor, mean,na.rm=T))
    rf_sum_yr <- with(data.frame(value=cellStats(rsub,mean),factor=yidx), tapply(value, factor, sum,na.rm=T))
    #
    # DAY of MIN/MAX
    sca_min_d <- with(data.frame(value=sca_weight_area,factor=yidx), tapply(value, factor, which.min))
    sca_min_day <- as.numeric(sapply(1:17,FUN = function(x) format(match_dates[yidx==x][sca_min_d[x]],"%j")))
    rf_max_d <- with(data.frame(value=rf_weight_area,factor=yidx), tapply(value, factor, which.max))
    rf_max_day <- as.numeric(sapply(1:17,FUN = function(x) format(match_dates[yidx==x][rf_max_d[x]],"%j")))
    #
    # days over total mean
    tm <- mean(cellStats(rsub,mean),na.rm=T)
    rf_ndays_overmean <- with(data.frame(value=cellStats(rsub,mean),factor=yidx), tapply(value, factor, function(x) sum(x > tm,na.rm=T)))
    tm <- mean(cellStats(ssub,mean),na.rm=T)
    sca_ndays_undermean <- with(data.frame(value=cellStats(ssub,mean),factor=yidx), tapply(value, factor, function(x) sum(x < tm,na.rm=T)))
    #
    # maxima?
    rf_maxima <- rep(1,17)
    # rf_maxima <- with(data.frame(value=cellStats(rsub,mean),factor=yidx), tapply(value, factor, function(x) length(localMaxima(x))))
    # 
    # QUICK TRENDS
    vals_list = list()
    ys = 1:17
    fit <- lm(sca_min_yr~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(sca_min_day~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(rf_mean_yr~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(rf_max_yr~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(rf_max_day~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(rf_sum_yr~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(rf_ndays_overmean~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(sca_ndays_undermean~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(rf_maxima~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    fit <- lm(rf_sca_frac~ys)
    vals_list = c(vals_list, coef(fit)['ys'],summary(fit)$r.squared)
    # 20 extra trend vals
    # 
    # add values to data frame
    df_glaciers[i,(ci+1)] <- mean(getValues(rsub),na.rm=T)
    df_glaciers[i,(ci+2)] <- sd(getValues(rsub),na.rm=T) 
    df_glaciers[i,(ci+3)] <- mean(getValues(ssub),na.rm=T) 
    df_glaciers[i,(ci+4)] <- sd(getValues(ssub),na.rm=T)
    df_glaciers[i,(ci+5)] <- cf
    df_glaciers[i,(ci+6)] <- lag
    df_glaciers[i,(ci+7)] <- max(cellStats(rsub,fun_count_na),na.rm=T) # max|mean n cell !na
    df_glaciers[i,(ci+8)] <- mean(cellStats(rsub,fun_count_na),na.rm=T)
    df_glaciers[i,(ci+9)] <- sd(cellStats(rsub,fun_count_na),na.rm=T)
    df_glaciers[i,(ci+10)] <- min(cellStats(ssub,fun_count_na),na.rm=T) 
    df_glaciers[i,(ci+11)] <- mean(cellStats(ssub,fun_count_na),na.rm=T)
    df_glaciers[i,(ci+12)] <- sd(cellStats(ssub,fun_count_na),na.rm=T)
    df_glaciers[i,(ci+13)] <- max(rf_mean_yr,na.rm=T) # RF max|mean area
    df_glaciers[i,(ci+14)] <- mean(rf_mean_yr,na.rm=T) # the yr mean of means
    df_glaciers[i,(ci+15)] <- sd(rf_mean_yr,na.rm=T)
    df_glaciers[i,(ci+16)] <- max(rf_max_yr,na.rm=T)
    df_glaciers[i,(ci+17)] <- mean(rf_max_yr,na.rm=T)
    df_glaciers[i,(ci+18)] <- sd(rf_max_yr,na.rm=T)
    df_glaciers[i,(ci+19)] <- min(sca_mean_yr,na.rm=T) # SCA min|mean area
    df_glaciers[i,(ci+20)] <- mean(sca_mean_yr,na.rm=T)
    df_glaciers[i,(ci+21)] <- sd(sca_mean_yr,na.rm=T)
    df_glaciers[i,(ci+22)] <- min(sca_min_yr,na.rm=T)
    df_glaciers[i,(ci+23)] <- mean(sca_min_yr,na.rm=T)
    df_glaciers[i,(ci+24)] <- sd(sca_min_yr,na.rm=T)
    df_glaciers[i,(ci+25)] <- mean(rf_max_day,na.rm=T) # RF day of max
    df_glaciers[i,(ci+26)] <- sd(rf_max_day,na.rm=T)
    df_glaciers[i,(ci+27)] <- sort(rf_max_day)[1] # earliest/latest
    df_glaciers[i,(ci+28)] <- tail(sort(rf_max_day),1)
    df_glaciers[i,(ci+29)] <- mean(sca_min_day,na.rm=T) # SCA day of min
    df_glaciers[i,(ci+30)] <- sd(sca_min_day,na.rm=T)
    df_glaciers[i,(ci+31)] <- sort(sca_min_day)[1] # earliest/latest
    df_glaciers[i,(ci+32)] <- tail(sort(sca_min_day),1)
    df_glaciers[i,(ci+33)] <- mean(rf_ndays_overmean,na.rm=T)
    df_glaciers[i,(ci+34)] <- sd(rf_ndays_overmean,na.rm=T)
    df_glaciers[i,(ci+35)] <- mean(sca_ndays_undermean,na.rm=T)
    df_glaciers[i,(ci+36)] <- sd(sca_ndays_undermean,na.rm=T)
    df_glaciers[i,(ci+37)] <- mean(rf_maxima,na.rm=T)
    df_glaciers[i,(ci+38)] <- mean(rf_sca_frac,na.rm=T)
    df_glaciers[i,(ci+39)] <- sd(rf_sca_frac,na.rm=T)
    # trends
    st= ci+40; en=st+19
    df_glaciers[i,(st:en)] <- unlist(vals_list)
    #
    print(i)
    #
    if (i%%20==0){
      end = Sys.time() - strt 
      print(end)
    }
    
  }
  # add names (move to end)
  names(df_glaciers) = c(o_names,c("rfmean","rfsd","scamean","scasd","lcorr","lag","n_rfmax","n_rfmean","n_rfsd",
                                   "n_scamax","n_scamean","n_scasd","rfarea_meanmax","rfarea_meanmean","rfarea_meansd",
                                   "rfarea_maxmean","rfarea_maxmax","rfarea_maxsd","scarea_meanmin",
                                   "scarea_meanmean","scarea_meansd","scarea_minmin","scarea_minmean","scarea_minsd",
                                   "rfmax_day_mean","rfmax_day_sd","rfmax_day_first","rfmax_day_last","scamin_day_mean",
                                   "scamin_day_sd","scamin_day_first","scamin_day_last",
                                   "rfover_mean","rfover_sd","sunder_mean","sunder_sd","rf_maxima","rfsca_f_mean","rfsca_f_sd",
                                   "tr_sarea_min","tr2_sarea_min","tr_sday_min","tr2_sday_min","tr_rarea_mean","tr2_rarea_mean",
                                   "tr_rarea_max","tr2_rarea_max","tr_rday_max","tr2_rday_max","tr_rsum","tr2_rsum",
                                   "tr_rover","tr2_rover","tr_sunder","tr2_sunder","tr_rmaxim","tr2_rmaxim",
                                   "tr_rfsca_f","tr2_rfsca_f"))

  end = Sys.time() - strt 
  print(end)
  return(df_glaciers)
}

## NEW 04/15/2020
#####################################################################################
glacier.month.trends <- function(df_glaciers, rgi,demL,noZero = TRUE,threshold = NULL,
                              rfpath,scapath,tfpath, t.func='mean'){
  # 
  # Data frame of cross correlation between RF and SCA 
  # 
  strt = Sys.time(); ci=ncol(df_glaciers)
  o_names <- names(df_glaciers)
  df_glaciers[,10:406] <- NA # make room for new variables
  ## LOOP THROUGH EACH GLACIER
  # !!!!
  # quick fix (last name has no SCA)
  for (i in 1:(length(df_glaciers$Object.ID)-1)){
    #
    g.match.rf.sca2(df_glaciers$Object.ID[i],rgi,demL,noZero = TRUE,threshold = threshold,
                   rfpath,scapath,tfpath)
    # ...creates: rsub | ssub | elv | match_dates
    
    #
    # Main Values
    rf_vals <- cellStats(rsub,mean)
    sca_vals <- cellStats(ssub,mean)
    #
    # split indices
    yid = format(match_dates,"%Y"); yidx = as.numeric(yid)-2000
    mid = format(match_dates,"%m"); midx = as.numeric(mid)
    ymidx = (yidx + midx/100)
    ym_sca <- stackApply(ssub,indices=ymidx, mean)
    ym_rf <- stackApply(rsub,indices=ymidx, mean)
    ymidx_u <- unique(ymidx)
    #
    ## YEAR_VALUES ------------->
    sca_mean_yr0 <- with(data.frame(value=sca_vals,factor=yidx), tapply(value, factor, mean,na.rm=T))
    sca_mean_yr = mean(sca_mean_yr0,na.rm=T)
    sca_sd_yr0 <- with(data.frame(value=sca_vals,factor=yidx), tapply(value, factor, sd,na.rm=T))
    sca_sd_yr = mean(sca_sd_yr0,na.rm=T)
    sca_min_yr0 <- with(data.frame(value=sca_vals,factor=yidx), tapply(value, factor, min,na.rm=T))
    sca_min_yr = mean(sca_min_yr0,na.rm=T)
    rf_mean_yr0 <- with(data.frame(value=rf_vals,factor=yidx), tapply(value, factor, mean,na.rm=T))
    rf_mean_yr = mean(rf_mean_yr0,na.rm=T)
    rf_sum_yr0 <- with(data.frame(value=rf_vals,factor=yidx), tapply(value, factor, sum,na.rm=T))
    rf_sum_yr = mean(rf_sum_yr0,na.rm=T)
    rf_sd_yr0 <- with(data.frame(value=rf_vals,factor=yidx), tapply(value, factor, sd,na.rm=T))
    rf_sd_yr = mean(rf_sd_yr0,na.rm=T)
    rf_max_yr0 <- with(data.frame(value=rf_vals,factor=yidx), tapply(value, factor, max,na.rm=T))
    rf_max_yr = mean(rf_max_yr0,na.rm=T)
    #
    ts = 1:16 # replace all (take unique years...)
    #
    ## YEAR TRENDS
    dfyv = data.frame(t = 1:16,sca_mean_yr0,sca_sd_yr0,sca_min_yr0,
                      rf_mean_yr0,rf_sum_yr0,rf_sd_yr0,rf_max_yr0)
    # run trends
    dtr1 = dfyv %>% melt(.,id.vars = 't') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~t, data=.)) %>% 
      tidy(dfit) %>% filter(term=="t")
    dtr2 = dfyv %>% melt(.,id.vars = 't') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~t, data=.)) %>% glance(dfit)
    # save values
    tr_est_yr = dtr1$estimate
    tr_pval_yr = dtr1$p.value
    tr_ste_yr = dtr1$std.error
    tr_sig_yr = dtr2$sigma
    tr_r2_yr = dtr2$r.squared
    #
    ## YEAR SCA~RF
    cor_scarf = dfyv %>%  
      dplyr::summarize(cor_mean = round(cor(sca_mean_yr0,rf_mean_yr0),2),
                       cor_smean_rsum = round(cor(sca_mean_yr0,rf_sum_yr0),2),
                       cor_smin_rsum = round(cor(sca_min_yr0,rf_sum_yr0),2),
                       cor_smin_rmean = round(cor(sca_min_yr0,rf_mean_yr0),2),
                       cor_sd = round(cor(sca_sd_yr0,rf_sd_yr0),2))
    # **if cor is significant then rerun with more lm
    scarf_mean1 = dfyv %>% do(dfit = lm(sca_mean_yr0~rf_mean_yr0, data=.)) %>% 
      tidy(dfit) %>% filter(term=="rf_mean_yr0") %>% select(-term,-statistic)
    scarf_mean2 = dfyv %>% do(dfit = lm(sca_mean_yr0~rf_mean_yr0, data=.)) %>% glance(dfit)
    tr_scarf_mean = c(as.numeric(scarf_mean1),scarf_mean2$r.squared)
    # COUNT npix
    npix = max(cellStats(ssub,fun_count_na,na.rm=T))
    snpix = cellStats(ssub,fun_count_na,na.rm=T); snpix[snpix < (npix*0.25)] = npix # if less than 1/4 of glacier
    ss50 = ssub; ss50[ss50 < 50] = NA # exclude pixels less than X
    snpix50 = cellStats(ss50,fun_count_na,na.rm=T);snpix50[snpix50 < (npix*0.25)] = npix
    spixmin <- with(data.frame(value=snpix,factor=yidx), tapply(value, factor, min,na.rm=T))
    spixmin50 <- with(data.frame(value=snpix50,factor=yidx), tapply(value, factor, min,na.rm=T))
    # trend for pixels
    ppx = data.frame(y=1:16,spixmin) %>%
      do(dfit = lm(spixmin~y, data=.)) %>% 
      tidy(dfit) %>% filter(term=="y") %>% select(-term,-statistic,-std.error)
    ppx50 = data.frame(y=1:16,spixmin50) %>%
      do(dfit = lm(spixmin50~y, data=.)) %>% 
      tidy(dfit) %>% filter(term=="y")  %>% select(-term,-statistic,-std.error)
    # save values
    spx_rat = (mean(spixmin,na.rm=T)/npix)
    spx_rat50 = (mean(spixmin50,na.rm=T)/npix)
    pix_vals = c(npix,spx_rat,spx_rat50)
    tr_pix = c(as.numeric(ppx),as.numeric(ppx50))
    #
    ## MONTH DFs ------------->
    sca_month <- with(data.frame(value=sca_vals,factor=ymidx), tapply(value, factor, t.func,na.rm=T))
    dfmos = cbind(y=1:16,as.data.frame(matrix(sca_month,ncol=12,byrow=TRUE)))
    rf_month <- with(data.frame(value=rf_vals,factor=ymidx), tapply(value, factor, t.func,na.rm=T))
    dfmor = cbind(y=1:16,as.data.frame(matrix(rf_month,ncol=12,byrow=TRUE)))
    #
    # MAIN MONTH VALUES
    sca_mean_month <- colMeans(dfmos)[2:13]
    sca_sd_month <- apply(dfmos, 2, sd)[2:13]
    sca_min_month <- apply(dfmos, 2, min)[2:13]
    rf_mean_month <- colMeans(dfmor)[2:13]
    rf_sum_month <- colSums(dfmor)[2:13]
    rf_sd_month <- apply(dfmor, 2, sd)[2:13]
    rf_max_month <- apply(dfmor, 2, max)[2:13]
    # 1 value
    sca_whmin_month <- as.numeric(which.min(sca_min_month))
    rf_whmax_month <- as.numeric(which.max(rf_mean_month))
    rf_whmaxsum_month <- as.numeric(which.max(rf_sum_month))
    #
    ## MONTH TRENDS
    # SCA
    dm1 = dfmos %>% melt(.,id.vars = 'y') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~y, data=.)) %>% 
      tidy(dfit) %>% filter(term=="y") %>% select(-term,-statistic)
    dm2 = dfmos %>% melt(.,id.vars = 'y') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~y, data=.)) %>% glance(dfit)
    # save values
    trsca_est_month = dm1$estimate
    trsca_pval_month = dm1$p.value
    trsca_ste_month = dm1$std.error
    trsca_sig_month = dm2$sigma
    trsca_r2_month = dm2$r.squared
    # SCA
    dm1 = dfmor %>% melt(.,id.vars = 'y') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~y, data=.)) %>% 
      tidy(dfit) %>% filter(term=="y") %>% select(-term,-statistic)
    dm2 = dfmor %>% melt(.,id.vars = 'y') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~y, data=.)) %>% glance(dfit)
    # save values (12 values)
    trf_est_month = dm1$estimate
    trf_pval_month = dm1$p.value
    trf_ste_month = dm1$std.error
    trf_sig_month = dm2$sigma
    trf_r2_month = dm2$r.squared
  
    ## YEAR SCA~RF
    mrf = dfmor %>%
      melt(id.vars="y")
    srf = dfmos %>%
      melt(id.vars="y")
    cor_month = cbind(srf, rvalue = mrf$value) %>% na.omit() %>% 
      group_by(variable) %>% 
      dplyr::summarize(corr = round(cor(rvalue,value),2)) %>% select(corr)
    dmon1 = cbind(srf, rvalue = mrf$value) %>%
      group_by(variable) %>%
      do(dfit = lm(value~rvalue, data=.)) %>%
      tidy(dfit) %>% filter(term=="rvalue") %>% select(-term,-statistic)
    dmon2 = cbind(srf, rvalue = mrf$value) %>%
      group_by(variable) %>%
      do(dfit = lm(value~rvalue, data=.)) %>% glance(dfit)
    # save values
    tr_scarf_est_month = dmon1$estimate
    tr_scarf_pval_month = dmon1$p.value
    tr_scarf_ste_month = dmon1$std.error
    tr_scarf_sig_month = dmon2$sigma
    tr_scarf_r2_month = dmon2$r.squared
    
    # COUNT npix (CHANGE BY MONTH!!!!!!)
    snpixm = cellStats(ssub,fun_count_na,na.rm=T) 
    spixmonth <- with(data.frame(value=snpixm,factor=midx), tapply(value, factor, min,na.rm=T))
    snpixm[snpixm < (npix*0.25)] = npix # if less than 1/4 of glacier
    spixmonth_c <- with(data.frame(value=snpixm,factor=midx), tapply(value, factor, min,na.rm=T))
    snpix50m = cellStats(ss50,fun_count_na,na.rm=T)
    spixmonth50 <- with(data.frame(value=snpix50m,factor=midx), tapply(value, factor, min,na.rm=T))
    snpix50m[snpix50m < (npix*0.25)] = npix
    spixmonth50_c <- with(data.frame(value=snpix50m,factor=midx), tapply(value, factor, min,na.rm=T))
    
    spix_allmonth <- with(data.frame(value=snpixm,factor=ymidx), tapply(value, factor, min,na.rm=T))
    dfspix = cbind(y=1:16,as.data.frame(matrix(spix_allmonth,ncol=12,byrow=TRUE)))
    
    # trend for pixels
    ppxm = dfspix %>% melt(.,id.vars="y") %>% group_by(variable) %>% 
      do(dfit = lm(value~y, data=.)) %>%
      tidy(dfit) %>% filter(term=="y") %>% select(-term,-statistic,-std.error)
    ppx50m = dfspix %>% melt(.,id.vars="y") %>% group_by(variable) %>% 
      do(dfit = lm(value~y, data=.)) %>%
      tidy(dfit) %>% filter(term=="y")  %>% select(-term,-statistic,-std.error)
    # save values
    spxmo_rat = spixmonth/npix
    spxmo_rat_c = spixmonth_c/npix
    spxmo_rat50_c = spixmonth50_c/npix
    tr_pix_month_est = ppxm$estimate
    tr_pix_month_pval = ppxm$p.value
    # tr_pix50_month_est = ppxm$estimate
    

    
    # can be saved as a list of values???
    # SAVE VALUES
    var_yr_mo_ls = list(sca_mean_yr,sca_sd_yr,sca_min_yr,  ## YEAR VALUES!!!
      rf_mean_yr,rf_sum_yr,rf_sd_yr,rf_max_yr, # first 7 are values
      tr_est_yr, # 7 trends (smean,ssd,smin,rmean,rsum,rsd,rmax)
      tr_pval_yr,tr_ste_yr,tr_sig_yr,tr_r2_yr, # (all 7 values each)
      cor_scarf, # 5 values (means, smean_rsum, smin_rsum, smin_rmean, sd)
      tr_scarf_mean, # 4 (estimate, std.error, p.value, r2)
      npix,spx_rat,spx_rat50, #  3 (npix, ave_minpix, ave_minpix50)
      tr_pix,   # 4 (tr_minpix_estimate, tr_minpix_p.value, tr_minpix50_estimate, tr_minpix50_p.value)
      sca_mean_month,  # MONTH VALUES!!! (12 values)
      sca_sd_month,
      sca_min_month,
      rf_mean_month,
      rf_sum_month,
      rf_sd_month,
      rf_max_month,
      sca_whmin_month, # 1 variable
      rf_whmax_month, # 1 variable
      rf_whmaxsum_month, # 1 variable
      trsca_est_month, # MONTHLY TRENDS (12 Vals)
      trsca_pval_month,
      trsca_ste_month,
      trsca_sig_month,
      trsca_r2_month,
      trf_est_month,
      trf_pval_month,
      trf_ste_month,
      trf_sig_month,
      trf_r2_month,
      cor_month, # SCARF by MONTH (12 values)
      tr_scarf_est_month,
      tr_scarf_pval_month,
      tr_scarf_ste_month,
      tr_scarf_sig_month,
      tr_scarf_r2_month,
      spxmo_rat,spxmo_rat_c,spxmo_rat50_c, 
      tr_pix_month_est,
      tr_pix_month_pval
    )
    # add pixels by month
    nm_h1 = c("smean","ssd","smin","rmean","rsum","rsd","rmax")
    
    # NAMES
    nms = c("sca_mean_yr","sca_sd_yr","sca_min_yr",  ## YEAR VALUES!!!
         "rf_mean_yr","rf_sum_yr","rf_sd_yr","rf_max_yr", # 1 each
         paste0(rep("tr_est_yr",7),nm_h1), # 7 trends (smean,ssd,smin,rmean,rsum,rsd,rmax)
         paste0(rep("tr_pval_yr",7),nm_h1),
         paste0(rep("tr_ste_yr",7),nm_h1),
         paste0(rep("tr_sig_yr",7),nm_h1),
         paste0(rep("tr_r2_yr",7),nm_h1),
         paste0(rep("cor_scarf",5),c('means', 'smean_rsum', 'smin_rsum', 'smin_rmean', 'sd')),
         paste0(rep("tr_scarf_mean",4),c('est', 'ste', 'pval', 'r2')),
         'npix','spx_rat','spx_rat50', #  3 (npix, ave_minpix, ave_minpix50)
         paste0(rep("tr_pix",4),c('est', 'pval', 'est50', 'pval50')),
         paste0(rep("sca_mean_month",12),1:12), # MONTH VALUES!!! (12 values)
         paste0(rep("sca_sd_month",12),1:12),
         paste0(rep("sca_min_month",12),1:12),
         paste0(rep("rf_mean_month",12),1:12),
         paste0(rep("rf_sum_month",12),1:12),
         paste0(rep("rf_sd_month",12),1:12),
         paste0(rep("rf_max_month",12),1:12),
         "sca_whmin_month", # 1 variable
         "rf_whmax_month", # 1 variable
         "rf_whmaxsum_month", # 1 variable
         paste0(rep("trsca_est_month",12),1:12), # MONTHLY TRENDS (12 Vals)
         paste0(rep("trsca_pval_month",12),1:12),
         paste0(rep("trsca_ste_month",12),1:12),
         paste0(rep("trsca_sig_month",12),1:12),
         paste0(rep("trsca_r2_month",12),1:12),
         paste0(rep("trf_est_month",12),1:12),
         paste0(rep("trf_pval_month",12),1:12),
         paste0(rep("trf_ste_month",12),1:12),
         paste0(rep("trf_sig_month",12),1:12),
         paste0(rep("trf_r2_month",12),1:12),
         paste0(rep("cor_month",12),1:12), # SCARF by MONTH (12 values)
         paste0(rep("tr_scarf_est_month",12),1:12),
         paste0(rep("tr_scarf_pval_month",12),1:12),
         paste0(rep("tr_scarf_ste_month",12),1:12),
         paste0(rep("tr_scarf_sig_month",12),1:12),
         paste0(rep("tr_scarf_r2_month",12),1:12),
         paste0(rep("spxmo_rat",12),1:12),
         paste0(rep("spxmo_rat_c",12),1:12),
         paste0(rep("spxmo_rat50_c",12),1:12),
         paste0(rep("tr_pix_month_est",12),1:12),
         paste0(rep("tr_pix_month_pval",12),1:12)
    )
    
    # print(paste("varnames = ",length(unlist(var_yr_mo_ls)) == length(nms)))
    
    st = (ci+1); en = st+(length(nms)-1) # 10:406
    df_glaciers[i,(st:en)] <- unlist(var_yr_mo_ls)
    #
    print(i)
    #
    if (i%%20==0){
      end = Sys.time() - strt 
      print(end)
    }
    
   
  }
  # add names
  names(df_glaciers) = c(o_names,nms)
  
  end = Sys.time() - strt 
  print(end)
  return(df_glaciers)
}


## NEW 04/23/2020 - old method for temporal average (pixel-based)
#####################################################################################
glacier.month.trends.pix <- function(df_glaciers, rgi,demL,noZero = TRUE,threshold = NULL,
                                 rfpath,scapath,tfpath, t.func='mean',t.func_r='mean'){
  # 
  # Data frame of cross correlation between RF and SCA 
  # 
  strt = Sys.time(); ci=ncol(df_glaciers)
  o_names <- names(df_glaciers)
  df_glaciers[,10:406] <- NA # make room for new variables
  ## LOOP THROUGH EACH GLACIER
  # !!!!
  # quick fix (last name has no SCA)
  for (i in 1:(length(df_glaciers$Object.ID)-1)){
    #
    g.match.rf.sca2(df_glaciers$Object.ID[i],rgi,demL,noZero = TRUE,threshold = threshold,
                   rfpath,scapath,tfpath)
    # ...creates: rsub | ssub | elv | match_dates
    
    #### !!!! FIND WAY TO CHANGE DATE RANGE @@@@ g.match function...
    
    # !!! = change vstk/vdates to other variable.names(
      
  
    # index values 
    yid = format(match_dates,"%Y"); yidx = as.numeric(yid)-2000
    mid = format(match_dates,"%m"); midx = as.numeric(mid)
    ymidx = (yidx + midx/100)
    # function type (FOR MONTHS)
    func_time_sca = t.func    # "min"  # "mean"
    func_time_rf =  t.func_r  # "max"  # "mean"
    # monthly indices
    m1i = midx == 1; m1 = yidx[m1i] # March
    m2i = midx == 2; m2 = yidx[m2i] # March
    m3i = midx == 3; m3 = yidx[m3i] # March
    m4i = midx == 4; m4 = yidx[m4i] # March
    m5i = midx == 5; m5 = yidx[m5i] # March
    m6i = midx == 6; m6 = yidx[m6i] # March
    m7i = midx == 7; m7 = yidx[m7i] # March
    m8i = midx == 8; m8 = yidx[m8i] # March
    m9i = midx == 9; m9 = yidx[m9i] # March
    m10i = midx == 10; m10 = yidx[m10i] # March
    m11i = midx == 11; m11 = yidx[m11i] # March
    m12i = midx == 12; m12 = yidx[m12i] # March
    # SCA MONTH
    m1_means0 <- cellStats(stackApply(ssub[[which(m1i)]],indices=m1, func_time_sca), mean)
    m2_means0 <- cellStats(stackApply(ssub[[which(m2i)]],indices=m2, func_time_sca), mean)
    m3_means0 <- cellStats(stackApply(ssub[[which(m3i)]],indices=m3, func_time_sca), mean)
    m4_means0 <- cellStats(stackApply(ssub[[which(m4i)]],indices=m4, func_time_sca), mean)
    m5_means0 <- cellStats(stackApply(ssub[[which(m5i)]],indices=m5, func_time_sca), mean)
    m6_means0 <- cellStats(stackApply(ssub[[which(m6i)]],indices=m6, func_time_sca), mean)
    m7_means0 <- cellStats(stackApply(ssub[[which(m7i)]],indices=m7, func_time_sca), mean)
    m8_means0 <- cellStats(stackApply(ssub[[which(m8i)]],indices=m8, func_time_sca), mean)
    m9_means0 <- cellStats(stackApply(ssub[[which(m9i)]],indices=m9, func_time_sca), mean)
    m10_means0 <- cellStats(stackApply(ssub[[which(m10i)]],indices=m10, func_time_sca), mean)
    m11_means0 <- cellStats(stackApply(ssub[[which(m11i)]],indices=m11, func_time_sca), mean)
    m12_means0 <- cellStats(stackApply(ssub[[which(m12i)]],indices=m12, func_time_sca), mean)
    # RF MONTH
    m1_means0r <- cellStats(stackApply(rsub[[which(m1i)]],indices=m1, func_time_rf), mean)
    m2_means0r <- cellStats(stackApply(rsub[[which(m2i)]],indices=m2, func_time_rf), mean)
    m3_means0r <- cellStats(stackApply(rsub[[which(m3i)]],indices=m3, func_time_rf), mean)
    m4_means0r <- cellStats(stackApply(rsub[[which(m4i)]],indices=m4, func_time_rf), mean)
    m5_means0r <- cellStats(stackApply(rsub[[which(m5i)]],indices=m5, func_time_rf), mean)
    m6_means0r <- cellStats(stackApply(rsub[[which(m6i)]],indices=m6, func_time_rf), mean)
    m7_means0r <- cellStats(stackApply(rsub[[which(m7i)]],indices=m7, func_time_rf), mean)
    m8_means0r <- cellStats(stackApply(rsub[[which(m8i)]],indices=m8, func_time_rf), mean)
    m9_means0r <- cellStats(stackApply(rsub[[which(m9i)]],indices=m9, func_time_rf), mean)
    m10_means0r <- cellStats(stackApply(rsub[[which(m10i)]],indices=m10, func_time_rf), mean)
    m11_means0r <- cellStats(stackApply(rsub[[which(m11i)]],indices=m11, func_time_rf), mean)
    m12_means0r <- cellStats(stackApply(rsub[[which(m12i)]],indices=m12, func_time_rf), mean)
    #
    #
    ## YEAR_VALUES ------------->
    sca_mean_yr0 <- cellStats(stackApply(ssub,indices=yidx, mean), mean)
    sca_mean_yr = mean(sca_mean_yr0,na.rm=T)
    sca_sd_yr0 <- cellStats(stackApply(ssub,indices=yidx, mean), sd)
    sca_sd_yr = mean(sca_sd_yr0,na.rm=T)
    sca_min_yr0 <- cellStats(stackApply(ssub,indices=yidx, min), mean)
    sca_min_yr = mean(sca_min_yr0,na.rm=T)
    #
    rf_mean_yr0 <- cellStats(stackApply(rsub,indices=yidx, mean), mean)
    rf_mean_yr = mean(rf_mean_yr0,na.rm=T)
    rf_sum_yr0 <- cellStats(stackApply(rsub,indices=yidx, sum), mean)
    rf_sum_yr = sum(rf_sum_yr0,na.rm=T)
    rf_sd_yr0 <- cellStats(stackApply(rsub,indices=yidx, mean), sd)
    rf_sd_yr = mean(rf_sd_yr0,na.rm=T)
    rf_max_yr0 <- cellStats(stackApply(rsub,indices=yidx, max), mean)
    rf_max_yr = mean(rf_max_yr0,na.rm=T)
    #
    ts = 1:16 # replace all (take unique years...)
    #
    ## YEAR TRENDS
    dfyv = data.frame(t = 1:16,sca_mean_yr0,sca_sd_yr0,sca_min_yr0,
                      rf_mean_yr0,rf_sum_yr0,rf_sd_yr0,rf_max_yr0)
    # run trends
    dtr1 = dfyv %>% melt(.,id.vars = 't') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~t, data=.)) %>% 
      tidy(dfit) %>% filter(term=="t")
    dtr2 = dfyv %>% melt(.,id.vars = 't') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~t, data=.)) %>% glance(dfit)
    # save values
    tr_est_yr = dtr1$estimate
    tr_pval_yr = dtr1$p.value
    tr_ste_yr = dtr1$std.error
    tr_sig_yr = dtr2$sigma
    tr_r2_yr = dtr2$r.squared
    #
    ## YEAR SCA~RF
    cor_scarf = dfyv %>%  
      dplyr::summarize(cor_mean = round(cor(sca_mean_yr0,rf_mean_yr0),2),
                       cor_smean_rsum = round(cor(sca_mean_yr0,rf_sum_yr0),2),
                       cor_smin_rsum = round(cor(sca_min_yr0,rf_sum_yr0),2),
                       cor_smin_rmean = round(cor(sca_min_yr0,rf_mean_yr0),2),
                       cor_sd = round(cor(sca_sd_yr0,rf_sd_yr0),2))
    # **if cor is significant then rerun with more lm
    scarf_mean1 = dfyv %>% do(dfit = lm(sca_mean_yr0~rf_mean_yr0, data=.)) %>% 
      tidy(dfit) %>% filter(term=="rf_mean_yr0") %>% select(-term,-statistic)
    scarf_mean2 = dfyv %>% do(dfit = lm(sca_mean_yr0~rf_mean_yr0, data=.)) %>% glance(dfit)
    tr_scarf_mean = c(as.numeric(scarf_mean1),scarf_mean2$r.squared)
    # COUNT npix
    npix = max(cellStats(ssub,fun_count_na,na.rm=T))
    snpix = cellStats(ssub,fun_count_na,na.rm=T); snpix[snpix < (npix*0.25)] = npix # if less than 1/4 of glacier
    ss50 = ssub; ss50[ss50 < 50] = NA # exclude pixels less than X
    snpix50 = cellStats(ss50,fun_count_na,na.rm=T);snpix50[snpix50 < (npix*0.25)] = npix
    spixmin <- with(data.frame(value=snpix,factor=yidx), tapply(value, factor, min,na.rm=T))
    spixmin50 <- with(data.frame(value=snpix50,factor=yidx), tapply(value, factor, min,na.rm=T))
    # trend for pixels
    ppx = data.frame(y=1:16,spixmin) %>%
      do(dfit = lm(spixmin~y, data=.)) %>% 
      tidy(dfit) %>% filter(term=="y") %>% select(-term,-statistic,-std.error)
    ppx50 = data.frame(y=1:16,spixmin50) %>%
      do(dfit = lm(spixmin50~y, data=.)) %>% 
      tidy(dfit) %>% filter(term=="y")  %>% select(-term,-statistic,-std.error)
    # save values
    spx_rat = (mean(spixmin,na.rm=T)/npix)
    spx_rat50 = (mean(spixmin50,na.rm=T)/npix)
    pix_vals = c(npix,spx_rat,spx_rat50)
    tr_pix = c(as.numeric(ppx),as.numeric(ppx50))
    #

    ## MONTH DFs ------------->
    dfmos = data.frame(y=1:16,m1_means0,m2_means0,m3_means0,m4_means0,m5_means0,m6_means0,
                       m7_means0,m8_means0,m9_means0,m10_means0,m11_means0,m12_means0)
    names(dfmos) <- c("y",paste0("V",1:12))
    dfmor = data.frame(y=1:16,m1_means0r,m2_means0r,m3_means0r,m4_means0r,m5_means0r,m6_means0r,
                       m7_means0r,m8_means0r,m9_means0r,m10_means0r,m11_means0r,m12_means0r)
    names(dfmor) <- c("y",paste0("V",1:12))
    #
    # MAIN MONTH VALUES
    sca_mean_month <- colMeans(dfmos)[2:13]
    sca_sd_month <- apply(dfmos, 2, sd)[2:13]
    sca_min_month <- apply(dfmos, 2, min)[2:13]
    rf_mean_month <- colMeans(dfmor)[2:13]
    rf_sum_month <- colSums(dfmor)[2:13]
    rf_sd_month <- apply(dfmor, 2, sd)[2:13]
    rf_max_month <- apply(dfmor, 2, max)[2:13]
    # 1 value
    sca_whmin_month <- as.numeric(which.min(sca_min_month))
    rf_whmax_month <- as.numeric(which.max(rf_mean_month))
    rf_whmaxsum_month <- as.numeric(which.max(rf_sum_month))
    #
    ## MONTH TRENDS
    # SCA
    dm1 = dfmos %>% melt(.,id.vars = 'y') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~y, data=.)) %>% 
      tidy(dfit) %>% filter(term=="y") %>% select(-term,-statistic)
    dm2 = dfmos %>% melt(.,id.vars = 'y') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~y, data=.)) %>% glance(dfit)
    # save values
    trsca_est_month = dm1$estimate
    trsca_pval_month = dm1$p.value
    trsca_ste_month = dm1$std.error
    trsca_sig_month = dm2$sigma
    trsca_r2_month = dm2$r.squared
    # SCA
    dm1 = dfmor %>% melt(.,id.vars = 'y') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~y, data=.)) %>% 
      tidy(dfit) %>% filter(term=="y") %>% select(-term,-statistic)
    dm2 = dfmor %>% melt(.,id.vars = 'y') %>% 
      group_by(variable) %>%
      do(dfit = lm(value~y, data=.)) %>% glance(dfit)
    # save values (12 values)
    trf_est_month = dm1$estimate
    trf_pval_month = dm1$p.value
    trf_ste_month = dm1$std.error
    trf_sig_month = dm2$sigma
    trf_r2_month = dm2$r.squared
    
    ## YEAR SCA~RF
    mrf = dfmor %>%
      melt(id.vars="y")
    srf = dfmos %>%
      melt(id.vars="y")
    cor_month = cbind(srf, rvalue = mrf$value) %>% na.omit() %>% 
      group_by(variable) %>% 
      dplyr::summarize(corr = round(cor(rvalue,value),2)) %>% select(corr)
    dmon1 = cbind(srf, rvalue = mrf$value) %>%
      group_by(variable) %>%
      do(dfit = lm(value~rvalue, data=.)) %>%
      tidy(dfit) %>% filter(term=="rvalue") %>% select(-term,-statistic)
    dmon2 = cbind(srf, rvalue = mrf$value) %>%
      group_by(variable) %>%
      do(dfit = lm(value~rvalue, data=.)) %>% glance(dfit)
    # save values
    tr_scarf_est_month = dmon1$estimate
    tr_scarf_pval_month = dmon1$p.value
    tr_scarf_ste_month = dmon1$std.error
    tr_scarf_sig_month = dmon2$sigma
    tr_scarf_r2_month = dmon2$r.squared
    
    # COUNT npix (CHANGE BY MONTH!!!!!!)
    snpixm = cellStats(ssub,fun_count_na,na.rm=T) 
    spixmonth <- with(data.frame(value=snpixm,factor=midx), tapply(value, factor, min,na.rm=T))
    snpixm[snpixm < (npix*0.25)] = npix # if less than 1/4 of glacier
    spixmonth_c <- with(data.frame(value=snpixm,factor=midx), tapply(value, factor, min,na.rm=T))
    snpix50m = cellStats(ss50,fun_count_na,na.rm=T)
    spixmonth50 <- with(data.frame(value=snpix50m,factor=midx), tapply(value, factor, min,na.rm=T))
    snpix50m[snpix50m < (npix*0.25)] = npix
    spixmonth50_c <- with(data.frame(value=snpix50m,factor=midx), tapply(value, factor, min,na.rm=T))
    
    spix_allmonth <- with(data.frame(value=snpixm,factor=ymidx), tapply(value, factor, min,na.rm=T))
    dfspix = cbind(y=1:16,as.data.frame(matrix(spix_allmonth,ncol=12,byrow=TRUE)))
    
    # trend for pixels
    ppxm = dfspix %>% melt(.,id.vars="y") %>% group_by(variable) %>% 
      do(dfit = lm(value~y, data=.)) %>%
      tidy(dfit) %>% filter(term=="y") %>% select(-term,-statistic,-std.error)
    ppx50m = dfspix %>% melt(.,id.vars="y") %>% group_by(variable) %>% 
      do(dfit = lm(value~y, data=.)) %>%
      tidy(dfit) %>% filter(term=="y")  %>% select(-term,-statistic,-std.error)
    # save values
    spxmo_rat = spixmonth/npix
    spxmo_rat_c = spixmonth_c/npix
    spxmo_rat50_c = spixmonth50_c/npix
    tr_pix_month_est = ppxm$estimate
    tr_pix_month_pval = ppxm$p.value
    # tr_pix50_month_est = ppxm$estimate
    
    
    
    # can be saved as a list of values???
    # SAVE VALUES
    var_yr_mo_ls = list(sca_mean_yr,sca_sd_yr,sca_min_yr,  ## YEAR VALUES!!!
                        rf_mean_yr,rf_sum_yr,rf_sd_yr,rf_max_yr, # first 7 are values
                        tr_est_yr, # 7 trends (smean,ssd,smin,rmean,rsum,rsd,rmax)
                        tr_pval_yr,tr_ste_yr,tr_sig_yr,tr_r2_yr, # (all 7 values each)
                        cor_scarf, # 5 values (means, smean_rsum, smin_rsum, smin_rmean, sd)
                        tr_scarf_mean, # 4 (estimate, std.error, p.value, r2)
                        npix,spx_rat,spx_rat50, #  3 (npix, ave_minpix, ave_minpix50)
                        tr_pix,   # 4 (tr_minpix_estimate, tr_minpix_p.value, tr_minpix50_estimate, tr_minpix50_p.value)
                        sca_mean_month,  # MONTH VALUES!!! (12 values)
                        sca_sd_month,
                        sca_min_month,
                        rf_mean_month,
                        rf_sum_month,
                        rf_sd_month,
                        rf_max_month,
                        sca_whmin_month, # 1 variable
                        rf_whmax_month, # 1 variable
                        rf_whmaxsum_month, # 1 variable
                        trsca_est_month, # MONTHLY TRENDS (12 Vals)
                        trsca_pval_month,
                        trsca_ste_month,
                        trsca_sig_month,
                        trsca_r2_month,
                        trf_est_month,
                        trf_pval_month,
                        trf_ste_month,
                        trf_sig_month,
                        trf_r2_month,
                        cor_month, # SCARF by MONTH (12 values)
                        tr_scarf_est_month,
                        tr_scarf_pval_month,
                        tr_scarf_ste_month,
                        tr_scarf_sig_month,
                        tr_scarf_r2_month,
                        spxmo_rat,spxmo_rat_c,spxmo_rat50_c, 
                        tr_pix_month_est,
                        tr_pix_month_pval
    )
    # add pixels by month
    nm_h1 = c("smean","ssd","smin","rmean","rsum","rsd","rmax")
    
    # NAMES
    nms = c("sca_mean_yr","sca_sd_yr","sca_min_yr",  ## YEAR VALUES!!!
            "rf_mean_yr","rf_sum_yr","rf_sd_yr","rf_max_yr", # 1 each
            paste0(rep("tr_est_yr",7),nm_h1), # 7 trends (smean,ssd,smin,rmean,rsum,rsd,rmax)
            paste0(rep("tr_pval_yr",7),nm_h1),
            paste0(rep("tr_ste_yr",7),nm_h1),
            paste0(rep("tr_sig_yr",7),nm_h1),
            paste0(rep("tr_r2_yr",7),nm_h1),
            paste0(rep("cor_scarf",5),c('means', 'smean_rsum', 'smin_rsum', 'smin_rmean', 'sd')),
            paste0(rep("tr_scarf_mean",4),c('est', 'ste', 'pval', 'r2')),
            'npix','spx_rat','spx_rat50', #  3 (npix, ave_minpix, ave_minpix50)
            paste0(rep("tr_pix",4),c('est', 'pval', 'est50', 'pval50')),
            paste0(rep("sca_mean_month",12),1:12), # MONTH VALUES!!! (12 values)
            paste0(rep("sca_sd_month",12),1:12),
            paste0(rep("sca_min_month",12),1:12),
            paste0(rep("rf_mean_month",12),1:12),
            paste0(rep("rf_sum_month",12),1:12),
            paste0(rep("rf_sd_month",12),1:12),
            paste0(rep("rf_max_month",12),1:12),
            "sca_whmin_month", # 1 variable
            "rf_whmax_month", # 1 variable
            "rf_whmaxsum_month", # 1 variable
            paste0(rep("trsca_est_month",12),1:12), # MONTHLY TRENDS (12 Vals)
            paste0(rep("trsca_pval_month",12),1:12),
            paste0(rep("trsca_ste_month",12),1:12),
            paste0(rep("trsca_sig_month",12),1:12),
            paste0(rep("trsca_r2_month",12),1:12),
            paste0(rep("trf_est_month",12),1:12),
            paste0(rep("trf_pval_month",12),1:12),
            paste0(rep("trf_ste_month",12),1:12),
            paste0(rep("trf_sig_month",12),1:12),
            paste0(rep("trf_r2_month",12),1:12),
            paste0(rep("cor_month",12),1:12), # SCARF by MONTH (12 values)
            paste0(rep("tr_scarf_est_month",12),1:12),
            paste0(rep("tr_scarf_pval_month",12),1:12),
            paste0(rep("tr_scarf_ste_month",12),1:12),
            paste0(rep("tr_scarf_sig_month",12),1:12),
            paste0(rep("tr_scarf_r2_month",12),1:12),
            paste0(rep("spxmo_rat",12),1:12),
            paste0(rep("spxmo_rat_c",12),1:12),
            paste0(rep("spxmo_rat50_c",12),1:12),
            paste0(rep("tr_pix_month_est",12),1:12),
            paste0(rep("tr_pix_month_pval",12),1:12)
    )
    
    # print(paste("varnames = ",length(unlist(var_yr_mo_ls)) == length(nms)))
    
    st = (ci+1); en = st+(length(nms)-1) # 10:406
    df_glaciers[i,(st:en)] <- unlist(var_yr_mo_ls)
    #
    print(i)
    #
    if (i%%20==0){
      end = Sys.time() - strt 
      print(end)
    }
    
    
  }
  # add names
  names(df_glaciers) = c(o_names,nms)
  
  end = Sys.time() - strt 
  print(end)
  return(df_glaciers)
}



##############################################################################
############ SCA/RF By Normalized glacier elevation
##############################################################################

# DONE
#####################################################################################
glacier.elv <- function(df_glaciers,rgi,demL,noZero = TRUE,threshold = c(60,15),
                        rfpath,scapath,tfpath){
  
  addcols = (7*4*5)-1
  onames <- names(df_glaciers)
  oc = ncol(df_glaciers)+1
  org <- (oc:(oc+addcols))
  df_glaciers[,org] <- NA
  for (i in 1:(length(df_glaciers$Object.ID)-1)){
    print(i)
    g.match.rf.sca(df_glaciers$Object.ID[i],rgi,demL,noZero = noZero,threshold = threshold,
                   rfpath,scapath,tfpath)
    # rsub | ssub | match_dates | elv
    gvals <- glacier.elv.trends(rsub,ssub,match_dates,elv,
                                   t.func = 'mean')
    
    df_glaciers[i,org] <- gvals
  }
  names(df_glaciers) <- c(onames,paste0(rep(c("yvals","mamvals","amjvals","mjjvals","jjavals","jasvals",
             "djfvals","yvals_s","mamvals_s","amjvals_s","mjjvals_s","jjavals_s","jasvals_s",
             "djfvals_s","ytr","mamtr","amjtr","mjjtr","jjatr","jastr",
             "djftr","ytr_s","mamtr_s","amjtr_s","mjjtr_s","jjatr_s","jastr_s",
             "djftr_s"),each=5),seq(0.2,1,0.2)))
  return(df_glaciers)
}


# DONE
#####################################################################################
glacier.elv.trends <- function(rsub,ssub,match_dates,elv,
                               t.func = 'mean'){
  # normalized elevation
  elvn <- (elv - elv@data@min)/(elv@data@max - elv@data@min)
  
  # indexes
  yid = format(match_dates,"%Y"); yidx = as.numeric(yid)-2000
  mid = format(match_dates,"%m"); midx = as.numeric(mid)
  mam = midx == 3 | midx == 4 | midx == 5; mamx = yidx[mam]
  amj = midx == 4 | midx == 5 | midx == 6; amjx = yidx[amj]
  mjj = midx == 5 | midx == 6 | midx == 7; mjjx = yidx[mjj]
  jja = midx == 6 | midx == 7 | midx == 8; jjax = yidx[jja]
  jas = midx == 7 | midx == 8 | midx == 9; jasx = yidx[jas]
  
  djf = midx == 1 | midx == 2 | midx == 12; djfx = yidx[djf]

  # elv intervals
  eint <- seq(0,1,0.2)
  h = hist(elvn,breaks=eint,plot=FALSE)
  npe = h$counts # number of pixels at each elevation
  
  ## RF 
  # take mean and sd
  yvals <- elv.vals(rsub,id=NULL,idx=NULL,elvn,npe,eint)
  mamvals <- elv.vals(rsub,id=mam,idx=NULL,elvn,npe,eint)
  amjvals <- elv.vals(rsub,id=amj,idx=NULL,elvn,npe,eint)
  mjjvals <- elv.vals(rsub,id=mjj,idx=NULL,elvn,npe,eint)
  jjavals <- elv.vals(rsub,id=jja,idx=NULL,elvn,npe,eint)
  jasvals <- elv.vals(rsub,id=jas,idx=NULL,elvn,npe,eint)
  djfvals <- elv.vals(rsub,id=djf,idx=NULL,elvn,npe,eint)
  
  # RF trends
  ytr <- elv.vals(rsub,id=NULL,idx=yidx,elvn,npe,eint)
  mamtr <- elv.vals(rsub,id=mam,idx=mamx,elvn,npe,eint)
  amjtr <- elv.vals(rsub,id=amj,idx=amjx,elvn,npe,eint)
  mjjtr <- elv.vals(rsub,id=mjj,idx=mjjx,elvn,npe,eint)
  jjatr <- elv.vals(rsub,id=jja,idx=jjax,elvn,npe,eint)
  jastr <- elv.vals(rsub,id=jas,idx=jasx,elvn,npe,eint)
  djftr <- elv.vals(rsub,id=djf,idx=djfx,elvn,npe,eint)
  
  ## SCA
  # take mean and sd
  yvals2 <- elv.vals(ssub,id=NULL,idx=NULL,elvn,npe,eint)
  mamvals2 <- elv.vals(ssub,id=mam,idx=NULL,elvn,npe,eint)
  amjvals2 <- elv.vals(ssub,id=amj,idx=NULL,elvn,npe,eint)
  mjjvals2 <- elv.vals(ssub,id=mjj,idx=NULL,elvn,npe,eint)
  jjavals2 <- elv.vals(ssub,id=jja,idx=NULL,elvn,npe,eint)
  jasvals2 <- elv.vals(ssub,id=jas,idx=NULL,elvn,npe,eint)
  djfvals2 <- elv.vals(ssub,id=djf,idx=NULL,elvn,npe,eint)
  
  # SCA trends
  ytr2 <- elv.vals(ssub,id=NULL,idx=yidx,elvn,npe,eint)
  mamtr2 <- elv.vals(ssub,id=mam,idx=mamx,elvn,npe,eint)
  amjtr2 <- elv.vals(ssub,id=amj,idx=amjx,elvn,npe,eint)
  mjjtr2 <- elv.vals(ssub,id=mjj,idx=mjjx,elvn,npe,eint)
  jjatr2 <- elv.vals(ssub,id=jja,idx=jjax,elvn,npe,eint)
  jastr2 <- elv.vals(ssub,id=jas,idx=jasx,elvn,npe,eint)
  djftr2 <- elv.vals(ssub,id=djf,idx=djfx,elvn,npe,eint)
  
  avals <- c(yvals,mamvals,amjvals,mjjvals,jjavals,jasvals,
             djfvals,yvals2,mamvals2,amjvals2,mjjvals2,jjavals2,jasvals2,
             djfvals2,ytr,mamtr,amjtr,mjjtr,jjatr,jastr,
             djftr,ytr2,mamtr2,amjvals2,mjjtr2,jjatr2,jastr2,
             djftr2)
  
  return(avals)
  
}

# DONE
#####################################################################################
elv.vals <- function(stk,id=NULL,idx=NULL,elvn,npe,
                     eint=seq(0,1,0.2), fun.m = 'mean'){
  # determine whether trend or vals
  if(is.null(idx)){
    if(is.null(id)){
      # determine yr or subset
      yr <- stackApply(stk,1, fun.m)
    } else{
      yr <- stackApply(stk[[which(id)]],1, fun.m)
    }
  } else{
    if(is.null(id)){
      yr <- stackApply(stk,indices = idx, fun.m)
    } else{
      yr <- stackApply(stk[[which(id)]],indices = idx, fun.m)
    }
    time = 1:nlayers(yr)
    fun_lm_time = function(x) { if (is.na(x[1])){ NA } else {lm(x ~ time)$coefficients[2] }} 
    yr <- calc(yr,fun_lm_time)
  }
  # subset based on elevation
  # y1 <- yr[elvn >= eint[1] & elvn <= eint[2]]
  # y2 <- yr[elvn > eint[2] & elvn <= eint[3]]
  # y3 <- yr[elvn > eint[3] & elvn <= eint[4]]
  # y4 <- yr[elvn > eint[4] & elvn <= eint[5]]
  # y5 <- yr[elvn > eint[5] & elvn <= eint[6]]
  yls <- list()
  for (i in 1:length(npe)){
    if(i==1){yvs <- yr[elvn >= eint[i] & elvn <= eint[i+1]]}else{
      yvs <- yr[elvn > eint[i] & elvn <= eint[i+1]]
      
    }
    yls <- c(yls,c(mean(yvs,na.rm=T),sd(yvs,na.rm=T)))
  }
  #yls <- unlist(yls) # only keep means
  ymeans <- unlist(yls)[c(TRUE,FALSE)]
  return(ymeans)
}


# MONTHLY ELEVATION
# NEW 04/24
#####################################################################################
glacier.elv.month <- function(df_glaciers,rgi,demL,noZero = TRUE,threshold = c(60,15),
                        rfpath,scapath,tfpath, fun.t = 'mean'){
  
  addcols = (13*4*5)-1
  onames <- names(df_glaciers)
  oc = ncol(df_glaciers)+1
  org <- (oc:(oc+addcols))
  df_glaciers[,org] <- NA
  for (i in 1:(length(df_glaciers$Object.ID)-1)){
    print(i)
    g.match.rf.sca(df_glaciers$Object.ID[i],rgi,demL,noZero = noZero,threshold = threshold,
                   rfpath,scapath,tfpath)
    # rsub | ssub | match_dates | elv
    gvals <- glacier.elv.trends.month(rsub,ssub,match_dates,elv,
                                t.func = fun.t)
    
    df_glaciers[i,org] <- gvals
  }
  names(df_glaciers) <- c(onames,paste0(rep(c("yvals",paste0("m",1:12,"vr"),"yvals_s",paste0("m",1:12,"vs"),
                                              "ytr",paste0("m",1:12,"rtr"),
                                              "ytr_s",paste0("m",1:12,"str")),each=5),seq(0.2,1,0.2)))
  return(df_glaciers)
}


# NEW 04/24 for MONTHS (goes with above)
#####################################################################################
glacier.elv.trends.month <- function(rsub,ssub,match_dates,elv,
                               t.func = 'mean'){
  # normalized elevation
  elvn <- (elv - elv@data@min)/(elv@data@max - elv@data@min)
  
  # indexes
  yid = format(match_dates,"%Y"); yidx = as.numeric(yid)-2000
  mid = format(match_dates,"%m"); midx = as.numeric(mid)
  # monthly indices
  m1i = midx == 1; m1 = yidx[m1i] # March
  m2i = midx == 2; m2 = yidx[m2i] # March
  m3i = midx == 3; m3 = yidx[m3i] # March
  m4i = midx == 4; m4 = yidx[m4i] # March
  m5i = midx == 5; m5 = yidx[m5i] # March
  m6i = midx == 6; m6 = yidx[m6i] # March
  m7i = midx == 7; m7 = yidx[m7i] # March
  m8i = midx == 8; m8 = yidx[m8i] # March
  m9i = midx == 9; m9 = yidx[m9i] # March
  m10i = midx == 10; m10 = yidx[m10i] # March
  m11i = midx == 11; m11 = yidx[m11i] # March
  m12i = midx == 12; m12 = yidx[m12i] # March
 
  
  # elv intervals
  eint <- seq(0,1,0.2)
  h = hist(elvn,breaks=eint,plot=FALSE)
  npe = h$counts # number of pixels at each elevation
  
  ## RF 
  # take mean and sd
  yvals <- elv.vals(rsub,id=NULL,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m1v <- elv.vals(rsub,id=m1i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m2v <- elv.vals(rsub,id=m2i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m3v <- elv.vals(rsub,id=m3i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m4v <- elv.vals(rsub,id=m4i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m5v <- elv.vals(rsub,id=m5i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m6v <- elv.vals(rsub,id=m6i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m7v <- elv.vals(rsub,id=m7i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m8v <- elv.vals(rsub,id=m8i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m9v <- elv.vals(rsub,id=m9i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m10v <- elv.vals(rsub,id=m10i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m11v <- elv.vals(rsub,id=m11i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m12v <- elv.vals(rsub,id=m12i,idx=NULL,elvn,npe,eint, fun.m = t.func)

  # RF trends
  ytr <- elv.vals(rsub,id=NULL,idx=yidx,elvn,npe,eint, fun.m = t.func)
  m1tr <- elv.vals(rsub,id=m1i,idx=m1,elvn,npe,eint, fun.m = t.func)
  m2tr <- elv.vals(rsub,id=m2i,idx=m2,elvn,npe,eint, fun.m = t.func)
  m3tr <- elv.vals(rsub,id=m3i,idx=m3,elvn,npe,eint, fun.m = t.func)
  m4tr <- elv.vals(rsub,id=m4i,idx=m4,elvn,npe,eint, fun.m = t.func)
  m5tr <- elv.vals(rsub,id=m5i,idx=m5,elvn,npe,eint, fun.m = t.func)
  m6tr <- elv.vals(rsub,id=m6i,idx=m6,elvn,npe,eint, fun.m = t.func)
  m7tr <- elv.vals(rsub,id=m7i,idx=m7,elvn,npe,eint, fun.m = t.func)
  m8tr <- elv.vals(rsub,id=m8i,idx=m8,elvn,npe,eint, fun.m = t.func)
  m9tr <- elv.vals(rsub,id=m9i,idx=m9,elvn,npe,eint, fun.m = t.func)
  m10tr <- elv.vals(rsub,id=m10i,idx=m10,elvn,npe,eint, fun.m = t.func)
  m11tr <- elv.vals(rsub,id=m11i,idx=m11,elvn,npe,eint, fun.m = t.func)
  m12tr <- elv.vals(rsub,id=m12i,idx=m12,elvn,npe,eint, fun.m = t.func)
  
  ## SCA 
  # take mean and sd
  yvals2 <- elv.vals(ssub,id=NULL,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m1v2 <- elv.vals(ssub,id=m1i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m2v2 <- elv.vals(ssub,id=m2i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m3v2 <- elv.vals(ssub,id=m3i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m4v2 <- elv.vals(ssub,id=m4i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m5v2 <- elv.vals(ssub,id=m5i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m6v2 <- elv.vals(ssub,id=m6i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m7v2 <- elv.vals(ssub,id=m7i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m8v2 <- elv.vals(ssub,id=m8i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m9v2 <- elv.vals(ssub,id=m9i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m10v2 <- elv.vals(ssub,id=m10i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m11v2 <- elv.vals(ssub,id=m11i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  m12v2 <- elv.vals(ssub,id=m12i,idx=NULL,elvn,npe,eint, fun.m = t.func)
  
  # SCA trends
  ytr2 <- elv.vals(ssub,id=NULL,idx=yidx,elvn,npe,eint, fun.m = t.func)
  m1tr2 <- elv.vals(ssub,id=m1i,idx=m1,elvn,npe,eint, fun.m = t.func)
  m2tr2 <- elv.vals(ssub,id=m2i,idx=m2,elvn,npe,eint, fun.m = t.func)
  m3tr2 <- elv.vals(ssub,id=m3i,idx=m3,elvn,npe,eint, fun.m = t.func)
  m4tr2 <- elv.vals(ssub,id=m4i,idx=m4,elvn,npe,eint, fun.m = t.func)
  m5tr2 <- elv.vals(ssub,id=m5i,idx=m5,elvn,npe,eint, fun.m = t.func)
  m6tr2 <- elv.vals(ssub,id=m6i,idx=m6,elvn,npe,eint, fun.m = t.func)
  m7tr2 <- elv.vals(ssub,id=m7i,idx=m7,elvn,npe,eint, fun.m = t.func)
  m8tr2 <- elv.vals(ssub,id=m8i,idx=m8,elvn,npe,eint, fun.m = t.func)
  m9tr2 <- elv.vals(ssub,id=m9i,idx=m9,elvn,npe,eint, fun.m = t.func)
  m10tr2 <- elv.vals(ssub,id=m10i,idx=m10,elvn,npe,eint, fun.m = t.func)
  m11tr2 <- elv.vals(ssub,id=m11i,idx=m11,elvn,npe,eint, fun.m = t.func)
  m12tr2 <- elv.vals(ssub,id=m12i,idx=m12,elvn,npe,eint, fun.m = t.func)
  
  avals <- c(yvals,m1v,m2v,m3v,m4v,m5v,m6v,m7v,m8v,m9v,m10v,m11v,m12v,
             yvals2,m1v2,m2v2,m3v2,m4v2,m5v2,m6v2,m7v2,m8v2,m9v2,m10v2,m11v2,m12v2,
             ytr,m1tr,m2tr,m3tr,m4tr,m5tr,m6tr,m7tr,m8tr,m9tr,m10tr,m11tr,m12tr,
             ytr2,m1tr2,m2tr2,m3tr2,m4tr2,m5tr2,m6tr2,m7tr2,m8tr2,m9tr2,m10tr2,m11tr2,m12tr2
             )
  # 52 (13x4) vars x 5(elvn) = 260 values
  return(avals)
  
}



########################################################
# EXTRA FUNCTIONS
########################################################
fun_lm_time = function(x) { if (is.na(x[1])){ NA } else {time = 1:nlayers(x);lm(x ~ time)$coefficients[2] }} 
fun_count_na = function(i, ...) sum(!is.na(i))
fun_count = function(i, ...) sum((i) < mean(i))
fun_count_mean = function(x) {
  comp <- x[1]
  main <- x[2:length(x)]
  length(main[main >= comp])
  #sum(main >= comp,na.rm=T)
}
# Function to fill na (with average or nearest neighbor)
fill.na.mean <- function(dat, takeMean=TRUE) {
  N <- length(dat)
  na.pos <- which(is.na(dat))
  if (length(na.pos) %in% c(0, N)) {
    return(dat)
  }
  non.na.pos <- which(!is.na(dat))
  intervals  <- findInterval(na.pos, non.na.pos,all.inside = TRUE)
  left.pos   <- non.na.pos[pmax(1, intervals)]
  right.pos  <- non.na.pos[pmin(N, intervals+1)]
  if(takeMean){
    # mean
    dat[na.pos] <- apply(cbind(dat[left.pos],dat[right.pos]), 1, mean)
  }else{
    # do nearest value (before or after)
    left.dist  <- na.pos - left.pos
    right.dist <- right.pos - na.pos
    dat[na.pos] <- ifelse(left.dist <= right.dist,
                          dat[left.pos], dat[right.pos])
  }
  return(dat)
}
fill.9999 <- function(dat, valfill= -9999, takeMean=TRUE) {
  N <- length(dat)
  na.pos <- which(dat==valfill)
  non.na.pos <- which(!is.na(dat))
  intervals  <- findInterval(na.pos, non.na.pos,all.inside = TRUE)
  left.pos   <- non.na.pos[pmax(1, intervals)]
  right.pos  <- non.na.pos[pmin(N, intervals+1)]
  if(takeMean){
    # mean
    dat[na.pos] <- apply(cbind(dat[left.pos],dat[right.pos]), 1, mean,na.rm=T)
  }else{
    # do nearest value (before or after)
    left.dist  <- na.pos - left.pos
    right.dist <- right.pos - na.pos
    dat[na.pos] <- ifelse(left.dist <= right.dist,
                          dat[left.pos], dat[right.pos])
  }
  return(dat)
}
fill.NA <- function(x,na.rm = F){ 
  fill <- approxfun(x, y = NULL, method = 'linear', rule = 2); 
  new <- fill(1:length(x)); 
  rec <- which(is.na(x)); 
  x[rec] <- new[rec]; 
  return(x) 
}
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}


########################################################
# correlation of rf and sca pixel by pixel
cor_table_means <- function(dates, condition, rf_raster, sca_raster, onlySCA=FALSE){
  df_cor <- data.frame(array(NA, dim=c(length(dates),2)),
                       row.names=dates)
  
  # # fill in NAs (first preference to prior day)
  # condition <- na.locf(condition, na.rm=F)
  # # if(is.na(condition[1])){condition <- na.locf(condition,fromLast = TRUE)}
  # if(is.na(condition[1])){condition[is.na(condition)] <- (min(condition,na.rm=T) - 1)}
  condition2 = condition[!is.na(condition)]
  ry = rf_raster[[condition2]]
  sy = sca_raster[[condition2]]
  
  
  
  # df_tmp = data.frame(r = cellStats(ry, mean), s = cellStats(sy, mean),
  #                     row.names = ????)
  
  ####$%^&* DOESN"T WORK ROW NAMES ARE NOT SAME LENGTH (HOW TO MATCH DATES?)
  
  # fill all NA values (make regular time interval)
  df_cor[match(row.names(df_tmp), row.names(df_cor)),1:2] <- df_tmp[,1:2]
  
  if(isTRUE(fillNA)){df_cor <- na.locf(df_cor)}
  if(isTRUE(onlySCA)){df_cor <- df_cor[,2]}
  return(df_cor)
}

## NOT FINISHED 
#####################################################################################
glacier.df.all <- function(df_glaciers, rsub, ssub, match_dates, threshold=NULL, 
                           func = 'mean'){
  # 
  # Data frame of RF and SCA for ALL GLACIERS IN DF 
  # 
  strt = Sys.time(); ci=ncol(df_glaciers)
  o_names <- names(df_glaciers)
  df_glaciers[,10:50] <- NA # was 38 (making room for new variables)
  ## LOOP THROUGH EACH GLACIER
  # !!!!
  # quick fix (last name has no SCA)
  for (i in 1:length(df_glaciers$Object.ID)){
    #
    g.match.rf.sca(df_glaciers$Object.ID[i],rgi,demL,noZero = TRUE,threshold = threshold,
                   rfpath,scapath,tfpath)
    # ...creates: rsub | ssub | elv | match_dates
    #
    # CROSS-CORRELATION
    rf_vals <- cellStats(rsub,mean)
    sca_vals <- cellStats(ssub,mean)
    continuous_weeks <- as.Date(seq(match_dates[1], tail(match_dates,1), by="7 days"))
    df_corr <- data.frame(time = continuous_weeks)
    df_corr$sca <- df_corr$rf <- NA
    match_weeks_idx <- match(match_dates,continuous_weeks)
    df_corr$rf[match_weeks_idx] <- rf_vals; df_corr$rf <- fill.na.mean(df_corr$rf)
    df_corr$sca[match_weeks_idx] <- sca_vals; df_corr$sca <- fill.na.mean(df_corr$sca)
    
    if (i==1){
      allDF <- df_corr
      allDF[,4:(length(df_glaciers$Object.ID)*2)-2]
    } else{
      allDF
    }
    
    print(i)
  }
  end = Sys.time() - strt 
  print(end)
  return(df_glaciers)
}