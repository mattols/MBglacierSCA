#
# GLACIER trends in SCA/RF vs G. mass balance (Maurer et al., 2019)
# RUN MAIN CODE
# 11/12/2019
#
# # # # # # # # # # # #

################################
## LOAD ALL FUNCTIONS
source("glacierSCA_functions.R")

################################
## LOAD DATA
on_mac = TRUE
load.paths(macpaths=on_mac)
load.data(macpaths=on_mac)
mbg <- mb.csv.new(macpaths=on_mac)
# DEM, df of glacier IDs, rgi polygons, subregion shapefiles

############################################################################
### RUN
## Call function and calculate trends (threshold=c(RF, SCA) or NULL for no threshold)

# glacier-wide trends based on single variable
mbg_sca15 <- glacier.var.vals(mbg, scapath, rgi, demL, threshold=15)
mbg_rf0 <- glacier.var.vals(mbg, rfpath, rgi, demL, threshold=60)

#  trends based on glacier-wide averages by month (inlcudes SCA + RF)
mbg_months_6015 <- glacier.month.trends(mbg, rgi, demL, noZero = TRUE, threshold = c(60,15),
                                                   rfpath,scapath,tfpath)
mbg_months_60 <- glacier.month.trends(mbg, rgi, demL, noZero = FALSE, threshold = 60,
                                        rfpath,scapath,tfpath)

## monthly trends by pixel
mbg_months_p6015 <- glacier.month.trends.pix(mbg, rgi, demL, noZero = FALSE, threshold = c(60,15),
                                           rfpath,scapath,tfpath)

# trends in glacier SCA/RF by elevation and month
mbg_elv_months6015 <- glacier.elv.month(mbg,rgi,demL,noZero = FALSE,threshold = c(60,15),
                                      rfpath,scapath,tfpath)


# Calculate trend based on the min glcier SCA instead of mean 
mbg_months_p6015_min <- glacier.month.trends.pix(mbg, rgi, demL, noZero = FALSE, threshold = c(60,15),
                                             rfpath,scapath,tfpath,t.func='min',t.func_r='mean')

#
mbg_elv_months6015_min <- glacier.elv.month(mbg,rgi,demL,noZero = FALSE,threshold = c(60,15),
                                        rfpath,scapath,tfpath, fun.t = 'min')

############################################################################

############################################################################
## View results (simple)
mbg_sca15 %>% mutate( region = recode_factor(region, `12` = "West H.", `13` = "Central H.", `14` = "East H.")) %>%
  ggplot(., aes(x=y_t, y=geoMassBal)) +
  geom_point(aes(colour=meanElev,size=y_r2)) + geom_smooth(method='lm',se=F,linetype='dashed',colour='firebrick4') + 
  labs(y= bquote('Geodetic mass balance'~(m ~yr^-1)),x=bquote("Trend in glacier SCA" ~("%"~yr^-1 ~ 2000-2017))) +
  facet_wrap(~region) + theme(plot.margin = unit(c(3.1,1,3.1,1), "cm"),panel.spacing = unit(1, "lines")) +
  scale_color_continuous(name = "Elevation (m)") +
  scale_size_continuous(name=bquote("SCA trend" ~R^2 ))
# correlation
mbg_sca15 %>% mutate( region = recode_factor(region, `12` = "West H.", `13` = "Central H.", `14` = "East H.")) %>%
  group_by(region) %>% na.omit() %>% 
  summarize(cor_mb = cor(y_t,geoMassBal))
#
# View trends by month and elevation
mbg_elv_months6015 %>% mutate( region = recode_factor(region, `12` = "West H.", `13` = "Central H.", `14` = "East H.")) %>%
  select(region,geoMassBal,meanElev, e = contains("str0.4")) %>%
  melt(.,id.vars=c("region","geoMassBal","meanElev")) %>% 
  ggplot(aes(x=variable,y=value, colour=region)) + geom_boxplot() +
  labs(x= "Month of year",y=bquote("Trend in SCA" ~("%"~yr^-1))) +
  scale_x_discrete(labels=as.character(1:12)) +
  geom_hline(yintercept = 0,colour='black',linetype="dashed") + 
  theme(plot.margin = unit(c(3.2,1,3.2,1), "cm"),
        panel.spacing = unit(0.2, "lines")) +
  facet_grid(~region)