####################################################################
#                                                                  #
#  THIS SCRIPT USES OBSERVED SWE TO PREDICT OBSERVED STREAMFLOW    #
#  USING ALL POSSIBLE COMBINATIONS (e.g. x-choose-4) OF FOUR       #
#  OBSERVATIONAL LOCATIONS FOR A GIVEN WATERSHED.  EACH WATERSHED  #
#  USES A LEAVE-ONE-OUT APPROACH (i.e. using n-1 years to predict  #
#  the nth year) TO CALIBRATE/VALIDATE THE SWE PREDICTIVE MODEL.   #
#                                                                  #
#  THE METHOD IS REPEATED WITH EACH MODEL FOR EACH RESPECTIVE      #
#  WATERSHED.  SWE "OBSERVATIONS" ARE MODELED GRID CELL SWE FOR    #
#  GRID CELL CONTAINING RESPECTIVE OBSERVATION LOCATIONS, AND      #
#  "OBSERVED" STREAMFLOW COMES FROM SUMMING RUNOFF IN EACH         #
#  RESPECTIVE WATERSHED.                                           #
#                                                                  #
####################################################################

rm(list=ls())

## Load libraries
library(MASS)
library(ncdf4)
library(fields)
library(maps)
library(hydroGOF)
library(plotrix)
library(rworldmap)
library(RColorBrewer)
library(colorRamps)
library(matrixcalc)
library(dataRetrieval)

## Set directory paths for I/O
data.dir = "/Users/andrewbadger/Documents/SARP/SNOTEL_Selction/"
snotel.dir = "/Volumes/LabShare/SARP/BCSD_Point_Predictions/"
write.dir = "/Users/andrewbadger/Documents/SARP/BCSD_Predictors/Jacknife_Data/"

## Set basin information, NOTE THAT ALL BASINS ARE NOT USED FOR PUBLICATION
sitenames = c("Bear_Lake_Thomas", "Merced_Pohono_Bridge", "Cedar_River", "Thunder_Creek", "Sauk_River_Whitechuck", "Middle_Fork_Flathead", "Stehekin_River", "American_River", "Pacific_Creek", "Greys_River", "SF_Payette_River", "Johnson_Creek", "Lochsa_River", "Gore_Creek", "Hams_Fork", "Yellowstone_Corwin_Springs", "Shell_Creek", "South_Fork_Shoshone", "Little_Bighorn", "Lodge_Grass_Creek", "North_Fork_Powder_River", "Castle_Creek", "EF_Carson", "WF_Carson", "Bear", "Blacksmith", "Clear")
huc.ids = c(18040006, 18040008, 17110012, 17110005, 17110006, 17010207, 17020009, 17030002, 17040101, 17040103, 17050120, 17060208, 17060303, 14010003, 14040107, 10070002, 10080010, 10080013, 10080016, 10080016, 10090201, 10120110, 16050201, 16050201, 16010101, 16010203, 16030003)
gages = c(11230500, 11266500, 12115000, 12175500, 12186000, 12186000, 12451000, 12488500, 13011500, 13023000, 13235000, 13313000, 13337000, 09065500, 09223000, 06191500, 06278300, 06280300, 06289000, 06291500, 06311000, 06409000, 10309000, 10310000, 10011500, 10113500, 10194200)
domain = c(1, 2, 3, 2, 3, 3, 4, 5, 1, 2, 1, 2, 0, 5, 1, 0, 3, 3, 3, 3, 1, 4, 2, 4, 2, 2, 3)
nsitenames = length(sitenames)

## Loop over each site
for(s in 1:nsitenames){
	sitename = sitenames[s]
	chosen.hucs = huc.ids[s]
	gage.id = gages[s]
	buffer = 2
	snotel.thresh = 100
	time.thresh = 21
	
	## Get gage lat-lon information
	if (nchar(gage.id) < 8){
  		toAdd = as.character(integer(8-nchar(gage.id)))
  		toAdd = paste(toAdd, collapse="")
  		gage.id = paste(c(toAdd, gage.id), collapse = "")
		}

	data = readNWISsite(siteNumber=gage.id)
	darea = data[30] #mi2
	cdarea = data[31] #mi2
	glon = data[8] + 360
	glat = data[7] 
	
	glon = as.numeric(glon)
	glat = as.numeric(glat)
	
	## Open HUC Map to get location of HUC in the the BCSD data - HUC information has been mapped to BCSD grid
	tmp.file = paste(data.dir,"REEDS_HUC.nc",sep='')
	p.nc = nc_open(tmp.file)
	huc.lon = ncvar_get(p.nc,'longitude')
	huc.lat = ncvar_get(p.nc,'latitude')
	huc.map = ncvar_get(p.nc,'huc')
	nc_close(p.nc)

	tmp.file = paste(data.dir,"HUC_ID.txt",sep='')
	huc.id = read.table(tmp.file)
	huc.id = as.matrix(huc.id)

	chosen.id = huc.id[match(chosen.hucs,huc.id[,2]),1]
	tmp.map = huc.map/huc.map
	tmp.map[which(is.na(match(huc.map,chosen.id))==FALSE)] = 2
	chosen.area = which(tmp.map == 2, arr.ind=TRUE)
	old.area = chosen.area
	
	# The "domain" variable allows for the approximation of the drainage area when the gage is not at the outlet of the watershed
	if(domain[s] == 0) chosen.area = chosen.area
	if(domain[s] == 1) chosen.area = chosen.area[which(huc.lat[chosen.area[,2]] >= glat & huc.lon[chosen.area[,1]] >= glon),]
	if(domain[s] == 2) chosen.area = chosen.area[which(huc.lat[chosen.area[,2]] <= glat & huc.lon[chosen.area[,1]] >= glon),]
	if(domain[s] == 3) chosen.area = chosen.area[which(huc.lat[chosen.area[,2]] <= glat & huc.lon[chosen.area[,1]] <= glon),]
	if(domain[s] == 4) chosen.area = chosen.area[which(huc.lat[chosen.area[,2]] >= glat & huc.lon[chosen.area[,1]] <= glon),]
		
	## Draw box around the extent of the chosen HUC and apply buffer of 2 grid cells; buffer is set earlier
	min.lon = huc.lon[min(old.area[,1]) - buffer] - 360
	max.lon = huc.lon[max(old.area[,1]) + buffer] - 360
	min.lat = huc.lat[min(old.area[,2]) - buffer]
	max.lat = huc.lat[max(old.area[,2]) + buffer]

	## Calculate basin area in meters-squared
	if(length(chosen.area) > 2) basin.area = sum(cos(pi* huc.lat[chosen.area[,2]]/180)*((110.5*1000/8))^2)
	if(length(chosen.area) == 2) basin.area = sum(cos(pi* huc.lat[chosen.area[2]]/180)*((110.5*1000/8))^2)
	
	## Observed SWE and Precipitaiton at SNOTEL sites - data has been parsed from the observational datasets that are described in the paper
	swe.file = "/Users/andrewbadger/Documents/SARP/BCSD_Predictors/SNOTEL_Sites_APR1_SWE_PredMatched"	
	pre.file = "/Users/andrewbadger/Documents/SARP/BCSD_Predictors/SNOTEL_Sites_APR1_PR_PredMatched"

	swe.data = as.matrix(read.table(swe.file))
	sites.lat = as.numeric(swe.data[,1])
	sites.lon = as.numeric(swe.data[,2])
	sites.ele = as.numeric(swe.data[,3])
	sites.huc = as.numeric(swe.data[,4])
	swe.vals = swe.data[,5:(dim(swe.data)[2])] * 25.4 #mm

	pre.data = as.matrix(read.table(pre.file))
	pre.vals = pre.data[,5:(dim(pre.data)[2])] * 25.4 #mm

	## Observational period is from 1979 to 2017, remiving data before 2005 due to that being projected future data in the BCSD simulations
	swe.years = seq(1979,2017)
	nyears = length(swe.years)
	y2005 = which(swe.years == 2005)
	swe.vals = swe.vals[,1:y2005]
	pre.vals = pre.vals[,1:y2005]
	
	## Count how many valid years we have at each site
	swe.nvalid = rowSums(swe.vals/swe.vals,na.rm=TRUE)
	pre.nvalid = rowSums(pre.vals/pre.vals,na.rm=TRUE)

	## See how many sites we have in our initial bounding box
	chosen.sites = which(sites.lon >= min.lon & sites.lon <= max.lon & sites.lat >= min.lat & sites.lat <= max.lat & swe.nvalid >= time.thresh & pre.nvalid >= time.thresh)
	chosen.ele = sites.ele[chosen.sites]
	nsites = length(chosen.sites)
	
	## If number of sites is less than 5, increase bounding box size until we get 5 observational sites within bounding box
	while(nsites <= 4){
		buffer = buffer + 1
		min.lon = huc.lon[min(old.area[,1]) - buffer] - 360
		max.lon = huc.lon[max(old.area[,1]) + buffer] - 360
		min.lat = huc.lat[min(old.area[,2]) - buffer]
		max.lat = huc.lat[max(old.area[,2]) + buffer]
	
		chosen.sites = which(sites.lon >= min.lon & sites.lon <= max.lon & sites.lat >= min.lat & sites.lat <= max.lat & swe.nvalid >= time.thresh & pre.nvalid >= time.thresh)
		chosen.ele = sites.ele[chosen.sites]
		nsites = length(chosen.sites)
		}
	
	## Write out which observational location are chosen for each watershed and the number of gridcells buffered around the lat-lon bounding box of the watershed	
	new.file = paste(write.dir,sitename,"_CHOSEN_SITES",sep='')
	write.table(chosen.sites, new.file, sep="\t", col.names = F, row.names = F, quote=FALSE)
	
	new.file = paste(write.dir,sitename,"_BUFFER",sep='')
	write.table(buffer, new.file, sep="\t", col.names = F, row.names = F, quote=FALSE)
	
	## Prepare SWE and Precip at chosen observation locations
	chosen.swe = swe.vals[chosen.sites,]
	chosen.pre = pre.vals[chosen.sites,]
	all.swe.valid = colSums(chosen.swe/chosen.swe,na.rm=TRUE)
	all.pre.valid = colSums(chosen.pre/chosen.pre,na.rm=TRUE)

	t.start = max(which(all.swe.valid == nsites)[1],which(all.pre.valid == nsites)[1])
	obs.swe.final = chosen.swe[,t.start:y2005]
	obs.pre.final = chosen.pre[,t.start:y2005]
	years.final = swe.years[t.start:y2005]
	nyears = length(years.final)
	obs.years.final = years.final
	
	## Read in observered streamflow data and convert from daily to desired predictand streamflow
	if (nchar(gage.id) < 8){
  		toAdd = as.character(integer(8-nchar(gage.id)))
  		toAdd = paste(toAdd, collapse="")
  		gage.id = paste(c(toAdd, gage.id), collapse = "")
		}

	data = whatNWISdata(siteNumber=gage.id, service = "dv") #service = dv indicates you want daily values

	#"data" contains a list of all the parameters observed by this gage.
	params = data$parm_cd  #00060 is streamflow, for example.

	#Obtain streamflow data for this gage using the info obtained in "data"
	siteNo = gage.id #site no
	pCode = data$parm_cd #code for streamflow
	start.date = data$begin_date #date of start collection
	end.date = data$end_date #date of end collection
	dis.id = which(pCode == "00060")
	pCode = pCode[dis.id]
	start.date = start.date[dis.id]
	end.date = end.date[dis.id]
	pinto = readNWISdv(siteNumbers = siteNo, parameterCd = pCode, startDate = start.date, endDate = end.date)
	
	datelist.day = pinto[,3]
	ndays = length(datelist.day)
	dates = as.numeric(unlist(strsplit(as.character(datelist.day), "-")))
	dim(dates) = c(3,ndays)
	flow.mon.d = dates[2,]
	flow.year.d = dates[1,]
	
	datelist = seq.Date(as.Date(start.date), as.Date(end.date), "months")
	flow.mon = as.numeric(format(datelist,'%m'))
	flow.year = as.numeric(format(datelist,'%Y'))
	
	flow.values = pinto[,4]*60*60*24*2.29569e-05
	
	flow.final = array(NA, c(nyears))
	for(i in 1:nyears){
		tmp = which(flow.year.d == years.final[i] & flow.mon.d >= 4 & flow.mon.d <= 7)
		flow.final[i] = sum(flow.values[tmp])
		}

	obs.flow.final = flow.final

	## ALL OBSERVATIONAL DATA READY

	## GET MODEL DATA
	swe.dir = "/Volumes/LabShare/SARP/rcp85_swe/"
	swe.list = list.files(swe.dir)
	n.swe.files = length(swe.list)

	run.dir = "/Volumes/LabShare/SARP/rcp85_total_runoff/"
	run.list = list.files(run.dir)
	n.run.files = length(run.list)
	r.flow.list = run.list

	pr.dir = "/Volumes/LabShare/SARP/rcp85_pr/"
	pr.list = list.files(pr.dir)
	n.pr.files = length(pr.list)

	print.names = gsub("_rcp85_r1i1p1_SWE.nc", "", swe.list)
	nmodels = length(swe.list)

	all.swe = array(NA, c(nsites,nyears,nmodels))
	all.pre = array(NA, c(nsites,nyears,nmodels))
	all.flow = array(NA, c(nyears,nmodels))

	for(m in 1:nmodels){
	
		mod.years = seq(1950,2099)
		mod.years.chosen = match(years.final,mod.years)
		mod.years = mod.years[mod.years.chosen]
		
		# READ IN DATA - parsed SWE and Precip BCSD data into text files at each observation location, BCSD data described in manuscript
		tmp.file = paste(snotel.dir,print.names[m],"_SNOTEL_SITES_Apr1_SWE",sep='')
		swe.all = read.table(tmp.file)
		swe.data  = as.matrix(swe.all[,-c(1,2)])
		swe.chosen = swe.data[chosen.sites,]
		swe.chosen = swe.chosen[,mod.years.chosen]
		all.swe[,,m] = swe.chosen
	
		tmp.file = paste(snotel.dir,print.names[m],"_SNOTEL_SITES_ONDJFM_PR",sep='')
		swe.all = read.table(tmp.file)
		swe.data  = as.matrix(swe.all[,-c(1,2)])
		swe.chosen = swe.data[chosen.sites,]
		swe.chosen = swe.chosen[,mod.years.chosen]
		all.pre[,,m] = swe.chosen
	
		## Read in streamflow data from model and sum up across the basin, then convert data to match observer predictand
		file.flow = paste(run.dir,r.flow.list[m],sep='')
		p.nc = nc_open(file.flow)
		tmp.flow = ncvar_get(p.nc,'total_runoff')
		nc_close(p.nc)
	
		ntime = dim(tmp.flow)[3]
		flow.values = array(NA, c(ntime))
		for(i in 1:ntime){
			tmp.mon = tmp.flow[,,i]
			flow.values[i] = sum(tmp.mon[chosen.area])
			if(length(chosen.area) > 2) flow.values[i] = sum(tmp.mon[chosen.area])
			if(length(chosen.area) == 2) flow.values[i] = tmp.mon[chosen.area[1],chosen.area[2]]
			}
		
		#flow.values =  flow.values * basin.area * 0.001 * 0.000810708486253
		if(length(chosen.area) > 2) flow.values = flow.values * basin.area * 0.001 * 0.000810708486253 / dim(chosen.area)[1]
		if(length(chosen.area) == 2) flow.values = flow.values * basin.area * 0.001 * 0.000810708486253
		
		datelist = seq.Date(as.Date("1950/01/01"), as.Date("2099/12/01"), "months")
		mon = as.numeric(format(datelist,'%m'))
		year = as.numeric(format(datelist,'%Y'))
	
		flow.final = array(NA, c(nyears))
		for(i in 1:nyears){
			tmp = which(year == mod.years[i] & mon >= 4 & mon <= 7)
			flow.final[i] = sum(flow.values[tmp])
			}
	
		all.flow[,m] = flow.final
		}

	## Now have both observation and model data for SWE and Q in a format to use for SVD Method
	
	## ALLOCATE COMBINATIONAL RESULTS - create list of all possible combinations and only choose options that use 4 observational sites
	combinations = as.matrix(expand.grid(lapply(numeric(nsites), function(x) c(0,1))))
	goodcom = which(rowSums(combinations) == 4)
	combinations = combinations[goodcom,]
	ncom = dim(combinations)[1]

	## Allocate array for results to be written out
	tot.calc = (nmodels + 1)*ncom
	tot.col = nyears+2
	to.write.flow = array(NA, c(tot.calc,tot.col))
	count = 0
	
	## CALCULATE PREDICTED FLOWS - Looping over models, then combinations, then years (leave one out approach)
	for(i in 1:(nmodels+1)){
		## First iteration uses observations, subsequent iterations use model data
		if(i == 1) tmp.all.swe = obs.swe.final
		if(i == 1) tmp.all.pre = obs.pre.final
		if(i == 1) tmp.all.flow = obs.flow.final
	
		if(i >= 2) tmp.all.swe = all.swe[,,(i-1)]
		if(i >= 2) tmp.all.pre = all.pre[,,(i-1)]
		if(i >= 2) tmp.all.flow = all.flow[,(i-1)]
	
		for(j in 1:ncom){
			id = i - 1
			tmp = combinations[j,]
			chosen.com = which(tmp == 1)
		
			pred.flow = array(NA, c(nyears))
			
			for(k in 1:nyears){
				
				## Get SWE data at desired locations based on combination chosen - Split into calibration and validation sets
				swe.cal = t(tmp.all.swe[chosen.com,-k])
				pre.cal = t(tmp.all.pre[chosen.com,-k])
				flow.cal = tmp.all.flow[-k]
				swe.val = t(tmp.all.swe[chosen.com,k])
				pre.val = t(tmp.all.pre[chosen.com,k])
				flow.val = tmp.all.flow[k]
			
				## Apply column of 1s to the data, this provides the "b" in the y = mx + b equation
				swe.tmp = cbind(1,swe.cal)
				
				## Use SVD to determine the "betas" and then predict flow for year left out
				swe.tmp[which(is.na(swe.tmp)==TRUE)] = 0
				xtx = t(swe.tmp) %*% swe.tmp
				xty = t(swe.tmp) %*% flow.cal
				if(is.positive.definite(xtx) == TRUE){
					xtx.inv = chol2inv(chol(xtx))
					beta = xtx.inv %*% xty
					flow.pred.cal = swe.tmp %*% beta
					dim(flow.cal) = c(length(flow.cal),1)
				
					swe.tmp = cbind(1,swe.val)
					pred.flow[k] = swe.tmp %*% beta
					dim(flow.val) = c(length(flow.val),1)
					}
				}
			dim(pred.flow) = c(1,nyears)
			
			## Save information to determine obs/model, combination and predicted flow
			tmp.write = cbind(id, j, pred.flow)
			count = count + 1
			to.write.flow[count,] = tmp.write
			}
		}
	
	## Write out all useful information - Predicted flow, observerd flow, modeled flow, combinations used		
	new.file = paste(write.dir,sitename,"_Jacknife_Pred_SWE_PR_Flow",sep='')
	write.table(to.write.flow, new.file, sep="\t", col.names = F, row.names = F, quote=FALSE)

	dim(obs.flow.final) = c(1,nyears)
	to.write.obs = cbind(0,obs.flow.final)
	new.file = paste(write.dir,sitename,"_Jacknife_Obs_Flow",sep='')
	write.table(to.write.obs, new.file, sep="\t", col.names = F, row.names = F, quote=FALSE)

	to.write.mod = cbind(seq(1,29),t(all.flow))
	new.file = paste(write.dir,sitename,"_Jacknife_Mod_Flow",sep='')
	write.table(to.write.mod, new.file, sep="\t", col.names = F, row.names = F, quote=FALSE)

	to.write.com = cbind(seq(1,70),combinations)
	new.file = paste(write.dir,sitename,"_Jacknife_Combinations",sep='')
	write.table(to.write.com, new.file, sep="\t", col.names = F, row.names = F, quote=FALSE)
	
	print(paste(sitename," FINISHED! - ",s," of ",nsitenames,sep=''))
	}		
			