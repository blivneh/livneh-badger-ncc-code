####################################################################
#                                                                  #
#  THIS SCRIPT USES BCSD MODEL OUTPUT TO APPLY SVD AND GENERATE    #
#  PREDICTIVE EQUATIONS AT ALL GRID POINTS THAT MEET OUR SELECTION #
#  CRITERIA (see manuscript for details).                          #
#                                                                  #
#  THE PREDICTORS INCLUDE:                                         #
#  APRIL 1 SWE                                                     #
#  APRIL 1 SOIL MOISTURE                                           #
#  ONDJFM ACCUMULATED PRECIPITATION                                #
#  MEAN ONDJFM TEMPERATURE                                         #  
#                                                                  #
#  THE PREDICTAND IS AMJJ TOTAL RUNOFF                             #  
#                                                                  #
#  EACH MODEL PRODUCES THREE OUTPUT FILES:                         #
#  GENERAL OUTPUT                                                  #
#  ETS INFORMATION                                                 #
#  GOODNESS OF FIT METRICS                                         #
#                                                                  #
#  THE OUTPUT IS IN THE FORMAT OF A TEXT FILE THAT LISTS EACH      #
#  POINT WITH THEIR X-Y COMPONENTS FOR PLOTTING/GRIDDING PURPOSES  #
#                                                                  #
####################################################################

rm(list=ls())

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

lon.thresh = -100

## Read in list of model data
swe.dir = "/Volumes/LabShare/SARP/rcp85_swe/"
swe.list = list.files(swe.dir)
n.swe.files = length(swe.list)

run.dir = "/Volumes/LabShare/SARP/rcp85_total_runoff/"
run.list = list.files(run.dir)
n.run.files = length(run.list)

pr.dir = "/Volumes/LabShare/SARP/rcp85_pr/"
pr.list = list.files(pr.dir)
n.pr.files = length(pr.list)

sm.dir = "/Volumes/LabShare/SARP/rcp85_smc/"
sm.list = list.files(sm.dir)
n.sm.files = length(sm.list)

tas.dir = "/Volumes/LabShare/SARP/rcp85_tas/"
tas.list = list.files(tas.dir)
n.tas.files = length(tas.list)

print.names = gsub("_rcp85_r1i1p1_SWE.nc", "", swe.list)

## Set output directory
output.dir = "/Volumes/LabShare/SARP/BCSD_Point_Predictions/"

## Loop over each input file - lists of each file were verified to make sure the models alligned properly in the list for each variable
for(f in 1:n.swe.files){

	## Read in the data
	file.swe = paste(swe.dir,swe.list[f],sep='')
	p.nc = nc_open(file.swe)
	swe.data = ncvar_get(p.nc,'swe') #mm
	vic.lon = ncvar_get(p.nc,'longitude')
	vic.lat = ncvar_get(p.nc,'latitude')
	nc_close(p.nc)
	nlon = length(vic.lon)
	nlat = length(vic.lat)
	
	file.run = paste(run.dir,run.list[f],sep='')
	p.nc = nc_open(file.run)
	run.data = ncvar_get(p.nc,'total_runoff') #mm/month
	nc_close(p.nc)
	
	file.pr = paste(pr.dir,pr.list[f],sep='')
	p.nc = nc_open(file.pr)
	pr.data = ncvar_get(p.nc,'pr') #mm/month
	nc_close(p.nc)

	file.sm = paste(sm.dir,sm.list[f],sep='')
	p.nc = nc_open(file.sm)
	sm.data = ncvar_get(p.nc,'smc') #mm
	nc_close(p.nc)
	
	file.tas = paste(tas.dir,tas.list[f],sep='')
	p.nc = nc_open(file.tas)
	tas.data = ncvar_get(p.nc,'tas') #C
	nc_close(p.nc)
	
	## Create dates for all files
	datelist = seq.Date(as.Date("1950/1/1"), as.Date("2099/12/1"), "months")
	vic.mon = as.numeric(format(datelist,'%m'))
	vic.year = as.numeric(format(datelist,'%Y'))
	vic.years = seq(1950,2099)
	vic.nyears = length(vic.years)
	
	vic.oct = which(vic.mon == 10)
	vic.nov = which(vic.mon == 11)
	vic.dec = which(vic.mon == 12)
	vic.jan = which(vic.mon == 1)
	vic.feb = which(vic.mon == 2)
	vic.mar = which(vic.mon == 3)
	vic.apr = which(vic.mon == 4)
	vic.may = which(vic.mon == 5)
	vic.june = which(vic.mon == 6)
	vic.july = which(vic.mon == 7)
	vic.hist = which(vic.years >=1980 & vic.years <= 2009)
	
	## Get model data and parse into correct predictors and predictands
	swe.apr = swe.data[,,vic.apr]
	swe.aprhist = swe.apr[,,vic.hist]
	
	run.amjj = run.data[,,vic.apr] + run.data[,,vic.may] + run.data[,,vic.june] + run.data[,,vic.july]
	
	pr.ond = pr.data[,,vic.oct] + pr.data[,,vic.nov] + pr.data[,,vic.dec]
	pr.jfm = pr.data[,,vic.jan] + pr.data[,,vic.feb] + pr.data[,,vic.mar]
	pr.tot = pr.ond[,,(vic.hist-1)] + pr.jfm[,,vic.hist]
	pr.tot.all = pr.ond[,,1:(vic.nyears-1)] + pr.jfm[,,2:(vic.nyears)]
	
	tas.ond = tas.data[,,vic.oct] + tas.data[,,vic.nov] + tas.data[,,vic.dec]
	tas.jfm = tas.data[,,vic.jan] + tas.data[,,vic.feb] + tas.data[,,vic.mar]
	tas.ondjfm = (tas.ond[,,1:(vic.nyears-1)] + tas.jfm[,,2:(vic.nyears)]) / 6
	
	sm.apr = sm.data[,,vic.apr]
	
	## Find points that match criteria for analysis, see manuscript for details
	swe.clim = apply(swe.aprhist,c(1,2),mean)
	pr.clim = apply(pr.tot,c(1,2),mean)
	swe.rat = swe.clim/pr.clim
	swe.100 = which(swe.clim >= 10 & swe.rat >= 0.25, arr.ind=TRUE)
	swe.100 = swe.100[which(swe.100[,1] <= which.min(abs(vic.lon - lon.thresh))),]
	npoints = dim(swe.100)[1]
	
	## Create list of calibration/validation start/end times
	pred.start = seq(1950,2065, by=5)
	pred.end = seq(1979,2094, by=5)
	val.start = seq(1980,2095, by=5)
	val.end = seq(1984,2099, by=5)
	nslices = length(val.end)
	
	mod.flow = array(NA, c(vic.nyears))
	pred.flow = array(NA, c(vic.nyears))
	
	## Create files to write output too and delete old versions of the file
	output.file = paste(output.dir,print.names[f],"_POINT_STATS_KITCHEN_SINK",sep='')
	output.file2 = paste(output.dir,print.names[f],"_THREAT_SCORE_KITCHEN_SINK",sep='')
	output.file3 = paste(output.dir,print.names[f],"_POINT_STATS_KITCHEN_SINK_Goodness_Of_Fit",sep='')

	if(file.exists(output.file)==TRUE) file.remove(output.file)
	if(file.exists(output.file2)==TRUE) file.remove(output.file2)
	if(file.exists(output.file3)==TRUE) file.remove(output.file3)
	
	## Loop over all the points that met selection criteria
	for(i in 1:npoints){
	
		## Get predictor/predictand information from valid point
		tmp.swe = swe.apr[swe.100[i,1],swe.100[i,2],]
		tmp.run = run.amjj[swe.100[i,1],swe.100[i,2],]
		tmp.pr = array(NA, c(vic.nyears))
		tmp.pr[2:vic.nyears] = pr.tot.all[swe.100[i,1],swe.100[i,2],]
		tmp.tas = array(NA, c(vic.nyears))
		tmp.tas[2:vic.nyears] = tas.ondjfm[swe.100[i,1],swe.100[i,2],]
		tmp.sm = sm.apr[swe.100[i,1],swe.100[i,2],]
		
		## Loop over the time slices (30 year validation + 5 year calibration)
		for(k in 1:nslices){
			cal.s = which(vic.years == pred.start[k])
			cal.e = which(vic.years == pred.end[k])
			val.s = which(vic.years == val.start[k])
			val.e = which(vic.years == val.end[k])
			
			## Use SVD to determine the "betas" and then predict flow for the calibration period
			## Gather predictor information in x.tmp
			## Following SVD to create xtx and xty matricies
			## Invert xtx matrix
			## Multiply inverted matrix by xty to get betas for regression
			## Apply betas to calibration period predictors to get flow (predictand)
			## Save modeled flow over the same calibration period for comparison purposes

			x.tmp = cbind(1,tmp.swe[cal.s:cal.e],tmp.pr[cal.s:cal.e],tmp.tas[cal.s:cal.e],tmp.sm[cal.s:cal.e])	
			xtx = t(x.tmp) %*% x.tmp
			xty = t(x.tmp) %*% tmp.run[cal.s:cal.e]
			if(is.positive.definite(xtx) == TRUE){
				xtx.inv = chol2inv(chol(xtx))
				beta = xtx.inv %*% xty
				x.tmp = cbind(1,tmp.swe[val.s:val.e],tmp.pr[val.s:val.e],tmp.tas[val.s:val.e],tmp.sm[val.s:val.e])	
				pred.flow[val.s:val.e] = x.tmp %*% beta
				mod.flow[val.s:val.e] = tmp.run[val.s:val.e]
				}
			}
		
		## Time periods for analysis
		historical = which(vic.years >=1980 & vic.years <= 2009)
		midcent = which(vic.years >=2036 & vic.years <= 2065)
		latecent = which(vic.years >=2070 & vic.years <= 2099)
		allyears = which(vic.years >=1980 & vic.years <= 2099)
		
		## ETS calculation
		## Allocate arrays for hits, misses, false alarms, correct rejections and ETS value
		drought.ets = array(NA, c(4)) #Hist, mid, late, all
		drought.hits = array(NA, c(4)) #Hist, mid, late, all
		drought.misses = array(NA, c(4)) #Hist, mid, late, all
		drought.false = array(NA, c(4)) #Hist, mid, late, all
		drought.rejections = array(NA, c(4)) #Hist, mid, late, all
		
		## Loop over time periods to calculate ETS
		for(d in 1:4){
			if(d == 1) chosen = historical
			if(d == 2) chosen = midcent
			if(d == 3) chosen = latecent
			if(d == 4) chosen = allyears
			
			## Calculate drough threshold of 20th percentile for modeled flow
			drought.thresh = quantile(mod.flow[chosen], c(0.2),na.rm=TRUE) 
			
			## Fill out contingency table components
			tmp.hits = length(which(mod.flow[chosen] <= drought.thresh & pred.flow[chosen] <= drought.thresh))
			tmp.misses = length(which(mod.flow[chosen] <= drought.thresh & pred.flow[chosen] > drought.thresh))
			tmp.false = length(which(mod.flow[chosen] > drought.thresh & pred.flow[chosen] <= drought.thresh))
			tmp.rejections = length(which(mod.flow[chosen] > drought.thresh & pred.flow[chosen] > drought.thresh))
			tmp.tot = tmp.hits + tmp.misses + tmp.false + tmp.rejections
			
			## Calulate "Hits Random" for ETS equation
			hitsr = (tmp.hits + tmp.false)*(tmp.hits + tmp.misses)/tmp.tot
			
			## Determine ETS value
			drought.ets[d] = (tmp.hits - hitsr)/(tmp.hits + tmp.misses + tmp.false - hitsr)
			
			## Save ETS components
			drought.hits[d] = tmp.hits
			drought.misses[d] = tmp.misses
			drought.false[d] = tmp.false
			drought.rejections[d] = tmp.rejections
			}
		
		## Write out results in text files, appending each file as we loop over points
		## Output File - General Output
		nstats = 29
		to.write = array(NA, c(nstats))
		
		to.write[1] = round(swe.100[i,1],digits=0) #x-point, lon
		to.write[2] = round(swe.100[i,2],digits=0) #y-point, lat
		
		chosen = historical
		tmp.mod = mod.flow[chosen]
		tmp.pred = pred.flow[chosen]
		swe.mod = tmp.swe[chosen]
		dim(tmp.mod) = c(length(tmp.mod),1)
		dim(tmp.pred) = c(length(tmp.pred),1)
		to.write[3] = round(NSE(tmp.pred,tmp.mod),digits=3) # Historical - NSE
		to.write[4] = round(cor(tmp.pred,tmp.mod),digits=3) # Historical - Correlation
		to.write[5]	= round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3)  # Historical - Percent Bias (Predicted - Observed)
		to.write[6] = round(sd(tmp.pred)/sd(tmp.mod),digits=3) # Historical - Ratio of SD (pred/mod)
		to.write[7] = round(sqrt(mse(tmp.pred,tmp.mod)),digits=3) # Histortical - RMSE
		to.write[8] = round(drought.ets[1],digits=3) # Historical - Drought ETS
		to.write[9] = round(mean(tmp.mod),digits=3) # Historical - Mean Flow
		to.write[10] = round(mean(swe.mod),digits=3) # Historical - Mean SWE
		to.write[11] = round(length(which(swe.mod == 0)),digits=0) # Historical - Number of no snow years on April 1
		
		chosen = midcent
		tmp.mod = mod.flow[chosen]
		tmp.pred = pred.flow[chosen]
		swe.mod = tmp.swe[chosen]
		dim(tmp.mod) = c(length(tmp.mod),1)
		dim(tmp.pred) = c(length(tmp.pred),1)
		to.write[12] = round(NSE(tmp.pred,tmp.mod),digits=3) # Mid Century - NSE
		to.write[13] = round(cor(tmp.pred,tmp.mod),digits=3) # Mid Century - Correlation
		to.write[14]	= round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3)  # Mid Century - Percent Bias (Predicted - Observed)
		to.write[15] = round(sd(tmp.pred)/sd(tmp.mod),digits=3) # Mid Century - Ratio of SD (pred/mod)
		to.write[16] = round(sqrt(mse(tmp.pred,tmp.mod)),digits=3) # Mid Century - RMSE
		to.write[17] = round(drought.ets[2],digits=3) # Mid Century - Drought ETS
		to.write[18] = round(mean(tmp.mod),digits=3) # Mid Century - Mean Flow
		to.write[19] = round(mean(swe.mod),digits=3) # Mid Century - Mean SWE
		to.write[20] = round(length(which(swe.mod == 0)),digits=0) # Mid Century - Number of no snow years on April 1
		
		chosen = latecent
		tmp.mod = mod.flow[chosen]
		tmp.pred = pred.flow[chosen]
		swe.mod = tmp.swe[chosen]
		dim(tmp.mod) = c(length(tmp.mod),1)
		dim(tmp.pred) = c(length(tmp.pred),1)
		to.write[21] = round(NSE(tmp.pred,tmp.mod),digits=3) # End Century - NSE
		to.write[22] = round(cor(tmp.pred,tmp.mod),digits=3) # End Century - Correlation
		to.write[23]	= round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3)  # End Century - Percent Bias (Predicted - Observed)
		to.write[24] = round(sd(tmp.pred)/sd(tmp.mod),digits=3) # End Century - Ratio of SD (pred/mod)
		to.write[25] = round(sqrt(mse(tmp.pred,tmp.mod)),digits=3) # End Century - RMSE
		to.write[26] = round(drought.ets[3],digits=3) # End Century - Drought ETS
		to.write[27] = round(mean(tmp.mod),digits=3) # End Century - Mean Flow
		to.write[28] = round(mean(swe.mod),digits=3) # End Century - Mean SWE
		to.write[29] = round(length(which(swe.mod == 0)),digits=0) # End Century - Number of no snow years on April 1
		
		write(to.write,output.file, ncolumns=nstats, append=TRUE)
		
		## Output File 2 - ETS Information
		
		nstats = 17
		to.write = array(NA, c(nstats))
		
		to.write[1] = round(swe.100[i,1],digits=0) #x-point, lon
		to.write[2] = round(swe.100[i,2],digits=0) #y-point, lat
		
		to.write[3] = round(drought.hits[1],digits=3) # Historical - Drought ETS Hits
		to.write[4] = round(drought.misses[1],digits=3) # Historical - Drought ETS Misses
		to.write[5] = round(drought.false[1],digits=3) # Historical - Drought ETS False Alarms
		to.write[6] = round(drought.rejections[1],digits=3) # Historical - Drought ETS Correct Rejection
		to.write[7] = round(drought.ets[1],digits=3) # Historical - Drought ETS
		
		to.write[8] = round(drought.hits[2],digits=3) # Mid Century - Drought ETS Hits
		to.write[9] = round(drought.misses[2],digits=3) # Mid Century - Drought ETS Misses
		to.write[10] = round(drought.false[2],digits=3) # Mid Century - Drought ETS False Alarms
		to.write[11] = round(drought.rejections[2],digits=3) # Mid Century - Drought ETS
		to.write[12] = round(drought.ets[2],digits=3) # Mid Century - Drought ETS
		
		to.write[13] = round(drought.hits[3],digits=3) # End Century - Drought ETS Hits
		to.write[14] = round(drought.misses[3],digits=3) # End Century - Drought ETS Misses
		to.write[15] = round(drought.false[3],digits=3) # End Century - Drought ETS False Alarms
		to.write[16] = round(drought.rejections[3],digits=3) # End Century - Drought ETS Correct Rejection
		to.write[17] = round(drought.ets[3],digits=3) # End Century - Drought ETS 
		
		write(to.write,output.file2, ncolumns=nstats, append=TRUE)
		
		## Output File 3 - Goodness of Fit Metrics
		nstats = 14
		to.write = array(NA, c(nstats))
		
		to.write[1] = round(swe.100[i,1],digits=0) #x-point, lon
		to.write[2] = round(swe.100[i,2],digits=0) #y-point, lat
		
		## Historical Goodness of Fit
		chosen = historical
		tmp.mod = mod.flow[chosen]
		tmp.pred = pred.flow[chosen]
		swe.mod = tmp.swe[chosen]
		dim(tmp.mod) = c(length(tmp.mod),1)
		dim(tmp.pred) = c(length(tmp.pred),1)
		to.write[3] = round(sqrt(mse(tmp.pred,tmp.mod))/mean(tmp.mod),digits=3) # RRMSE (RMSE normalized by mean flow)
		to.write[4] = round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3) # Percent Bias
		to.write[5] = round((mean(tmp.pred) - mean(tmp.mod)),digits=3) # Bias
		to.write[6]	= round(cor(tmp.mod,tmp.pred)**2,digits=3)  # R-squared
		
		## Mid Century Goodness of Fit
		chosen = midcent
		tmp.mod = mod.flow[chosen]
		tmp.pred = pred.flow[chosen]
		swe.mod = tmp.swe[chosen]
		dim(tmp.mod) = c(length(tmp.mod),1)
		dim(tmp.pred) = c(length(tmp.pred),1)
		to.write[7] = round(sqrt(mse(tmp.pred,tmp.mod))/mean(tmp.mod),digits=3) # RRMSE (RMSE normalized by mean flow)
		to.write[8] = round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3) # Percent Bias
		to.write[9] = round((mean(tmp.pred) - mean(tmp.mod)),digits=3) # Bias
		to.write[10] = round(cor(tmp.mod,tmp.pred)**2,digits=3)  # R-squared
		
		## End of Century Goodness of Fit
		chosen = latecent
		tmp.mod = mod.flow[chosen]
		tmp.pred = pred.flow[chosen]
		swe.mod = tmp.swe[chosen]
		dim(tmp.mod) = c(length(tmp.mod),1)
		dim(tmp.pred) = c(length(tmp.pred),1)
		to.write[11] = round(sqrt(mse(tmp.pred,tmp.mod))/mean(tmp.mod),digits=3) # RRMSE (RMSE normalized by mean flow)
		to.write[12] = round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3) # Percent Bias
		to.write[13] = round((mean(tmp.pred) - mean(tmp.mod)),digits=3) # Bias
		to.write[14] = round(cor(tmp.mod,tmp.pred)**2,digits=3)  # R-squared
		
		write(to.write,output.file3, ncolumns=nstats, append=TRUE)
		}
	print(print.names[f])
	}