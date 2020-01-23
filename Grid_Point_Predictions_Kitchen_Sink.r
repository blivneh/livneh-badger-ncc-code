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

swe.dir = "/Volumes/LabShare/SARP/rcp85_swe/"
swe.list = list.files(swe.dir)
swe.list = swe.list[-11]
swe.list = swe.list[-11]
n.swe.files = length(swe.list)

run.dir = "/Volumes/LabShare/SARP/rcp85_total_runoff/"
run.list = list.files(run.dir)
run.list = run.list[-11]
n.run.files = length(run.list)

pr.dir = "/Volumes/LabShare/SARP/rcp85_pr/"
pr.list = list.files(pr.dir)
pr.list = pr.list[-11]
pr.list = pr.list[-11]
n.pr.files = length(pr.list)

sm.dir = "/Volumes/LabShare/SARP/rcp85_smc/"
sm.list = list.files(sm.dir)
sm.list = sm.list[-11]
n.sm.files = length(sm.list)

tas.dir = "/Volumes/LabShare/SARP/rcp85_tas/"
tas.list = list.files(tas.dir)
tas.list = tas.list[-11]
tas.list = tas.list[-11]
n.tas.files = length(tas.list)

tasmax.dir = "/Volumes/LabShare/SARP/rcp85_tasmax/"
tasmax.list = list.files(tasmax.dir)
tasmax.list = tasmax.list[-11]
n.tasmax.files = length(tasmax.list)

tasmin.dir = "/Volumes/LabShare/SARP/rcp85_tasmin/"
tasmin.list = list.files(tasmin.dir)
tasmin.list = tasmin.list[-11]
n.tasmin.files = length(tasmin.list)

output.dir = "/Volumes/LabShare/SARP/BCSD_Point_Predictions/"

print.names = gsub("_rcp85_r1i1p1_SWE.nc", "", swe.list)

for(f in 1:n.swe.files){
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
	
	file.tasmax = paste(tasmax.dir,tasmax.list[f],sep='')
	p.nc = nc_open(file.tasmax)
	tasmax.data = ncvar_get(p.nc,'tasmax') #C
	nc_close(p.nc)
	
	file.tasmin = paste(tasmin.dir,tasmin.list[f],sep='')
	p.nc = nc_open(file.tasmin)
	tasmin.data = ncvar_get(p.nc,'tasmin') #C
	nc_close(p.nc)
	
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
	#vic.hist = which(vic.years <= 2006)
	vic.hist = which(vic.years >=1980 & vic.years <= 2009)
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
	
	tasmax.ond = tasmax.data[,,vic.oct] + tasmax.data[,,vic.nov] + tasmax.data[,,vic.dec]
	tasmax.jfm = tasmax.data[,,vic.jan] + tasmax.data[,,vic.feb] + tasmax.data[,,vic.mar]
	tasmax.ondjfm = (tasmax.ond[,,1:(vic.nyears-1)] + tasmax.jfm[,,2:(vic.nyears)]) / 6
	
	tasmin.ond = tasmin.data[,,vic.oct] + tasmin.data[,,vic.nov] + tasmin.data[,,vic.dec]
	tasmin.jfm = tasmin.data[,,vic.jan] + tasmin.data[,,vic.feb] + tasmin.data[,,vic.mar]
	tasmin.ondjfm = (tasmin.ond[,,1:(vic.nyears-1)] + tasmin.jfm[,,2:(vic.nyears)]) / 6
	
	sm.apr = sm.data[,,vic.apr]
	
	swe.clim = apply(swe.aprhist,c(1,2),mean)
	pr.clim = apply(pr.tot,c(1,2),mean)
	swe.rat = swe.clim/pr.clim
	swe.100 = which(swe.clim >= 10 & swe.rat >= 0.25, arr.ind=TRUE)
	swe.100 = swe.100[which(swe.100[,1] <= which.min(abs(vic.lon - lon.thresh))),]
	npoints = dim(swe.100)[1]
	
	pred.start = seq(1950,2065, by=5)
	pred.end = seq(1979,2094, by=5)
	val.start = seq(1980,2095, by=5)
	val.end = seq(1984,2099, by=5)
	nslices = length(val.end)
	
	mod.flow = array(NA, c(vic.nyears))
	pred.flow = array(NA, c(vic.nyears))
		
	output.file = paste(output.dir,print.names[f],"_POINT_STATS_KITCHEN_SINK",sep='')
	output.file2 = paste(output.dir,print.names[f],"_THREAT_SCORE_KITCHEN_SINK",sep='')
	if(file.exists(output.file)==TRUE) file.remove(output.file)
	if(file.exists(output.file2)==TRUE) file.remove(output.file2)
	for(i in 1:npoints){
		tmp.swe = swe.apr[swe.100[i,1],swe.100[i,2],]
		tmp.run = run.amjj[swe.100[i,1],swe.100[i,2],]
		tmp.pr = array(NA, c(vic.nyears))
		tmp.pr[2:vic.nyears] = pr.tot.all[swe.100[i,1],swe.100[i,2],]
		tmp.tas = array(NA, c(vic.nyears))
		tmp.tas[2:vic.nyears] = tas.ondjfm[swe.100[i,1],swe.100[i,2],]
		tmp.tasmax = array(NA, c(vic.nyears))
		tmp.tasmax[2:vic.nyears] = tasmax.ondjfm[swe.100[i,1],swe.100[i,2],]
		tmp.tasmin = array(NA, c(vic.nyears))
		tmp.tasmin[2:vic.nyears] = tasmin.ondjfm[swe.100[i,1],swe.100[i,2],]
		tmp.sm = sm.apr[swe.100[i,1],swe.100[i,2],]
		
		print(paste("f = ",f,"; i = ",i,sep=''))
		
		for(k in 1:nslices){
			cal.s = which(vic.years == pred.start[k])
			if(k==1) cal.s = which(vic.years == (pred.start[k]+1))
			cal.e = which(vic.years == pred.end[k])
			val.s = which(vic.years == val.start[k])
			val.e = which(vic.years == val.end[k])
			
			x.tmp = cbind(1,tmp.swe[cal.s:cal.e],tmp.pr[cal.s:cal.e],tmp.tas[cal.s:cal.e],tmp.sm[cal.s:cal.e])	
			xtx = t(x.tmp) %*% x.tmp
			xty = t(x.tmp) %*% tmp.run[cal.s:cal.e]
			if(length(which(diag(xtx) == 0)) >= 1) xtx[which(diag(xtx) == 0), which(diag(xtx) == 0)] = 0.0000001
			if(is.positive.definite(xtx) == TRUE){
				if(length(which(diag(xtx) == 0)) >= 1) xtx[which(diag(xtx) == 0), which(diag(xtx) == 0)] = 0.0000001
				xtx.inv = chol2inv(chol(xtx))
				beta = xtx.inv %*% xty
				x.tmp = cbind(1,tmp.swe[val.s:val.e],tmp.pr[val.s:val.e],tmp.tas[val.s:val.e],tmp.sm[val.s:val.e])	
				pred.flow[val.s:val.e] = x.tmp %*% beta
				mod.flow[val.s:val.e] = tmp.run[val.s:val.e]
				}
			}
		historical = which(vic.years >=1980 & vic.years <= 2009)
		midcent = which(vic.years >=2036 & vic.years <= 2065)
		latecent = which(vic.years >=2070 & vic.years <= 2099)
		allyears = which(vic.years >=1980 & vic.years <= 2099)
		
		drought.ets = array(NA, c(4)) #Hist, mid, late, all
		drought.hits = array(NA, c(4)) #Hist, mid, late, all
		drought.misses = array(NA, c(4)) #Hist, mid, late, all
		drought.false = array(NA, c(4)) #Hist, mid, late, all
		drought.rejections = array(NA, c(4)) #Hist, mid, late, all
		for(d in 1:4){
			if(d == 1) chosen = historical
			if(d == 2) chosen = midcent
			if(d == 3) chosen = latecent
			if(d == 4) chosen = allyears
			
			drought.thresh = quantile(mod.flow[chosen], c(0.2),na.rm=TRUE) 
			tmp.hits = length(which(mod.flow[chosen] <= drought.thresh & pred.flow[chosen] <= drought.thresh))
			tmp.misses = length(which(mod.flow[chosen] <= drought.thresh & pred.flow[chosen] > drought.thresh))
			tmp.false = length(which(mod.flow[chosen] > drought.thresh & pred.flow[chosen] <= drought.thresh))
			tmp.rejections = length(which(mod.flow[chosen] > drought.thresh & pred.flow[chosen] > drought.thresh))
			tmp.tot = tmp.hits + tmp.misses + tmp.false + tmp.rejections
			hitsr = (tmp.hits + tmp.false)*(tmp.hits + tmp.misses)/tmp.tot
			drought.ets[d] = (tmp.hits - hitsr)/(tmp.hits + tmp.misses + tmp.false - hitsr)
			drought.hits[d] = tmp.hits
			drought.misses[d] = tmp.misses
			drought.false[d] = tmp.false
			drought.rejections[d] = tmp.rejections
			
			}
		nstats = 29
		# nstats = 38 # If Writing All Yeats Stats
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
		
		# chosen = allyears
		# tmp.mod = mod.flow[chosen]
		# tmp.pred = pred.flow[chosen]
		# swe.mod = tmp.swe[chosen]
		# dim(tmp.mod) = c(length(tmp.mod),1)
		# dim(tmp.pred) = c(length(tmp.pred),1)
		# to.write[30] = round(NSE(tmp.pred,tmp.mod),digits=3) # All Years - NSE
		# to.write[31] = round(cor(tmp.pred,tmp.mod),digits=3) # All Years - Correlation
		# to.write[32]	= round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3)  # All Years - Percent Bias (Predicted - Observed)
		# to.write[33] = round(sd(tmp.pred)/sd(tmp.mod),digits=3) # All Years - Ratio of SD (pred/mod)
		# to.write[34] = round(sqrt(mse(tmp.pred,tmp.mod)),digits=3) # All Years - RMSE
		# to.write[35] = round(drought.ets[4],digits=3) # All Years - Drought ETS
		# to.write[36] = round(mean(tmp.mod),digits=3) # All Years - Mean Flow
		# to.write[37] = round(mean(swe.mod),digits=3) # All Years - Mean SWE
		# to.write[38] = round(length(which(swe.mod == 0)),digits=0) # All Years - Number of no snow years on April 1
		
		write(to.write,output.file, ncolumns=nstats, append=TRUE)
		
		nstats = 17
		# nstats = 38 # If Writing All Yeats Stats
		to.write = array(NA, c(nstats))
		
		to.write[1] = round(swe.100[i,1],digits=0) #x-point, lon
		to.write[2] = round(swe.100[i,2],digits=0) #y-point, lat
		
		to.write[3] = round(drought.hits[1],digits=3) # Historical - Drought ETS
		to.write[4] = round(drought.misses[1],digits=3) # Historical - Drought ETS
		to.write[5] = round(drought.false[1],digits=3) # Historical - Drought ETS
		to.write[6] = round(drought.rejections[1],digits=3) # Historical - Drought ETS
		to.write[7] = round(drought.ets[1],digits=3) # Historical - Drought ETS
		
		to.write[8] = round(drought.hits[2],digits=3) # Historical - Drought ETS
		to.write[9] = round(drought.misses[2],digits=3) # Historical - Drought ETS
		to.write[10] = round(drought.false[2],digits=3) # Historical - Drought ETS
		to.write[11] = round(drought.rejections[2],digits=3) # Historical - Drought ETS
		to.write[12] = round(drought.ets[2],digits=3) # Historical - Drought ETS
		
		to.write[13] = round(drought.hits[3],digits=3) # Historical - Drought ETS
		to.write[14] = round(drought.misses[3],digits=3) # Historical - Drought ETS
		to.write[15] = round(drought.false[3],digits=3) # Historical - Drought ETS
		to.write[16] = round(drought.rejections[3],digits=3) # Historical - Drought ETS
		to.write[17] = round(drought.ets[3],digits=3) # Historical - Drought ETS
		
		
		write(to.write,output.file2, ncolumns=nstats, append=TRUE)
		
		}
	print(print.names[f])
	}