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
run.list = run.list[-11]
n.run.files = length(run.list)

pr.dir = "/Volumes/LabShare/SARP/rcp85_pr/"
pr.list = list.files(pr.dir)
pr.list = pr.list[-11]
pr.list = pr.list[-11]
n.pr.files = length(pr.list)

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
		
	output.file = paste(output.dir,print.names[f],"_POINT_STATS_revisions",sep='')
	if(file.exists(output.file)==TRUE) file.remove(output.file)
	for(i in 1:npoints){
		tmp.swe = swe.apr[swe.100[i,1],swe.100[i,2],]
		tmp.run = run.amjj[swe.100[i,1],swe.100[i,2],]
		
		for(k in 1:nslices){
			cal.s = which(vic.years == pred.start[k])
			cal.e = which(vic.years == pred.end[k])
			val.s = which(vic.years == val.start[k])
			val.e = which(vic.years == val.end[k])
			
			x.tmp = cbind(1,tmp.swe[cal.s:cal.e])	
			xtx = t(x.tmp) %*% x.tmp
			xty = t(x.tmp) %*% tmp.run[cal.s:cal.e]
			if(is.positive.definite(xtx) == TRUE){
				xtx.inv = chol2inv(chol(xtx))
				beta = xtx.inv %*% xty
				x.tmp = cbind(1,tmp.swe[val.s:val.e])
				pred.flow[val.s:val.e] = x.tmp %*% beta
				mod.flow[val.s:val.e] = tmp.run[val.s:val.e]
				}
			}
		historical = which(vic.years >=1980 & vic.years <= 2009)
		midcent = which(vic.years >=2036 & vic.years <= 2065)
		latecent = which(vic.years >=2070 & vic.years <= 2099)
		allyears = which(vic.years >=1980 & vic.years <= 2099)
		
		nstats = 14
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
		to.write[3] = round(sqrt(mse(tmp.pred,tmp.mod))/mean(tmp.mod),digits=3) # RRMSE
		to.write[4] = round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3) # Percent Bias
		to.write[5] = round((mean(tmp.pred) - mean(tmp.mod)),digits=3) # Bias
		to.write[6]	= round(cor(tmp.mod,tmp.pred)**2,digits=3)  # R-squared
		
		chosen = midcent
		tmp.mod = mod.flow[chosen]
		tmp.pred = pred.flow[chosen]
		swe.mod = tmp.swe[chosen]
		dim(tmp.mod) = c(length(tmp.mod),1)
		dim(tmp.pred) = c(length(tmp.pred),1)
		to.write[7] = round(sqrt(mse(tmp.pred,tmp.mod))/mean(tmp.mod),digits=3) # RRMSE
		to.write[8] = round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3) # Percent Bias
		to.write[9] = round((mean(tmp.pred) - mean(tmp.mod)),digits=3) # Bias
		to.write[10] = round(cor(tmp.mod,tmp.pred)**2,digits=3)  # R-squared
		
		chosen = latecent
		tmp.mod = mod.flow[chosen]
		tmp.pred = pred.flow[chosen]
		swe.mod = tmp.swe[chosen]
		dim(tmp.mod) = c(length(tmp.mod),1)
		dim(tmp.pred) = c(length(tmp.pred),1)
		to.write[11] = round(sqrt(mse(tmp.pred,tmp.mod))/mean(tmp.mod),digits=3) # RRMSE
		to.write[12] = round(100*(mean(tmp.pred) - mean(tmp.mod))/mean(tmp.mod),digits=3) # Percent Bias
		to.write[13] = round((mean(tmp.pred) - mean(tmp.mod)),digits=3) # Bias
		to.write[14] = round(cor(tmp.mod,tmp.pred)**2,digits=3)  # R-squared
		
		write(to.write,output.file, ncolumns=nstats, append=TRUE)
		
		}
	print(print.names[f])
	}