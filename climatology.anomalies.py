#!/usr/bin/python
import numpy as np
from scipy.stats import norm
import netCDF4 as nc
import math
from datetime import timedelta, date
#=======================================================================================
"""
Program that calculates anomalies by removing the smoothed mean annual cycle and 
trends of 3-D [time, latitude, longitude] data

The mean annual cycle is smoothed using the first three harmonics of the mean annual
cycle

The input variables below need to be defined by the user. Either retrieved from the
dataset or created by the user.

TRENDS: if you want to analyze trends, comment the block under 'Removing Trends'

The program saves anomalies as NetCDF files.

Input:
  tot    : total number of data in a single year (can be retrieved from dataset)
  ntot   : total number of data points in time (can be retrieved from dataset)
  nlat   : total number of data points in latitude (can be retrieved from dataset)
  nlon   : total number of data points in longitude (can be retrieved from dataset)
  lats   : array with latitude values (can be retrieved from dataset)
  lons   : array with longitude values (can be retrieved from dataset)
  missval: value for missing values (can be retrieved from dataset)
  dper   : minimum percentage [0.,100.] of missing data that can be tolerated
  pathout: path of the directory to where the output will be saved
  prefix : The first part of the name of output files.
Output:
  A set of NetCDF files containing gridded values of anomalies.

"""
tot=365
missval = -999.00     # can be defined after reading the data
dper    =    5.       # set this parameter to zero if you don't want to mask any regions
pathout='/home/login/subsets/' # set this path to your directory of choice
prefix='anom.variable.'        # choice of name. For example: anom.Temp.

#====================================== Reading Data ===================================
"""
This block needs to be edited by the user. Read your input dataset and retrieve (or create)
the Input variables described above.
"""

ntot=
nlat=
nlon=

""" Example of how to create latitude and longitude arrays"""
lats=np.arange(-90,90,0.5)
lons=np.arange(0.,360.,0.5)

""" Example of how to create the dataset array"""
prec=np.zeros((ntot,nlat,nlon))

"""Example of how to find the value for missing data. Usually a very large nagative or positive number. For a very large positive number, change .min() to .max() """
missval=prec.min()   # missing value

#=======================================================================================
#---------------------------------------------------------------------------------------
"""    No Further Editing is Required from this point on (see comment about trends)  """
#---------------------------------------------------------------------------------------
#=======================================================================================
print("Data read.")

def julian(dd,mm,yy):
    """
    Function that calculates julian days [Day of Year] from day,month, and year imput
    Imput:
       yy:     year [integer]
       mm:     month [integer]
       dd:     day [integer]
    Output:
       jday:   COrresponding Julian day or Day of Year [1,366]
    Example:
    --------
      >>> jday = julian(1,1,1980)
    """
    if yy % 4 != 0 or yy % 100 == 0:
       mon=[31,28,31,30,31,30,31,31,30,31,30,31]
    if yy % 4 == 0 and yy % 100 != 0 or yy % 400 == 0:
       mon=[31,29,31,30,31,30,31,31,30,31,30,31]
    if mm == 1:
       jday=dd
    if mm > 1:
       jday=sum(mon[0:mm-1])+dd
    return jday

#=======================================================================================
def dates_daily(y0,m0,d0,mtot,noleap):
    """
    Function that calculates arrays of dates (hour, day,month, year)

    Imput:
       y0:     initial year [integer]
       m0:     initial month [integer]
       d0:     initial day [integer]
       mtot:   total number of elements in the time series [integer]
       noleap: flag to indicate whether or not leap years should be considered.
               noleap = 0 --> time searies contain Feb 29
               noleap = 1 --> time series does not contain Feb 29
    Output:
       day:    array with values of days
       month:  array with values of months
       year:   array with values of years
    Example:
    --------
      >>> day,month,year=dates_daily(1980,1,1,365,0)
    """
    day=np.zeros((mtot))
    month=np.zeros((mtot))
    year=np.zeros((mtot))
    start_date=date(y0,m0,d0)
    deltad = timedelta(days=1)
    single_date=start_date
    dt=0
    while dt < mtot:
          day[dt]=single_date.day
          month[dt]=single_date.month
          year[dt]=single_date.year
          if noleap == 0:
             single_date=single_date+deltad
          if noleap == 1:
             if year[dt] % 4 != 0:
                single_date=single_date+deltad
             if year[dt] % 4 == 0:
                if month[dt] != 2:
                   single_date=single_date+deltad
                if month[dt] == 2:
                   if day[dt] < 28:
                      single_date=single_date+deltad
                   if day[dt] == 28:
                      single_date=single_date+2*deltad
          dt=dt+1
    return day, month, year

#=======================================================================================
def mk_test_sens_slope(tmp, alpha=0.05):
    """
    Copyright (c) 2017 Michael Schramm

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    This function calculates the Mann-Kendall Test and the Sen's slope calculation. It was
    partially inspired by the function created by Michael Schramm
    https://github.com/mps9506/Mann-Kendall-Trend/blob/master/mk_test.py.

    Input:
        tmp:   a vector of data
        alpha: significance level (0.05 default)
    Output:
        slope: the Sen's slope value
        h:     True (if trend is present) or False (if trend is absence)
        pval:  p value of the significance test
        zmk:   normalized test statistics
    Examples
    --------
      >>> tmp = np.random.rand(100)
      >>> slope,h,pval,zmk = mk_test_sens_slope(tmp,0.05)
    """
    num = len(tmp)
    index=np.arange(0,num,1)
    #    Creating arrays with dimension size equal to the number of all possible differences
    sgn = np.zeros((int(num*(num-1)/2)))
    slo = np.zeros((int(num*(num-1)/2)))
    #    Calculaing the indicator function (sgn) and the Sen's slope (slo) using only
    # one loop.
    beg=0                       # first position of sgn vector for the kth iteration
    for kt in range(1,num):     # loop in time
        ned=beg+num-kt          # last position of sgn vector for the kth iteration
        sgn[beg:ned]=tmp[kt:num]-tmp[kt-1]
        slo[beg:ned]=(tmp[kt:num]-tmp[kt-1])/(index[kt:num]-index[kt-1])
        beg=ned
    #    Calculaing the mean
    es = np.sum(np.sign(sgn))
    # ---- Calculating the variance -----
    # calculating the correction term for tied observations (qt = number of data points in each tie group)
    corr=0.
    un = np.unique(tmp)
    for kt in un:
        id=np.where(tmp == kt)
        corr += len(id[0])*(len(id[0])-1.)*(2*len(id[0])+5)
    var = (num*(num - 1.)*(2.*num+5.) - corr)/18.
    if es > 0.:
       zmk = (es-1.)/np.sqrt(var)
    if es == 0.:
       zmk = 0.
    if es < 0.:
       zmk = (es+1.)/np.sqrt(var)
    # calculate the p_value
    pval = 2*(1-norm.cdf(abs(zmk)))  # two tail test
    h = abs(zmk) > norm.ppf(1-alpha/2)
    slope = np.median(slo)
    return slope, h, pval, zmk

#=======================================================================================
"""
Funtion that calculates the Fourier coefficients and the explained variance of the Nth
first harmonics of a time series

Input:
   tseries: input time series
   nmodes : number of harmonics to retain (N)
   coefa  : Array with N (or 'nmodes') elements
   coefb  : Array with N (or 'nmodes') elements
   hvar   : Array with N (or 'nmodes') elements
   missval: Falg value for missing data
Output:
   coefa: Array of A coefficients of the Nth first harmonics
   coefb: Array of B coefficients of the Nth first harmonics
   hvar : Array of explained variance of the Nth first harmonics
"""
def Harmonics(coeafa,coefb,hvar,tseries,nmodes,missval):
    mtot=len(tseries)
    time=np.arange(1,mtot+1,1.)
    newdim=len(tseries[tseries!=missval])  # removing missing data
    time=time[tseries!=missval]
    tdata=tseries[tseries!=missval]
    svar=sum((tdata[:]-np.mean(tdata))**2)/(newdim-1)
    nm=nmodes
    if 2*nm > newdim:
       nm=newdim/2
    coefa=np.zeros((nm))
    coefb=np.zeros((nm))
    hvar=np.zeros((nm))
    for tt in range(0,nm):
        Ak=np.sum(tdata[:]*np.cos(2.*math.pi*(tt+1)*time[:]/float(newdim)))
        Bk=np.sum(tdata[:]*np.sin(2.*math.pi*(tt+1)*time[:]/float(newdim)))
        coefa[tt]=Ak*2./float(newdim)
        coefb[tt]=Bk*2./float(newdim)
        hvar[tt]=newdim*(coefa[tt]**2+coefb[tt]**2)/(2.*(newdim-1)*svar)
    return coefa,coefb,hvar

#================================= Formatting Data ======================================
print("Formatting Data...")


"""
Removing Feb 29th. This block averages Feb 28 and 29 in leap years. 
"""
day,month,year=dates_daily(1979,1,1,ntot,0)
id=np.where((month == 2.) & (day == 29.))
id2=[x-1 for x in id]         # creating a list equal to (id-1) 
prec2=0.5*(prec[id2,:,:]+prec[id,:,:]) 
prec[id2,:,:]=prec2[0,:,:,:]  # python adds a dimention to the result
prec=np.delete(prec,id,axis=0)
year=np.delete(year,id,axis=0)
month=np.delete(month,id,axis=0)
day=np.delete(day,id,axis=0)
ntot=len(year)
prec[prec<0.]=missval
#---------------------------------------------------------------------------------------
"""
Masking Missing values. Sometimes datasets have significant amounts of
missing data (e.g. land only data, ocean only data, complex topography).
In such cases it is sometimes useful to mask those regions. Masking can
prevent code errors and improve code efficiency. The minimum percentage
of missing data allowed is determined by namelist variable "dper". If a
grid pint has more missing data then the minimum percentage, that grid
point will be masked at all times.
"""
thres=(ntot-0.01*dper*ntot) # minimum threshold of non-missing data
mask=np.zeros((nlat,nlon))
for it in range(0,nlat):
    for jt in range(0,nlon):
        id=np.where(prec[:,it,jt]!=missval)
        if len(id[0]) >= thres:
           mask[it,jt]=1.
        if len(id[0]) < thres:
           mask[it,jt]=0.

print("Data Formatted.")
#================================= Removing Trends =====================================
print("Removing Trends... This might take a while.")
"""
Removing trends using the Mann-Kendall test and and Sen's slope methods.
Only statistically significant (at 5% level) trends are removed.

Comment this block if you are interested in analyzing trends
"""
time=np.arange(0,ntot,1)
for it in range(0,nlat):
    for jt in range(0,nlon):
        if mask[it,jt] == 1.:
           tmp=prec[:,it,jt]
           id=np.where(tmp != missval)
           slope,h,pval,zmk=mk_test_sens_slope(tmp[id[0]],alpha=0.05)
           print(slope)
           if h:
              prec[id[0],it,jt]=prec[id[0],it,jt]-time[id[0]]*slope

print("Trends Removed.")
#========================== Calculating the Mean Annual Cycle ==========================
print("Calculating the mean annual cycle... This might take a while.")
"""
This block will calculate the mean annual cycle for the whole time series. It will need
to be adapted if the period of interest is only a portion of the total time series.
For example, climatologies of a 30 year period of reference: 1981-2010.
"""

""" calculating Julian days """
jday=np.zeros((len(prec[:,0,0])))
for tt in range(0,ntot):
    jday[tt]=julian(int(day[tt]),int(month[tt]),int(year[tt]))

""" calculating the mean annual cycle for each grid point """
cycle=np.zeros((tot,nlat,nlon))
for tt in range(0,tot):
    id=np.where(jday[:] == jday[tt])
    for it in range(0,nlat):
        for jt in range(0,nlon):
            if mask[it,jt] == 1.:
               tmp=prec[id[0],it,jt]
               id2=np.where(tmp >= 0.)
               if len(id2[0]) > 1:      # have to specify id2[0] because id2 is a tuple
                  cycle[tt,it,jt]=np.mean(tmp[id2[0]])


print("Mean annual cycle calculated.")
#=========================== Smoothing the Mean Annual Cycle ===========================
print("Smoothing the mean annual cycle... This might take a while.")
"""
This block will calculate the smoothed mean annual cycle. This is anessential part of
calculating anomalies (often overlooked).

For details, see: Bombardi RJ and Carvalho LMV (2017). Simple Practices in Climatological Analyses: A Review. Revista Brasileira de Meteorologia. 32 (3), 311-320
"""
time=None
time=np.arange(1,tot+1,1.)
smoothed=np.zeros((tot,nlat,nlon))
for it in range(0,nlat):
    for jt in range(0,nlon):
        if mask[it,jt] == 1.:
           coefa=np.zeros((3))
           coefb=np.zeros((3))
           hvar=np.zeros((3))
           tseries=cycle[:,it,jt]
           coefa,coefb,hvar=Harmonics(coefa,coefb,hvar,tseries,3,missval)
           smoothed[:,it,jt]=np.mean(cycle[:,it,jt])
           for pp in range(0,3):
               smoothed[:,it,jt]=smoothed[:,it,jt]+coefa[pp]*np.cos(2.*math.pi*time[:]*(pp+1)/float(tot))+coefb[pp]*np.sin(2.*math.pi*time[:]*(pp+1)/float(tot))

print("Mean annual cycle smoothed.")
#=================================== Calculating Anomalies =============================
print("Calculating anomalies...")
anomalies=np.zeros((ntot,nlat,nlon))
for tt in range(0,ntot,tot):
    beg=tt
    if ned <= ntot:
       ned=beg+tot
    if ned > ntot:
       ned=ntot
    anomalies[beg:ned,:,:]=prec[beg:ned,:,:,]-smoothed[0:ned-beg,:,:]
    ned-beg,tot

"""
Making sure missing values in the original input data are transmitted to the anomalies
"""
for it in range(0,nlat):
    for jt in range(0,nlon):
        if mask[it,jt] == 1.:
           id=np.where(prec[:,it,jt]==missval)
           anomalies[id,it,jt]=missval
        if mask[it,jt] == 0.:
           anomalies[:,it,jt]=missval

print("Anomalies Calculated.")
#======================================= Saving results ================================
print("Saving Results...")
"""
This block will save NetCDF files of anomalies separated into years. That is, each year
of data will be saved in a separate NetCDF file to avoid memory issues printing or
reading the data.

These NetCDF files can be directly opened in GrADS. Change the value of 'var.long_name'
to acurately describe the data
"""

#----- saving anomalies as individual years -----
for tt in range(0,ntot,tot):
    beg=tt
    if ned <= ntot:
       ned=beg+tot
    if ned > ntot:
       ned=ntot
    outfile = pathout+prefix+str(year[tt])[0:4]+".nc"
    from netCDF4 import Dataset
    rootgrp = Dataset(outfile, "w", format="NETCDF4")
    rootgrp.close()
    rootgrp = Dataset(outfile, "a")
    # Creating dimensions
    lon = rootgrp.createDimension("lon", len(lons[0:nlon]))
    lat = rootgrp.createDimension("lat", len(lats[0:nlat]))
    lev = rootgrp.createDimension("lev", 1)
    time = rootgrp.createDimension("time",tot)
    #Creating coordinates
    times = rootgrp.createVariable("time","i4",("time",))
    levels = rootgrp.createVariable("lev","i4",("lev",))
    latitudes = rootgrp.createVariable("lat","f8",("lat",))
    longitudes = rootgrp.createVariable("lon","f8",("lon",))
    # Filling coordinates
    longitudes.units='degrees_east'
    longitudes.long_name='Longitude'
    longitudes[:]=lons[0:nlon]
    latitudes.units='degrees_north'
    latitudes.long_name='Latitude'
    latitudes[:]=lats[0:nlat]
    levels.units='millibar'
    levels.long_name='Level'
    levels[:]=1
    times.long_name='Time'
    times.units='days since '+str(year[tt])[0:4]+'-01-01 00:00'
    times[:]=jday[0:ned-beg]-1.
    # Creating variables 
    var = rootgrp.createVariable("anom","f4", ("time","lev","lat","lon",),fill_value=missval)
    # Filling Variable
    var.long_name = 'variable anomalies [units]'
    var[:,0,:,:] = anomalies[tt:tt+ned-beg,:,:]

#=======================================================================================
#---------------------------------------------------------------------------------------
"""                                    END OF PROGRAM                                """
#---------------------------------------------------------------------------------------
#=======================================================================================

