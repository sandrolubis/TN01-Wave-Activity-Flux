; Sandro Lubis, July 2014 (updated Nov, 2016).
; PhD Student, GEOMAR.
; 3-D wave-activity flux derived by Takaya and Nakamura (1999, 2001)
; See (38) of Takaya and Nakamura (2001, JAS)
; Dr. Sandro Lubis M.Sc - GEOMAR

;==============================================================================

; Input data: 
;   Monthly-mean data 
;   geopotential height (hgt: m) 
;
;  Monthly climatology of U, V, and T
;   air tempeature (air: degC)   
;   zonal wind (uwnd: m/s)     
;   meridional wind (vwnd: m/s) 
;
; Notes:
; The unit of level is hPa or Pa
; Basic state: monthly climatology
; Perturbation: deviation from climatology
; Calculating the flux for each month of input data.
; Note that 
; (i) "p" (pressure) in Eq. (5.7) is dimensionless (pressure/1000hPa)

; Output files are
; x-component:  TN2001-Fx
; y-component:  TN2001-Fy
; z-component:  TN2001-Fz
; QG stream function anomaly

;==============================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin


;  Some Constants
;===============================================================================
gc=290.             				   ;  Gas constant
ga=9.80665         				   ;  Gravitational acceleration
re=6378388.        				   ;  Radius of the earth
sclhgt=8000.                                       ;  Scale height
pi = atan(1.0)*4.                                  ;  Pi
;===============================================================================

; monthly climatology data

tfile = addfile("./input/climo/t200.seas.nc","r")
ufile = addfile("./input/climo/u200.seas.nc","r")
vfile = addfile("./input/climo/v200.seas.nc","r")
zfile = addfile("./input/climo/z200.seas.nc","r")

ctvar = tfile->var130(2,0,:,:) ; jja
cuvar = ufile->var131(2,0,:,:) ; jja
cvvar = vfile->var132(2,0,:,:) ; jja
czvar = zfile->var129(2,0,:,:) ; jja

printVarSummary(ctvar)
printVarSummary(cuvar)
printVarSummary(cvvar)

;===============================================================================

  do year=1979,1979
     print(year)

zfile = addfile("./input/anom/erai_Z_"+year+".anom.nc","r") ; already in geopotential daily anom
zvar = zfile->var129
printVarSummary(zvar)

; pressure correction for hPa input levels
  if (max(zvar&lev).gt.1000) then
    p_corr = 0.01d
  else
    p_corr = 1.d
  end if


time  = zfile->time
lat   = zfile->lat
lon   = zfile->lon
level = zfile->lev*p_corr

ntime = dimsizes(time)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlevel = dimsizes(level)


; Coriolis parameter
f      =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
f!0    = "lat"
f&lat  = lat
f@_FillValue = -1.e+21

; missing for 10S - 10N
do ilat = 0, nlat-1
if (abs(lat(ilat) ).lt. 10. ) then
f(ilat)= f@_FillValue
end if
end do

; cosine
coslat = cos(lat(:)*pi/180.)

print(level)

; 1-D -> 4-D
leveltmp = conform_dims(dimsizes(zvar),level,1)
coslattmp = conform_dims(dimsizes(zvar),coslat,2)
ftmp = conform_dims(dimsizes(zvar),f,2)

; vertical gradient of potential temperature  (K/m)
dthetadz = center_finite_diff_n(ctvar*(1000./leveltmp)^0.286,-sclhgt*log(level/1000),False,0,1)

; Brunt Vaisala frequency
NN = (gc*(leveltmp/1000.)^0.286)/sclhgt * dthetadz
NN@_FillValue = -1.e+21
NN = where(NN .gt. 0, NN, NN@_FillValue)

; magnitude of climatological wind
cumag = sqrt(cuvar^2 + cvvar^2)
cumag@_FillValue = -1.e+21
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

; QG steam function anomaly
psidev = zvar /ftmp

;dpsidev/dlon
dpsidevdlon =  center_finite_diff_n(psidev,lon*pi/180.,True,0,3)

;ddpsidev/dlonlon
ddpsidevdlonlon =  center_finite_diff_n(dpsidevdlon,lon*pi/180.,True,0,3)

;dpsidev/dlat
dpsidevdlat = center_finite_diff_n(psidev, lat*pi/180., False,0,2)

;ddpsidev/dlonlat
ddpsidevdlonlat =  center_finite_diff_n(dpsidevdlon,lat*pi/180.,False,0,2)

;ddpsidev/dlatdlat
ddpsidevdlatlat = center_finite_diff_n(dpsidevdlat, lat*pi/180.,False,0,2)

;dpsidev/dz
dpsidevdz = center_finite_diff_n(psidev, -sclhgt*log(level/1000),False,0,1)

;ddpsidev/dlondz
ddpsidevdlonz = center_finite_diff_n(dpsidevdlon, -sclhgt*log(level/1000),False,0,1)

;ddpsidev/dlatdz
ddpsidevdlatz = center_finite_diff_n(dpsidevdlat, -sclhgt*log(level/1000),False,0,1)


xuterm = (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
xvterm = (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)
yvterm = (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)
zuterm = (dpsidevdlon*dpsidevdz - psidev*ddpsidevdlonz)
zvterm = (dpsidevdlat*dpsidevdz - psidev*ddpsidevdlatz)

;  Mask out where westerlies is small or negative (less than 5 m/s).
;  by using mask

;x-component of (38)
Fx = mask( leveltmp/1000./(2.*cumag*re*re)*( cuvar/coslattmp * xuterm + cvvar * xvterm), cuvar.lt.5,False)

;y-component 
Fy = mask(leveltmp/1000./(2.*cumag*re*re)*( cuvar*xvterm + coslattmp*cvvar*yvterm),cuvar.lt.5,False)

;z-component 
Fz = mask(leveltmp/1000.*ftmp*ftmp/(2.*cumag*NN*re)*(cuvar*zuterm + coslattmp*cvvar*zvterm), cuvar.lt.5, False)


Fx@units = "m^2/s^2"
Fx@units = "m^2/s^2"
Fz@units = "m^2/s^2"
psidev@units = "m^2/s"


 copy_VarCoords(zvar,Fx)
 copy_VarCoords(zvar,Fy)
 copy_VarCoords(zvar,Fz)
 copy_VarCoords(zvar,psidev)

printVarSummary(Fx)
printVarSummary(Fy)
printVarSummary(Fz)
printVarSummary(psidev)

; Save the output in netcdf
;==================================
setfileoption("nc","Format","LargeFile")
system( "rm " + "TN01.mon."+year+".nc" )
output = addfile("TN01.mon."+year+".nc","c")
output->Fx = Fx
output->Fy = Fy
output->Fz = Fz
output->psidev = psidev
output->zvar  = zvar

; Delete temporary files.
;==================================
 
delete(zfile)
delete(zvar)
delete(time)
delete(lat)
delete(lon)
delete(level)
delete(ntime)
delete(nlat)
delete(nlon)
delete(nlevel)
delete(f)
delete(coslat)
delete(leveltmp)
delete(coslattmp)
delete(ftmp)
delete(dthetadz)
delete(NN)
delete(cumag)
delete(psidev)
delete(dpsidevdlon)
delete(ddpsidevdlonlon)
delete(dpsidevdlat)
delete(ddpsidevdlonlat)
delete(ddpsidevdlatlat)
delete(dpsidevdz)
delete(ddpsidevdlonz)
delete(ddpsidevdlatz)
delete(xuterm)
delete(xvterm)
delete(yvterm)
delete(zuterm)
delete(zvterm)
delete(Fx)
delete(Fy)
delete(Fz)
delete(output)

  end do
end
