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

ufile = addfile("./input/climo/u200.seas.nc","r")
vfile = addfile("./input/climo/v200.seas.nc","r")
zfile = addfile("./input/climo/z200.seas.nc","r")

cuvar = ufile->var131(2,0,:,:) ; jja
cvvar = vfile->var132(2,0,:,:) ; jja
czvar = zfile->var129(2,0,:,:) ; jja

printVarSummary(czvar)
printVarSummary(cuvar)
printVarSummary(cvvar)

;===============================================================================

  do year=1979,2021
     print(year)

zfile = addfile("/glade/p/univ/uchi0007/ERA5_Data/Z/daily/z200/anom/z3."+year+".nc","r") ; already in geopotential daily anom
zvar = zfile->var129(:,{20000:20000},:,:)
printVarSummary(zvar)

; pressure correction for hPa input levels
  if (max(zvar&plev).gt.1000) then
    p_corr = 0.01d
  else
    p_corr = 1.d
  end if


time  = zfile->time
lat   = zfile->lat
lon   = zfile->lon
level = zfile->plev({20000:20000})*p_corr

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

; 2-D -> 4-D

cuvar_tmp = conform_dims(dimsizes(zvar),cuvar,(/2,3/))
cvvar_tmp = conform_dims(dimsizes(zvar),cvvar,(/2,3/))
czvar_tmp = conform_dims(dimsizes(zvar),czvar,(/2,3/))



; magnitude of climatological wind
cumag = sqrt(cuvar_tmp^2 + cvvar_tmp^2)
cumag@_FillValue = -1.e+21
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

; QG steam function anomaly
;psidev = (zvar-czvar_tmp) /ftmp
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


xuterm = (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
xvterm = (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)
yvterm = (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)


;  Mask out where westerlies is small or negative (less than 5 m/s).
;  by using mask

;x-component of (38)
;Fx = mask( leveltmp/1000./(2.*cumag*re*re)*( cuvar_tmp/coslattmp * xuterm + cvvar_tmp * xvterm), cuvar_tmp.lt.5,False)

;y-component 
;Fy = mask(leveltmp/1000./(2.*cumag*re*re)*( cuvar_tmp*xvterm + coslattmp*cvvar_tmp*yvterm),cuvar_tmp.lt.5,False)


;No Mask
;x-component of (38)
Fx = leveltmp/1000./(2.*cumag*re*re)*( cuvar_tmp/coslattmp * xuterm + cvvar_tmp * xvterm)

;y-component 
Fy = leveltmp/1000./(2.*cumag*re*re)*( cuvar_tmp*xvterm + coslattmp*cvvar_tmp*yvterm)


Fx@units = "m^2/s^2"
Fx@units = "m^2/s^2"
psidev@units = "m^2/s"


 copy_VarCoords(zvar,Fx)
 copy_VarCoords(zvar,Fy)
 copy_VarCoords(zvar,psidev)

printVarSummary(Fx)
printVarSummary(Fy)
printVarSummary(psidev)

; Save the output in netcdf
;==================================
setfileoption("nc","Format","LargeFile")
system( "rm " + "waf.jja."+year+".nc" )
output = addfile("waf.jja."+year+".nc","c")
output->Fx = Fx
output->Fy = Fy
output->psidev = psidev

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
delete(cuvar_tmp)
delete(cvvar_tmp)
delete(czvar_tmp)
delete(cumag)
delete(psidev)
delete(dpsidevdlon)
delete(ddpsidevdlonlon)
delete(dpsidevdlat)
delete(ddpsidevdlonlat)
delete(ddpsidevdlatlat)
delete(xuterm)
delete(xvterm)
delete(yvterm)
delete(Fx)
delete(Fy)
delete(output)

  end do
end
