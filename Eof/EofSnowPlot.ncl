

begin
  vname  = "snow"               ; variable name
  fname  = "EofSnowData.nc"   ; CDC: NCEP Reanalysis




  optEOF = False      
  optETS = False

; ==============================================================
; Open the file: Read only the user specified period
; ==============================================================
  f      = addfile (fname, "r")



  VAR    = f->$vname$
                        ; variable overview                  

; ==============================================================
; compute desired global seasonal mean: month_to_season (contributed.ncl) 
; ==============================================================
  var    = VAR
  nyrs   = dimsizes(var&time)

; =================================================================
; create weights:  sqrt(cos(lat))   [or sqrt(gw) ] for covariance
; =================================================================
  rad    = 4.0*atan(1.0)/180.0
  clat   = cos(rad*var&lat)           
  clat   = where(clat.lt.0, 0.0, clat)  ; avoid a potential numerical issue at pole
  clat   = sqrt( clat )  ; avoid a potential numerical issue at pole

; =================================================================
; weight all observations 
; =================================================================
  wvar   = var                                   ; copy meta data
  wvar   = var*conform(var, clat, 1)

; =================================================================
; EOF
; Reorder (lat,lon,time) the *weighted* input data
; Access the area of interest via coordinate subscripting
; =================================================================
  x      = wvar(lat|:,lon|:,time|:)

  eof    = eofunc_Wrap(x, 3, optEOF)  
if (abs(max(eof)).ge.abs(min(eof))) then
   temp=abs(max(eof)) ; write branch
else
    temp=abs(min(eof)); write else branch
end if
do neof=1,3
wks=gsn_open_wks("png", "SnowEof"+neof)
res=True
res@gsnRightString="EOF "+neof
res@pmTickMarkDisplayMode="Always"
res@mpMinLatF=min(eof&lat)
res@mpMaxLatF=max(eof&lat)
res@mpMinLonF=min(eof&lon)
res@mpMaxLonF=max(eof&lon)
res@mpFillOn=False
res@cnFillOn=True
res@cnLinesOn=False
res@cnLevelSelectionMode="ExplicitLevels"
res@cnLevels=fspan(-temp,temp, 17)
res@gsnAddCyclic=False
plot=gsn_csm_contour_map(wks, eof(neof-1,:,:), res)
end do



end