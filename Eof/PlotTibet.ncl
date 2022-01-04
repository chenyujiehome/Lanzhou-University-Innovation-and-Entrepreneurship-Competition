 begin
  f      = addfile ("EofSnowData.nc", "r")



  snowF    = f->snow
  snow=dim_avg_n(snowF, 0)
  snow!0="lat"
  snow&lat=snowF&lat
  snow!1="lon"
  snow&lon=snowF&lon
  
  if (abs(max(snow)).ge.abs(min(snow))) then
   temp=abs(max(snow)) ; write branch
else
    temp=abs(min(snow)); write else branch
end if
wks=gsn_open_wks("png", "Snow")
res=True
res@gsnRightString="snow "
res@pmTickMarkDisplayMode="Always"
res@mpMinLatF=min(snow&lat)
res@mpMaxLatF=max(snow&lat)
res@mpMinLonF=min(snow&lon)
res@mpMaxLonF=max(snow&lon)

res@mpFillOn=False
res@mpGridSpacingF=10.
res@cnFillOn=True
res@cnLinesOn=False
res@cnLevelSelectionMode="ExplicitLevels"
res@cnLevels=fspan(-temp,temp, 17)
res@gsnAddCyclic=False
plot=gsn_csm_contour_map(wks, snow, res)

end