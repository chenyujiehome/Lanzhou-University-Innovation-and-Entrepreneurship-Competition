begin
do nsvd = 1, 3
    ; write loop content


 f    = addfile("RightSvd.nc", "r")
 snow=f->snow

wks=gsn_open_wks("png", "RhotoSvd"+nsvd)
res=True
  res@cnFillPalette        = "BlWhRe"     ; choose colormap
res@gsnRightString="SVD "+nsvd
res@pmTickMarkDisplayMode="Always"
res@mpMinLatF=min(snow&lat)
res@mpMaxLatF=max(snow&lat)
res@mpMinLonF=min(snow&lon)
res@mpMaxLonF=max(snow&lon)

res@mpFillOn=False
res@mpGridAndLimbOn=True

res@mpGridLineThicknessF=0.5
res@mpGridLineDashPattern=2
res@mpGridSpacingF=10.
res@cnFillOn=True
res@cnLinesOn=False
res@cnLevelSelectionMode="ExplicitLevels"
res@cnLevels=fspan(min(snow(nsvd-1,:,:)), max(snow(nsvd-1,:,:)), 10)
res@gsnAddCyclic=False
plot=gsn_csm_contour_map(wks, snow(nsvd-1,:,:), res)

end do
 end
