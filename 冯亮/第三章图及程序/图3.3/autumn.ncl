f=addfile("1961-2018-month-pr.nc","r")
shptp="/home/heyl/tianwenlong/fengliang/nc_coutour/add_TP/add_TP/DBATP_Line.shp"
a = f ->pr
lat=f->lat
lon=f->lon
printVarSummary(a)
p = reshape(a,(/58,12,31,61/)) 
printVarSummary(p)
m=p(:,8:10,:,:)
pp=dim_avg_n(m,1)
copyatt(pp,a)
year=ispan(1961,2018,1)
pp!0="year"
pp&year=year
;printVarSummary(pp)
rc=regCoef_n(year,pp,0,0)
copyatt(rc,pp(0,:,:))
printMinMax(rc,1)
tval=onedtond(rc@tval,dimsizes(rc))
 df=onedtond(rc@nptxy,dimsizes(rc))-2
 b=tval
 b=0.5
 prob=betainc(df/(df+tval^2),df/2.0,b)
 copyatt(prob,rc)
wks = gsn_open_wks("png", "autumn_trend")
res2=True
  res2@gsnDraw=False
   res2@gsnFrame=False
    res2@cnFillOn=True
     res2@cnLinesOn=False
      res2@cnLineLabelsOn=False
       res2@cnInfoLabelOn=False
     res2@lbLabelBarOn=False
         res2@cnMonoFillPattern=False
           res2@cnLevelSelectionMode="ExplicitLevels"
             res2@cnLevels=(/0.05/)
             res2@cnFillPatterns=(/3,-1/)
                 res2@cnMonoFillScale=False
                   res2@cnFillScales=(/0.5,2.5/)
                     res2@cnFillColors=(/1,0/)
                       res2@gsnRightString=""
                         res2@gsnLeftString=""
                           plot1=gsn_csm_contour(wks,prob,res2)
res = True   
res@mpMinLatF           = 25
res@mpMaxLatF           = 40
res@mpMinLonF           = 75
res@mpMaxLonF           = 105
res@mpFillOn                     = True
res@mpOutlineOn                  = True

res@cnLinesOn                    =True
res@cnLineDashPattern            =0
res@cnLineColor                  = "black"
res@cnLineThicknessF             = 2
res@cnLevelSelectionMode         = "ManualLevels"
res@cnMinLevelValF               = -1
res@cnMaxLevelValF               = 1
res@cnLevelSpacingF              = 0.1


;res@mpGridAndLimbOn       =  True
;res@mpGridLatSpacingF     =  1
;res@mpGridLonSpacingF     =  1  
;res@mpGridLineColor       = "gray" 
;res@mpGridLineThicknessF  =  1 


res@mpDataSetName         = "Earth..4"
res@mpDataBaseVersion     = "MediumRes"
res@mpOutlineBoundarySets = "AllBoundaries"

res@tmYLLabelStride       =  1  
res@gsnMajorLatSpacing    =  5  
res@gsnMinorLatSpacing    =  1
res@tmYLMajorLineColor    = "black"
res@tmYLMinorLineColor    = "grey20"
res@tmYLLabelFontColor    = "black" 

res@tmXBLabelStride       =  1
res@gsnMajorLonSpacing    =  5
res@gsnMinorLonSpacing    =  1
res@tmXBMajorLineColor    = "black"
res@tmXBMinorLineColor    = "grey20"
res@tmXBLabelFontColor    = "black"

res@mpGridLineDashPattern =  3
res@pmTickMarkDisplayMode = "Always"

res@gsnMaximize           = True
res@cnFillOn              = True
res@cnLinesOn             =  False
res@cnLineLabelsOn        = False
res@cnLevelSelectionMode  = "ManualLevels"

res@tfDoNDCOverlay        = True
res@mpLimitMode           = "Corners"
res@mpLeftCornerLatF      = 25.
res@mpLeftCornerLonF      = 75.
res@mpRightCornerLatF     = 40.
res@mpRightCornerLonF     = 105.

;res@mpDataBaseVersion     = "MediumRes"
;res@mpOutlineBoundarySets = "National"
;res@mpDataBaseVersion     = "MediumRes"
;res@mpOutlineBoundarySets = "National"
plot = gsn_csm_contour_map(wks,rc,res)
res@cnFillPalette = "cmp_b2r"
tpres        = True
tpres@gsLineColor = "black"     ;设置廓线颜色
tpres@gsLineThicknessF = 8              ;设置廓线宽度
   shp_plot2     = gsn_add_shapefile_polylines(wks,plot,shptp,tpres)    
overlay(plot,plot1)                                                                                               
 draw(plot)
  frame(wks)
