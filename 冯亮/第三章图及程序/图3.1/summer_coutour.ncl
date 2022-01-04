f=addfile("1961-2018-month-pr.nc","r")
fn = addfile("mask_data.nc","r")
shptp="/home/heyl/tianwenlong/fengliang/nc_coutour/add_TP/add_TP/DBATP_Line.shp"
data = fn ->data
a = f ->pr
lat=f->lat
lon=f->lon
;==================================
m = new((/696,31,61/),float)               ;mask青藏高原
copyatt(m,a)
m = m@_FillValue
do i= 0,695
 do j = 0,30
  do k= 0,60
m(i,j,k) = where (ismissing(data(j,k)),a@_FillValue,a(i,j,k))
end do
end do 
end do
;===============================
p = reshape(m,(/58,12,31,61/))            ;分为年，月，经，维四维数据
;printVarSummary(p)
n=p(:,5:7,:,:)                       ;提取某月数据
pp=dim_avg_n(n,1)
;printVarSummary(pp)
ppp  = dim_avg_n(pp,0)
copyatt(ppp,m(0,:,:))
;printVarSummary(ppp)
printMinMax(ppp,1)
;copyatt(pp,m)
;year=ispan(1961,2018,1)
;pp!0="year"
;pp&year=year
;printVarSummary(pp)
;rc=regCoef_n(year,pp,0,0)
;copyatt(rc,pp(0,:,:))
;printMinMax(rc,1)
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
res@cnMinLevelValF               = 0
res@cnMaxLevelValF               = 200
res@cnLevelSpacingF              = 10


res@mpGridAndLimbOn       =  True
res@mpGridLatSpacingF     =  1
res@mpGridLonSpacingF     =  1  
res@mpGridLineColor       = "gray" 
res@mpGridLineThicknessF  =  1 

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
res@mpLeftCornerLonF      = 70.
res@mpRightCornerLatF     = 40.
res@mpRightCornerLonF     = 105.

res@mpDataSetName         = "Earth..4"
res@mpDataBaseVersion     = "MediumRes"
res@mpOutlineBoundarySets = "AllBoundaries"

;res@mpDataBaseVersion     = "MediumRes"
;res@mpOutlineBoundarySets = "National"
;res@mpDataBaseVersion     = "MediumRes"
;res@mpOutlineBoundarySets = "National"

;res@tiMainString   = "Year Average Precipitation Distribution"
res@gsnLeftString  = "Precipitation"
res@gsnRightString = "mm/a"
wks = gsn_open_wks("png", "summer")
plot = gsn_csm_contour_map(wks,ppp,res)
tpres        = True
 tpres@gsLineColor = "black"     ;设置廓线颜色
 tpres@gsLineThicknessF = 3              ;设置廓线宽度
 shp_plot2     = gsn_add_shapefile_polylines(wks,plot,shptp,tpres)        ;base为底图
