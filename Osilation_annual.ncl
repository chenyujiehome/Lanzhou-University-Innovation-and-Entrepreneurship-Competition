load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"



begin

  dir        = "/mnt/zfm/ERA5-monthly/"
  f_obs      = addfile (dir + "ERA5-mon_geopt.grib", "r")

;-------------------------------------------------------------
; Read the data at 200 hPa
;-------------------------------------------------------------
  time_julian      = f_obs->initial_time0_hours
  time_utc         = calendar_decode2(time_julian,2)
  print(time_utc)

  geopt_Jun        = f_obs->Z_GDS0_ISBL_S123(5:485:12,{200},:,:)  ; June   of 1979~2019
  geopt_Jul        = f_obs->Z_GDS0_ISBL_S123(6:486:12,{200},:,:)  ; July   of 1979~2019
  geopt_Aug        = f_obs->Z_GDS0_ISBL_S123(7:487:12,{200},:,:)  ; August of 1979~2019

  geoheight_Jun    = geopt_Jun/9.8
  geoheight_Jul    = geopt_Jul/9.8
  geoheight_Aug    = geopt_Aug/9.8
  geoheight_Jun@unit   = "m"
  geoheight_Jul@unit   = "m"
  geoheight_Aug@unit   = "m"
  copy_VarCoords(geopt_Jun,geoheight_Jun)
  copy_VarCoords(geopt_Jun,geoheight_Jul)
  copy_VarCoords(geopt_Jun,geoheight_Aug)
  printVarSummary(geoheight_Jun)


;-------------------------------------------------------------
; Compute the Summer mean results in each year
;-------------------------------------------------------------
  geoheight_JJA     = (geoheight_Jun + geoheight_Jul + geoheight_Aug)/3.0
  copy_VarCoords(geoheight_Jun,geoheight_JJA)
  printVarSummary(geoheight_JJA)


;-------------------------------------------------------------
; Select the regional data and compute its average
;-------------------------------------------------------------
  geoheight1_JJA       = dim_avg_n(geoheight_JJA(:,{20:27.5},{85:115}),(/1,2/))    ; Southeast
  geoheight2_JJA       = dim_avg_n(geoheight_JJA(:,{27.5:35},{50:80}),(/1,2/))     ; Northwest
  printVarSummary(geoheight1_JJA)
  printVarSummary(geoheight2_JJA)


;-------------------------------------------------------------
; Compute the difference (Ise-Inw) and normalize it
;-------------------------------------------------------------
  geoheight_diff_JJA   = geoheight1_JJA - geoheight2_JJA
  Index                = dim_standardize(geoheight_diff_JJA,1) 

  printVarSummary(geoheight_diff_JJA)
  printVarSummary(Index)


  time_label = ispan(0,dimsizes(Index)-1,1)
  printVarSummary(time_label)


;-------------------------------------------------------------
; Start to plot the index
;-------------------------------------------------------------
  wks = gsn_open_wks("ps","index_se-nw")
  gsn_define_colormap(wks,"wh-bl-gr-ye-re")

  res = True
  res@vpXF          = 0.15
;  res@vpYF          = 0.8
  res@vpHeightF     = 0.3               ; change aspect ratio of plot
  res@vpWidthF      = 0.8

  res@gsnDraw          = False
  res@gsnFrame         = False

  res@tiMainString  = "Oscillation index of SAH"    ; title
  res@tiYAxisString = "Oscillation index"                               ; y axis title
  res@tiXAxisString = "Year"                                 ; x axis title

  res@xyLineColor       = "Black"
  res@xyLineThicknessF  = 6.0
  res@xyDashPattern     = 0.0

;  res@xyMarkLineModes   = (/"Lines","Lines","Lines","Lines","Lines","Lines"/)      ; whether each line on the plot is
;  res@xyMarkers         = (/16,    16,     16,     16,     16,     16/)            ; represented soly by a line, by markers, or both.
;  res@xyMarkerSizes     = (/0.012,  0.012,  0.012,  0.012,  0.012,  0.012/)        ; represented soly by a line, by markers, or both.
;
;  res@xyMarkerColors    = (/"Red","Blue","Green","Orange","Yellow","Purple"/)     ; color

;-----------------------------------------------
; Set the resources of x/y axis
;-----------------------------------------------
  res@tmXBMode             = "Explicit"
  res@tmXBValues           = (/0,4,8,12,16,20,\
                               24,28,32,36,40/)
  res@tmXBLabels           = (/"1979","1983","1987",\
                               "1991","1995","1999",\
                               "2003","2007","2011",\
                               "2015","2019"/)
  res@tmXBMinorValues      = (/1,2,3,5,6,7,9,10,11,13,14,15,17,18,19,\
                               21,22,23,25,26,27,29,30,31,33,34,35,37,38,39/)
;  res@tmXBMinorLengthF     = 0.01     ; It's a dynamic value, which is controled
                                       ; by NCL itself, you can modify too!
;  res@tmXBMinorOutwardLengthF  = 0.0

  res@tmXBLabelAngleF       = 0.0

;  res@trXMinF              = 0
;  res@trXMaxF              = dimsizes(time_fcst)-1

  res@trYMinF              = -3
  res@trYMaxF              = 3

;-------------------------------------------------------------
; Set the resources of legends
;-------------------------------------------------------------
;  res@pmLegendDisplayMode    = "Always"       ; turn on legend
;  res@pmLegendSide           = "Bottom"       ; set the location
;  res@pmLegendParallelPosF   = .85            ; adjust the location of legend
;  res@pmLegendOrthogonalPosF = -0.55          ; adjust the location of legend
;  res@pmLegendWidthF         = 0.10           ; resize the legend
;  res@pmLegendHeightF        = 0.12           ; resize the legend
;  res@lgPerimOn              = False          ; turn off box around
;  res@lgLabelFontHeightF     = 0.015           ; label font height
;  res@xyExplicitLegendLabels = (/"SLAB","NOAH","SLAB_SM","NOAH_SM"/)

;---------------------------------------------------------------
; Set the resources of tickmarks
;---------------------------------------------------------------
  res@tmYROn = False     ; Remove the tickmarks of right y axis
  res@tmXTOn = False     ; Remove the tickmarks of top x axis
  res@tmYLOn = True      ; Retain the tickmarks of left y axis
  res@tmXBOn = True      ; Retain the tickmarks of bottom x axis


;-----------------------------------------------
; Add a referece line to represent fog period
;-----------------------------------------------
  res@gsnYRefLine              = (/-0.5,0.0,0.5/)
  res@gsnYRefLineDashPatterns  = (/16,0,16/)
  res@gsnYRefLineThicknessF    = 3
  res@gsnYRefLineColor         = "grey"


  plot = gsn_csm_xy(wks,time_label,Index,res)


;---------------------------------------------
; Add text to represent figure number
;---------------------------------------------
;  txres_number                         = True                     ; text mods desired
;  txres_number@txFontHeightF           = 0.03                     ; font smaller. default big
;  txres_number@txBackgroundFillColor   = "transparent"
;
;  dum_number   = gsn_add_text(wks,plot,"(a)",1.0,900,txres_number)


;--------------------------------
; These two step are required!
;--------------------------------
  draw(plot)                     ; note we are drawing the first one!
  frame(wks)  ; Now that we are done drawing, draw the frame



end




