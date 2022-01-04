;--------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "shapefile_utils.ncl"
;--------------------------------------------------------------------------------------

begin

  Year    = 1981

  dir     = "./"


;---------------------------------------------------------------------------
; Read the data in spring and winter
;---------------------------------------------------------------------------
  f            = addfile(dir + "snow_" + Year + ".nc", "r")

  snow_spring  = f->snow_spring({45:25},{65:105})
  snow_winter  = f->snow_winter({45:25},{65:105})

  printVarSummary(snow_spring)
  printVarSummary(snow_winter)


;-------------------------------------------
; Mask the data outside of Tibetan Plateau
;-------------------------------------------
  Tibet_shp       = "Tibet.shp"
  mask_filename   = "Mask_Tibet.nc"

  opt             = True
  opt@return_mask = True    ; This forces the return of a 0s and 1s mask array
  mask_array1     = shapefile_mask_data(snow_winter,Tibet_shp,opt)
  mask_array2     = shapefile_mask_data(snow_spring,Tibet_shp,opt)

  snow_winter_mask  = where(mask_array1.eq.1,snow_winter,-1.0)   ; Keep the data over TP only
  snow_spring_mask  = where(mask_array2.eq.1,snow_spring,-1.0)   ; Keep the data over TP only

  copy_VarAtts(snow_winter,snow_winter_mask)
  copy_VarCoords(snow_winter,snow_winter_mask)
  snow_winter_mask@_FillValue   = -1.0

  copy_VarAtts(snow_spring,snow_spring_mask)
  copy_VarCoords(snow_spring,snow_spring_mask)
  snow_spring_mask@_FillValue   = -1.0

  printVarSummary(snow_spring_mask)


;---------------------------------------
; Start to plot:
;---------------------------------------
  plot = new(2,graphic)                       ; create a plot array
  wks = gsn_open_wks("x11","snow")
  gsn_define_colormap(wks,"precip3_16lev")

  res                      = True
  res@gsnDraw             = False           ; don't draw
  res@gsnFrame            = False           ; don't advance frame

;  res@vpXF          = 0.17
;  res@vpYF          = 0.8
;  res@vpHeightF     = 0.4               ; change aspect ratio of plot
;  res@vpWidthF      = 0.5

  res@tiMainString    = ""
  res@gsnCenterString = ""

  res@tiMainFontHeightF           = 0.020

  res@tiYAxisFontHeightF  = 0.030
  res@tiXAxisFontHeightF  = 0.030
  res@gsnLeftStringFontHeightF    = 0.025
  res@gsnCenterStringFontHeightF  = 0.025
  res@gsnRightStringFontHeightF   = 0.025

  res@tmXBLabelFontHeightF     = 0.025
  res@tmXTLabelFontHeightF     = 0.025
  res@tmYLLabelFontHeightF     = 0.025
  res@tmYRLabelFontHeightF     = 0.025

  res@cnFillOn             = True                ; turn on color
  res@gsnSpreadColors      = True                ; use full colormap
  res@cnLinesOn            = False               ; turn off contour lines
  res@cnLineLabelsOn       = False
  res@cnInfoLabelOn        = False
  res@lbLabelAutoStride    = True                ; automatic lb label stride
  res@lbLabelBarOn         = False

  res@cnLevelSelectionMode = "ManualLevels"      ; manually set cn levels
  res@cnMinLevelValF       = 2                 ; min level
  res@cnMaxLevelValF       = 22                ; max level
  res@cnLevelSpacingF      = 2                 ; contour level spacing

  res@mpLimitMode          = "LatLon"            ;Default is LatLon
  res@mpMinLonF            = 65
  res@mpMaxLonF            = 105
  res@mpMinLatF            = 25
  res@mpMaxLatF            = 45
  res@gsnAddCyclic         = False

;  res@mpFillDrawOrder             = "PostDraw"
  res@mpOutlineOn                 = True         ; True or False
                                                 ;produce different plot!!!
  res@mpOutlineBoundarySets       = "National"   ; turn on country boundaries
  res@mpNationalLineColor         = "Black"
  res@mpGeophysicalLineColor      = "Black"
  res@mpGridAndLimbOn             = False
  res@mpGridLineColor             = "Black"
  res@mpLimbLineColor             = "Black"
  res@mpNationalLineColor         = "Black"
  res@mpPerimLineColor            = "Black"
;  res@mpUSStateLineColor          = "Black"
;  res@mpDataSetName               = "Earth..4"   ; divisions for other countries.
;  res@mpDataBaseVersion           = "MediumRes"  ; Medium resolution database
;  res@mpOutlineSpecifiers      = (/"China:states","Taiwan","United States:states"/)
;  res@mpOutlineBoundarySets    = "National"
;  res@mpFillAreaSpecifiers        = (/"United States","Canada","Mexico"/)
;  res@mpSpecifiedFillColors       = (/"white","white","white"/)
;  res@mpUSStateLineColor          = "Black"
;  res@mpGeophysicalLineThicknessF = 1.5   ; for Offshore Line thickness
;  res@mpNationalLineThicknessF    = 1.5   ; for National Line thickness
;  res@mpUSStateLineThicknessF     = 1.2   ; for Province Line thickness

;------------------------------------------------------
; These 4 lines are the newly learnt usefull map info
;------------------------------------------------------
  res@mpGeophysicalLineColor = "Black"     ; color of continental outlines
;  res@mpUSStateLineColor     = "Black"
  res@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed
  res@mpOutlineBoundarySets  = "GeophysicalAndUSStates"
  res@mpUSStateLineColor     = "Black"
  res@mpGridAndLimbOn        = False
  res@mpGridLatSpacingF =  2.            ; change latitude  line spacing
  res@mpGridLonSpacingF =  2.            ; change longitude line spacing
  res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
  res@tmXBTickSpacingF  = 5.0            ; change tickmark interval
  res@tmYLTickSpacingF  = 5.0            ; change tickmark interval
;------------------------------------------------------

;  res@mpFillOn                    = True
;  res@mpFillAreaSpecifiers        = (/"Water","China","Taiwan","Land"/)
;  res@mpSpecifiedFillColors       = (/"white","transparent","transparent","white"/)
;  res@mpOutlineSpecifiers         = (/"China:states"/)
;  res@mpMaskOutlineSpecifiers     = (/"water"/)

;---------------------------------------------------------------
; Set the resources of tickmarks
;---------------------------------------------------------------
  res@tmYROn = False     ; Remove the tickmarks of right y axis
  res@tmXTOn = False     ; Remove the tickmarks of top x axis
  res@tmYLOn = True      ; Retain the tickmarks of left y axis
  res@tmXBOn = True      ; Retain the tickmarks of bottom x axis


  res@lbTitleOn                 = False
  res@lbTitleString             = "PBL height (m)"
  res@lbTitleFontHeightF        = 0.015
  res@pmLabelBarOrthogonalPosF  = 0.2
;  res@pmLabelBarParallelPosF    = 0.6
;  res@pmLabelBarWidthF          = 0.1
;  res@pmLabelBarHeightF         = 0.6


;------------------------------------------------------
; For totrain heat flux
;------------------------------------------------------
  opts_obs1                 = res
  opts_obs1@gsnLeftString   = "Winter"
  opts_obs1@gsnRightString  = ""
  plot(0) = gsn_csm_contour_map_ce(wks,snow_winter_mask,opts_obs1)

  opts_ctl1                 = res
  opts_ctl1@gsnLeftString   = "Spring"
  opts_ctl1@gsnRightString  = ""
  plot(1) = gsn_csm_contour_map_ce(wks,snow_spring_mask,opts_ctl1)


;-------------------------------------------
; Attach the shapefile of Tibetan Plateau
;-------------------------------------------
  Tibet_shp              = "/mnt/zfm/TPV_postprocess/Run_info/shp_tp/Tibet.shp"

  lnres                  = True
  lnres@gsLineColor      = "Black"
  lnres@gsLineThicknessF = 3

  Tibet_id0  = gsn_add_shapefile_polylines(wks,plot(0),Tibet_shp,lnres)
  Tibet_id1  = gsn_add_shapefile_polylines(wks,plot(1),Tibet_shp,lnres)

;------------------------------------------------
; Create panel
;------------------------------------------------
  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance panel plot
  resP@gsnPanelLabelBar = True                  ; do not use common colorbar
  resP@pmLabelBarHeightF      = 0.08
  resP@pmLabelBarWidthF       = 0.60
  resP@lbLabelFontHeightF     = 0.015
  resP@txString         = "   "    ; "Bill total precipitation"
  resP@gsnMaximize      = True                   ; use full page
  resP@gsnPanelBottom   = 0.05                   ; save space at bottom
  resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)","(e)","(f)"/)
  resP@amJust               = "TopLeft"
  resP@amOrthogonalPosF     = -0.5
  resP@amParallelPosF       = -0.5
  resP@gsnPanelFigureStringsFontHeightF  = 0.015
  resP@gsnPanelFigureStringsPerimOn = False
  resP@gsnPanelFigureStringsBackgroundFillColor = "white"
  resPgsnPanelYWhiteSpacePercent = 4
  resPgsnPanelXWhiteSpacePercent = 2
  resP@gsnPanelRowSpec           = True
  gsn_panel(wks,plot,(/1,2/),resP)


frame(wks) ; Advance the frame.


end



