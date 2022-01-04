
begin
   diri = "/home/heyl/tianwenlong/fengliang/shuju/"
  shptp="/home/heyl/tianwenlong/fengliang/nc_coutour/add_TP/add_TP/DBATP_Line.shp"
  tistrat = 196101
   tiend = 201812

   f1 = addfile(diri+"uwnd.mon.mean.nc", "r")
   ;level = f1->level
   ;print(level)(9)     200
   ;exit
   time = f1->time
   yyyymm = cd_calendar(time,-1)
   u = f1->uwnd(ind(yyyymm.eq.tistrat):ind(yyyymm.eq.tiend),9,{0:180},{0:180})
   u := month_to_season(u, "JJA")
   printVarSummary(u)

   f2 = addfile(diri+"vwnd.mon.mean.nc","r")
   v = f2->vwnd(ind(yyyymm.eq.tistrat):ind(yyyymm.eq.tiend),9,{0:180},{0:180})
   v := month_to_season(v, "JJA")
   printVarSummary(v)

   eof_ts0 = asciiread(diri+"eof_ts1.txt", -1, "float")
   printVarSummary(eof_ts0)

   wind = wind_speed(u, v)
   printVarSummary(wind)
   ;exit
   ;计算相关
   rc = escorc(eof_ts0,wind(lat|:,lon|:,time|:))
   copy_VarCoords(wind(0,:,:),rc)
   printVarSummary(rc)

    df   = 58-2
    t    = rc*sqrt((58-2)/(1-rc^2))      
    prob = student_t(t, df)
    prob = where(prob.eq.0,prob@_FillValue,prob)
    copy_VarCoords(rc,prob)
    printVarSummary(prob)
    print(max(prob))
    print(min(prob))
    ;exit

    wks =  gsn_open_wks("ps","relation_JJA")
    gsn_define_colormap(wks, "cmp_b2r")

 ;************************************************
 ; plotting parameters（画相关系数rc） 
 ;************************************************
   
   res                       = True 
   res@gsnDraw               = False
   res@gsnFrame              = False    
   
   res@gsnMaximize           = True             ; make large
   res@gsnAddCyclic          = False
   
   res@cnFillOn              = True             ; turn on color填色等值线图
   res@cnLinesOn             = False            ; turn off contour lines关闭等值线
   res@cnLineLabelsOn        = False            ; turn off contour line labels关闭等值线标签
   ;;res@cnFillMode            = "RasterFill"
   res@cnLevelSelectionMode  = "ManualLevels"   ; set manual contour levels
   res@cnMinLevelValF        =  -0.5             ; set min contour level
   res@cnMaxLevelValF        =   0.5              ; set max contour level
   res@cnLevelSpacingF       =   0.05            ; set contour interval
  
   res@mpOutlineBoundarySets       = "National" ;地理边界线 
    res@mpCenterLonF        = 180.
    res@mpMinLatF           = 0
    res@mpMaxLatF           = 90
    res@mpMinLonF           = 0
    res@mpMaxLonF           = 180

   ;res@mpCenterLonF                = 180.
   res@mpGeophysicalLineThicknessF = 2.         ; double the thickness of geophysical boundaries地图边界线宽度
   res@mpNationalLineThicknessF    = 2.        ; double the thickness of national boundaries国家边界线宽度
   res@mpOutlineOn                 = True       ; turn off outline
   res@mpFillOn                    = False    ; turn off map fill 关闭地图填色   

   res@tmXTOn   = False
   res@tmYROn  = False
   res@tmXBOn = True
   base = gsn_csm_contour_map(wks,rc,res)  
 ;============================================<
 ;add 画通过显著性检验的区域
 ;============================================<
   res2 = True
   res2@gsnDraw = False;不画这个图
   res2@gsnFrame = False;不画这个图的框架，如果只改了上面那句话没有改这句话，会出现一个有经纬度外框的空白图
   
   res2@cnFillOn = True
   res2@cnLinesOn = False;在这里不能是pattern，否则下面的109，110，111这三句话会出错
   res2@cnLineLabelsOn = False
   res2@cnInfoLabelOn = False;关闭等值线信息标签,即右下角显示的contour from x to x by x
   res2@lbLabelBarOn = False;关闭色标
   res2@cnMonoFillPattern = False
   res2@cnLevelSelectionMode = "ExplicitLevels"
   res2@cnLevels = (/0.1/);set to significance level设置显著性水平
   res2@cnFillPatterns = (/17,-1/);-1表示透明,res2@cnLevels = (/0.05/)把区域分成了>0.05（用7形状填充）和<0.05（用-1形状填充）两个区域
   res2@cnFillColors = (/1,-1/);-1表示透明,<0.05的区域用1填色，>0.05的区域用-1填色
   res2@gsnLeftString = ""
   res2@cnFillScaleF = 0.8
   res2@cnFillDotSizeF = 0.003
   plot = gsn_csm_contour(wks,prob,res2)
tpres        = True
tpres@gsLineColor = "black"     ;设置廓线颜色
tpres@gsLineThicknessF = 8              ;设置廓线宽度
   shp_plot2     = gsn_add_shapefile_polylines(wks,plot,shptp,tpres)    
   overlay(base, plot) 
   draw(base)
   frame(wks)
   end
