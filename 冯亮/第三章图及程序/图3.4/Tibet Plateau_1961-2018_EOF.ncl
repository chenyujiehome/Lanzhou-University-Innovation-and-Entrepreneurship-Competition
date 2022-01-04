begin

    f = addfile("1961-2018-month-pr.nc", "r")
    pr = f->pr;(48:647,:,:)
    printVarSummary(pr)

;*********EOF**************
    neof = 3
    ;pr_avg = clmMonTLL(pr);计算逐月平均值
    ;pr := calcMonAnomTLL(pr, pr_avg);计算距平
    pr := month_to_season(pr, "JJA")
    printVarSummary(pr)
    pr := dim_rmvmean_n_Wrap(pr,0)
    eof = eofunc_n_Wrap(pr,neof,False,0);这两个函数官网有详细说明 经验正交数据并保留元数据
    eof_ts = eofunc_ts_n_Wrap(pr,eof,False,0);这两个函数官网有详细说明 计算eof的时间序列并保留元数据
    eof_ts = dim_standardize_n(eof_ts,0,1);时间序列标准化
    printVarSummary(eof_ts)

;*********plot*******************
    output = "JJA_pr_eof"
    wks = gsn_open_wks("png",output)        ; send graphics to PNG file
    gsn_define_colormap(wks, "BlWhRe")

    plot = new(neof,graphic)                ; create graphic array
                                          ; only needed if paneling
    res                      = True         
    res@gsnDraw              = False        ; don't draw yet
    res@gsnFrame             = False        ; don't advance frame yet
    res@gsnAddCyclic         = False        ; data not cyclic
    
    res@cnFillPalette        = "BlWhRe"   ; choose colormap
    res@cnFillOn             = True         ; turn on color fill
    res@cnLinesOn            = True         ; True is default
    res@lbLabelBarOn         = False        ; turn off individual lb's

    res@mpCenterLonF        = 180.
    res@mpMinLatF           = min(pr&lat)
    res@mpMaxLatF           = max(pr&lat)
    res@mpMinLonF           = min(pr&lon)
    res@mpMaxLonF           = max(pr&lon)


    res@mpDataBaseVersion       = "MediumRes"
    res@mpDataSetName           = "Earth..4"
    res@mpFillDrawOrder      = "PreDraw"
    res@mpGeophysicalLineColor = "white"

    symMinMaxPlt(eof, 30, False, res); contributed.ncl

; panel plot only resources
  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar

  ;resP@gsnPanelMainString  = "1961-2018-EOF"
  do n=0,neof-1
     res@gsnLeftString  = "EOF "+(n+1)
     res@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
     plot(n) = gsn_csm_contour_map(wks,eof(n,:,:),res)
  end do
  gsn_panel(wks,plot,(/neof,1/),resP)     ; draw all 'neof' as one plot

;*******************************************
; time series (principal component) plot
;*******************************************
  eof_ts@long_name = "Amplitude"

  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
  rts@gsnScale  = True        ; force text scaling               

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

  rts@vpHeightF = 0.40        ; Changes the aspect ratio
  rts@vpWidthF  = 0.85
  rts@vpXF      = 0.10        ; change start locations
  rts@vpYF      = 0.75        ; the plot

  rts@gsnYRefLine           = 0.              ; reference line   
  ;rts@gsnXYBarChart         = True  
  rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
  rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue

; panel plot only resources
  rtsP                     = True             ; modify the panel plot
  rtsP@gsnMaximize         = True             ; large format
      
  ;resP@gsnPanelMainString  = "SCRIPPS: Pierce"
  ;这里时间设置要注意
  yyyymm = ispan(1961,2018,1);横坐标时间设置
  ;rts@tmXBMode = "Explicit"
  ;rts@tmXBValues = (/0,120,240,360,480,600/)
  ;rts@tmXBLabels = (/196501,197501,198501,199501,200501,201501/)

  do n=0,neof-1
     rts@gsnLeftString  = "EOF "+(n+1)
     rts@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
     plot(n) = gsn_csm_xy (wks,yyyymm,eof_ts(n,:),rts)
  end do                                     
  gsn_panel(wks,plot,(/neof,1/),rtsP)        ; draw all 'neof' as one plot

  ;输出eof_ts(0)
;  asciiwrite("eof_ts0.txt", (eof_ts(0,:)))
end
