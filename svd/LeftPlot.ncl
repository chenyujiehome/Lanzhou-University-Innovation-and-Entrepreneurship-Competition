    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

nsvd=3
 f    = addfile("LeftSvd.nc", "r")
 sic=f->sic


 
  wks = gsn_open_wks("png","LhomoSvd")         ; send graphics to PNG file
  plot = new(nsvd,graphic)                ; create graphic array
                                          ; only needed if paneling
; svd patterns

  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@gsnPolar             = "NH"

  res@mpFillOn             = False        ; turn off map fill
  res@mpMinLatF            = 66.5         ; zoom in on map
  res@mpMaxLatF            = 90.0

  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "BlWhRe"     ; choose colormap
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False        ; turn off individual lb's

                                          ; set symmetric plot min/max
  symMinMaxPlt(sic, 16, False, res)       ; contributed.ncl

; panel plot only resources
  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar

;*******************************************
; Plot standard patterns
;*******************************************
  do n=0,nsvd-1
     res@gsnLeftString  = "svd "+(n+1)
     plot(n)=gsn_csm_contour_map_polar(wks,sic(n,:,:),res)
  end do
  gsn_panel(wks,plot(0:2),(/3,1/),resP)     ; only plot the 1st four
 ;gsn_panel(wks,plot,(/nsvd/2,2/),resP)     ; now draw as one plot
 





end
