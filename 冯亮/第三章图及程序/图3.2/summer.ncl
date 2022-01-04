   ncol  = 2
   ntim  = numAsciiRow("summer.txt")
   data  = asciiread("summer.txt", (/ntim,ncol/), "float")
   data@_FillValue = -9999.0 

   x     = data(:,0)              ; model time: units are days              
   y     = data(:,1)              ; model value: units degK

;************************************************
; calculate the regression coefficient (slope)    
;************************************************
   rc    = regline(x, y)          ; slope           
;   rc@units = "degK/day"         
   print(rc)

;************************************************
; create an array to hold both the original data
; and the calculated regression line   
;     ---------
;     y = mx+b 
;     m is the slope:       rc      returned from regline
;     b is the y intercept: rc@yave attribute of rc returned from regline
;************************************************
   pltarry   = new ( (/2,ntim/), typeof(data), data@_FillValue)

   pltarry(0,:) = y                                 ; use markers
   pltarry(1,:) = rc*x + rc@yintercept              ; use solid line

;************************************************
; plotting parameters
; This illustrates one approach. Overlays could also be used.
;************************************************
   wks  = gsn_open_wks("ps","summer")             ; send graphics to PNG file
   
   res                     = True                   ; plot mods desired
;   res@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
;   res@xyMarkers           = 16                     ; choose type of marker 
;  res@xyMarkerColor       = "red"                  ; Marker color
;   res@xyMarkerSizeF       = 0.005                  ; Marker size (default 0.01)
;   res@xyDashPatterns      = 1                      ; solid line 
;   res@xyLineThicknesses   = (/1,2/)                ; set second line to 2
;   res@tmYLFormat          = "f"                    ; not necessary but nicer labels 

res@tiMainString        = "Summer Precipitation"  ; title
res@tiXAxisString="year"
res@tiYAxisString="precipitation(mm/a)"
res@xyLineThicknessF= 5.0
res@xyLineColors      = (/"black","red"/)
   plot  = gsn_csm_xy (wks,x,pltarry,res)           ; create plot

