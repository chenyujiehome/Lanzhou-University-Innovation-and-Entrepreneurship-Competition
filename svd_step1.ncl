;--------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "/mnt/zfm/TPV_postprocess/Run_info/shp_tp/shapefile_utils.ncl"
;--------------------------------------------------------------------------------------

begin

  Year    = 1981

  dir     = "./"


;---------------------------------------------------------------------------
; Read the data in spring and winter
;---------------------------------------------------------------------------
  f            = addfile(dir + "snow_" + Year + ".nc", "r")

  snow_spring  = f->snow_spring({25:45},{65:105})
  snow_winter  = f->snow_winter({25:45},{65:105})

  printVarSummary(snow_spring)
  printVarSummary(snow_winter)


;-------------------------------------------
; Mask the data outside of Tibetan Plateau
;-------------------------------------------
  Tibet_shp       = "/mnt/zfm/TPV_postprocess/Run_info/shp_tp/Tibet.shp"
  mask_filename   = "Mask_Tibet.nc"

  opt             = True
  opt@return_mask = True    ; This forces the return of a 0s and 1s mask array
  mask_array1     = shapefile_mask_data(snow_winter,Tibet_shp,opt)
  mask_array2     = shapefile_mask_data(snow_spring,Tibet_shp,opt)

  snow_winter_mask  = where(mask_array1.eq.1,snow_winter,-999.0)   ; Keep the data over TP only
  snow_spring_mask  = where(mask_array2.eq.1,snow_spring,-999.0)   ; Keep the data over TP only

  copy_VarAtts(snow_winter,snow_winter_mask)
  copy_VarCoords(snow_winter,snow_winter_mask)
  snow_winter_mask@_FillValue   = -999.0

  copy_VarAtts(snow_spring,snow_spring_mask)
  copy_VarCoords(snow_spring,snow_spring_mask)
  snow_spring_mask@_FillValue   = -999.0
do z=0,160
    do j = 0, 80

            ; write loop content
            print(snow_winter_mask(j,z))
        ; write loop content
    end do
    
    
    ; write loop content
end do

;---------------------------------------

end



