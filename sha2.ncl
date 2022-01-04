;--------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "/mnt/zfm/TPV_postprocess/Run_info/shp_tp/shapefile_utils.ncl"
;--------------------------------------------------------------------------------------

begin
do Year =1981,2018
  dir     = "./"


;---------------------------------------------------------------------------
; Read the data in spring and winter
;---------------------------------------------------------------------------
  f            = addfile(dir + "snow_" + Year + ".nc", "r")

  snow_spring  = f->snow_spring({45:25},{65:105})
  snow_winter  = f->snow_winter({45:25},{65:105})


;-------------------------------------------
; Mask the data outside of Tibetan Plateau
;-------------------------------------------
  Tibet_shp       = "/mnt/zfm/TPV_postprocess/Run_info/shp_tp/Tibet.shp"
  mask_filename   = "Mask_Tibet.nc"

  opt             = True
  opt@return_mask = True    ; This forces the return of a 0s and 1s mask array
  mask_array1     = shapefile_mask_data(snow_winter,Tibet_shp,opt)
  mask_array2     = shapefile_mask_data(snow_spring,Tibet_shp,opt)

  snow_winter_mask  = where(mask_array1.eq.1,snow_winter,0.0)   ; Keep the data over TP only
  snow_spring_mask  = where(mask_array2.eq.1,snow_spring,0.0)   ; Keep the data over TP only
 snow_mask  = where(mask_array2.eq.1,1,0)   ; Keep the data over TP only
  
copy_VarAtts(snow_winter,snow_winter_mask)
  copy_VarCoords(snow_winter,snow_winter_mask)
  snow_winter_mask@_FillValue   = -1.0

copy_VarAtts(snow_winter,snow_mask)
  copy_VarCoords(snow_winter,snow_mask)
  snow_winter_mask@_FillValue   = -1.0

  copy_VarAtts(snow_spring,snow_spring_mask)
  copy_VarCoords(snow_spring,snow_spring_mask)
  snow_spring_mask@_FillValue   = -1.0
s_wi=0.0
s_sp=0.0
ma=0.0
do lo =0,160
do la =0,80
s_wi=s_wi+snow_winter_mask(la,lo)
s_sp=s_sp+snow_spring_mask(la,lo)
ma=ma+snow_mask(la,lo)
end do
end do
av_wi=s_wi/(ma)
av_sp=s_sp/(ma)
print(av_wi)
print(av_sp)
end do



end



