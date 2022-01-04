;--------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;--------------------------------------------------------------------------------------

begin

  Year          = 1981

  dir_spring    = "/home/zhangfm/Wansh/Snow_depth/snowdepth-" + Year + "/spring/"
  dir_winter    = "/home/zhangfm/Wansh/Snow_depth/snowdepth-" + Year + "/winter/"

  input_spring  = asciiread(dir_spring + "file_list.dat",-1,"string")
  input_winter  = asciiread(dir_winter + "file_list.dat",-1,"string")


;---------------------------------------------------------------------------
; Read the data in spring and winter
;---------------------------------------------------------------------------
  nrows        = 161
  ncols        = 321
  snow_spring  = new((/dimsizes(input_spring),nrows,ncols/),"float")
  snow_winter  = new((/dimsizes(input_winter),nrows,ncols/),"float")


  do mm        = 0, dimsizes(input_spring)-1

    f_spring   = dir_spring + input_spring(mm)

    snow_spring(mm,:,:)   = readAsciiTable(f_spring,ncols,"float",(/6,0/))

  end do


  do nn        = 0, dimsizes(input_winter)-1

    f_winter   = dir_winter + input_winter(nn)

    snow_winter(nn,:,:)   = readAsciiTable(f_winter,ncols,"float",(/6,0/))

  end do

  snow_spring@_FillValue  = -1
  snow_winter@_FillValue  = -1
  printVarSummary(snow_spring)


;---------------------------------------------------------------------------
; Do the average of spring and winter
;---------------------------------------------------------------------------
  snow_avg_spring      = dim_avg_n(snow_spring,0)
  snow_avg_winter      = dim_avg_n(snow_winter,0)

  Latitude             = fspan(55,15,161)
  Latitude@long_name   = "latitude"
  Longitude            = fspan(60,140,321)
  Longitude@long_name  = "longitude"

  snow_avg_spring!0    = "lat"
  snow_avg_spring&lat  = Latitude
  snow_avg_spring&lat@units = "degrees_north"  
  snow_avg_spring!1    = "lon"
  snow_avg_spring&lon  = Longitude
  snow_avg_spring&lon@units = "degrees_east"  

  copy_VarCoords(snow_avg_spring,snow_avg_winter)
  copy_VarAtts(snow_avg_spring,snow_avg_winter)

  printVarSummary(snow_avg_spring)


;---------------------------------------------------------------
;  Rewrite to NetCDF-3 file.
;---------------------------------------------------------------
  system("rm -f ./snow_" + Year + ".nc")         ; remove any pre-existing file
  ncdf = addfile("./snow_" + Year + ".nc" ,"c")  ; open output netCDF file

;-----------------------------
; Creat global attributes
;-----------------------------
  fAtt               = True              ; assign file attributes
  fAtt@title         = "Spring and Winter mean snow depth"
  fAtt@Conventions   = "None"
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )               ; copy file attributes
;-------------------------------------------------
; make time an UNLIMITED dimension (recommended)
; You must commnet this at this time !!!
;-------------------------------------------------
  filedimdef(ncdf,"time",-1,True)

;-----------------------------------------------
; Make outputs (some variables are not output)
;-----------------------------------------------
  ncdf->snow_spring       = snow_avg_spring
  ncdf->snow_winter       = snow_avg_spring

  ncdf->lat               = Latitude
  ncdf->lon               = Longitude

end



