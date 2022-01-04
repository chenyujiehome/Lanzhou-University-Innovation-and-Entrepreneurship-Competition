;--------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;--------------------------------------------------------------------------------------

begin




  dir_winter    = "./"
 

  input_winter  = asciiread(dir_winter + "file_list.dat",-1,"string")


;---------------------------------------------------------------------------
; Read the data in spring and winter
;---------------------------------------------------------------------------
  nrows        = 161
  ncols        = 321
  snow_winter  = new((/dimsizes(input_winter),nrows,ncols/),"float")
do m=0,dimsizes(input_winter)-1

snow_winter(m,:,:)=readAsciiTable(dir_winter + input_winter(m),ncols,"float",(/6,0/))
end do
end



