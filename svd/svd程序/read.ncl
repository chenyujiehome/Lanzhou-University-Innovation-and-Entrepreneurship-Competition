;  Grads control file
; ---------------
; dset snow.grd
; title snow depth  winter
; options sequential big_endian
; undef -9999.00               
; xdef 29 linear 20. 5.
; ydef 8 linear 35.  5.
; zdef 1 linear 1 1
; tdef 60 linear jan1936 1yr
; vars 1
; p 0 0 precipitation 
; endvars
; ---------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
begin
  
  diri  = "./"                            ; input directory
  fili  = "h.grd"                      ; DSET
  fName = diri+fili                       ; path

  nlat  = 13                               ; YDEF
  nlon  = 17                             ; XDEF
  nlev  = 1
  year  = 1951                          ; TDEF
  ntim  = 60                              ; time steps

  
                                          ; create an array to contain data
  UNDEF = -9999.                          ; UNDEF
  x     = new ( (/nlon,nlat,ntim/), float, UNDEF)

  x@long_name = "Snow Depth"              ; VARS
  x@units     = "mm"                 
                                    
  setfileoption("bin","ReadByteOrder","BigEndian")
                                         
  do nt=0,ntim-1                          ; read each record: store in x
     x(:,:,nt) = fbinrecread(fName, nt, (/nlon,nlat/), "float")
  end do

  printVarSummary(x)
  print ("min(x)="+min(x))
  print ("max(x)="+max(x))
end