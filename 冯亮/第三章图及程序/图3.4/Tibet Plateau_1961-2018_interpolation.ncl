begin
;*********读入数据*************
    f1 = "pre_month_tibet_plateau_66stas_1961-2018.txt" 
    nrow=numAsciiRow(f1)
    ncol=numAsciiCol(f1)
    data = asciiread(f1,(/nrow,ncol/),"float")
    ;printVarSummary(data);66*696 58yrs

    f2 = "66_stas_num_lat_lon.txt"
    nrow:=numAsciiRow(f2)
    ncol:=numAsciiCol(f2)
    station = asciiread(f2,(/nrow,ncol/),"float")
    printVarSummary((station))
    lat = station(:,1)
    lon = station(:,2)

;********插值成0.5分辨率格点数据*********
    minlat      =   25
    maxlat      =   40
    minlon      =   75
    maxlon      =   105
    olon = new(61,"float")
    olat = new(31,"float")
    do i=0,60
        olon(i) = minlon+i*0.5;75-105经度差30°
    end do
    do l=0,30
        olat(l) = minlat+l*0.5;25-40维度差15°
    end do

    olon!0          = "lon"
    olon@long_name  = "lon"
    olon@units      = "degrees-east"
    olon&lon        = olon
    olat!0          = "lat"
    olat@long_name  = "lat"
    olat@units      = "degrees_north"
    olat&lat        = olat

    lon@units       = olon@units
    lat@units       = olat@units

    data@_FillValue = -9999.0
    rscan = (/10,5,3/)

    pr= new((/696,31,61/),"float")
    pr!0 = "time"
    do i=0,695
        pr(i,:,:)= obj_anal_ic_Wrap(lon,lat,data(:,i),olon,olat,rscan, False)
    end do
    printVarSummary(pr)

    fout = "1961-2018-month-pr.nc";输出nc文件
    ncdf = addfile(fout, "c")
    ncdf->pr = pr
  
end
