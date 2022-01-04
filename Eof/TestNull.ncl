
  f      = addfile ("EofSnowData.nc", "r")



  snow    = f->snow
  print(snow(1,{30},{85}))