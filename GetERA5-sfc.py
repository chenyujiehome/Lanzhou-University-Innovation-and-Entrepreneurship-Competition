import cdsapi
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type':'reanalysis',
        'grid': "1.0/1.0",
        'format':'grib',
        'variable':[
            '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
            '2m_temperature','land_sea_mask','mean_sea_level_pressure',
            'sea_ice_cover','orography','sea_surface_temperature','skin_temperature',
            'snow_depth','snow_density','soil_temperature_level_1','soil_temperature_level_2',
            'soil_temperature_level_3','soil_temperature_level_4','surface_pressure',
            'volumetric_soil_water_layer_1','volumetric_soil_water_layer_2','volumetric_soil_water_layer_3',
            'volumetric_soil_water_layer_4'
        ],
        'date':'20060810/20060816',
        'area':'global',
        'time':[
            '00:00',
            '06:00',
            '12:00',
            '18:00'
        ]
    },
    'ERA5-20060810-20060816-sfc.grib')

