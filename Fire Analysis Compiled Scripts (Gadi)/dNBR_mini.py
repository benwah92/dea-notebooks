import sys
import datacube
from datacube.helpers import write_geotiff
from datetime import datetime
from datetime import timedelta
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("/scratch/wj97/ab4513/dea-notebooks/Scripts")
import dea_datahandling
from dea_datahandling import load_ard
from dea_plotting import rgb
from dea_plotting import display_map
from dea_bandindices import calculate_indices

print(datacube.__version__)

dc = datacube.Datacube(app="dNBR")

central_lat = -35.783333
central_lon = 148.016667
crs = 'EPSG:32755'

buffer = 0.2

study_area_lat = (central_lat - buffer, central_lat + buffer)
study_area_lon = (central_lon - buffer, central_lon + buffer)

prefire_start = '2019-11-01'
prefire_end = '2020-01-06'
postfire_start = '2020-01-07'
postfire_end = '2020-05-01'

query_1 = {"x": (central_lon - buffer, central_lon + buffer),
         "y": (central_lat - buffer, central_lat + buffer),
         "time": (prefire_start, prefire_end),
         "output_crs": "EPSG:32755",
         "resolution": (-10, 10)}

prefire_data = load_ard(dc=dc,
    products=['s2a_ard_granule', 's2b_ard_granule'],
    measurements=['nbart_nir_1', 'nbart_swir_3'],
    min_gooddata=0,
    # dask_chunks={'x': 'auto', 'y': 'auto'},

    group_by='solar_day',**query_1)

prefire_image = prefire_data.median(dim='time')
prefire_image = calculate_indices(prefire_image,
                                  index='NBR', 
                                  collection='ga_s2_1', 
                                  drop=False)
prefire_burnratio = prefire_image.NBR
prefire_burnratio.data

query_2 = {"x": (central_lon - buffer, central_lon + buffer),
         "y": (central_lat - buffer, central_lat + buffer),
         "time": (postfire_start, postfire_end),
         "output_crs": "EPSG:32755",
         "resolution": (-10, 10)}

postfire_data = load_ard(dc=dc,
    products=['s2a_ard_granule', 's2b_ard_granule'],
    measurements=['nbart_nir_1', 'nbart_swir_3'],
    min_gooddata=0,
    # dask_chunks={'x': 'auto', 'y': 'auto'},

    group_by='solar_day',**query_2)

postfire_image = postfire_data.median(dim='time')
postfire_image = calculate_indices(postfire_image,
                                  index='NBR', 
                                  collection='ga_s2_1', 
                                  drop=False)
postfire_burnratio = postfire_image.NBR
postfire_burnratio.data

delta_NBR = prefire_burnratio - postfire_burnratio
# delta_NBR.compute()
dnbr_dataset = delta_NBR.to_dataset(name='delta_NBR')
write_geotiff(f'./scratch/wj97/ab4513/tumbarumba_dNBR_mini.tif', dnbr_dataset)