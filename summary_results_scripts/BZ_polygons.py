# -*- coding: utf-8 -*-
"""
Generate buffer zone polygons from surveillance results to create graphs

@author: Martina Cendoya
"""

import model_functions as fn
import pandas as pd
import numpy as np
from shapely.ops import unary_union
import geopandas as gpd
# enable shapely.speedups which makes some of the spatial queries running faster
import shapely.speedups
shapely.speedups.enable()


xy = pd.read_csv("coords.csv")
coords = xy[['x', 'y']].values
data = gpd.GeoDataFrame(xy, geometry=gpd.points_from_xy(xy.x, xy.y))

survey = ["one_step_00", "one_step_01", "two_step_00", "two_step_01"]
r_zt = [2500, 5000]
r_er = [50, 100]
r_b = [0, 50, 90]

for surv_t in survey:
    for zt in r_zt:
        for er in r_er:
            for b in r_b:
                surv = np.loadtxt("./results/"+str(surv_t)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/surveys_median.txt")
                surv = np.transpose(surv)
                b_zt, bound_zt, bound_e = gpd.GeoSeries(), gpd.GeoSeries(), gpd.GeoSeries()
                for i in range(len(surv)):
                    positivos_detectados = np.where(surv[i]>0)[0]
                    b_zt_i = fn.bound_delim(centroid=positivos_detectados, coords=coords, r=zt+50)
                    b_zt = pd.concat([b_zt, b_zt_i])
                    b_zt = gpd.GeoSeries(unary_union(b_zt))
                    bound_zt = pd.concat([bound_zt, b_zt])
                bound_zt.to_file("./results/"+str(surv_t)+"/zt"+str(zt)+"_e"+str(er)+"_b"+str(b)+"/bound_zt.shp")

