"""
helper functions for icepyx
"""

import icepyx
import xopr
import geopandas as gpd

opr = xopr.opr_access.OPRConnection(cache_dir="/tmp")


antarctic_link = 'https://storage.googleapis.com/opr_stac/reference_geometry/measures_boundaries_4326.geojson'


gdf = gpd.read_file(antarctic_link)
# icesat_gdf = gpd.read_parquet('../data/antarctic_merged_icesat2_tracks.parquet')
# gdf = gdf.dissolve()


region = xopr.geometry.get_antarctic_regions(name='Thwaites', merge_regions=True)

stac_items = opr.query_frames(geometry=region)








def icepyx_query(x_points, roi=None):
    
     
    
    
    return



