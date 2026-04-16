# iterate over all points in intersections_success
# to extract the bed echo power
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import xopr
import holoviews as hv

def extract_layer_peak_power(radar_ds, layer_twtt, margin_twtt):
    """Extract peak power (dB) and its TWTT within a margin around a layer pick."""
    t_start = np.minimum(radar_ds.slow_time.min(), layer_twtt.slow_time.min())
    t_end = np.maximum(radar_ds.slow_time.max(), layer_twtt.slow_time.max())
    layer_twtt = layer_twtt.sel(slow_time=slice(t_start, t_end))
    radar_ds = radar_ds.sel(slow_time=slice(t_start, t_end))
    layer_twtt = layer_twtt.reindex(
        slow_time=radar_ds.slow_time,
        method="nearest",
        tolerance=pd.Timedelta(seconds=1),
        fill_value=np.nan,
    )
    # print("layer_twtt, max", layer_twtt.max().values)
    # print("layer_twtt, min", layer_twtt.min().values)

    start_twtt = layer_twtt - margin_twtt
    end_twtt = layer_twtt + margin_twtt
    data_within_margin = radar_ds.where(
        (radar_ds.twtt >= start_twtt) & (radar_ds.twtt <= end_twtt),
        drop=True,
    )
    if data_within_margin.twtt.size == 0:
        nan_array = xr.full_like(layer_twtt, fill_value=np.nan)
        return nan_array, nan_array
    
    power_dB = 10 * np.log10(np.abs(data_within_margin.Data))
    peak_twtt_index = power_dB.argmax(dim="twtt")
    peak_twtt = power_dB.twtt[peak_twtt_index]
    peak_power = power_dB.isel(twtt=peak_twtt_index)

    peak_twtt = peak_twtt.drop_vars("twtt")
    peak_power = peak_power.drop_vars("twtt")

    return peak_twtt, peak_power

#     )
def extract_basal_echo_power(stac_item, layer, coord, margin_twtt):
    frame = opr.load_frame(stac_item)
    frame = xopr.radar_util.add_along_track(frame)

    x_int, y_int = coord
    frame_proj = xopr.geometry.project_dataset(frame, "EPSG:3031")
    dist = np.sqrt((frame_proj['x'] - x_int)**2 + (frame_proj['y'] - y_int)**2)
    idx_closest = int(np.argmin(dist.data))

    slow_time_closest = frame_proj.slow_time.isel(slow_time=idx_closest).values

    layer_twtt = layer['twtt']
    peak_twtt, peak_power = extract_layer_peak_power(
        frame_proj, layer_twtt, margin_twtt
    )

    pt = float(peak_twtt.sel(slow_time=slow_time_closest, method='nearest').values)
    pp = float(peak_power.sel(slow_time=slow_time_closest, method='nearest').values)

    return pt, pp

def crossover_echo_power(intersections_radio, stac_items_df):
    print(f"Processing {len(intersections_radio)} intersections")

    margin_twtt = 1e-6
    opr = xopr.OPRConnection(cache_dir="radar_cache")


    for idx, intersect_row in intersections_radio.iterrows():
        stac_1 = stac_items_df.loc[intersect_row['id_1']].to_dict()
        stac_2 = stac_items_df.loc[intersect_row['id_2']].to_dict()

        frame_1 = opr.load_frame(stac_1)
        frame_2 = opr.load_frame(stac_2)

        layer_1 = opr.get_layers(frame_1)
        layer_2 = opr.get_layers(frame_2)

        # intersect coordinates
        x_int, y_int, *_ = intersect_row.geometry.coords[0]

        # get layer object
        peak_twtt_surf_1, peak_power_surf_1 = extract_basal_echo_power(stac_1, layer_1["standard:surface"], (x_int, y_int), margin_twtt=margin_twtt)
        peak_twtt_base_1, peak_power_base_1 = extract_basal_echo_power(stac_1, layer_1["standard:bottom"], (x_int, y_int), margin_twtt=margin_twtt)
        peak_twtt_surf_2, peak_power_surf_2 = extract_basal_echo_power(stac_2, layer_2["standard:surface"], (x_int, y_int), margin_twtt=margin_twtt)
        peak_twtt_base_2, peak_power_base_2 = extract_basal_echo_power(stac_2, layer_2["standard:bottom"], (x_int, y_int), margin_twtt=margin_twtt)

        # delta
        dtwtt_1 = peak_twtt_base_1 - peak_twtt_surf_1
        dtwtt_2 = peak_twtt_base_2 - peak_twtt_surf_2
        dtpower_1 = peak_power_base_1 - peak_power_surf_1
        dtpower_2 = peak_power_base_2 - peak_power_surf_2

        if idx % 10 == 0:
            print("Processed {} intersections".format(idx))
            # print("TWTT delta 1:", dtwtt_1)
            # print("TWTT delta 2:", dtwtt_2)
            # print("Power delta 1:", dtpower_1)
            # print("Power delta 2:", dtpower_2)

        # recombine into a geopandas dataframe
        intersections_radio.loc[idx, "dtwtt_1"] = dtwtt_1
        intersections_radio.loc[idx, "dtwtt_2"] = dtwtt_2
        intersections_radio.loc[idx, "dtpower_1"] = dtpower_1
        intersections_radio.loc[idx, "dtpower_2"] = dtpower_2

    return intersections_radio
        
