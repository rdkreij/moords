"""Module to organise mooring simulation data"""

import numpy as np
import pandas as pd
import xarray as xr

from moords.mooring_design import Mooring


def convert_ds_to_regular_mooring(ds_interpolated: xr.Dataset):
    """
    Convert simulation results dataset from interpolated to original mooring.
    This conversion works on both a dataset from an instance and a time series.
    """

    weight_under_anchor = ds_interpolated.weight_under_anchor

    # Set tension above inline
    ds_interpolated = ds_interpolated.rename_vars({"tension": "tension_top"})

    # Set tension below inline
    idx_inline = ~ds_interpolated.bool_clampon
    tension_bottom = ds_interpolated.tension_top.isel(idx=idx_inline).shift(
        idx=1, fill_value=0
    )
    ds_interpolated["tension_bottom"] = tension_bottom

    # Find id of all clamps and inlines
    id_all = ds_interpolated.id.values
    id_unique = np.unique(id_all)

    ds_list = []  # Allocate list of all elements
    for id_i in id_unique:
        idx = np.where(id_all == id_i)[0]
        ds_segment = ds_interpolated.isel(idx=idx)

        if len(idx) > 1:  # Element consists out of segments in ds_interpolated
            # Combine values of segments
            bool_clampon = ds_segment["bool_clampon"].isel(idx=0)
            id = ds_segment["id"].isel(idx=0)
            X = ds_segment["X"].mean(dim="idx")
            Y = ds_segment["Y"].mean(dim="idx")
            Z = ds_segment["Z"].mean(dim="idx")
            psi = ds_segment["psi"].mean(dim="idx")
            dZ = ds_segment["dZ"].mean(dim="idx")
            tension_bottom = (
                ds_segment["tension_bottom"]
                .isel(idx=ds_segment["id_sub"].argmin(dim="idx"))
                .drop_vars("idx")
            )
            tension_top = (
                ds_segment["tension_top"]
                .isel(idx=ds_segment["id_sub"].argmax(dim="idx"))
                .drop_vars("idx")
            )
            ds_segment = xr.merge(
                [bool_clampon, id, X, Y, Z, psi, dZ, tension_bottom, tension_top]
            )
            ds_segment = ds_segment.expand_dims(idx=(ds_segment.id.size,))
            ds_list.append(ds_segment)
        else:  # Element is single item in ds_interpolated
            ds_segment = ds_segment.drop_vars("id_sub")
            ds_segment = ds_segment.drop_vars("weight_under_anchor")
            ds_list.append(ds_segment)
    # Combine all elements into single dataset
    ds_regular = xr.concat(ds_list, dim="idx")
    ds_regular = ds_regular.swap_dims({"idx": "id"}).drop_vars("idx")
    ds_regular["weight_under_anchor"] = weight_under_anchor
    return ds_regular


def add_details_to_ds_sim(ds: xr.Dataset, mooring: Mooring):
    """
    Add mooring basic element details to existing mooring simulation dataset.
    The function can be applied to an interpolated or regular mooring and
    to a single instance or a time series.
    """
    df = mooring.make_df_combined()

    idx = [np.where(df["id"] == id)[0][0] for id in ds.id.values]
    df_c = df.loc[idx]
    ds_updated = ds.assign_coords(
        {
            "type": ("id", df_c["type"].values),
            "name": ("id", df_c["name"].values),
            "serial": ("id", df_c["serial"].values),
            "bottom_height": ("id", df_c["bottom_height"].values),
            "height": ("id", df_c["height"].values),
        }
    )
    return ds_updated


def make_df_summary_series(ds_sim_series: xr.Dataset):
    """
    Make a dataframe with summary statistic of each element for a time series.
    Includes a summary of psi (vertical tilt) and dZ (vertical pulldown).
    """
    item_id = ds_sim_series.id.values
    item_type = ds_sim_series.type.values
    name = ds_sim_series.name.values
    serial = ds_sim_series.serial.values

    psi_mean = ds_sim_series["psi"].mean(dim="time").values
    psi_std = ds_sim_series["psi"].std(dim="time").values
    psi_01 = ds_sim_series["psi"].quantile(0.01, dim="time").values
    psi_99 = ds_sim_series["psi"].quantile(0.99, dim="time").values
    psi_min = ds_sim_series["psi"].min(dim="time").values
    psi_max = ds_sim_series["psi"].max(dim="time").values

    dZ_mean = ds_sim_series["dZ"].mean(dim="time")
    dZ_std = ds_sim_series["dZ"].std(dim="time")
    dZ_01 = ds_sim_series["dZ"].quantile(0.01, dim="time")
    dZ_99 = ds_sim_series["dZ"].quantile(0.99, dim="time")
    dZ_min = ds_sim_series["dZ"].min(dim="time")
    dZ_max = ds_sim_series["dZ"].max(dim="time")

    df_sim_series_summary = pd.DataFrame(
        {
            "id": item_id,
            "type": item_type,
            "name": name,
            "serial": serial,
            "psi_mean": psi_mean,
            "psi_std": psi_std,
            "psi_01": psi_01,
            "psi_99": psi_99,
            "psi_min": psi_min,
            "psi_max": psi_max,
            "dZ_mean": dZ_mean,
            "dZ_std": dZ_std,
            "dZ_01": dZ_01,
            "dZ_99": dZ_99,
            "dZ_min": dZ_min,
            "dZ_max": dZ_max,
        }
    )
    df_sim_series_summary = df_sim_series_summary.set_index("id")
    return df_sim_series_summary
