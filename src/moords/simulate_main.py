"""Classes to simulate mooring in flow"""

import numpy as np
import xarray as xr

import moords.simulation_data as simulate_data
import moords.simulation_dynamics as simulate_dynamics
import moords.simulation_plot as simulation_plot
from moords.mooring_design import Mooring


class SimulateInstance:
    """ "Class to simulate mooring in single instance / timestep of flow"""

    def __init__(self, mooring: Mooring, ds_flow_instance: xr.Dataset):
        self.mooring = mooring
        self.ds_flow_instance = ds_flow_instance
        self.ds_sim_instance_interpolated = None
        self.ds_sim_instance = None

    def run(self, print_anchor_summary: bool = True):
        """Run mooring simulation"""
        self.ds_sim_instance_interpolated = simulate_dynamics.simulate_mooring_instance(
            self.mooring,
            self.ds_flow_instance,
            print_anchor_summary=print_anchor_summary,
        )
        self.ds_sim_instance = simulate_data.convert_ds_to_regular_mooring(
            self.ds_sim_instance_interpolated
        )
        self.ds_sim_instance = simulate_data.add_details_to_ds_sim(
            self.ds_sim_instance, self.mooring
        )

    def plot_instance(self):
        """Plot mooring in flow instance"""
        return simulation_plot.plot_simulated_instance(
            self.ds_sim_instance_interpolated, self.ds_flow_instance
        )


class SimulateSeries:
    """ "Class to simulate mooring in time series of flow"""

    def __init__(self, mooring: Mooring, ds_flow_series: xr.Dataset):
        self.mooring = mooring
        self.ds_flow_series = ds_flow_series
        self.ds_sim_series_interpolated = None
        self.ds_sim_series = None
        self.df_summary_series = None

    def run(self):
        """Run mooring simulation"""
        self.ds_sim_series_interpolated = simulate_dynamics.simulate_mooring_series(
            self.mooring, self.ds_flow_series
        )
        self.ds_sim_series = simulate_data.convert_ds_to_regular_mooring(
            self.ds_sim_series_interpolated
        )
        self.ds_sim_series = simulate_data.add_details_to_ds_sim(
            self.ds_sim_series, self.mooring
        )
        self.df_summary_series = simulate_data.make_df_summary_series(
            self.ds_sim_series
        )

    def plot_flow_series(self):
        """Plot flow profile as a function of time"""
        return simulation_plot.plot_flow_series(self.ds_flow_series)

    def plot_instance(self, time_idx: int):
        """Plot mooring in flow instance"""
        ds_sim_instance = self.ds_sim_series_interpolated.isel(time=time_idx)
        ds_flow_instance = self.ds_flow_series.isel(time=time_idx)

        fig, ax = simulation_plot.plot_simulated_instance(
            ds_sim_instance, ds_flow_instance
        )
        timestamp = ds_sim_instance.time.values
        fig.suptitle(np.datetime_as_string(timestamp, unit="s"))
        return fig, ax

    def plot_anchor_stats_series(self):
        """Plot anchor statistics as a function of time"""
        return simulation_plot.plot_anchor_stats_series(self.ds_sim_series_interpolated)

    def plot_simulated_item_series(self, item_id: int):
        """Plot statistics of an item (clampon or inline) as a function of time"""
        ds_sim_series_item = self.ds_sim_series.sel(id=item_id)
        return simulation_plot.plot_simulated_item_series(ds_sim_series_item)
