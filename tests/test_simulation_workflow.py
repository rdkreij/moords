import unittest

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import moords


def construct_mooring():
    path_mooring_database = "data/mooring_elements/mooring_elements_data.csv"
    df_database = moords.load_database_from_csv(path_mooring_database)

    mooring = moords.Mooring(df_database, "test")
    mooring.add_inline("1 Railway Wheel")
    mooring.add_inline("16mmBS")
    mooring.add_inline("3/8 wire rope", line_length=50)
    mooring.add_inline("16mmBS")
    mooring.add_inline("30in float")
    mooring.add_clampon_by_height("FLNTUSB (UWA)", serial="2927", height=20)
    mooring.add_clampon_by_height("SBE39 T (UWA)", serial="3921", height=30)
    mooring.add_clampon_by_height("SBE39 T (UWA)", serial="4263", height=35)
    return mooring


class TestSimulationInstanceWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Initialize test moorings and set up temporary directories.
        """
        cls.mooring = construct_mooring()

        bottom_depth = 245
        z = np.linspace(bottom_depth, 0, 100)
        cls.ds_flow_instance = xr.Dataset(
            {
                "U": (["z"], z / bottom_depth),
                "V": (["z"], np.zeros_like(z)),
                "W": (["z"], np.zeros_like(z)),
                "rho": (["z"], 1025 * np.ones_like(z)),
                "bottom_depth": bottom_depth,
            },
            coords={
                "z": z,
            },
        )

        cls.instance = moords.SimulateInstance(cls.mooring, cls.ds_flow_instance)
        cls.instance.run()

    def test_simulation_instance(self):
        try:
            fig, ax = self.instance.plot_instance()
            plt.close(fig)
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")


class TestSimulationSeriesWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Initialize test moorings and set up temporary directories.
        """
        cls.mooring = construct_mooring()
        ds_flow_series = xr.open_dataset("data/flow_series/flow_series_S245.nc")
        ds_flow_series = ds_flow_series.isel(time=slice(0, 4, 1))
        cls.series = moords.SimulateSeries(cls.mooring, ds_flow_series)
        cls.series.run()

    def test_series_flow_plot(self):
        try:
            fig, ax = self.series.plot_flow_series()
            plt.close(fig)
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_simulation_series_anchor_plot(self):
        try:
            fig, ax = self.series.plot_anchor_stats_series()
            plt.close(fig)
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_simulation_series_instance_plot(self):
        try:
            fig, ax = self.series.plot_instance(time_idx=2)
            plt.close(fig)
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_simulation_series_item_plot(self):
        try:
            fig, ax = self.series.plot_simulated_item_series(item_id=2)
            plt.close(fig)
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_simulation_series_summary_ds(self):
        try:
            self.series.df_summary_series
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")


if __name__ == "__main__":
    unittest.main()
