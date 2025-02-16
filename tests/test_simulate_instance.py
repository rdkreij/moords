import unittest

import numpy as np
import xarray as xr

import moords


class TestMooringSimulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Set up common test data (e.g., mooring database) for all tests.
        """
        cls.path_mooring_database = "data/mooring_elements/mooring_elements_data.csv"
        cls.df_database = moords.load_database_from_csv(cls.path_mooring_database)

    def create_simple_mooring(self):
        """
        Helper function to create simple mooring.
        """
        mooring = moords.Mooring(self.df_database, "test")
        mooring.add_inline("1 Railway Wheel")
        mooring.add_inline("16mmBS")
        mooring.add_inline("3/8 wire rope", line_length=30)
        mooring.add_inline("16mmBS")
        mooring.add_inline("30in float")
        return mooring

    def create_flow_instance(
        self, start_height=0, stop_height=100, n_z=100, theta_deg=0, speed=0
    ):
        """
        Helper function to create a flow instance with the specified parameters.
        """
        theta_rad = np.deg2rad(theta_deg)
        z = np.linspace(start_height, stop_height, n_z)
        return xr.Dataset(
            {
                "U": (["z"], speed * np.ones(n_z) * np.cos(theta_rad)),
                "V": (["z"], speed * np.ones(n_z) * np.sin(theta_rad)),
                "W": (["z"], np.zeros(n_z)),
                "rho": (["z"], 1025 * np.ones(n_z)),
                "bottom_depth": stop_height,
            },
            coords={"z": z},
        )

    def test_mooring_simulation_anchor_too_light_warning(self):
        """
        Tries to simulate the mooring while the anchor is too light.
        """
        ds_flow_instance = self.create_flow_instance(theta_deg=20, speed=1)
        mooring = self.create_simple_mooring()
        mooring.add_inline("61in ORE")
        instance = moords.SimulateInstance(mooring, ds_flow_instance)
        with self.assertRaises(Exception) as context:
            instance.run(print_anchor_summary=False)
        print(f"Test passed: Caught expected error - {context.exception}")

    def test_clampon_effect_top(self):
        """
        Test that attaching the same top float as clampon doubles the tension.
        """
        ds_flow_instance = self.create_flow_instance(theta_deg=20, speed=1)

        list_tension = []
        for add_float in [False, True]:
            mooring = self.create_simple_mooring()

            if add_float:
                mooring.add_clampon_along_inline("30in float", 0.2)

            instance = moords.SimulateInstance(mooring, ds_flow_instance)
            instance.run(print_anchor_summary=False)
            df_sim_instance = instance.ds_sim_instance.to_dataframe()
            list_tension.append(
                df_sim_instance[~df_sim_instance.bool_clampon].iloc[-1].tension_bottom
            )

        self.assertTrue(
            np.isclose(list_tension[0] * 2, list_tension[1]),
            (
                f"Test failed: Expected tension on inline {list_tension[0] * 2} "
                f"!= simulated tension with clampon {list_tension[1]}"
            ),
        )
        print("Test passed: Adding clampon gives expected tension.")

    def test_mooring_flow_direction(self):
        """
        Simulates the mooring system and validates the tension under the second float.
        """
        theta_deg = 20
        ds_flow_instance = self.create_flow_instance(theta_deg=theta_deg, speed=1)
        mooring = self.create_simple_mooring()
        instance = moords.SimulateInstance(mooring, ds_flow_instance)
        instance.run(print_anchor_summary=False)
        df_sim_instance = instance.ds_sim_instance.to_dataframe()

        # Get horizontal rotation of top float
        inline_dict = df_sim_instance.iloc[-1]
        X = inline_dict["X"]
        Y = inline_dict["Y"]
        simulated_theta_deg = np.rad2deg(np.arctan(Y / X))

        self.assertTrue(
            np.isclose(simulated_theta_deg, theta_deg),
            (
                f"Test failed: Flow angle {theta_deg} "
                f"!= mooring angle {simulated_theta_deg}"
            ),
        )
        print("Test passed: Simulated mooring matches flow angle.")

    def test_mooring_simulation_incomplete_flow_profile_top_warning(self):
        """
        Tries to simulate the mooring while the flow profile is incomplete
        at the top, which should raise an error.
        """
        ds_flow_instance = self.create_flow_instance(stop_height=20)
        mooring = self.create_simple_mooring()
        instance = moords.SimulateInstance(mooring, ds_flow_instance)
        with self.assertRaises(Exception) as context:
            instance.run(print_anchor_summary=False)

        print(f"Test passed: Caught expected error - {context.exception}")

    def test_mooring_simulation_incomplete_flow_profile_bottom_warning(self):
        """
        Tries to simulate the mooring while the flow profile is incomplete
        at the bottom, which should raise an error.
        """
        ds_flow_instance = self.create_flow_instance(start_height=100)
        mooring = self.create_simple_mooring()
        instance = moords.SimulateInstance(mooring, ds_flow_instance)
        with self.assertRaises(Exception) as context:
            instance.run(print_anchor_summary=False)

        print(f"Test passed: Caught expected error - {context.exception}")

    def test_mooring_line_tension_in_flow(self):
        """
        Simulates the mooring system and validates the tension under the second float.
        """
        theta_deg = 20
        speed = 2
        ds_flow_instance = self.create_flow_instance(theta_deg=theta_deg, speed=speed)
        mooring = self.create_simple_mooring()
        instance = moords.SimulateInstance(mooring, ds_flow_instance)
        instance.run(print_anchor_summary=False)
        df_sim_instance = instance.ds_sim_instance.to_dataframe()

        simulated_tension = df_sim_instance.iloc[-1].tension_bottom
        name_last_inline = mooring.inline[-1].name
        inline_dict = self.df_database.loc[name_last_inline]
        buoyancy = inline_dict["buoyancy_kg"] * 9.81
        drag_coef = inline_dict["drag"]
        diameter = inline_dict["diameter_m"]
        rho = 1025
        area = np.pi * (diameter / 2) ** 2
        Qhorizontal = 0.5 * rho * drag_coef * area * speed**2
        expected_tension = np.sqrt(Qhorizontal**2 + buoyancy**2) / 9.81

        # Test result
        self.assertAlmostEqual(
            expected_tension,
            simulated_tension,
            places=5,
            msg=f"Test failed: Simulated tension {simulated_tension} "
            f"!= Expected {expected_tension}",
        )

        print("Test passed: Tension matches expectation.")

    def test_mooring_line_tension_in_rest(self):
        """
        Simulates the mooring system and validates the tension under the second float.
        """
        ds_flow_instance = self.create_flow_instance(speed=0)
        mooring = self.create_simple_mooring()
        mooring.add_inline("30in float")
        instance = moords.SimulateInstance(mooring, ds_flow_instance)
        instance.run(print_anchor_summary=False)
        df_sim_instance = instance.ds_sim_instance.to_dataframe()

        # Get tension under the 2nd float from the top
        dict_element = df_sim_instance.iloc[-2]
        tension_bottom = dict_element.tension_bottom

        # Calculate tension based on buoyancy
        buoyancy_element = self.df_database.loc["30in float"].buoyancy_kg
        expected_tension = 2 * buoyancy_element

        self.assertAlmostEqual(
            tension_bottom,
            expected_tension,
            places=5,
            msg=f"Test failed: Tension bottom {tension_bottom} "
            f"!= Expected {expected_tension}",
        )

        print("Test passed: Tension matches expected buoyancy.")

    def test_mooring_full_line_stretching(self):
        """
        Simulates the mooring system and validates the line stretching.
        """
        ds_flow_instance = self.create_flow_instance(speed=0)
        mooring = self.create_simple_mooring()
        mooring.inline[2].buoyancy = 0
        mooring.inline[2].width = 0.01
        mooring.inline[2].material = 5  # Polyethy
        instance = moords.SimulateInstance(mooring, ds_flow_instance)
        instance.run(print_anchor_summary=False)
        z_design_top = mooring.inline[4].bottom_height + mooring.inline[4].length / 2
        z_sim_top = instance.ds_sim_instance.Z.values[4]
        line_stretch_sim = z_sim_top - z_design_top

        # Material property
        mod_elast = 6.9 * 10**8  # Modulus of elasticity for Polyethy
        length = mooring.inline[2].length
        width = mooring.inline[2].width
        buoyancy_above_line_kg = mooring.inline[3].buoyancy + mooring.inline[4].buoyancy
        tension = buoyancy_above_line_kg * 9.81
        line_stretch_expected = length * (tension * 4 / (np.pi * width**2 * mod_elast))

        # Test result
        self.assertAlmostEqual(line_stretch_expected, line_stretch_sim, places=4)

        print("Test passed: Line stretch matches expectation.")

    def test_clampon_effect_bottom(self):
        """
        Test that attaching the same anchor as clampon increases the weight under the
        anchor with the right amount.
        """
        ds_flow_instance = self.create_flow_instance(theta_deg=20, speed=1)
        list_woa = []
        for add_clampon_anchor in [False, True]:
            mooring = self.create_simple_mooring()
            if add_clampon_anchor:
                mooring.add_clampon_by_height("1 Railway Wheel", 0.1)
            instance = moords.SimulateInstance(mooring, ds_flow_instance)
            instance.run(print_anchor_summary=False)
            list_woa.append(instance.ds_sim_instance.weight_under_anchor.item())

        woa_no_added_anchor = list_woa[0]
        woa_added_anchor = list_woa[1]
        buoyancy_element = self.df_database.loc["1 Railway Wheel"].buoyancy_kg
        expected_woa_added_anchor = woa_no_added_anchor + buoyancy_element

        self.assertTrue(
            np.isclose(woa_added_anchor, expected_woa_added_anchor),
            (
                f"Test failed: Expected woa {expected_woa_added_anchor} "
                f"!= simulated woa {woa_added_anchor}"
            ),
        )
        print("Test passed: Adding clampon gives expected tension.")


if __name__ == "__main__":
    unittest.main()
