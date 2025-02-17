import os
import shutil
import unittest

import matplotlib.pyplot as plt

import moords


class TestMooringWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Initialize test moorings and set up temporary directories.
        """
        cls.path_mooring_database = "data/mooring_elements/mooring_elements_data.csv"
        cls.df_database = moords.load_database_from_csv(cls.path_mooring_database)

        cls.temp_dir = "tests/temp"
        os.makedirs(cls.temp_dir, exist_ok=True)

        mooring_1 = moords.Mooring(cls.df_database, "test")
        mooring_1.add_inline("1 Railway Wheel")
        mooring_1.add_inline("16mmBS")
        mooring_1.add_inline("3/8 wire rope", line_length=50)
        mooring_1.add_inline("16mmBS")
        mooring_1.add_inline("30in float")
        mooring_1.add_clampon_by_height("FLNTUSB (UWA)", serial="2927", height=20)
        mooring_1.add_clampon_by_height("SBE39 T (UWA)", serial="3921", height=30)
        mooring_1.add_clampon_by_height("SBE39 T (UWA)", serial="4263", height=35)
        cls.mooring_1 = mooring_1

        mooring_2 = moords.Mooring(cls.df_database, "test")
        mooring_2.add_inline("1 Railway Wheel")
        mooring_2.add_inline("16mmBS")
        mooring_2.add_inline("3/8 wire rope", line_length=70)
        mooring_2.add_inline("16mmBS")
        mooring_2.add_inline("30in float")
        mooring_1.add_clampon_by_height("SBE39 T (UWA)", serial="2672", height=30)
        mooring_1.add_clampon_by_height("SBE39 T (UWA)", serial="4263", height=35)
        cls.mooring_2 = mooring_2

    @classmethod
    def tearDownClass(cls):
        """
        Remove temporary test directories after tests are complete.
        """
        if os.path.exists(cls.temp_dir):
            shutil.rmtree(cls.temp_dir)

    def test_create_table_mooring(self):
        """
        Test generating mooring component tables without errors.
        """
        try:
            self.mooring_1.make_df_clampon()
            self.mooring_1.make_df_inline()
            self.mooring_1.make_df_combined()
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_create_plot_mooring(self):
        """
        Test generating a mooring design plot without errors.
        """
        try:
            fig, ax = self.mooring_1.plot_design(
                label_rigging=True, drop_strings=["(UWA)"], line_ratio_plot=0.8
            )
            plt.close(fig)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_create_combined_table_mooring(self):
        """
        Test creating a combined table from multiple moorings.
        """
        try:
            overview = moords.Overview([self.mooring_1, self.mooring_2])
            overview.df_combined_moorings
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_create_print_combined_table_mooring(self):
        """
        Test generating a PDF of the combined mooring tables.
        """
        try:
            overview = moords.Overview([self.mooring_1, self.mooring_2])
            overview.generate_overview_pdf(
                file_path=f"{self.temp_dir}/example_tables.pdf",
                header=r"V1, \today",
                replacements=[("_", r"\_"), ("(UWA)", ""), ("(passive)", "")],
            )
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_create_print_dubplicated_serial(self):
        """
        Test checking for duplicated serial numbers in mooring components.
        """
        try:
            overview = moords.Overview([self.mooring_1, self.mooring_2])
            overview.print_duplicated_serial()
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_create_table_missing_serial(self):
        """
        Test generating a table of mooring components with missing serial numbers.
        """
        try:
            overview = moords.Overview([self.mooring_1, self.mooring_2])
            overview.df_missing_serial
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")

    def test_create_export_csv_mooring(self):
        """
        Test exporting mooring data to a CSV file.
        """
        try:
            self.mooring_1.export_csv(f"{self.temp_dir}/test_export.csv")
            self.assertTrue(True)
        except Exception as e:
            self.fail(f"Code raised an exception: {e}")


if __name__ == "__main__":
    unittest.main()
