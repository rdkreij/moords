"Class to generate overview documents and check serial numbers"

import pandas as pd

import moords.overview_to_pdf as overview_to_pdf
import moords.overview_tools as overview_tools
from moords.mooring_design import Mooring


class Overview:
    """Class to generate overview documents and check serial number"""

    def __init__(self, moorings: Mooring | list[Mooring]):
        if not isinstance(moorings, list):
            moorings = [moorings]
        self.moorings = moorings
        self.df_combined_moorings = overview_tools.make_df_combined_moorings(
            self.moorings
        )

    def generate_overview_pdf(
        self,
        file_path: str = "",
        header: str = None,
        replacements: list[tuple[str]] = None,
    ) -> None:
        """Generate mooring design summary in pdf using latex."""
        overview_to_pdf.generate_overview_pdf(
            df=self.df_combined_moorings,
            file_path=file_path,
            header=header,
            replacements=replacements,
        )

    def print_duplicated_serial(self) -> None:
        """Print a summary of the duplicated serial numbers"""
        overview_tools.print_duplicated_serial(self.df_combined_moorings)

    @property
    def df_missing_serial(self) -> pd.DataFrame:
        """Make dataframe of all elements with a missing serial number"""
        return overview_tools.make_df_missing_serial(self.df_combined_moorings)
