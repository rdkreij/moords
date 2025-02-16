"""Classes to design a mooring"""

import os

import numpy as np
import pandas as pd

from moords.plot_design import plot_design


def load_database_from_csv(path_csv: str) -> pd.DataFrame:
    """load mooring database from csv file"""
    dtype_dict = {
        "type": "str",
        "name": "str",
        "buoyancy_kg": "float64",
        "length_m": "float64",
        "width_m": "float64",
        "diameter_m": "float64",
        "drag": "float64",
        "material": "float64",
        "comment": "str",
    }
    df_database = pd.read_csv(path_csv, dtype=dtype_dict, na_values=None)
    df_database["material"] = df_database["material"].astype("Int64")
    df_database["comment"] = df_database["comment"].replace(float("NaN"), "")
    df_database.set_index("name", inplace=True)
    return df_database


def get_property_database(
    df: pd.DataFrame,
    name: str,
    bool_line: bool,
    passive_weight: bool = False,
    passive_drag: bool = False,
) -> dict:
    """Request properties given the name of an element."""
    property_dict = df.loc[name].to_dict()
    if passive_weight:
        property_dict["buoyancy_kg"] = 0
    if passive_drag:
        property_dict["diameter_m"] = 0
        property_dict["width_m"] = 0
        property_dict["drag"] = 0

    # Raise error missing properties
    check_property = [
        "buoyancy_kg",
        "length_m",
        "width_m",
        "diameter_m",
        "drag",
        "material",
    ]
    for prop in check_property:
        value = property_dict[prop]
        try:
            bool_material = prop == "material"
            bool_missing = pd.isna(value)
            if bool_missing and (not bool_line and not bool_material or bool_line):
                set_value = int(bool_material)
                raise ValueError(
                    f"Element {name} has missing property "
                    f"{prop} in database; {prop} set to {set_value}."
                )
        except ValueError as e:
            print(f"Error: {e}")
            property_dict[prop] = set_value
    return property_dict


class ClampOn:
    """Class of a clamp-on object."""

    attributes = [
        "type",
        "name",
        "serial",
        "section",
        "height",
        "height_along_inline",
        "clamped_to_name",
        "clamped_to_serial",
        "length",
        "bool_fit",
        "buoyancy",
        "width",
        "diameter",
        "drag",
        "material",
        "comment",
        "bool_clampon",
        "bool_line",
        "inline_bottom",
        "inline_top",
        "num_inline",
        "id",
        "passive_weight",
        "passive_drag",
    ]

    def __init__(
        self,
        df_database: pd.DataFrame,
        name: str,
        serial: str = None,
        passive_weight=False,
        passive_drag=False,
    ) -> None:
        self.name = name
        self.serial = serial
        self.passive_weight = passive_weight
        self.passive_drag = passive_drag

        # Retrieve properties from the mooring database
        prop = get_property_database(
            df_database, name, False, self.passive_weight, self.passive_drag
        )

        # Set attributes based on the properties from the database
        self.type = prop.get("type")
        self.length = prop.get("length_m")
        self.width = prop.get("width_m")
        self.material = prop.get("material")
        self.comment = prop.get("comment")
        self.buoyancy = prop.get("buoyancy_kg")
        self.diameter = prop.get("diameter_m")
        self.drag = prop.get("drag")

        # Default attributes for clamp-on configuration
        self.height = None
        self.height_along_inline = None
        self.clamped_to_name = None
        self.clamped_to_serial = None
        self.bool_fit = None
        self.section = None
        self.bool_clampon = True
        self.bool_line = False
        self.inline_bottom = None
        self.inline_top = None
        self.num_inline = None
        self.id = None

    def get_dict(self) -> dict:
        """Turn object into a dict."""

        properties = {}
        for attr in self.attributes:
            properties[attr] = getattr(self, attr)

        return properties

    def __repr__(self):
        return (
            f"ClampOn({self.type},"
            f"{self.name},"
            f"key={self.serial},"
            f"l={self.length},"
            f"h={self.height},"
            f"clamped_to={self.clamped_to_name})"
        )


class InLine:
    """Class of an in-line object."""

    attributes = [
        "type",
        "name",
        "serial",
        "section",
        "bottom_height",
        "length",
        "bool_line",
        "buoyancy",
        "width",
        "diameter",
        "drag",
        "material",
        "comment",
        "top_height",
        "bool_clampon",
        "total_buoyancy",
        "num_inline",
        "id",
    ]

    def __init__(
        self,
        df_database: pd.DataFrame,
        name: str,
        line_length: float = None,
        serial: str = None,
        section: str = None,
    ) -> None:
        self.name = name
        self.serial = serial
        self.section = section
        self.bool_line = line_length is not None

        # Retrieve properties from the mooring database
        prop = get_property_database(df_database, name, self.bool_line)

        # Set length and buoyancy based on whether line_length is provided
        self.length = line_length if line_length is not None else prop.get("length_m")
        self.buoyancy = prop.get("buoyancy_kg") * (self.length if line_length else 1)

        # Set other attributes from the database properties
        self.type = prop.get("type")
        self.width = prop.get("width_m")
        self.diameter = prop.get("diameter_m")
        self.drag = prop.get("drag")
        self.material = prop.get("material")
        self.comment = prop.get("comment")

        # Default attributes for inline configuration
        self.bool_clampon = False
        self.bottom_height = None
        self.num_inline = None
        self.total_buoyancy = self.buoyancy
        self.clamp_ons = []
        self.id = None

    @property
    def top_height(self):
        """Inline top height"""
        return self.bottom_height + self.length

    @property
    def center_height(self):
        """Ïnline center height"""
        return self.bottom_height + 0.5 * self.length

    def get_dict(self) -> dict:
        """Turn object into a dict."""
        properties = {}
        for attr in self.attributes:
            properties[attr] = getattr(self, attr)

        return properties

    # @property
    # def interpolate(self) -> List["InLine"]:
    #     """Interpolate into segments"""
    #     inlines = []
    #     if self.bool_line:
    #         length = self.length
    #         height = self.bottom_height
    #         clamp_ons = self.clamp_ons

    #         # segment line
    #         if length < 0.2:
    #             segment_target_length = length
    #         if (length > 0.2) and (length <= 5):
    #             segment_target_length = 0.2
    #         elif (length > 5) and (length <= 50):
    #             segment_target_length = 0.5
    #         elif (length > 50) & (length <= 100):
    #             segment_target_length = 1
    #         elif (length > 100) & (length <= 500):
    #             segment_target_length = 2
    #         else:
    #             segment_target_length = 5
    #         num_segment = int(np.round(length / segment_target_length))
    #         segment_fraction = 1 / num_segment
    #         segment_length = segment_fraction * length
    #         segment_buoyancy = segment_fraction * self.buoyancy

    #         # collect segments
    #         inlines = []
    #         for i in range(num_segment):
    #             segment_inline = copy.deepcopy(self)
    #             segment_inline.clamp_ons = []

    #             segment_inline.bottom_height = height + i * segment_length
    #             segment_inline.length = segment_length
    #             segment_inline.buoyancy = segment_buoyancy
    #             segment_inline.total_buoyancy = segment_buoyancy

    #             # re-attach clampons
    #             for clamp_on in clamp_ons:
    #                 height_clamp_on = clamp_on.height
    #                 if (height_clamp_on >= segment_inline.bottom_height) and (
    #                     height_clamp_on < segment_inline.top_height
    #                 ):
    #                     segment_inline.clamp_ons.append(clamp_on)
    #                     segment_inline.total_buoyancy += clamp_on.buoyancy
    #             inlines.append(segment_inline)
    #     else:
    #         inlines.append(self)

    #     return inlines

    def __repr__(self):
        return (
            f"InLine({self.type},"
            f"{self.name},"
            f"key={self.serial},"
            f"l={self.length},"
            f"bottom_h={self.bottom_height})"
        )


class Mooring:
    """Class of mooring object."""

    def __init__(self, df_database: pd.DataFrame, name: str = ""):
        self.df_database = df_database
        self.name = name

        # Default attributes for mooring configuration
        self.inline = []
        self.id_track = -1  # id assigning to inline and clampons

    @property
    def num_inline(self):
        """Total number of inline elements"""
        return len(self.inline)

    def __repr__(self):
        string = []
        for elem in self.inline:
            string.append(f"{elem}")
            for clamp in elem.clamp_ons:
                string.append(f"{clamp}")
        return "\n".join(string)

    def add_inline(
        self,
        name: str,
        line_length: float = None,
        serial: str = None,
        section: str = None,
    ) -> None:
        """Add an inline element to mooring."""
        new_inline = InLine(self.df_database, name, line_length, serial, section)

        if self.num_inline == 0:
            new_inline.bottom_height = 0
            new_inline.num_inline = 1
        else:
            previous_inline = self.inline[-1]
            new_inline.bottom_height = previous_inline.top_height
            new_inline.num_inline = previous_inline.num_inline + 1

        self.id_track = self.id_track + 1
        new_inline.id = self.id_track

        self.inline.append(new_inline)

    def add_clampon_along_inline(
        self,
        name: str,
        height_along_inline: float,
        serial: str = None,
        passive_weight: bool = False,
        passive_drag: bool = False,
    ) -> None:
        """Add clampon along in-line element."""
        new_clampon = ClampOn(
            self.df_database,
            name,
            serial,
            passive_weight,
            passive_drag,
        )
        new_clampon.height_along_inline = height_along_inline
        target_inline = self.inline[-1]
        self.attach_clampon(target_inline, new_clampon)

    def check_total_buoyancy(self):
        """Check total buoyancy of mooring and raise warning when needed"""
        total_buoyancy = np.array([elem.total_buoyancy for elem in self.inline])
        cumsum_buoyancy = np.cumsum(total_buoyancy[::-1])[::-1]
        idx_negative = cumsum_buoyancy < 0
        any_negative = np.any(idx_negative[1:])
        positive_anchor = cumsum_buoyancy[0] > 0
        flag = True

        try:
            if any_negative:
                idx_list = np.where(idx_negative[1:])[0] + 1
                for idx in idx_list:
                    elem = self.inline[idx]
                    raise ValueError(
                        f"Not enough buoyancy, {elem.name} "
                        f"starting at {elem.bottom_height} lacks "
                        f"tension {cumsum_buoyancy[idx]:.2f} kg"
                    )
        except ValueError as e:
            print(f"Error: {e}")
            flag = False

        try:
            if positive_anchor:
                raise ValueError(
                    f"Anchor too light,"
                    f"netto mooring buoyancy "
                    f"= +{cumsum_buoyancy[0]:.2f}."
                )
        except ValueError as e:
            print(f"Error: {e}")
            flag = False

        if flag:
            print("Buoyancy check succesful")

    def add_clampon_by_height(
        self,
        name: str,
        height: float,
        serial: str = None,
        passive_weight: bool = False,
        passive_drag: bool = False,
    ) -> None:
        """Add clampon by height."""
        new_clampon = ClampOn(
            self.df_database,
            name,
            serial,
            passive_weight,
            passive_drag,
        )
        new_clampon.height = height

        top_heights = np.array([elem.top_height for elem in self.inline])
        bottom_heights = np.array([elem.bottom_height for elem in self.inline])

        if height >= top_heights[-1]:
            raise ValueError(
                f"Clamp-on {name} height {height:.2f} is "
                f"set above mooring height {top_heights[-1]:.2f}."
            )

        idx = (height >= bottom_heights) & (height < top_heights)
        idx = int(np.where(idx)[0])
        target_inline = self.inline[idx]
        self.attach_clampon(target_inline, new_clampon)

    def attach_clampon(self, elem: InLine, clamp: ClampOn) -> None:
        """Attach clamp-on to in-line object."""
        # Attach to inline
        elem.clamp_ons.append(clamp)

        # Add to buoyancy
        elem.total_buoyancy += clamp.buoyancy

        # Fill in clamp targets
        clamp.clamped_to_name = elem.name
        clamp.clamped_to_serial = elem.serial

        # Fill in heights
        if clamp.height is None:
            clamp.height = elem.bottom_height + clamp.height_along_inline
        if clamp.height_along_inline is None:
            clamp.height_along_inline = clamp.height - elem.bottom_height

        # Fill in elem details
        clamp.inline_bottom = elem.bottom_height
        clamp.inline_top = elem.top_height
        clamp.num_inline = elem.num_inline

        self.id_track = self.id_track + 1
        clamp.id = self.id_track

        # Add section to clamp
        if elem.section is not None:
            clamp.section = elem.section

        # Check if element fits
        bool_fit_height = (clamp.height >= elem.bottom_height) & (
            clamp.height <= elem.top_height
        )
        bool_fit_below = clamp.height - clamp.length >= elem.bottom_height
        bool_fit_above = clamp.height + clamp.length <= elem.top_height
        bool_fit = bool_fit_height & (bool_fit_below | bool_fit_above)
        try:
            if not bool_fit:
                raise ValueError(
                    (
                        f"Clamp-on {clamp.name} may not fit at {clamp.height:.2f} "
                        f"with length {clamp.length:.2f}; target in-line "
                        f"reaches from {elem.bottom_height:.2f}"
                        f" to {elem.top_height:.2f}."
                    )
                )
        except ValueError as e:
            print(f"Error: {e}")
        clamp.bool_fit = bool_fit

        # Sort clampons by height
        clamp_heights = [clamp.height for clamp in elem.clamp_ons]
        sorted_idx = np.argsort(clamp_heights)
        elem.clamp_ons = [elem.clamp_ons[i] for i in sorted_idx]

    # @property
    # def interpolate(self) -> "Mooring":
    #     interpolated_inline = []
    #     for inline in self.inline:
    #         interpolated_inline = interpolated_inline + inline.interpolate
    #     interpolated_mooring = copy.deepcopy(self)
    #     interpolated_mooring.inline = interpolated_inline
    #     return interpolated_mooring

    def make_df_inline(self):
        """Create a dataframe for in-line elements."""
        inline_data = []
        for elem in self.inline:
            inline_data.append(elem.get_dict())
        return pd.DataFrame(inline_data)

    def make_df_clampon(self):
        """Create a dataframe for clamp-on elements."""
        # Gather all clampon elements from inline
        clampon_list = []
        for elem in self.inline:
            for clamp in elem.clamp_ons:
                clampon_list.append(clamp)

        # Create a DataFrame for clampon attributes
        num_clampon = len(clampon_list)
        if num_clampon == 0:
            return pd.DataFrame(columns=ClampOn.attributes)
        else:
            clampon_data = []
            for clamp in clampon_list:
                clampon_data.append(clamp.get_dict())
            return pd.DataFrame(clampon_data)

    def make_df_combined(self) -> pd.DataFrame:
        """Make combined dataframe of clamp-on and in-line elements."""
        num_clampon = 0
        combined_list = []
        for elem in self.inline:
            combined_list.append(elem.get_dict())
            for clamp in elem.clamp_ons:
                combined_list.append(clamp.get_dict())
                num_clampon += 1
        df = pd.DataFrame(combined_list)
        if num_clampon == 0:
            df = df.reindex(columns=df.columns.union(ClampOn.attributes))
        df["mooring"] = self.name

        # Sort columns
        desired_order = [
            "type",
            "name",
            "serial",
            "bottom_height",
            "height",
            "length",
            "section",
        ]
        new_order = desired_order + [
            col for col in df.columns if col not in desired_order
        ]
        df_reordered = df[new_order]
        return df_reordered

    def export_csv(self, file_path: str) -> None:
        """Export combined dataframe to .csv"""
        # Make dir if not existing
        file_dir = os.path.dirname(file_path)
        os.makedirs(file_dir, exist_ok=True)

        df = self.make_df_combined()
        df.to_csv(file_path, index=False)

    def plot_design(
        self,
        label_rigging: bool = True,
        drop_strings: bool = None,
        show_serial: bool = True,
        show_length: bool = True,
        fontsize: int = 6,
        line_ratio_plot: float = 0.5,
    ):
        """
        Plot the structure of a mooring system including inline elements and clamp-ons.
        The function adjusts the plot dimensions based on specified ratios and minimum
        sizes. Elements are color-coded by type based on a predefined color dictionary.


        Parameters
        ----------
        mooring : object
            A mooring object containing attributes and metadata for inline elements
            and clamp-ons.
        label_rigging : bool, optional
            Whether to label rigging components (default is True).
        drop_strings : None, str or list(str), optional
            Strings to be removed from element names for display (default is None).
        show_serial : bool, optional
            Whether to display serial numbers of components (default is True).
        show_length : bool, optional
            Whether to display the lengths or vertical positions of components (default
            is True).
        fontsize : int, optional
            Font size for labels in the plot (default is 6).
        line_ratio_plot : float, optional
            Proportion of the plot height allocated to inline line elements (default is
            0.5).

        Returns
        -------
        fig : matplotlib.figure.Figure
            The generated matplotlib figure object.
        ax : matplotlib.axes._axes.Axes
            The matplotlib axes object for the plot.
        """
        return plot_design(
            self,
            label_rigging=label_rigging,
            drop_strings=drop_strings,
            show_serial=show_serial,
            show_length=show_length,
            fontsize=fontsize,
            line_ratio_plot=line_ratio_plot,
        )
