"""Module with tools to generate an overview of the mooring design (incl. serial number
checks)."""

import numpy as np
import pandas as pd

import moords.format_data as format_data
from moords.mooring_design import Mooring


def make_df_combined_moorings(moorings: list[Mooring]) -> pd.DataFrame:
    """Combine dataframes of multiple moorings."""
    list_df = [mooring.make_df_combined() for mooring in moorings]
    return pd.concat(list_df).reset_index(drop=True)


def make_df_element_summary(df: pd.DataFrame, mooring: str) -> tuple[pd.DataFrame, str]:
    """Make a summary dataframe of all elements in the specified mooring."""
    df_mooring = df.loc[df["mooring"] == mooring]

    value_counts = df_mooring.name.value_counts().sort_index()
    names = format_data.format_array(value_counts.index)
    counts = format_data.format_array(value_counts.values, "{:.0f}")

    line_length = []
    for item in names:
        total_length = 0
        idx = np.where(df_mooring["name"] == item)[0]
        for i in idx:
            data = df_mooring.iloc[i]
            if data["bool_line"]:
                total_length += data["length"]
        if total_length == 0:
            line_length.append("")
        else:
            line_length.append(f"{total_length:.1f}")

    buoyancy = []
    for item in names:
        total_buoyancy = 0
        idx = np.where(df_mooring["name"] == item)[0]
        for i in idx:
            data = df_mooring.iloc[i]
            total_buoyancy += data["buoyancy"]
        if data["bool_line"]:
            buoyancy.append(f"{total_buoyancy:.1f}")
        else:
            buoyancy.append(f"{total_buoyancy / len(idx):.1f}")

    df_element_summary = pd.DataFrame(
        {
            "Element": names,
            "Qty": counts,
            "Length [m]": line_length,
            "Buoyancy [kg]": buoyancy,
        }
    )
    info = f"{mooring}: Element summary"
    return (df_element_summary, info)


def make_df_clamp_with_summary(
    df: pd.DataFrame, mooring: str
) -> tuple[pd.DataFrame, str]:
    df_mooring = df.loc[df["mooring"] == mooring]

    clampon_names = df_mooring.loc[df_mooring["bool_clampon"], "name"].unique()
    clampon_names = sorted(clampon_names)

    clamp_with_info = {}
    for clampon_name in clampon_names:
        idx_clampon = df_mooring["name"] == clampon_name
        df_clampon = df_mooring.loc[idx_clampon]
        for i in range(len(df_clampon)):
            data = df_clampon.iloc[i]
            clamped_with = data["clamped_with"]
            if (
                clamped_with is not None
                and clamped_with != ""
                and str(clamped_with) != "nan"
            ):
                if clampon_name not in clamp_with_info:
                    clamp_with_info[clampon_name] = {}
                if clamped_with not in clamp_with_info[clampon_name]:
                    clamp_with_info[clampon_name][clamped_with] = 0
                clamp_with_info[clampon_name][clamped_with] += 1

    for key in clamp_with_info:
        clamp_with_info[key] = dict(
            sorted(clamp_with_info[key].items(), key=lambda item: item[1], reverse=True)
        )

    clamp_names = []
    clamped_with_names = []
    qtys = []
    for clamp_name, clamped_with_dict in clamp_with_info.items():
        for i, (clamped_with_name, qty) in enumerate(clamped_with_dict.items()):
            if i == 0:
                clamp_names.append(f"{clamp_name}")
            else:
                clamp_names.append("")
            clamped_with_names.append(f"{clamped_with_name}")
            qtys.append(f"{qty:.0f}")

    df_clamp_with_summary = pd.DataFrame(
        {"Clamp-on": clamp_names, "Clamped with": clamped_with_names, "Qty": qtys}
    )
    info = f"{mooring}: Clamped with summary"
    return (df_clamp_with_summary, info)


def make_df_clampon_summary(df: pd.DataFrame, mooring: str) -> tuple[pd.DataFrame, str]:
    """Make a summary of all the clamp-on elements in the specified mooring."""
    df_mooring = df.loc[df["mooring"] == mooring]

    idx = df_mooring["bool_clampon"]
    df_crop = df_mooring.loc[idx]
    name = format_data.format_array(df_crop["name"])
    serial = format_data.format_array(df_crop["serial"])
    height = format_data.format_array(df_crop["height"], "{:.1f}")
    clamped_to_name = format_data.format_array(df_crop["clamped_to_name"])
    section = format_data.format_array(df_crop["section"])
    height_along_inline = format_data.format_array(
        df_crop["height_along_inline"], "{:.1f}"
    )
    clamped_with = format_data.format_array(df_crop["clamped_with"])

    df_clampon_summary = pd.DataFrame(
        {
            "Name": name,
            "Serial": serial,
            "H[m]": height,
            "Clamped to": clamped_to_name,
            "Clamped with": clamped_with,
            "Section": section,
            "Inline H[m]": height_along_inline,
        }
    )
    info = f"{mooring}: Clamp-on summary"
    return (df_clampon_summary, info)


def make_df_assambly(df: pd.DataFrame, mooring: str) -> tuple[pd.DataFrame, str]:
    """Make an assambly summary for the specified mooring."""
    df_mooring = df.loc[df["mooring"] == mooring]

    name = format_data.format_array(df_mooring["name"])
    serial = format_data.format_array(df_mooring["serial"])

    line_length = []
    for i in range(len(df_mooring)):
        data = df_mooring.iloc[i]
        if data["bool_line"]:
            line_length.append(f"{data['length']:.1f}")
        else:
            line_length.append("")
    section = format_data.format_array(df_mooring["section"])
    height = []
    for i in range(len(df_mooring)):
        data = df_mooring.iloc[i]
        if data["bool_clampon"]:
            height.append(
                f"{data['height']:.1f} [{data['height_along_inline']:.1f} m AE]"
            )
        else:
            height.append(f"{data['bottom_height']:.1f}")
    clamped_with = format_data.format_array(df_mooring["clamped_with"])

    df_assembly = pd.DataFrame(
        {
            "Element": name,
            "Serial": serial,
            "Length [m]": line_length,
            "Section": section,
            "Clamped with": clamped_with,
            "Height [in ASB]": height,
        }
    )
    info = f"{mooring}: Assembly summary"
    return (df_assembly, info)


def make_df_count_all(df) -> tuple[pd.DataFrame, str]:
    """Make summary of the total count of all elements across all moorings."""
    count_df = df.groupby(["name", "mooring"]).size().unstack(fill_value=0)
    count_df["Total"] = count_df.sum(axis=1)
    df_count_all = count_df.reset_index()
    # turn into string with no decimals
    for col in df_count_all.columns[1:]:
        df_count_all[col] = df_count_all[col].apply(lambda x: f"{x:.0f}")
    df_count_all.rename(columns={"name": "Name"}, inplace=True)
    info = "All : Total count of each element"
    return (df_count_all, info)


def make_df_count_clamped_with_all(df) -> tuple[pd.DataFrame, str]:
    clampon_names = df.loc[df["bool_clampon"], "name"].unique()
    clampon_names = sorted(clampon_names)
    clamped_with_info = {}

    mooring_names = df["mooring"].unique()
    # keep track of each count per mooring
    for clampon_name in clampon_names:
        for mooring in mooring_names:
            idx_clampon = (
                (df["name"] == clampon_name)
                & (df["mooring"] == mooring)
                & df["bool_clampon"]
            )
            df_clampon = df.loc[idx_clampon]
            for i in range(len(df_clampon)):
                data = df_clampon.iloc[i]
                clamped_with = data["clamped_with"]
                if (
                    clamped_with is not None
                    and clamped_with != ""
                    and str(clamped_with) != "nan"
                ):
                    if clampon_name not in clamped_with_info:
                        clamped_with_info[clampon_name] = {}
                    if clamped_with not in clamped_with_info[clampon_name]:
                        clamped_with_info[clampon_name][clamped_with] = {}
                    if mooring not in clamped_with_info[clampon_name][clamped_with]:
                        clamped_with_info[clampon_name][clamped_with][mooring] = 0
                    clamped_with_info[clampon_name][clamped_with][mooring] += 1

    # sort
    for clampon in clamped_with_info:
        clamped_with_info[clampon] = dict(
            sorted(
                clamped_with_info[clampon].items(),
                key=lambda item: sum(item[1].values()),
                reverse=True,
            )
        )

    # make dataframe
    clampon_names = []
    clamped_with_names = []
    qtys_per_mooring = {mooring: [] for mooring in mooring_names}
    total = []

    for clampon_name, clamped_with_dict in clamped_with_info.items():
        for i, (clamped_with_name, mooring_dict) in enumerate(
            clamped_with_dict.items()
        ):
            if i == 0:
                clampon_names.append(clampon_name)
            else:
                clampon_names.append("")

            clamped_with_names.append(clamped_with_name)

            for mooring in mooring_names:
                qtys_per_mooring[mooring].append(f"{mooring_dict.get(mooring, 0):.0f}")

            total.append(f"{sum(mooring_dict.values()):.0f}")

    df_count_clamped_with_all = pd.DataFrame(
        {
            "Clamp-on": clampon_names,
            "Clamped with": clamped_with_names,
            **qtys_per_mooring,
            "Total": total,
        }
    )

    info = "All : Total count clamp types"
    return (df_count_clamped_with_all, info)


def make_df_simple_section_sum(df) -> tuple[pd.DataFrame, str]:
    """Make a section summary across all moorings."""
    idx_line = df["bool_line"]
    df_crop = df.loc[idx_line]
    mooring = format_data.format_array(df_crop["mooring"])
    section = format_data.format_array(df_crop["section"])
    name = format_data.format_array(df_crop["name"])
    length = format_data.format_array(df_crop["length"], "{:.1f}")
    df_simple_section_sum = pd.DataFrame(
        {"Mooring": mooring, "Section": section, "Material": name, "Length": length}
    )
    info = "All: simple section summary"
    return (df_simple_section_sum, info)


def make_df_summary_element_all(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """Make summary of instruments across all moorings."""
    idx_instruments = (
        (df["type"] == "miscs") | (df["type"] == "arcels") | (df["type"] == "cms")
    )
    dfc = df.loc[idx_instruments]

    unique_elem = dfc["name"].unique()
    unique_elem = sorted(unique_elem)

    list_collect = []
    for elem in unique_elem:
        idx_elem = np.where(dfc["name"] == elem)[0]
        n = len(idx_elem)

        dfcc = dfc.iloc[idx_elem]
        unique_mooring = dfcc["mooring"].unique()

        serial_list = []
        clamped_with_list = []
        mooring_list = []

        for mooring in unique_mooring:
            idx_mooring = dfcc["mooring"] == mooring
            dfccc = dfcc.loc[idx_mooring]
            dfccc = dfccc.sort_values(by="serial")
            serial_list += list(dfccc["serial"])
            clamped_with_list += list(dfccc["clamped_with"])
            mooring_list += list(dfccc["mooring"])
        serial_list = format_data.format_array(serial_list)
        clamped_with_list = format_data.format_array(clamped_with_list)
        mooring_list = format_data.format_array(mooring_list)

        for i in range(n):
            serial = serial_list[i]
            clamped_with = clamped_with_list[i]
            mooring = mooring_list[i]
            if i == 0:
                list_collect.append(
                    {
                        "Element": f"{elem} [n={n:.0f}]",
                        "Serial number": serial,
                        "Clamped with": clamped_with,
                        "Mooring": mooring,
                    }
                )
            else:
                list_collect.append(
                    {
                        "Element": "",
                        "Serial number": serial,
                        "Clamped with": clamped_with,
                        "Mooring": mooring,
                    }
                )
    df_summary_element_all = pd.DataFrame(list_collect)
    info = "All: summary of instruments"
    return df_summary_element_all, info


def make_df_list_section_all(df: pd.DataFrame) -> tuple[list[pd.DataFrame], list[str]]:
    """Make summary of each section across all moorings"""

    df_list = []
    info_list = []

    moorings = df["mooring"].unique()
    for mooring in moorings:
        idx_mooring = df["mooring"] == mooring
        df_crop = df.loc[idx_mooring]

        sections = df_crop["section"].unique()
        sections = [x for x in sections if x is not None]  # drop None

        for section in sections:
            idx_section = df_crop["section"] == section
            idx_clampon = df_crop["bool_clampon"]
            data = df_crop.loc[idx_section & ~idx_clampon].iloc[0].to_dict()
            info_list.append(
                (
                    f"Summary of {data['mooring']}, "
                    f"{section}, "
                    f"length: {data['length']:.1f} m, "
                    f"material: {data['name']}"
                )
            )

            dfcc = df_crop.loc[idx_section & idx_clampon]
            name = format_data.format_array(dfcc["name"])
            serial = format_data.format_array(dfcc["serial"])
            height = format_data.format_array(dfcc["height"], "{:.1f}")
            height_along_inline = format_data.format_array(
                dfcc["height_along_inline"], "{:.1f}"
            )

            df_list.append(
                pd.DataFrame(
                    {
                        "Name": name,
                        "Serial number": serial,
                        "Height [m]": height,
                        "Along element": height_along_inline,
                    }
                )
            )

    return df_list, info_list


def print_duplicated_serial(df: pd.DataFrame) -> None:
    """Print a summary of the duplicated serial numbers"""
    mask = df["serial"].notna()
    duplicated_idx = df.loc[mask, "serial"].duplicated(keep=False)
    duplicated_list = df.loc[mask & duplicated_idx, "serial"].unique()

    print("DUPLICATED SERIAL NUMBERS")
    if len(duplicated_list) == 0:
        print(" no duplicates")
    else:
        for serial in duplicated_list:
            idx = df["serial"] == serial
            df_crop = df.loc[idx]

            print(f"{serial}")
            for i in range(len(df_crop)):
                data = df_crop.iloc[i]
                if data["bool_clampon"]:
                    print(
                        (
                            f" {data['mooring']}, "
                            f"{data['name']}, "
                            f"height {data['height']:.1f}m, "
                            f"clamped to: {data['clamped_to_name']}"
                        )
                    )
                else:
                    print(
                        (
                            f" {data['mooring']}, "
                            f"{data['name']}, "
                            f"bottom height {data['bottom_height']:.1f}m"
                        )
                    )


def make_df_missing_serial(df: pd.DataFrame) -> pd.DataFrame:
    """Make dataframe of all elements with a missing serial number"""
    mask = ~df["serial"].notna() | (df["serial"] == "")
    return df.loc[mask]
