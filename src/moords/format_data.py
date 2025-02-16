"""Module to format data in an array or dataframe."""

import pandas as pd


def replace_in_df(
    df: pd.DataFrame, replacements: list[tuple[str, str]]
) -> pd.DataFrame:
    """Replace strings in a dataframe."""
    df_copy = df.copy()

    # Define a function to apply replacements to strings
    def escape_and_replace(s):
        if isinstance(s, str):
            for old, new in replacements:
                s = s.replace(old, new)
            return s
        return s  # Leave non-strings unchanged

    # Apply the function to each column
    return df_copy.apply(
        lambda col: col.map(escape_and_replace) if col.dtype == "object" else col
    )


def remove_trailing_spaces_in_df(df: pd.DataFrame) -> pd.DataFrame:
    """Remove trailing spaces from string columns."""
    return df.apply(
        lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x)
        if col.dtype == "object"
        else col
    )


def format_array(arr, format_str="{}") -> list:
    """Format an array to a string given the format."""
    formatted = []
    for item in arr:
        if item is None:
            formatted.append("")
        else:
            formatted.append(format_str.format(item))
    return formatted
