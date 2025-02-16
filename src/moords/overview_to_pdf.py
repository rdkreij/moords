"""Module to generate mooring design summary in pdf."""

import math
import os

import pandas as pd

import moords.format_data as format_data
import moords.overview_tools as overview_tools


def df_to_latex_table(df: pd.DataFrame, caption: str, rows_max: int = 55) -> str:
    """Dataframe to latex table."""

    # Determine the number of parts
    num_parts = math.ceil(len(df) / rows_max)

    latex_str = ""

    for part in range(num_parts):
        # Start building the LaTeX string for each part
        latex_str += r"\begin{table}[!htbp]" + "\n"
        latex_str += r"\centering" + "\n"

        # Adjust the caption based on the number of parts
        if num_parts == 1:
            latex_str += r"\caption{" + caption + "}" + "\n"
        else:
            latex_str += (
                r"\caption{"
                + caption
                + f" (Part {part + 1} of {num_parts})"
                + "}"
                + "\n"
            )

        # Make table
        latex_str += r"\begin{tabular}{" + "l" * len(df.columns) + "}" + "\n"
        latex_str += r"\toprule" + "\n"
        latex_str += " & ".join(df.columns) + r" \\" + "\n"
        latex_str += r"\midrule" + "\n"

        # Get the specific rows for this part
        start_row = part * rows_max
        end_row = min((part + 1) * rows_max, len(df))
        df_part = df.iloc[start_row:end_row]

        for _, row in df_part.iterrows():
            latex_str += " & ".join(map(str, row.values)) + r" \\" + "\n"

        latex_str += r"\bottomrule" + "\n"
        latex_str += r"\end{tabular}" + "\n"
        latex_str += r"\end{table}" + "\n\n"

    return latex_str


def generate_latex_summary(df: pd.DataFrame, header: str = None) -> str:
    """Genrate latex summary of mooring design in design."""
    rows_max = 55
    moorings = df["mooring"].unique()

    latex_str = ""

    # Configure latex document
    latex_str += r"\documentclass{article}" + "\n"
    latex_str += r"\usepackage{tabularx}" + "\n"
    latex_str += r"\usepackage{booktabs}" + "\n"
    latex_str += r"\usepackage[margin=.5in]{geometry}" + "\n"
    if header is not None:
        latex_str += r"\usepackage{fancyhdr}" + "\n"
        latex_str += r"\pagestyle{fancy}" + "\n"
        latex_str += r"\renewcommand{\headrulewidth}{0pt}" + "\n"
        latex_str += r"\fancyhead[C]{" + f"{header}" + r"}" + "\n"
    latex_str += r"\begin{document}" + "\n\n"

    for mooring in moorings:
        df_element_summary, caption = overview_tools.make_df_element_summary(
            df, mooring
        )
        latex_str += df_to_latex_table(df_element_summary, caption, rows_max)

        df_clampon_summary, caption = overview_tools.make_df_clampon_summary(
            df, mooring
        )
        latex_str += df_to_latex_table(df_clampon_summary, caption, rows_max)

        df_assembly, caption = overview_tools.make_df_assambly(df, mooring)
        latex_str += df_to_latex_table(df_assembly, caption, rows_max)

    df_summary_element_all, caption = overview_tools.make_df_summary_element_all(df)
    latex_str += df_to_latex_table(df_summary_element_all, caption, rows_max)

    df_count_all, caption = overview_tools.make_df_count_all(df)
    latex_str += df_to_latex_table(df_count_all, caption, rows_max)

    df_simple_section_sum, caption = overview_tools.make_df_simple_section_sum(df)
    latex_str += df_to_latex_table(df_simple_section_sum, caption, rows_max)

    df_list, info_list = overview_tools.make_df_list_section_all(df)
    for idx, df_i in enumerate(df_list):
        latex_str += df_to_latex_table(df_i, info_list[idx], rows_max)

    latex_str += r"\end{document}"
    return latex_str


def generate_overview_pdf(
    df: pd.DataFrame,
    file_path: str = "",
    header: str = None,
    replacements: str | list[str] = None,
):
    """Generate mooring design summary in pdf using latex."""

    # Make dir if not existing
    file_dir = os.path.dirname(file_path)
    os.makedirs(file_dir, exist_ok=True)

    file_name = os.path.splitext(os.path.basename(file_path))[0]

    if replacements is not None:
        df = format_data.replace_in_df(df, replacements)
    df = format_data.remove_trailing_spaces_in_df(df)

    latex_str = generate_latex_summary(df, header)

    tex_file_path = f"{file_dir}/{file_name}.tex"

    # Write to .tex file
    with open(tex_file_path, "w", encoding="utf-8") as tex_file:
        tex_file.write(latex_str)

    # Convert to pdf
    os.system(f"pdflatex -output-directory={file_dir} {tex_file_path}")
