"""Module to generate mooring design summary in pdf."""

import os

import pandas as pd

import moords.format_data as format_data
import moords.overview_tools as overview_tools


def df_to_latex_table(df: pd.DataFrame, caption: str, rows_max: int = 55) -> str:
    """Dataframe to latex table."""
    max_len_array = [df[col].astype(str).str.len().max() for col in df.columns]
    sum_max_len = sum(max_len_array)

    num_rows, _ = df.shape

    col_specs = ["|"]
    if sum_max_len <= 70:
        col_specs.extend(["l|"] * len(df.columns))
    else:
        for col in df.columns:
            max_len = df[col].astype(str).str.len().max()
            if max_len > 15:
                col_specs.append("X|")
            else:
                col_specs.append("l|")

    latex_str = ""
    latex_str += r"\Needspace{" + str(int(num_rows + 5)) + r"\baselineskip}" + "\n"
    latex_str += r"\begin{tabularx}{\linewidth}{" + "".join(col_specs) + "}" + "\n"
    latex_str += r"\caption{" + caption + r"}\\" + "\n"
    latex_str += r"\hline" + "\n"
    latex_str += " & ".join([rf"\textbf{{{i}}}" for i in df.columns]) + r" \\" + "\n"
    latex_str += r"\hline" + "\n"
    for _, row in df.iterrows():
        latex_str += " & ".join(map(str, row.values)) + r" \\ \hline" + "\n"
    latex_str += r"\end{tabularx}" + "\n"

    return latex_str


def generate_latex_summary(df: pd.DataFrame, header: str | None = None) -> str:
    """Genrate latex summary of mooring design in design."""
    rows_max = 55
    moorings = df["mooring"].unique()

    latex_str = ""

    # Configure latex document
    latex_str += r"\documentclass{article}" + "\n"
    latex_str += r"\usepackage{ltablex}" + "\n"
    latex_str += r"\keepXColumns" + "\n"
    latex_str += r"\usepackage{needspace}" + "\n"
    latex_str += r"\usepackage{booktabs}" + "\n"
    latex_str += r"\usepackage[table]{xcolor}" + "\n"
    latex_str += r"\usepackage[margin=.5in]{geometry}" + "\n"
    if header is not None:
        latex_str += r"\usepackage{fancyhdr}" + "\n"
        latex_str += r"\pagestyle{fancy}" + "\n"
        latex_str += r"\renewcommand{\headrulewidth}{0pt}" + "\n"
        latex_str += r"\fancyhead[C]{" + f"{header}" + r"}" + "\n"
    latex_str += r"\begin{document}" + "\n\n"
    latex_str += r"\arrayrulecolor[gray]{0.7}" + "\n"

    for mooring in moorings:
        df_element_summary, caption = overview_tools.make_df_element_summary(
            df, mooring
        )
        latex_str += df_to_latex_table(df_element_summary, caption, rows_max)

        df_clampon_summary, caption = overview_tools.make_df_clampon_summary(
            df, mooring
        )
        latex_str += df_to_latex_table(df_clampon_summary, caption, rows_max)

        df_clamp_with_summary, caption = overview_tools.make_df_clamp_with_summary(
            df, mooring
        )
        latex_str += df_to_latex_table(df_clamp_with_summary, caption, rows_max)

        df_assembly, caption = overview_tools.make_df_assambly(df, mooring)
        latex_str += df_to_latex_table(df_assembly, caption, rows_max)

    df_summary_element_all, caption = overview_tools.make_df_summary_element_all(df)
    latex_str += df_to_latex_table(df_summary_element_all, caption, rows_max)

    df_count_all, caption = overview_tools.make_df_count_all(df)
    latex_str += df_to_latex_table(df_count_all, caption, rows_max)

    df_count_clamped_with_all, caption = overview_tools.make_df_count_clamped_with_all(
        df
    )
    latex_str += df_to_latex_table(df_count_clamped_with_all, caption, rows_max)

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
    header: str | None = None,
    replacements: list[tuple[str, str]] | None = None,
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
