"""Module for plotting mooring design"""

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Rectangle


def plot_design(
    mooring,
    label_rigging=True,
    drop_strings=None,
    show_serial=True,
    show_length=True,
    fontsize=6,
    line_ratio_plot=0.5,
):
    """
    Plot the structure of a mooring system including inline elements and clamp-ons.
    The function adjusts the plot dimensions based on specified ratios and minimum
    sizes. Elements are color-coded by type based on a predefined color dictionary.


    Parameters
    ----------
    mooring : object
        A mooring object containing attributes and metadata for inline elements and
        clamp-ons.
    label_rigging : bool, optional
        Whether to label rigging components (default is True).
    drop_strings : None, str or list(str), optional
        Strings to be removed from element names for display (default is None).
    show_serial : bool, optional
        Whether to display serial numbers of components (default is True).
    show_length : bool, optional
        Whether to display the lengths or vertical positions of components (default is
        True).
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
    # Set properties
    minimum_element_width = 0.008
    minimum_element_length = minimum_element_width
    minimum_line_length = 0.05
    fig_width_in, fig_height_in = 3, 10
    text_height_in = fontsize / 72  # Approximate text height in inches
    relative_height = (
        text_height_in / fig_height_in
    )  # Text height relative to figure height
    x_offset_inline_text = -0.10
    y_offset_text = relative_height * 1.6
    x_buffer_inline_text = 0.005
    x_buffer_sub_text = 0.25
    x_line_margin = 0.015
    x_offset_clamp = 0.1
    x_offset_clamp_space = 0.01
    x_offset_clamp_text = 0.1
    x_buffer_clamp_text = x_buffer_inline_text
    color_dict = {
        "wires": "k",
        "chains": "k",
        "floats": "orange",
        "anchors": "grey",
        "rigging": "brown",
        "miscs": "g",
        "arcels": "b",
        "cms": "r",
    }

    # Inline properties
    bottom_height = np.array([inline.bottom_height for inline in mooring.inline])
    top_height = np.array([inline.top_height for inline in mooring.inline])
    width = np.array([inline.width for inline in mooring.inline])
    diameter = np.array([inline.diameter for inline in mooring.inline])
    bool_line = np.array([inline.bool_line for inline in mooring.inline])
    type_inline = np.array([inline.type for inline in mooring.inline])
    name = np.array([inline.name for inline in mooring.inline])
    serial = np.array([inline.serial for inline in mooring.inline])

    # Compute additional variables
    n_inline = len(mooring.inline)
    length = top_height - bottom_height
    total_length = np.sum(length)
    total_line_length = np.sum(length * bool_line)
    line_ratio = total_line_length / total_length
    line_factor = line_ratio_plot / line_ratio / total_length
    non_line_ratio_plot = 1 - line_ratio_plot
    non_line_ratio = 1 - line_ratio
    non_line_factor = non_line_ratio_plot / non_line_ratio / total_length
    length_plot = (
        bool_line * line_factor + np.logical_not(bool_line) * non_line_factor
    ) * length
    bool_sphere = diameter > 0

    # Modify inline plot length according to ratio
    line_length_plot = length_plot[bool_line]
    partitions = line_length_plot / np.sum(line_length_plot)
    adjusted_partitions = calculate_partition_with_min(partitions, minimum_line_length)
    length_plot[bool_line] = adjusted_partitions * np.sum(line_length_plot)
    top_height_plot = np.cumsum(length_plot)
    bottom_height_plot = np.cumsum(np.hstack([[0], length_plot[:-1]]))

    # Modify inline plot width according minimum width
    width_plot = (
        bool_line * line_factor + np.logical_not(bool_line) * non_line_factor
    ) * width
    width_plot[width_plot < minimum_element_width] = minimum_element_width

    # Clampon properties
    heightco = []
    heightco_along_inline = []
    heightco_plot = []
    lengthco_plot = []
    widthco_plot = []
    diameterco = []
    typeco = []
    nameco = []
    serialco = []
    for i in range(n_inline):
        clampons = mooring.inline[i].clamp_ons
        for clamp in clampons:
            ratio_clamp = clamp.height_along_inline / length[i]
            heightco_plot.append(bottom_height_plot[i] + ratio_clamp * length_plot[i])
            lengthco_plot.append(non_line_factor * clamp.length)
            widthco_plot.append(non_line_factor * clamp.width)
            heightco.append(clamp.height)
            diameterco.append(clamp.diameter)
            heightco_along_inline.append(clamp.height_along_inline)
            typeco.append(clamp.type)
            nameco.append(clamp.name)
            serialco.append(clamp.serial)
    heightco_plot = np.array(heightco_plot)
    heightco = np.array(heightco)
    heightco_along_inline = np.array(heightco_along_inline)
    lengthco_plot = np.array(lengthco_plot)
    widthco_plot = np.array(widthco_plot)
    diameterco = np.array(diameterco)
    typeco = np.array(typeco)
    nameco = np.array(nameco)
    n_clampon = len(heightco_plot)
    bool_sphereco = diameterco > 0

    # Modify clampon plot dimensions
    widthco_plot[widthco_plot < minimum_element_width] = minimum_element_width
    lengthco_plot[lengthco_plot < minimum_element_length] = minimum_element_length

    x_positions_clamp = update_positions_no_overlap(
        heightco_plot, widthco_plot, lengthco_plot, x_offset_clamp, x_offset_clamp_space
    )

    # Remove unwanted strings
    if drop_strings is not None:
        if not isinstance(drop_strings, list):
            drop_strings = [drop_strings]

        for drop in drop_strings:
            name = np.array([i.replace(drop, "").strip() for i in name])
            nameco = np.array([i.replace(drop, "").strip() for i in nameco])

    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in))
    ax.set_aspect("equal")
    ax.text(0, 1.02, mooring.name, va="bottom", ha="center")

    # Plot inline
    for i in range(n_inline):
        # Plot shape
        if bool_line[i]:
            plt.plot(
                [0, 0],
                [bottom_height_plot[i], top_height_plot[i]],
                "-",
                color=color_dict[type_inline[i]],
                lw=2,
                zorder=-2,
            )
        elif bool_sphere[i]:
            radius = length_plot[i] / 2
            circle = Circle(
                (0, bottom_height_plot[i] + radius),
                radius,
                ec="k",
                fc=color_dict[type_inline[i]],
                lw=1,
            )
            ax.add_patch(circle)
        else:
            rectangle = Rectangle(
                [-width_plot[i] / 2, bottom_height_plot[i]],
                width_plot[i],
                length_plot[i],
                ec="k",
                fc=color_dict[type_inline[i]],
                lw=1,
            )
            ax.add_patch(rectangle)

        # Add label
        if label_rigging or (type_inline[i] != "rigging"):
            if i == 0:
                y_text = bottom_height_plot[i] + length_plot[i] / 2
            else:
                if (y_text + y_offset_text) >= (
                    bottom_height_plot[i] + length_plot[i] / 2
                ):
                    y_text = y_text + y_offset_text
                else:
                    y_text = bottom_height_plot[i] + length_plot[i] / 2

            text = name[i]
            if show_serial:
                if serial[i] is not None:
                    text = f"#{serial[i]} - {text}"
            ax.text(
                x_offset_inline_text,
                y_text,
                text,
                va="center",
                ha="right",
                fontsize=fontsize,
            )

            if show_length:
                if bool_line[i]:
                    sub_text = f" [{length[i]:.1f}m sect.]"
                else:
                    sub_text = f" [{bottom_height[i]:.1f}m to {top_height[i]:.1f}m ASB]"
                ax.text(
                    x_offset_inline_text - x_buffer_sub_text,
                    y_text,
                    sub_text,
                    va="center",
                    ha="right",
                    fontsize=fontsize,
                )

            ax.plot(
                [
                    x_offset_inline_text + x_buffer_inline_text,
                    x_offset_inline_text + x_line_margin,
                    -x_line_margin,
                    0,
                ],
                [
                    y_text,
                    y_text,
                    bottom_height_plot[i] + length_plot[i] / 2,
                    bottom_height_plot[i] + length_plot[i] / 2,
                ],
                "k-",
                lw=0.5,
                zorder=-5,
            )

    # Plot clampon
    for j in range(n_clampon):
        if bool_sphereco[j]:
            radius = lengthco_plot[j] / 2
            circle = Circle(
                (x_positions_clamp[j], heightco_plot[j] + radius),
                radius,
                ec="k",
                fc=color_dict[typeco[j]],
                lw=1,
            )
            ax.add_patch(circle)
        else:
            rect = Rectangle(
                (
                    x_positions_clamp[j] - widthco_plot[j] / 2,
                    heightco_plot[j] - lengthco_plot[j] / 2,
                ),
                widthco_plot[j],
                lengthco_plot[j],
                ec="k",
                fc=color_dict[typeco[j]],
                lw=1,
            )
            ax.add_patch(rect)

        # Add label
        if j == 0:
            y_text = heightco_plot[j]
        else:
            if (y_text + y_offset_text) >= (heightco_plot[j]):
                y_text = y_text + y_offset_text
            else:
                y_text = heightco_plot[j]
        text = nameco[j]
        if show_serial:
            if serialco[j] is not None:
                text = f"{text} #{serialco[j]}"
        ax.text(
            x_offset_clamp + x_offset_clamp_text,
            y_text,
            text,
            va="center",
            ha="left",
            fontsize=fontsize,
        )

        if show_length:
            sub_text = f"[{heightco[j]:.1f}m ASB | {heightco_along_inline[j]:.1f}m AE]"
            ax.text(
                x_offset_clamp + x_offset_clamp_text + x_buffer_sub_text,
                y_text,
                sub_text,
                va="center",
                ha="left",
                fontsize=fontsize,
            )
        ax.plot(
            [
                x_offset_clamp + x_offset_clamp_text - x_buffer_clamp_text,
                x_offset_clamp + x_offset_clamp_text - x_line_margin,
                x_offset_clamp + 2 * x_line_margin,
                x_positions_clamp[j],
            ],
            [y_text, y_text, heightco_plot[j], heightco_plot[j]],
            "k-",
            lw=0.5,
            zorder=-5,
        )

    ax.axis("off")

    return fig, ax


def calculate_partition_with_min(partitions, min_size):
    """Adjust partitions given minimum size"""
    adjustments_needed = partitions < min_size
    adjusted_partitions = partitions.copy()
    adjusted_partitions[adjustments_needed] = min_size
    increase_needed = np.sum(
        adjusted_partitions[adjustments_needed] - partitions[adjustments_needed]
    )
    larger_partitions = ~adjustments_needed
    scaling_factor = 1 - increase_needed
    adjusted_partitions[larger_partitions] *= scaling_factor / np.sum(
        partitions[larger_partitions]
    )
    adjusted_partitions /= np.sum(adjusted_partitions)  # esure they sum to 1
    return adjusted_partitions


def update_positions_no_overlap(heights, widths, length, initial_x, spacing_x):
    """Update positions of shapes to prevent overlap"""
    x_positions = [initial_x] * len(heights)
    n = len(heights)
    for i in range(n):
        for j in range(i):
            overlap_x = (
                abs(x_positions[i] - x_positions[j])
                < (widths[i] + widths[j]) / 2 + spacing_x
            )
            overlap_y = abs(heights[i] - heights[j]) < (length[i] + length[j]) / 2
            if overlap_x and overlap_y:
                x_positions[i] = (
                    x_positions[j] + (widths[j] + widths[i]) / 2 + spacing_x
                )
    return x_positions


def generate_design_pdf(file_path, list_fig):
    """Save generated mooring figures"""
    # Ensure list
    list_fig = list_fig if isinstance(list_fig, list) else [list_fig]

    # Make dir if not existing
    file_dir = os.path.dirname(file_path)
    os.makedirs(file_dir, exist_ok=True)

    with PdfPages(file_path) as pdf:
        for fig in list_fig:
            pdf.savefig(fig, bbox_inches="tight")
