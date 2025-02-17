"""Module to plot mooring simulation results"""

import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def plot_simulated_instance(
    ds_sim_instance: xr.Dataset, ds_flow_instance: xr.Dataset, title: str = ""
):
    """Plot simulated mooring instance"""
    fig, ax = plt.subplots(
        1, 6, figsize=(14, 5), sharey=True, gridspec_kw={"wspace": 0.2}
    )
    fig.suptitle(title)
    tension = ds_sim_instance.tension
    ax[0].plot(ds_flow_instance.U, ds_flow_instance.z, "k--", label="$U$")
    ax[0].plot(ds_flow_instance.V, ds_flow_instance.z, "r", label="$V$", zorder=-1)
    im = ax[1].scatter(ds_sim_instance.X, ds_sim_instance.Z, 5, c=tension)
    ax[2].scatter(ds_sim_instance.Y, ds_sim_instance.Z, 5, c=tension)
    ax[3].scatter(ds_sim_instance.dZ, ds_sim_instance.Z, 5, c=tension)
    ax[4].scatter(ds_sim_instance.psi, ds_sim_instance.Z, 5, c=tension)
    ax[5].scatter(tension, ds_sim_instance.Z, 5, c=tension)

    ax[0].set(ylabel="$z$ (m)", xlabel="vel (m/s)")
    ax[1].set(xlabel="$x$ (m)")
    ax[2].set(xlabel="$y$ (m)")
    ax[3].set(xlabel="$\\Delta z$ (m)")
    ax[4].set(xlabel="$\\psi$ ($\\degree$)")
    ax[5].set(xlabel="$T$ (kg)")

    for axi in ax:
        axi.axhline(ds_flow_instance.bottom_depth, ls=":", c="b", lw=1)
        axi.set_ylim([0, 1.1 * ds_flow_instance.bottom_depth])

    cbar = fig.colorbar(im, ax=ax, location="right", shrink=1, aspect=40, pad=0.01)
    cbar.set_label("$T$ (kg)")
    ax[0].legend(frameon=False, loc="best")
    return fig, ax


def plot_flow_series(ds_flow_series):
    """
    Plots flow series data including U, V, W velocity components and density (rho)
    """
    # Flip the data arrays for proper orientation
    U = np.flip(ds_flow_series.U.values.T, axis=0)
    V = np.flip(ds_flow_series.V.values.T, axis=0)
    W = np.flip(ds_flow_series.W.values.T, axis=0)
    rho = np.flip(ds_flow_series.rho.values.T, axis=0)

    z = ds_flow_series.z.values
    time = ds_flow_series.time.values
    time_numeric = mdates.date2num(time)
    bottom_depth = ds_flow_series.bottom_depth.item()

    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    # Define velocity and density limits for color scaling
    vmin, vmax = -np.max(np.abs([U, V])), np.max(np.abs([U, V]))
    vminW, vmaxW = -np.max(np.abs(W)), np.max(np.abs(W))

    # Plot each variable
    variables = [U, V, W, rho]
    labels = ["$U$ (m/s)", "$V$ (m/s)", "$W$ (m/s)", "$\\rho$ (kg/m$^3$)"]
    cmaps = [
        matplotlib.cm.seismic,
        matplotlib.cm.seismic,
        matplotlib.cm.seismic,
        matplotlib.cm.viridis,
    ]
    vmins = [vmin, vmin, vminW, None]
    vmaxs = [vmax, vmax, vmaxW, None]

    for i, (var, label, cmap, vmin, vmax) in enumerate(
        zip(variables, labels, cmaps, vmins, vmaxs)
    ):
        im = axes[i].imshow(
            var,
            aspect="auto",
            cmap=cmap,
            extent=[time_numeric[0], time_numeric[-1], z[0], z[-1]],
            vmin=vmin,
            vmax=vmax,
        )
        cbar = fig.colorbar(
            im, ax=axes[i], location="right", shrink=1, aspect=20, pad=0.01
        )
        cbar.set_label(label)

        axes[i].set_ylabel("$z$ (m)")
        axes[i].invert_yaxis()
        axes[i].axhline(bottom_depth, ls="--", c="b", label="Bottom Depth")
        axes[i].set_ylim([0, bottom_depth * 1.1])

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    axes[-1].set_xlabel("Time")
    return fig, axes


def plot_anchor_stats_series(ds_sim_series_interpolated):
    """Plot simulated anchor statistics as a function of time"""
    time = ds_sim_series_interpolated.time
    tension = ds_sim_series_interpolated.tension.isel(idx=0)
    psi = ds_sim_series_interpolated.psi.isel(idx=0)

    tension_vertical = tension * np.cos(np.deg2rad(psi))
    tension_horizontal = tension * np.sin(np.deg2rad(psi))
    weight_under_anchor = ds_sim_series_interpolated.weight_under_anchor

    fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
    ax[0].plot(time, tension, "k--", zorder=10, label="total")
    ax[0].plot(time, tension_vertical, "r-", label="vertical")
    ax[0].plot(time, tension_horizontal, "b-", label="horizontal")
    ax[1].plot(time, weight_under_anchor, "k-")

    ax[1].xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    ax[0].set(ylabel="load on\nanchor (kg)")
    ax[0].legend(frameon=False)
    ax[1].set(ylabel="weight under\nanchor (kg)", xlabel="$t$")
    return fig, ax


def plot_simulated_item_series(ds_sim_series_item: xr.Dataset):
    """Plot simulation results of an item (clampon or inline) as a function of time"""
    fig, ax = plt.subplots(5, 1, sharex=True, figsize=(10, 10))
    time = ds_sim_series_item.time
    ax[0].plot(time, ds_sim_series_item.X, "k-")
    ax[1].plot(time, ds_sim_series_item.Y, "k-")
    ax[2].plot(time, ds_sim_series_item.dZ, "k-")
    ax[3].plot(time, ds_sim_series_item.psi, "k-")
    ax[4].plot(time, ds_sim_series_item.tension_bottom, "k--", label="bottom")
    ax[4].plot(time, ds_sim_series_item.tension_top, "k-", label="top")

    ax[0].set(ylabel="$X$ (m)")
    ax[1].set(ylabel="$Y$ (m)")
    ax[2].set(ylabel="$dZ$ (m)")
    ax[3].set(ylabel="$\\psi$ ($\\degree$)")
    ax[4].set(ylabel="$T$ (kg)", xlabel="$t$")
    ax[4].legend(frameon=False)

    name = ds_sim_series_item.name.item()
    serial = ds_sim_series_item.serial.item()
    if serial is None:
        serial_str = ""
    else:
        serial_str = f" [{serial}]"
    bool_clampon = ds_sim_series_item.bool_clampon.item()
    if bool_clampon:
        height = ds_sim_series_item.height.values
        height_str = f"height={height:.1f} m"
    else:
        bottom_height = ds_sim_series_item.bottom_height.item()
        height_str = f"bottom={bottom_height:.1f} m"

    ax[0].set(title=f"{name}{serial_str}, {height_str}")
    return fig, ax
