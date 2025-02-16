"""Module to simulate mooring dynamics"""

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d

from moords.mooring_design import Mooring


def simulate_mooring_instance(
    mooring: Mooring, ds_flow_profile: xr.Dataset, print_anchor_summary: bool = False
):
    """
    Function to simulate mooring in a single instance
    and compute the difference in rest
    """
    df0, dfco0, _ = simulate_mooring_in_rest(mooring)
    df, dfco, WoB = simulate_mooring_in_flow(mooring, ds_flow_profile)

    df["dZ"] = df["Z"] - df0["Z"]
    dfco["dZ"] = dfco["Z"] - dfco0["Z"]

    anchor_dict = df.iloc[-1]
    anchor_tension = anchor_dict.tension
    anchor_vert_load = anchor_tension * np.cos(anchor_dict.psi)
    anchor_hor_load = anchor_tension * np.sin(anchor_dict.psi)
    weight_under_anchor = WoB

    # print summary
    if print_anchor_summary:
        print(f"Total Tension on Anchor = {anchor_tension:.1f} kg")
        print(f"Vertical load = {anchor_vert_load:.1f} kg")
        print(f"Horizontal load = {anchor_hor_load:.1f} kg")
        print(f"Weight under anchor = {WoB:.1f} kg  (negative is down)")

    df_sim_inline = df
    df_sim_clampon = dfco

    df_sim = pd.merge(df_sim_inline, df_sim_clampon, how="outer")
    df_sim.index.name = "idx"

    ds = xr.Dataset(df_sim)
    ds["weight_under_anchor"] = weight_under_anchor
    return ds


def simulate_mooring_series(
    mooring: Mooring, ds_flow_series: xr.Dataset, print_progress: bool = True
):
    """
    Function to simulate mooring in a time series
    and compute the difference in rest
    """
    df0, dfco0, _ = simulate_mooring_in_rest(mooring)

    time = ds_flow_series.time.values
    n_time = len(time)

    for i in np.arange(n_time):
        print(f"\rProgress: {i + 1}/{n_time}", end="")

        ds_flow_instance = ds_flow_series.isel(time=i)
        df, dfco, WoB = simulate_mooring_in_flow(mooring, ds_flow_instance)

        df["dZ"] = df["Z"] - df0["Z"]
        dfco["dZ"] = dfco["Z"] - dfco0["Z"]

        df_sim_inline = df
        df_sim_clampon = dfco
        df_sim = pd.merge(df_sim_inline, df_sim_clampon, how="outer")
        df_sim.index.name = "idx"
        ds = xr.Dataset(df_sim)
        ds["weight_under_anchor"] = WoB

        list_expand_dims = [
            "X",
            "Y",
            "Z",
            "dZ",
            "psi",
            "tension",
            "weight_under_anchor",
        ]
        ds = ds.assign_coords(time=time[i])
        ds[list_expand_dims] = ds[list_expand_dims].expand_dims("time")
        if i == 0:
            ds_merge = ds
        else:
            ds_merge = xr.merge([ds_merge, ds[list_expand_dims]])
    return ds_merge


def simulate_mooring_in_rest(mooring: Mooring):
    """Function to simulate flow in rest"""
    mooring_length = np.sum([elem.length for elem in mooring.inline])
    bottom_depth = mooring_length * 2
    n = 100
    ds_flow_profile_rest = xr.Dataset(
        {
            "U": (["z"], np.zeros(n)),
            "V": (["z"], np.zeros(n)),
            "W": (["z"], np.zeros(n)),
            "rho": (["z"], 1025 * np.ones(n)),
            "bottom_depth": bottom_depth,
        },
        coords={
            "z": np.linspace(bottom_depth, 0, n),
        },
    )
    return simulate_mooring_in_flow(mooring, ds_flow_profile_rest)


def simulate_mooring_in_flow(mooring: Mooring, ds_flow_profile: xr.Dataset):
    """
    Simulate mooring in flow instance. This is a Python conversion of
        https://web.uvic.ca/~rdewey/mooring/moordyn.php
        Mooring Design & Dynamics
        A Matlab Package for Designing and Testing
        Oceanographic Moorings And Towed Bodies
        by Richard K. Dewey,
        University of Victoria
        RDewey@UVic.ca
    Note that this version only works for submerged moorings and no surface moorings.
    Also no wind is included at the surface
    """

    # Variables flow profile
    U = ds_flow_profile.U.values
    V = ds_flow_profile.V.values
    W = ds_flow_profile.W.values
    rho = ds_flow_profile.rho.values
    z = ds_flow_profile.z.values
    Zw = ds_flow_profile.bottom_depth.item()  # water depth

    # Variables mooing - inline
    id = np.flip(mooring.make_df_inline()["id"].to_numpy())
    B = np.flip(mooring.make_df_inline()["buoyancy"].to_numpy())
    L = np.flip(mooring.make_df_inline()["length"].to_numpy())
    Width = np.flip(mooring.make_df_inline()["width"].to_numpy())
    Diam = np.flip(mooring.make_df_inline()["diameter"].to_numpy())
    Type = np.flip(
        mooring.make_df_inline()["bool_line"].to_numpy().astype(int)
    )  # 1 = line/chain
    Cd = np.flip(mooring.make_df_inline()["drag"].to_numpy())
    ME_index = np.flip(mooring.make_df_inline()["material"].to_numpy())

    # Variabes mooring - clampon
    idCO = np.flip(mooring.make_df_clampon()["id"].to_numpy())
    id_subCO = np.zeros_like(idCO, dtype=int)
    BCO = np.flip(mooring.make_df_clampon()["buoyancy"].to_numpy())
    LCO = np.flip(mooring.make_df_clampon()["length"].to_numpy())
    ZCO = np.flip(mooring.make_df_clampon()["height"].to_numpy())
    WidthCO = np.flip(mooring.make_df_clampon()["width"].to_numpy())
    DiamCO = np.flip(mooring.make_df_clampon()["diameter"].to_numpy())
    TypeCO = np.flip(
        mooring.make_df_clampon()["bool_line"].to_numpy().astype(int)
    )  # 1 = line/chain
    CdCO = np.flip(mooring.make_df_clampon()["drag"].to_numpy())

    # Inline properties
    mm = len(B)  # Initial number of elements (before interpolation)
    H = np.vstack([L, Width, Diam, Type])  # Property matrix inline elements
    B[B == 0] = -1e-4  # No neutral buoyancy to prevent errors
    S = np.sum(H[0, :])  # Mooring length
    material_properties = {
        1: 1.38e11,  # Steel
        2: 3.45e8,  # Nylon
        3: 6.9e8,  # Dacron
        4: 3.45e8,  # Polyprop
        5: 6.9e8,  # Polyethy
        6: 6.9e10,  # Kevlar
        7: 7.6e10,  # Aluminum
        8: 1.0e11,  # Dyneema
    }  # Modulus of elasticity array
    ME = np.array([material_properties.get(index, np.nan) for index in ME_index])

    if S > Zw:
        raise ValueError(
            f"Mooring length ({S:.1f}m) exceeds water depth ({Zw:.1f}m)"
            f" - surface mooring is not yet implemented"
        )
    if S > np.max(z):
        raise ValueError(
            f"Mooring length ({S:.1f}m) exceeds flow profile height ({max(z):.1f}m)"
            f" - consider adding a boundary condition"
        )
    if np.all(z > 0):
        raise ValueError("Flow profile is missing boundary condition at bottom (z = 0)")

    # Clampon properties
    mmco = len(BCO)  # Number of clampond
    flag_clampons = mmco > 0
    HCO = np.vstack([LCO, WidthCO, DiamCO, TypeCO])  # Property matrix clampon elements
    Jobj = np.zeros(mmco).astype(int)  # Initialize clampon index
    Pobj = np.zeros(mmco)  # Initialize clampon ratio along inline
    heights = np.flip(np.append(0, np.cumsum(H[0, ::-1])))
    for ico in range(mmco):
        if ZCO[ico] == S:
            i = 0
        else:
            i = np.where((ZCO[ico] < heights[:-1]) & (ZCO[ico] >= heights[1:]))[0][
                0
            ]  # Find inline index where clampon is attached
        Jobj[ico] = i  # Store index
        bottom_height = heights[i + 1]  # Inline element bottom height
        dz = ZCO[ico] - bottom_height
        Pobj[ico] = dz / H[0, i]  # Fraction of length along inline

    # Masses/buoyancies into forces
    g = 9.81
    Bw = g * B
    BwCO = g * BCO
    # Bmax = Bw[0]

    # Interpolate mooring (from top to bottom)
    idi = np.array(id[mm - 1])  # Interpolated id
    idi_sub = np.array(0)  # Submerged interpolated id
    Zi = np.array(H[0, mm - 1])  # Height of the top of the anchor, start of mooring
    Hi = np.array([H[:, mm - 1]]).T  # Setup interpolated H, B, Cd variables
    Bi = np.array(Bw[mm - 1])
    Cdi = np.array(Cd[mm - 1])
    MEi = np.array(ME[mm - 1])
    Elindx = np.zeros([mm, 2]).astype(int)  # Element indexing interpolation
    z0 = H[0, mm - 1]  # Height of top of the first element (anchor)
    j = 1  # Interpolated element index
    for i in range(mm - 2, -1, -1):
        idi_sub_count = 0
        if H[3, i] == 1:  # This section is wire/chain
            Hw = H[0, i]
            if Hw <= 0.2:
                dz
            elif 0.2 < Hw <= 5:
                dz = 0.2
            elif 5 < Hw <= 50:
                dz = 0.5
            elif 50 < Hw <= 100:
                dz = 1
            elif 100 < Hw <= 500:
                dz = 2
            else:
                dz = 5

            n = round(Hw / dz)
            dz = Hw / n
            Elindx[i, 0] = j

            for jj in range(j, j + n):
                Zi = np.append(Zi, z0 + dz / 2)
                z0 += dz
                Hjj = np.array([[dz, H[1, i], H[2, i], H[3, i]]]).T
                Hi = np.hstack([Hi, Hjj])
                Bi = np.append(Bi, Bw[i] * dz / H[0, i])
                Cdi = np.append(Cdi, Cd[i])
                MEi = np.append(MEi, ME[i])

                idi = np.append(idi, id[i])
                idi_sub = np.append(idi_sub, idi_sub_count)
                idi_sub_count += 1

            j += n - 1
            Elindx[i, 1] = j

        else:
            Elindx[i, 0:2] = [j, j]
            Zi = np.append(Zi, z0 + H[0, i] / 2)
            z0 += H[0, i]
            Hi = np.hstack([Hi, np.array([H[:, i]]).T])
            Bi = np.append(Bi, Bw[i])
            Cdi = np.append(Cdi, Cd[i])
            MEi = np.append(MEi, ME[i])

            idi = np.append(idi, id[i])
            idi_sub = np.append(idi_sub, idi_sub_count)
        j += 1

    J = len(Zi) - 1  # Initialize J

    # Find clampon locations interpolated mooring
    if flag_clampons:
        Iobj = np.zeros(mmco).astype(int)
        PIobj = np.zeros(mmco)

        # Reverse ZCO, HCO, CdCO, Pobj, and Jobj
        ZiCO = ZCO[::-1]
        HiCO = HCO[:, ::-1]
        CdiCO = CdCO[::-1]
        Piobj = Pobj[::-1]
        Jiobj = Jobj[::-1]

        for jco in range(mmco):
            # Compute interpolated element index
            range_element = Elindx[Jiobj[jco], 1] - Elindx[Jiobj[jco], 0]
            if range_element == 0:
                Iobj[jco] = Elindx[Jiobj[jco], 0]
            else:
                Iobj[jco] = int(range_element * Piobj[jco]) + Elindx[Jiobj[jco], 0]
            PIobj[jco] = (ZiCO[jco] - Zi[Iobj[jco]] + Hi[0, Iobj[jco]] / 2) / Hi[
                0, Iobj[jco]
            ]
    else:
        Iobj = np.array([])
        PIobj = np.array([])

    # Modify Elindx: effectively flip this, now top to bottom
    Elindx = J - Elindx

    # Initialize variables
    dz = 1
    dz0 = np.mean(np.abs(np.diff(z)))
    # maxz = np.sum(H[0, :])

    # Flow properties
    # num_grid = len(z)

    if dz0 < 1:  # If the velocity profile is already 1m or less
        Ui = U
        Vi = V
        Wi = W
        rhoi = rho
        zi = z
    else:
        # Adjust dz based on the direction of z
        if z[0] > z[1]:
            dz = -1
        if abs(z[-1] - z[0]) < 10:
            dz = np.sign(dz) * 0.1

        # Create the new depth grid
        zi = np.arange(z[0], z[-1] + dz, dz)

        # Interpolate U, V, W, and rho to the new depth grid
        Ui = interp1d(z, U, kind="linear", fill_value="extrapolate")(zi)
        Vi = interp1d(z, V, kind="linear", fill_value="extrapolate")(zi)
        Wi = interp1d(z, W, kind="linear", fill_value="extrapolate")(zi)
        rhoi = interp1d(z, rho, kind="linear", fill_value="extrapolate")(zi)

    # Compute the total current speed
    Umag = np.sqrt(Ui**2 + Vi**2 + Wi**2)

    # Calculate the new number of "interpolated" in-line mooring elements
    N = len(Bi)

    # Initialize drag forces
    Qx = np.zeros(N)
    Qy = np.zeros(N)
    Qz = np.zeros(N)

    # Compute drag from bottom-to-top
    for j in range(N):
        if flag_clampons:  # If there are clamp-on devices
            ico = np.where(Iobj == j)[0]  # Find the indices of clamp-on devices
        else:
            ico = []

        # Find the index in the depth grid corresponding to this segment
        i = np.where((zi >= (Zi[j] - 0.5)) & (zi <= (Zi[j] + 0.5)))[0]
        i = i[0]  # Take the first match

        # Compute the exposed area
        if Hi[2, j] == 0:  # Cylinder/wire/chain section
            A = Hi[0, j] * Hi[1, j]  # Exposed area of cylinder
        else:
            A = np.pi * (Hi[2, j] / 2) ** 2  # Exposed area of sphere

        # Drag in X and Y directions
        Qx[j] = 0.5 * rhoi[i] * Cdi[j] * A * Umag[i] * Ui[i]
        Qy[j] = 0.5 * rhoi[i] * Cdi[j] * A * Umag[i] * Vi[i]

        # Handle clamp-on devices
        Qxco = 0
        Qyco = 0
        if len(ico) > 0:
            for icoc in ico:
                if HiCO[2, icoc] == 0:  # Cylinder device
                    Axco = HiCO[0, icoc] * HiCO[1, icoc]
                    Ayco = Axco
                    Cdjxco = CdiCO[icoc]
                    Cdjyco = CdiCO[icoc]
                else:  # Sphere device
                    Axco = np.pi * (HiCO[2, icoc] / 2) ** 2
                    Ayco = Axco
                    Cdjxco = CdiCO[icoc]
                    Cdjyco = Cdjxco
                Qxco += 0.5 * rhoi[i] * Cdjxco * Axco * Umag[i] * Ui[i]
                Qyco += 0.5 * rhoi[i] * Cdjyco * Ayco * Umag[i] * Vi[i]

        Qx[j] += Qxco
        Qy[j] += Qyco

        # Drag in Z direction
        if Hi[2, j] == 0:  # Cylinder/wire/chain section
            A = np.pi * (Hi[1, j] / 2) ** 2  # Area of bottom of cylinder
            if Hi[3, j] == 1:  # If wire, vertical area = 0
                A = 0

        Qz[j] = 0.5 * rhoi[i] * Cdi[j] * A * Umag[i] * Wi[i]

        # Handle vertical drag for clamp-on devices
        Qzco = 0
        if len(ico) > 0:
            for icoc in ico:
                if HiCO[2, icoc] == 0:  # Cylinder section
                    Azco = np.pi * (HiCO[1, icoc] / 2) ** 2
                    Cdjzco = CdiCO[icoc]
                else:  # Sphere
                    Azco = np.pi * (HiCO[2, icoc] / 2) ** 2
                    Cdjzco = CdiCO[icoc]
                Qzco += 0.5 * rhoi[i] * Cdjzco * Azco * Umag[i] * Wi[i]

        Qz[j] += Qzco

    # Flip mooring right side up (indices now start at top)
    idi = np.flip(idi)
    idi_sub = np.flip(idi_sub)
    Qx = np.flip(Qx)
    Qy = np.flip(Qy)
    Qz = np.flip(Qz)
    Hi = np.flip(Hi, axis=1)
    Bi = np.flip(Bi)
    Cdi = np.flip(Cdi)
    MEi = np.flip(MEi)
    Iobj = (N - 1) - np.flip(Iobj)
    PIobj = 1 - np.flip(PIobj)

    # First Pass
    # Now we solve for first order wire angles, starting at the top of the mooring,
    #  where there is no tension from above, Ti=0.
    # Then there are three equations and three unknowns at each element.
    # 1) Qx(i) + T(i)*cos(theta(i))*sin(psi(i)) = T(i+1)*cos(theta(i+1))*sin(psi(i+1))
    # 2) Qy(i) + T(i)*sin(theta(i))*sin(psi(i)) = T(i+1)*sin(theta(i+1))*sin(psi(i+1))
    # 3) Wz(i) + Qz(i) + T(i)*cos(psi(i)) = T(i+1)*cos(psi(i+1))
    # where the Q's are the drags, Wz is the weight/buoyancy,
    #    T(i) is the tension from above, T(i+1) is the tension from below,
    #    psi(i) is the wire angle from z, theta(i) angle in x-y plane
    #    Calculate T(i+1), phi(i+1) and theta(i+1) at element i,
    #    working from the top T(1)=0, to the bottom.
    #  All cylinders have a tangential drag coefficient of 0.01
    #  Here it is assumed that the top of the mooring is a float, or
    #  at least has positive buoyancy and will "lift" the first few elements.

    # Initialize variables
    Ti = np.zeros(N)
    theta = np.zeros(N)
    psi = np.zeros(N)

    # First element initialization
    Ti[0] = 0
    theta[0] = 0
    psi[0] = 0

    # Second element initialization
    gamma = 1
    gamma_float = 1
    b = gamma * (Bi[0] + Qz[0])
    if Iobj.size > 0:
        ico = np.where(Iobj == 0)[0]
    if len(ico) > 0:
        b += np.sum(BwCO[ico])  # Add buoyancy of clamp-on devices to top
    theta[1] = np.arctan2(Qy[0], Qx[0])
    Ti[1] = np.sqrt(Qx[0] ** 2 + Qy[0] ** 2 + b**2)
    psi[1] = np.real(np.arccos(b / Ti[1]))

    # Solve from top to bottom
    for i in range(1, N - 1):
        ico = []
        if Iobj.size > 0:
            ico = np.where(Iobj == i)[0]

        ip1 = i + 1
        xx = Qx[i] + Ti[i] * np.cos(theta[i]) * np.sin(psi[i])  # Force in x direction
        yy = Qy[i] + Ti[i] * np.sin(theta[i]) * np.sin(psi[i])  # Force in y direction
        zz = Bi[i] + Qz[i] + Ti[i] * np.cos(psi[i])  # Vertical force in z direction

        if len(ico) > 0:
            zz += np.sum(BwCO[ico])  # Add buoyancy of clamp-on devices

        theta[ip1] = np.arctan2(yy, xx)
        Ti[ip1] = np.sqrt(xx**2 + yy**2 + zz**2)
        if Ti[ip1] != 0:
            psi[ip1] = np.real(np.arccos(zz / Ti[ip1]))
        else:
            psi[ip1] = psi[i]

    # Integrate from bottom to top to compute first-order positions [x, y, z]
    X = np.zeros(N)
    Y = np.zeros(N)
    Z = np.zeros(N)

    # Now integrate from the bottom to the top to get the first order [x,y,z]
    # Allow wire/rope sections to stretch under tension

    for _ in range(2):  # stretch in 2 steps
        X[N - 1] = 0
        Y[N - 1] = 0
        Z[N - 1] = Hi[0, N - 1]  # Reference to top of anchor
        dx0 = dy0 = dz0 = 0

        for i in range(N - 2, -1, -1):
            if (Hi[1, i] != 0) and (Hi[3, i]):
                dL = 1 + (Ti[i] * 4 / (np.pi * Hi[1, i] ** 2 * MEi[i]))
            else:
                dL = 1

            LpdL = Hi[0, i] * dL
            X[i] = X[i + 1] + LpdL * np.cos(theta[i]) * np.sin(psi[i]) / 2 + dx0
            Y[i] = Y[i + 1] + LpdL * np.sin(theta[i]) * np.sin(psi[i]) / 2 + dy0
            Z[i] = Z[i + 1] + LpdL * np.cos(psi[i]) / 2 + dz0

            dx0 = LpdL * np.cos(theta[i]) * np.sin(psi[i]) / 2
            dy0 = LpdL * np.sin(theta[i]) * np.sin(psi[i]) / 2
            dz0 = LpdL * np.cos(psi[i]) / 2

            # Adjust for surface or bottom constraints
            if Z[i] > Zw and Hi[3, i] == 1 and Bi[i] > 0:  # Surface line
                Z[i] = Zw
                psi[i] = np.pi / 2
            if Z[i] <= Z[N - 1]:  # Bottom chain
                Z[i] = Z[N - 1]
                psi[i] = np.pi / 2

    # Now with the first order positions, we must re-estimate the new
    # drags at the new heights (Zi) and for cylinders tilted by psi in flow.

    # Re-estimate drags at new heights (Zi) and adjust for surface float moorings
    breaknow = False
    # iconv = 0
    # icnt = 0
    iavg = 0
    isave = -1
    dg = 0.1
    # gf = 2
    # dgf = 0
    # dgc = 0

    # Solver
    # iprt = 100  #  If solution isn't converging, set this to 50-100 and watch to see
    # what's happening
    Zsave = []  # Allocate list for top float z-position
    deltaz = 0.01  # Higher precision for convergence (1 cm)

    # Initialization for the iterative loop
    # ilines = 1
    ico = []
    # iiprt = 0
    # dgci = 10

    # Main iteration/convergence loop

    while not breaknow:  # Loop until convergence
        # Update the iteration counter
        isave += 1

        if (
            isave >= 20
        ):  # If having problems converging, start a running average when close
            iavg += 1
            if iavg == 1:
                Tiavg = Ti.copy()
                psiavg = psi.copy()
                Zavg = Z.copy()
                Z1 = [Z[0]]
                Xavg = X.copy()
                Yavg = Y.copy()
                # Uio = Ui.copy()
            else:
                Tiavg += Ti
                psiavg += psi
                Zavg += Z
                Z1.append(Z[0])
                Xavg += X
                Yavg += Y
                # Z1std = np.std(Z1)

        if (
            iavg > 20
        ):  # After 20 averaging iterations, use averages to assist convergence
            X = Xavg / iavg
            Y = Yavg / iavg
            Z = Zavg / iavg
            Ti = Tiavg / iavg
            psi = psiavg / iavg

        # if iavg > iprt:
        #     print([Z[0], (Z[0] - Z1[isave - 1])])

        # Calculate angles and current magnitude
        # phix = np.arctan2((np.cos(theta) * np.sin(psi)), np.cos(psi))
        # phiy = np.arctan2((np.sin(theta) * np.sin(psi)), np.cos(psi))
        Umag = np.sqrt(
            Ui**2 + Vi**2 + Wi**2
        )  # Current magnitude for interpolated profile

        # Reset all drag forces to zero
        Qx = np.zeros(N)
        Qy = np.zeros(N)
        Qz = np.zeros(N)

        for j in range(
            N
        ):  # Loop through the interpolated in-line segments and any clamp-on devices
            if flag_clampons:
                ico = [ic for ic, val in enumerate(Iobj) if val == j]
            else:
                ico = []

            if j == N - 1:  # at anchor
                psi_j = 0  # anchor always stays up right
            else:
                psi_j = psi[j]

            # Find indices of velocity values near the current segment
            i = [k for k in range(len(zi)) if (Z[j] - 1.0) <= zi[k] <= (Z[j] + 1.0)]

            if j == 0:
                i = [
                    k
                    for k in range(len(zi))
                    if (Z[j] - Hi[0, 0]) <= zi[k] <= (Z[j] + Hi[0, 0])
                ]
                if len(i) == 0:
                    i = [0]  # Take the top velocity value if none found

            if len(i) == 0:
                print(f"Check this configuration: {j} {Z[0]} {Z[j]}")
                raise ValueError("Can't find the velocity at this element!")

            i = i[0]  # Use the first index in case multiple values found

            # Calculate angles and velocity magnitudes
            theta2 = np.arctan2(Vi[i], Ui[i])
            UVLmag = np.sqrt(Ui[i] ** 2 + Vi[i] ** 2) * np.cos(theta[j] - theta2)
            UL = UVLmag * np.cos(theta[j])
            VL = UVLmag * np.sin(theta[j])
            Up = Ui[i] - UL
            Vp = Vi[i] - VL
            theta3 = np.arctan2(VL, UL)
            thetap = np.arctan2(Vp, Up)

            if Hi[2, j] == 0:  # Cylinder/wire/chain section
                A = Hi[0, j] * Hi[1, j]
                Cdjxy = Cdi[j]
                Qh = 0.5 * rhoi[i] * Cdjxy * A * (Up**2 + Vp**2)
                Qx[j] = Qh * np.cos(thetap)
                Qy[j] = Qh * np.sin(thetap)
            else:  # Sphere section
                A = np.pi * (Hi[2, j] / 2) ** 2
                Cdj = Cdi[j]
                Qh = 0.5 * rhoi[i] * Cdj * A * (Ui[i] ** 2 + Vi[i] ** 2)
                Qx[j] = Qh * np.cos(theta2)
                Qy[j] = Qh * np.sin(theta2)
                Qz[j] = 0.5 * rhoi[i] * Cdj * A * abs(Wi[i]) * Wi[i]

            Qxco = 0
            Qyco = 0
            Qzco = 0

            if ico:  # Loop through clamp-on devices
                for icoc in ico:
                    if HCO[2, icoc] == 0:  # Cylinder device
                        A = HCO[0, icoc] * HCO[1, icoc]
                        Cdjco = CdCO[icoc] + HCO[1, icoc] * np.pi * 0.01 * (
                            1 - ((np.pi / 2) - psi_j) / (np.pi / 2)
                        )
                        Qhco = 0.5 * rhoi[i] * Cdjco * A * (Up**2 + Vp**2)
                        Qxco += Qhco * np.cos(thetap)
                        Qyco += Qhco * np.sin(thetap)
                    else:  # Sphere device
                        A = np.pi * (HCO[2, icoc] / 2) ** 2
                        Cdjco = CdCO[icoc]
                        Qh = 0.5 * rhoi[i] * Cdjco * A * (Ui[i] ** 2 + Vi[i] ** 2)
                        Qxco += Qh * np.cos(theta2)
                        Qyco += Qh * np.sin(theta2)
                        Qzco += 0.5 * rhoi[i] * Cdjco * A * abs(Wi[i]) * Wi[i]

            Qx[j] += Qxco
            Qy[j] += Qyco
            Qz[j] += Qzco

            # The next section represents the lift/normal drag due to tilted
            # cylinder/wire segments, associate only with the flow in the theta plane
            # Hoerner (1965): drag coeficients are reduced by cos(psi)^3=sin(pi/2-psi)^3
            # Chapter III, equations (22) & (23)
            # There are two additional horizontal forces associated with the component
            # of u & v along theta and a lift component coming from w.
            # There are two additional vertical components associated with w along theta
            # and a lift component from u & v along theta.
            # see Dewey's mooreleang.m rountine to plot these components.

            # Additional forces for tilted cylinders
            psi2 = psi_j - np.pi / 2
            if Hi[2, j] == 0:  # Tilted cylinder device/wire segment
                A = Hi[0, j] * Hi[1, j]
                CdUV = Cdi[j] * np.cos(psi_j) ** 3 + Hi[1, j] * np.pi * 0.01 * (
                    1 - ((np.pi / 2) - psi_j) / (np.pi / 2)
                )
                CdW = Cdi[j] * np.cos(psi2) ** 3 + Hi[1, j] * np.pi * 0.01 * (
                    1 - ((np.pi / 2) - psi2) / (np.pi / 2)
                )
                sl = np.sign(np.sin(theta[j])) * np.sign(np.sin(theta3))
                sl = 1 if sl == 0 else sl
                CdLUV = -Cdi[j] * np.cos(psi_j) ** 2 * np.sin(psi_j)
                CdLW = Cdi[j] * np.cos(psi2) ** 2 * np.sin(psi2)

                QhUV = 0.5 * rhoi[i] * A * CdUV * UVLmag**2
                Qx[j] += QhUV * np.cos(theta3)
                Qy[j] += QhUV * np.sin(theta3)

                QhLW = 0.5 * rhoi[i] * A * CdLW * abs(Wi[i]) * Wi[i]
                Qx[j] += QhLW * np.cos(theta[j])
                Qy[j] += QhLW * np.sin(theta[j])

                Qz[j] += 0.5 * rhoi[i] * A * CdLUV * UVLmag**2 * sl
                Qz[j] += 0.5 * rhoi[i] * A * CdW * abs(Wi[i]) * Wi[i]

            if ico:  # Lift terms for tilted clamp-on devices
                Qxco = 0
                Qyco = 0
                Qzco = 0
                for icoc in ico:
                    if HCO[2, icoc] == 0:  # Tilted clamp-on cylinder
                        Aco = HCO[0, icoc] * HCO[1, icoc]
                        Aeco = np.pi * (HCO[0, icoc] / 2) ** 2
                        CdUV = CdCO[icoc] * np.cos(psi_j) ** 3 + HCO[
                            1, icoc
                        ] * np.pi * 0.01 * (1 - ((np.pi / 2) - psi_j) / (np.pi / 2))
                        CdW = CdCO[icoc] * np.cos(psi2) ** 3 + HCO[
                            1, icoc
                        ] * np.pi * 0.01 * (1 - ((np.pi / 2) - psi2) / (np.pi / 2))
                        sl = np.sign(np.sin(theta[j])) * np.sign(np.sin(theta3))
                        sl = 1 if sl == 0 else sl
                        CdLUV = -CdCO[icoc] * np.cos(psi_j) ** 2 * np.sin(psi_j)
                        CdLW = CdCO[icoc] * np.cos(psi2) ** 2 * np.sin(psi2)

                        QhUV = 0.5 * rhoi[i] * CdUV * Aco * UVLmag**2
                        Qxco += QhUV * np.cos(theta3)
                        Qyco += QhUV * np.sin(theta3)
                        Qzco += 0.5 * rhoi[i] * CdLUV * Aco * UVLmag**2 * sl

                        QhLW = 0.5 * rhoi[i] * CdLW * Aco * abs(Wi[i]) * Wi[i]
                        Qxco += QhLW * np.cos(theta[j])
                        Qyco += QhLW * np.sin(theta[j])
                        Qzco += 0.5 * rhoi[i] * CdW * Aco * abs(Wi[i]) * Wi[i]

                        Qhe = (
                            0.5 * rhoi[i] * 0.65 * abs(np.sin(psi_j)) * Aeco * UVLmag**2
                        )
                        Qxco += Qhe * np.cos(theta3)
                        Qyco += Qhe * np.sin(theta3)
                        Qzco += (
                            0.5
                            * rhoi[i]
                            * 0.65
                            * abs(np.cos(psi_j))
                            * Aeco
                            * abs(Wi[i])
                            * Wi[i]
                        )

                Qx[j] += Qxco
                Qy[j] += Qyco
                Qz[j] += Qzco

        # Now re-solve for displacements with new positions/drags

        # Initialize arrays
        Ti = np.zeros(N)
        thetaNew = np.zeros(N)
        psiNew = np.zeros(N)

        # No tension above top element
        Ti[0] = 0

        # Top element (float)
        b = Bi[0] + Qz[0]
        if Iobj.size > 0:
            ico = np.where(Iobj == 0)[0]
        if len(ico) > 0:
            b += np.sum(BwCO[ico])  # Add buoyancy of clamp-on devices to top
        thetaNew[1] = np.arctan2(Qy[0], Qx[0])
        Ti[1] = np.sqrt(Qx[0] ** 2 + Qy[0] ** 2 + b**2)

        T_h = np.sqrt((gamma_float * Qx[0]) ** 2 + (gamma_float * Qy[0]) ** 2)
        psiNew[1] = np.real(np.arcsin(T_h / Ti[1]))
        psiNew[0] = 0  # Top float is inclined upwards
        thetaNew[0] = thetaNew[1]  # Top float is tilted in X/Y as a solid object

        # Solve from top to bottom
        for i in range(1, N - 1):
            ico = np.where(Iobj == i)[0] if len(Iobj) > 0 else []
            ip1 = i + 1

            xx = Qx[i] + Ti[i] * np.cos(thetaNew[i]) * np.sin(psiNew[i])
            yy = Qy[i] + Ti[i] * np.sin(thetaNew[i]) * np.sin(psiNew[i])
            zz = Bi[i] + Qz[i] + Ti[i] * np.cos(psiNew[i])

            if len(ico) > 0:
                zz += np.sum(BwCO[ico])

            thetaNew[ip1] = np.arctan2(yy, xx)
            Ti[ip1] = np.sqrt(xx**2 + yy**2 + zz**2)

            if Ti[ip1] != 0:
                psiNew[ip1] = np.real(np.arccos(zz / Ti[ip1]))
            else:
                psiNew[ip1] = psiNew[i]

        thetaNew = np.real(thetaNew)
        psiNew = np.real(psiNew)

        # Now integrate/sum positions from the bottom to the
        # top to get the second order [x,y,z]
        # Allow wire/rope to stretch under tension
        X = np.zeros(N)
        Y = np.zeros(N)
        Z = np.zeros(N)

        X[N - 1] = 0
        Y[N - 1] = 0
        Z[N - 1] = Hi[0, N - 1]

        Zii = True
        iint = 0

        while Zii:
            Zii = False
            dx0 = dy0 = dz0 = 0
            iint += 1

            for i in range(N - 2, -1, -1):
                if Hi[1, i] != 0 and MEi[i] < np.inf and Hi[3, i]:
                    dL = 1 + (Ti[i] * 4 / (np.pi * Hi[1, i] ** 2 * MEi[i]))
                else:
                    dL = 1
                LpdL = Hi[0, i] * dL

                dX = LpdL * np.cos(thetaNew[i]) * np.sin(psiNew[i])
                dY = LpdL * np.sin(thetaNew[i]) * np.sin(psiNew[i])
                dZ = LpdL * np.cos(psiNew[i])

                X[i] = X[i + 1] + dX / 2 + dx0 / 2
                Y[i] = Y[i + 1] + dY / 2 + dy0 / 2
                Z[i] = Z[i + 1] + dZ / 2 + dz0 / 2

                if (
                    Z[i] > Zw and Hi[3, i] == 1 and Bi[i] >= 0
                ):  # force line to lie on surface
                    Zii = True
                    Z[i] = Zw
                    psi[i] = np.pi / 2

                if Z[i] < 0:  # force line to lie on anchor/bottom
                    Zii = True
                    Z[i] = 0
                    psi[i] = np.pi / 2

                dx0, dy0, dz0 = dX, dY, dZ

            if iint > 4:
                Zii = False

        psi[N - 1] = psi[N - 2]

        Z = np.real(Z)
        X = np.real(X)
        Y = np.real(Y)
        theta = np.real(theta)
        psi = np.real(psi)

        # Scale the correction to get faster convergence
        scale = 0.5
        psi = psi + (psiNew - psi) * scale
        theta = theta + (thetaNew - theta) * scale

        # Recompute the positions
        dx0 = dy0 = dz0 = 0

        for i in range(N - 2, -1, -1):  # i decreases from N-1 (anchor) to 0 (float)
            if Hi[1, i] != 0 and MEi[i] < np.inf and Hi[1, i]:
                dL = 1 + (
                    Ti[i] * 4 / (np.pi * Hi[1, i] ** 2 * MEi[i])
                )  # Stretching of wire/rope
            else:
                dL = 1

            LpdL = Hi[0, i] * dL  # Length of this element

            dX = (
                LpdL * np.cos(theta[i]) * np.sin(psi[i])
            )  # Length in x direction (spherical coordinates)
            dY = LpdL * np.sin(theta[i]) * np.sin(psi[i])  # Length in y direction
            dZ = LpdL * np.cos(psi[i])  # Length in z direction

            X[i] = X[i + 1] + dX / 2 + dx0 / 2
            Y[i] = Y[i + 1] + dY / 2 + dy0 / 2
            Z[i] = Z[i + 1] + dZ / 2 + dz0 / 2

            if (
                Z[i] > Zw and Hi[3, i] == 1 and Bi[i] >= 0
            ):  # Force line to lie on surface
                Z[i] = Zw
                psi[i] = np.pi / 2

            if Z[i] < 0:  # Force line to lie on anchor/bottom
                Z[i] = 0
                psi[i] = np.pi / 2

            dx0, dy0, dz0 = dX, dY, dZ

        Z = np.real(Z)
        X = np.real(X)
        Y = np.real(Y)
        psi = np.real(psi)

        # Calculate final float position
        Zf = Z[0] - Hi[0, 0] / 2
        if np.max(Z) > Zw:
            raise ("Water depth < mooring length")

        # Convergence checks
        if isave > 1:  # Must do at least three iterations
            if (
                abs(Zsave[isave - 1] - Z[0]) < deltaz
                and abs(Zsave[isave - 2] - Zsave[isave - 1]) < deltaz
            ):  # Check close calls
                if Zw > (Zf + Hi[0, 0]):  # Sub-surface solution converged
                    breaknow = True

            if iavg == 120 or (
                iavg > 100 and dg < 1e-10
            ):  # Force convergence after many iterations
                X = Xavg / iavg
                Y = Yavg / iavg
                Z = Zavg / iavg
                Ti = Tiavg / iavg
                psi = psiavg / iavg
                breaknow = True
                # iconv = 1

        Zsave.append(Z[0])

        # Slacken the threshold if not converging
        if isave % 100 == 0:
            deltaz *= 2

    Xfco = np.empty(mmco)
    Yfco = np.empty(mmco)
    Zfco = np.empty(mmco)
    psifco = np.empty(mmco)

    # If there are clamp-on devices, figure out their final positions
    if flag_clampons:  # Check if there are any clamp-on devices
        # Loop through each clamp-on device
        for jco in range(mmco):
            idx_obj = Iobj[jco]
            pidx_obj = PIobj[jco]

            if idx_obj == N - 1:  # attached to anchor
                Xfco[jco] = 0
                Yfco[jco] = 0
                Zfco[jco] = (1.5 - pidx_obj) * Hi[0, idx_obj]
                psifco[jco] = 0

            else:
                Xfco[jco] = (
                    X[idx_obj]
                    + np.cos(theta[idx_obj])
                    * np.sin(psi[idx_obj])
                    * (0.5 - pidx_obj)
                    * Hi[0, idx_obj]
                )
                Yfco[jco] = (
                    Y[idx_obj]
                    + np.sin(theta[idx_obj])
                    * np.sin(psi[idx_obj])
                    * (0.5 - pidx_obj)
                    * Hi[0, idx_obj]
                )
                Zfco[jco] = (
                    Z[idx_obj]
                    + np.cos(psi[idx_obj]) * (0.5 - pidx_obj) * Hi[0, idx_obj]
                )
                psifco[jco] = psi[Iobj[jco] + 1]

    # Calculations related to weight and forces
    # ba = psi[N - 1]
    # Wa = Ti[N - 1] / 9.81
    # VWa = Wa * np.cos(ba)
    # HWa = Wa * np.sin(ba)
    WoB = (Bi[N - 1] + Qz[N - 1] + Ti[N - 1]) / 9.81  # weight under anchor
    if Iobj.size > 0:
        ico = np.where(Iobj == N - 1)[0]
    if len(ico) > 0:
        WoB += np.sum(BwCO[ico]) / 9.81  # Add buoyancy of clamp-on devices to top
        # print(f'{np.sum(BwCO[ico])/9.81}\n{ico}')

    if WoB >= 0:
        raise ValueError(
            f"Simulation stopped - Anchor too light: "
            f"weight under anchor = {WoB:.1f} [kg] (negative is down)"
        )

    # # Display the calculated forces and weights
    # print(f"Total Tension on Anchor [kg] = {Wa:.1f}")
    # print(f"Vertical load [kg] = {VWa:.1f}  Horizontal load [kg] = {HWa:.1f}")

    # # Calculate and display safe wet anchor mass
    # TWa = 1.5 * (VWa + HWa / 0.6)
    # print(f"Safe wet anchor mass = {TWa:.1f} [kg] = {TWa * 2.2:.1f} [lb]")
    # print(
    #     f"Safe dry steel anchor mass = {(TWa / 0.87):.1f} [kg]"
    #     f" = {(TWa * 2.2 / 0.87):.1f} [lb]"
    # )
    # print(
    #     f"Safe dry concrete anchor mass = {(TWa / 0.65):.1f} [kg] "
    #     f"= {(TWa * 2.2 / 0.65):.1f} [lb]"
    # )
    # print(f"Weight under anchor = {WoB:.1f} [kg]  (negative is down)")

    # Check if the anchor is too light
    # if abs(B[-1]) < TWa:
    #     print("*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*")
    #     print("*!*!*!*  Warning. Anchor is likely TOO light!   *!*!*!*")
    #     print("*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*")

    numi = np.empty(len(Zi), dtype=object)
    for i, idx_range in enumerate(Elindx):
        for j in range(idx_range[1], idx_range[0] + 1):
            numi[j] = mm - i

    psi = np.hstack([psi[1:], 0])

    df = pd.DataFrame(
        {
            "id": idi,
            "id_sub": idi_sub,
            "bool_clampon": np.zeros_like(idi, bool),
            "X": X,
            "Y": Y,
            "Z": Z,
            "psi": np.rad2deg(psi),
            "tension": Ti / 9.81,
        }
    )

    dfco = pd.DataFrame(
        {
            "id": idCO,
            "id_sub": id_subCO,
            "bool_clampon": np.ones_like(idCO, bool),
            "X": Xfco,
            "Y": Yfco,
            "Z": Zfco,
            "psi": np.rad2deg(psifco),
        }
    )

    return df, dfco, WoB
