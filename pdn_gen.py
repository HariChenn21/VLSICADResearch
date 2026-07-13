#!/usr/bin/env python3
"""
PDN Mesh Generator for LTSpice
Generates a parametric two-layer backside PDN mesh netlist for any m×n array.

Topology (per cell at intersection i,j):
  - n_i_j_h : upper layer node (BSM3, horizontal laterals)
  - n_i_j   : lower layer node (BSM1, vertical laterals + current source)
  - Via Rvia_cell connects the two

Pads (voltage supply points):
  - Configurable via pad_spacing + pad_start, or explicit pad_locations
  - Each pad connects Vout -> upper-layer node through Rsup
  - Rsup scales with pad count so parallel R stays at tile-level

Voltage excitation:
  - Vsupply Vout 0 <vout_dc> is auto-generated (no manual edits)
  - .op analysis and .end included so the .cir runs standalone

Auto-calibration:
  - MNA solve computes exact Iload for target worst-case droop.

Usage (CLI):
    python pdn_gen.py m n [-o FILE]
                          [--pad-spacing S | --pad-spacing SX SY]
                          [--pad-start I J]
                          [--target-droop MV]
                          [--vout V]
Examples:
    python pdn_gen.py 20 20                                  # 4 corners default
    python pdn_gen.py 20 20 --pad-spacing 3 --pad-start 2 2  # ~7x7 pads
    python pdn_gen.py 20 20 --pad-spacing 3 --vout 1.0 --target-droop 100
    python pdn_gen.py 32 32 --pad-spacing 4 -o big.cir
"""

import sys
import numpy as np


# ---- Pad generation ----
def pad_grid(m, n, pad_spacing, start=(1, 1)):
    """Generate (i, j) pad locations on a regular grid.

    pad_spacing: int (same both directions) or (sx, sy) tuple
    start:       (i, j) first pad position, 1-indexed
    Returns:     list of (i, j) tuples clipped to the m×n mesh
    """
    if isinstance(pad_spacing, (tuple, list)):
        sx, sy = pad_spacing
    else:
        sx = sy = pad_spacing
    si, sj = start
    pads = []
    i = si
    while i <= m:
        j = sj
        while j <= n:
            pads.append((i, j))
            j += sy
        i += sx
    return pads


# ---- Numeric mirror of the .param resistance formulas (for calibration) ----
def _mesh_resistances():
    """Evaluate the per-layer resistances numerically. Must stay in sync with
    the .param block emitted below. Returns (Rlat, Rtap, Rsup_tile) in ohms.

    Rsup_tile is the *tile-level* (parallel-equivalent) supply-side resistance,
    matching the original 4-corner '*4' formula. Per-pad Rsup is scaled by
    Npads at generation time so the parallel combination stays constant."""
    tile = 400e-6; PnTSV = 4e-6
    Pb = 57e-6; alpha = 3; Pp = 150e-6
    hVia = 10e-6; Dvia = 50e-6
    PI = 3.14159

    RtapOne = 0.5 * (2e-6 / PnTSV) ** 2
    Ntap = (tile / (2 * PnTSV)) ** 2
    Rtap = RtapOne / Ntap

    Nbump = (tile / Pb) ** 2
    Rbump = 10e-3 / Nbump
    RviaOne = 1.68e-8 * hVia / (PI * (Dvia / 2) ** 2)
    Rvia = RviaOne / Nbump

    wBus = Pb * (alpha + 1) / alpha
    Nbus = tile / Pp
    Rrdl3 = 2e-3 * tile / wBus / 3 / Nbus

    Npillar = (tile / Pp) ** 2
    Rrdl2 = 2e-3 * 0.25 / Npillar
    Rrdl1 = 2e-3 * 0.25 / Npillar
    RpillarOne = 1.68e-8 * 120e-6 / (PI * (75e-6 / 2) ** 2)
    Rpillar = RpillarOne / Npillar

    Rlat = 40e-3
    # Tile-level (matches original *4 corner assumption)
    Rsup_tile = (Rpillar + Rrdl1 + 2 * Rvia + Rrdl2 + Rrdl3 + Rbump) * 4
    return Rlat, Rtap, Rsup_tile


def solve_droop(m, n, Iload, pads):
    """Nodal (MNA) solve of the linear mesh. Returns droop (Vout - V) at every
    lower/current-source node as an m x n array, in volts.

    pads: list of (i, j) supply positions on the upper layer."""
    Rlat, Rtap, Rsup_tile = _mesh_resistances()
    Ncells = m * n
    Rvia_cell = Rtap * Ncells
    Npads = len(pads)
    # Scale so parallel combination of Npads supplies equals tile-level R.
    # Reduces to the original *4 formula when Npads == 4.
    Rsup = (Rsup_tile / 4) * Npads
    Icell = Iload / Ncells

    def lo(i, j): return (i - 1) * n + (j - 1)
    def hi(i, j): return Ncells + (i - 1) * n + (j - 1)

    Nn = 2 * Ncells
    G = np.zeros((Nn, Nn))
    Ivec = np.zeros(Nn)

    def add(a, b, g):  # conductance a<->b; b<0 => Vout reference (0 V)
        G[a, a] += g
        if b >= 0:
            G[b, b] += g
            G[a, b] -= g
            G[b, a] -= g

    glat, gvia, gsup = 1.0 / Rlat, 1.0 / Rvia_cell, 1.0 / Rsup
    for i in range(1, m + 1):
        for j in range(1, n):
            add(hi(i, j), hi(i, j + 1), glat)          # horizontal laterals
    for i in range(1, m):
        for j in range(1, n + 1):
            add(lo(i, j), lo(i + 1, j), glat)          # vertical laterals
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            add(hi(i, j), lo(i, j), gvia)              # vias
            Ivec[lo(i, j)] -= Icell                    # current out of load node
    for (i, j) in pads:
        add(hi(i, j), -1, gsup)                        # supply stacks -> Vout

    v = np.linalg.solve(G, Ivec)
    return -v[:Ncells].reshape(m, n)


def calibrate_iload(m, n, target_droop_mv, pads):
    """Total Iload placing the worst-case node droop at target_droop_mv.
    Droop is linear in Iload, so one solve + one scale is exact."""
    droop_per_amp = solve_droop(m, n, Iload=1.0, pads=pads)
    Iload = (target_droop_mv * 1e-3) / droop_per_amp.max()
    return Iload, droop_per_amp * Iload


def generate_pdn(m, n, output_file,
                 pad_spacing=None, pad_start=(1, 1), pad_locations=None,
                 target_droop_mv=75, vout_dc=0.6):
    """Generate the .cir file for an m×n mesh with configurable pad grid + Vsupply.

    Args:
        m, n:            mesh rows / columns
        output_file:     path to output .cir file
        pad_spacing:     int or (sx, sy) tuple. Cells between adjacent pads.
                         If None (and pad_locations is None), defaults to 4 corners.
        pad_start:       (i, j) first-pad position on the grid (1-indexed).
        pad_locations:   explicit list of (i, j) pad positions.
                         Overrides pad_spacing/pad_start if given.
        target_droop_mv: worst-case droop target (auto-calibrates Iload).
        vout_dc:         DC voltage forced onto Vout by Vsupply (V).
    """
    if m < 2 or n < 2:
        sys.exit("Error: m and n must be >= 2")

    # ---- Resolve pad list: explicit > grid > default corners ----
    if pad_locations is not None:
        pads = [tuple(p) for p in pad_locations]
    elif pad_spacing is not None:
        pads = pad_grid(m, n, pad_spacing, pad_start)
    else:
        pads = [(1, 1), (1, n), (m, 1), (m, n)]

    if not pads:
        sys.exit("Error: pad list is empty; no supplies would be generated.")

    N_cells = m * n
    N_pads = len(pads)
    Iload, droop = calibrate_iload(m, n, target_droop_mv, pads)

    out = []
    out.append(f"* {m}x{n} Two-Layer Backside PDN Model")
    out.append(f"* Total cells: {N_cells}   Total pads: {N_pads}")
    preview = pads[:8]
    pad_str = ", ".join(f"({i},{j})" for i, j in preview)
    if N_pads > 8:
        pad_str += ", ..."
    out.append(f"* Pad locations (first 8): {pad_str}")
    out.append(f"* Interface: 'Vout' (driven by Vsupply at {vout_dc} V)")
    out.append("")

    # ---- Design knobs ----
    out.append("*=========== Design Parameters ===========")
    # Auto-calibrate Iload so the worst-case node droop lands at target_droop_mv.
    # Solved exactly from the linear mesh (see calibrate_iload).
    out.append(f".param Iload={Iload:.4g} PnTSV=4u tile=400u")
    out.append(f".param Vout_dc={vout_dc}")
    out.append(".param Pb=57u alpha=3 Pp=150u")
    out.append(".param hVia=10u Dvia=50u")
    out.append("")

    # ---- Per-layer resistance derivations ----
    out.append("*=========== Per-Layer Resistance ===========")
    out.append(".param Nbpr={tile/0.105u}")
    out.append(".param Nseg={tile/PnTSV}")
    out.append(".param Rbpr={50*PnTSV/12/(Nbpr*Nseg)}")
    out.append(".param NnTSV={(tile/PnTSV)*(tile/(2*0.105u))}")
    out.append(".param RnTSV={2/NnTSV}")
    out.append(".param Lbsm1={2*PnTSV}")
    out.append(".param RsegBSM1={80m*Lbsm1/0.5u}")
    out.append(".param NsegBSM1={(tile/PnTSV)*(tile/(2*PnTSV))}")
    out.append(".param Rbsm1={RsegBSM1/12/NsegBSM1}")
    out.append(".param RtapOne={0.5*(2u/PnTSV)**2}")
    out.append(".param Ntap={(tile/(2*PnTSV))**2}")
    out.append(".param Rtap={RtapOne/Ntap}")
    out.append(".param Nbump={(tile/Pb)**2}")
    out.append(".param Rbsm3={40m/4/Nbump}")
    out.append(".param Rbump={10m/Nbump}")
    out.append(".param wBus={Pb*(alpha+1)/alpha}")
    out.append(".param Nbus={tile/Pp}")
    out.append(".param Rrdl3={2m*tile/wBus/3/Nbus}")
    out.append(".param RviaOne={1.68e-8*hVia/(3.14159*(Dvia/2)**2)}")
    out.append(".param Rvia={RviaOne/Nbump}")
    out.append(".param Npillar={(tile/Pp)**2}")
    out.append(".param Rrdl2={2m*0.25/Npillar}")
    out.append(".param Rrdl1={2m*0.25/Npillar}")
    out.append(".param RpillarOne={1.68e-8*120u/(3.14159*(75u/2)**2)}")
    out.append(".param Rpillar={RpillarOne/Npillar}")
    out.append("")

    # ---- Mesh aggregate values (Rsup now scales with Npads, not fixed *4) ----
    out.append(f"*=========== Mesh Aggregate ({m}x{n}, {N_pads} pads) ===========")
    out.append(f".param Ncells={N_cells}")
    out.append(".param Rlat=40m")
    out.append(f".param Rsup={{(Rpillar+Rrdl1+2*Rvia+Rrdl2+Rrdl3+Rbump)*{N_pads}}}")
    out.append(".param Rvia_cell={Rtap*Ncells}")
    out.append("")

    # ---- Horizontal laterals (upper layer) ----
    out.append("*=========== Upper Layer Horizontal Laterals ===========")
    for i in range(1, m+1):
        for j in range(1, n):
            out.append(f"RH_{i}_{j} n_{i}_{j}_h n_{i}_{j+1}_h {{Rlat}}")
    out.append("")

    # ---- Vertical laterals (lower layer) ----
    out.append("*=========== Lower Layer Vertical Laterals ===========")
    for i in range(1, m):
        for j in range(1, n+1):
            out.append(f"RV_{i}_{j} n_{i}_{j} n_{i+1}_{j} {{Rlat}}")
    out.append("")

    # ---- Vias ----
    out.append("*=========== Vias (Upper <-> Lower) ===========")
    for i in range(1, m+1):
        for j in range(1, n+1):
            out.append(f"RVia_{i}_{j} n_{i}_{j}_h n_{i}_{j} {{Rvia_cell}}")
    out.append("")

    # ---- Current sources ----
    out.append("*=========== Current Source Loads ===========")
    for i in range(1, m+1):
        for j in range(1, n+1):
            out.append(f"Ic_{i}_{j} n_{i}_{j} 0 {{Iload/Ncells}}")
    out.append("")

    # ---- Supply stacks (one per pad; count matches Npads) ----
    out.append(f"*=========== Supply Stacks ({N_pads} pads) ===========")
    for idx, (i, j) in enumerate(pads, 1):
        out.append(f"RSup_{idx} Vout n_{i}_{j}_h {{Rsup}}")
    out.append("")

    # ---- Voltage excitation + analysis (auto-generated, no manual add) ----
    out.append("*=========== Voltage Excitation + Analysis ===========")
    out.append("Vsupply Vout 0 {Vout_dc}")
    out.append(".op")
    out.append(".end")
    out.append("")

    with open(output_file, 'w') as f:
        f.write('\n'.join(out) + '\n')

    total_components = m*(n-1) + (m-1)*n + 2*N_cells + N_pads + 1  # +1 for Vsupply
    print(f"\nGenerated {m}x{n} mesh -> {output_file}")
    print(f"  Pads:                {N_pads} (spacing={pad_spacing}, start={pad_start})")
    print(f"  Horizontal laterals: {m*(n-1)}")
    print(f"  Vertical laterals:   {(m-1)*n}")
    print(f"  Vias:                {N_cells}")
    print(f"  Current sources:     {N_cells}")
    print(f"  Supply stacks:       {N_pads}")
    print(f"  Total components:    {total_components}")
    print(f"  Calibrated Iload:    {Iload:.4g} A  (Icell = {Iload/N_cells*1e3:.3f} mA)")
    print(f"  Droop @ load nodes:  {droop.min()*1e3:.1f} mV (best) -> "
          f"{droop.max()*1e3:.1f} mV (worst), target {target_droop_mv} mV")
    print(f"  Vout DC:             {vout_dc} V")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(
        description="Generate a parametric two-layer backside PDN netlist for LTSpice.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("m", type=int, help="mesh rows")
    p.add_argument("n", type=int, help="mesh cols")
    p.add_argument("-o", "--output", default=None,
                   help="output .cir path (default pdn_mesh_MxN.cir)")
    p.add_argument("--pad-spacing", type=int, nargs="+", default=None,
                   metavar="S",
                   help="cells between pads (one int for uniform, two for SX SY)")
    p.add_argument("--pad-start", type=int, nargs=2, default=(1, 1),
                   metavar=("I", "J"), help="pad-grid start cell (1-indexed)")
    p.add_argument("--target-droop", type=float, default=75,
                   help="target worst-case droop in mV")
    p.add_argument("--vout", type=float, default=0.6,
                   help="Vsupply DC value in volts")
    args = p.parse_args()

    # Normalize pad-spacing: allow one int (uniform) or two (sx sy)
    ps = args.pad_spacing
    if ps is None:
        pad_spacing = None
    elif len(ps) == 1:
        pad_spacing = ps[0]
    elif len(ps) == 2:
        pad_spacing = tuple(ps)
    else:
        sys.exit("Error: --pad-spacing takes 1 or 2 integers")

    out_file = args.output or f"pdn_mesh_{args.m}x{args.n}.cir"
    generate_pdn(args.m, args.n, out_file,
                 pad_spacing=pad_spacing,
                 pad_start=tuple(args.pad_start),
                 target_droop_mv=args.target_droop,
                 vout_dc=args.vout)
