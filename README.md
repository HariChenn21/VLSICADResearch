# VLSI CAD Research — Backside Power Delivery Network (PDN) Modeling

Research code and simulation models for analyzing **backside power delivery networks** in
advanced 2.5D/3D VLSI, together with the **buck converter** that supplies them. The core
tool is a parametric netlist generator that builds arbitrarily large two-layer PDN meshes,
auto-calibrated so the worst-case IR droop lands on a target value.

---

## Contents

| File / group | Description |
|---|---|
| `pdn_gen.py` | Parametric two-layer backside PDN mesh netlist generator (with a built-in MNA solver for auto-calibration). |
| `pdn_mesh_*.cir` | Generated ngspice/LTspice netlists for meshes from `4x4` up to `160x160`. |
| `pdn_mesh_10x10.raw` | A small worked simulation output, kept as an example. |
| `BuckConverter*.asc` | LTspice buck converter schematics that drive the PDN (`BuckConverter`, `-WithLayers`, `-WithLayers2D`). |
| `*.plt` | LTspice plot settings (droop / `V(vout)-V(vload)` and per-corner node views). |
| `*.log` | LTspice run logs. |

> Large binary simulation outputs (`*.raw`) are intentionally **not** tracked in git — they
> are regenerable and can be hundreds of MB. Only the small `10x10` example is included.

---

## The PDN model

Each mesh cell sits at an intersection `(i, j)` of a two-layer backside network:

- `n_i_j_h` — **upper layer** node (BSM3, horizontal laterals)
- `n_i_j` — **lower layer** node (BSM1, vertical laterals + the current-source load)
- a **via** (`RVia_i_j`) connects the two layers at every cell

**Supply pads** inject `Vout` into the upper layer through a supply-stack resistance `Rsup`
(pillar + RDL + vias + bump). Pad placement is configurable — four corners by default, or a
regular grid via `--pad-spacing` / `--pad-start`. `Rsup` scales with the number of pads so the
parallel (tile-level) supply resistance stays constant regardless of pad count.

All per-layer resistances (`Rlat`, `Rvia_cell`, `Rsup`, taps, bumps, RDL, pillars) are derived
from physical parameters — tile pitch, TSV pitch `PnTSV`, bump pitch `Pb`, via geometry, etc. —
emitted as `.param` expressions so the `.cir` files stay self-documenting and runnable standalone.

### Auto-calibration

Droop is linear in load current, so the generator does one exact **Modified Nodal Analysis (MNA)**
solve of the mesh (`solve_droop`) to find the total `Iload` that places the *worst-case* node
droop at a target (default 75 mV). No iteration, no manual tuning — the emitted netlist already
carries the calibrated `Iload` and a `.op` analysis.

---

## Usage

Requires Python 3 with NumPy:

```bash
pip install numpy
```

Generate a netlist:

```bash
# 20x20 mesh, default 4-corner pads
python pdn_gen.py 20 20

# ~7x7 pad grid, spacing 3, starting at cell (2,2)
python pdn_gen.py 20 20 --pad-spacing 3 --pad-start 2 2

# Target 100 mV worst-case droop at Vout = 1.0 V
python pdn_gen.py 20 20 --pad-spacing 3 --vout 1.0 --target-droop 100

# Custom output path
python pdn_gen.py 32 32 --pad-spacing 4 -o big.cir
```

**Options**

| Flag | Meaning | Default |
|---|---|---|
| `m n` | mesh rows / columns (positional, ≥ 2) | — |
| `-o, --output` | output `.cir` path | `pdn_mesh_MxN.cir` |
| `--pad-spacing S [S2]` | cells between pads (one value = uniform, two = `SX SY`) | 4 corners |
| `--pad-start I J` | first-pad position, 1-indexed | `1 1` |
| `--target-droop MV` | worst-case droop target (mV); auto-calibrates `Iload` | `75` |
| `--vout V` | `Vsupply` DC value | `0.6` |

Each run prints a summary: pad count, component counts (laterals / vias / current sources /
supply stacks), the calibrated `Iload` and per-cell current, and the best→worst droop spread.

Run the generated netlist in **LTspice** or **ngspice**:

```bash
ngspice pdn_mesh_20x20.cir
```

---

## Buck converter models

The `BuckConverter*.asc` LTspice schematics model the DC-DC converter feeding the PDN. Variants
add progressively more physical detail (`-WithLayers`, `-WithLayers2D`) and sweep the TSV pitch
`pntsv` across `2u / 4u / 8u / 16u`. The `.plt` files capture the droop view
(`V(vout)-V(vload)`) and per-corner node voltages used to evaluate supply quality under load.

---

## Repository notes

- `*.raw` outputs are gitignored (regenerate by re-running the simulation); the `10x10` example
  is force-included for reference.
- Netlists are plain text and diff cleanly, so mesh changes are easy to review in git.
