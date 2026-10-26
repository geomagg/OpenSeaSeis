# Depth focusing velocity analysis: 2D synthetic

Complete velocity-analysis loop on a 2D OBN synthetic with a known velocity model:
modelling, migration with focusing analysis, tau panels to velocity coherence,
picks to interval velocity, new model, and iterations. Starting from 90% of the
sediment velocity, the model converges to within 2% of the true one in two updates.

The step-by-step description (in Portuguese) is in
`doc/SPK_depth_focusing_velocity_analysis_pt.pdf`.

## Requirements

- SeaSeis built with the SPK modules, `seaseis` in the PATH
- Python 3 with numpy, scipy and matplotlib (`pip install --user numpy scipy matplotlib`)

## Run

```
cd example_flows/spk_focus_synthetic
./run_all.sh 4        # 4 = modelling runs in parallel (default 2)
```

About 25 minutes with 4 CPUs, most of it in the modelling of the 17 nodes.
If `node*.cseis` already exist, the modelling is skipped.

| Step | What | Outputs |
| --- | --- | --- |
| 1 | Velocity models (`mkmod.py nx nz dz factor file`): true, and 90% in the sediments | `VTRUE10`, `VTRUE25`, `V90_25` |
| 2 | Modelling of 17 nodes, SPK_XDVSSYN2D (`runfd.sh`, `fd.flow`) | `node1.cseis` ... `node17.cseis` |
| 3 | Migration with focusing analysis, SPK_XVSPMIG1 (`MIG.flow`), true (control) and 90% models; SPK_TAU2VFIX (`T2V.flow`) | `TRUE_*`, `V90_*` |
| 4 | Picks (T, Vavg) -> Vint -> new model (`update.py`), then 3 iterations of steps 3-4 | `VIT1_25` ... `VIT4_25`, `IT1_*` ... `IT3_*` |
| 5 | Tau deviation and reflector depths per migration (`stats.py`), figures (`figs.py`) | table, `fig_*.png` |

Expected output of step 5 (DT = tau deviation in ms at the 4 reflectors, z = image depth in m;
true depths 1500, 2100, 2800, 3500 m):

```
V90 DT [-172, -276, -372, -404] z [1445, 1975, 2580, 3200]
IT1 DT [44, -52, -132, -292] z [1505, 2085, 2735, 3355]
IT2 DT [-12, 4, -28, -20] z [1495, 2100, 2785, 3480]
IT3 DT [44, -36, 28, -20] z [1505, 2090, 2800, 3490]
TRUE DT [-12, -12, -20, -20] z [1495, 2095, 2795, 3485]
```

IT2 is at the level of the run with the true velocity: the loop has converged.

## Single steps

```
seaseis -D V=VTRUE25 -D P=TRUE -f MIG.flow    # migration + focusing analysis
seaseis -D P=TRUE -f T2V.flow                 # tau panels -> velocity coherence
python3 update.py TRUE_vfix.sgy 9 VNEW_25     # picks -> Vint -> model
```

Note: SPK_TAU2VFIX is used with `tau_dt 8` (tau spacing of the decimated panels)
and `factor 2`. With the factor built into the migration's own picks (`_vtps.txt`,
f = 0.8 of xvspmig1) the correction is about 3 times too large for this geometry.
