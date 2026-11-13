# A Newton/Picard semi-implicit algorithm for finite-strain plasticity: code of the revised manuscript

P. Areias, T. Rabczuk, N. Silvestre, M.F. Vaz, *Computational Mechanics* (revised version).

Fortran 2003 sources of the constitutive update of the paper and of the programs behind its material-point results.
Section, figure and table numbers are those of the revised manuscript.

## What is here

| File | Contents |
|---|---|
| `npicard_plasticity.f90` | The module. `np_update_alg1`: Algorithm 1 (scalar iteration on the multiplier at a frozen flow direction, Picard update of the direction, consistent tangent of Section 5). `np_update_sub`: the adaptive substepping of Section 6.5 with the exact tangent of the substepped increment. `np_update_base`: the staggered-hardening Newton baseline of Section 8.2. Elastic laws and yield functions of Section 8.1 (Hosford, the polycarbonate criterion, Yld2004-18p). The header of the file states the model, the conventions and the status codes. |
| `npicard_guo_kernel.f90` | Energy, stress and tangent of the PLA model of Section 8.1.2 (code generated with AceGen). |
| `acegen/yld_ghorbel.nb` | Mathematica/AceGen notebook of the polycarbonate criterion of Section 8.1.1; its output is the generated part of the criterion in `npicard_plasticity.f90`. |
| `acegen/guo2008voidfromforce.nb`, `acegen/guo2008drivingforce.nb`, `acegen/guo2008modified.nb` | Mathematica/AceGen notebooks of the three routines of `npicard_guo_kernel.f90`: void fraction from the driving force, driving force (largest transverse stretch), and energy with its stress and tangent. |
| `pctests.f90` | Tension and shear tests of the polycarbonate: Figures 2 and 3, the numbers of Section 8.1.1. |
| `wkmap.f90` | Stored energy and hardening variable over the strain plane: Figures 4 and 7. |
| `isomap.f90` | Iso-error maps: Table 3, Figures 9 to 13. |
| `basemap.f90` | The baseline on the PLA map: Figure 10(f) and the baseline numbers of Section 8.2. |
| `limitscan.f90` | Largest trial ratio of a single increment along rays: the ratios 1.18, 1.29, 1.28 and 1.03 of Section 8.2. |
| `raydiag.f90` | The first ray of a single increment (no root, root passed by the doubling): the shares of Section 8.2. |
| `countmap.f90` | Iteration counts of Section 6.4. |
| `dencheck.f90`, `symcheck.f90` | Denominator of the multiplier increment and symmetry of the tangent: the numbers of Section 5.2. |
| `nnsens.f90` | A fixed number of Picard updates beyond the contraction range: the amplification quoted in Section 8.3. |
| `Makefile`, `run_paper.sh` | Build and runs. |
| `check/` | Reference values of the quick check. |

Not here: the framed beam of Section 8.3 was computed with this module inside the authors' finite element code,
which is not part of this package, and the figures were drawn from the output files with the authors' plotting
scripts. Table 2 and Figure 8 are evaluations of the Yld2004-18p function (`np_yield` with `np_barlat_set`).

## Build and run

Requirements: gfortran, GNU make, bash and awk. No library is needed. Mathematica with AceGen is needed only to
regenerate the routines of `acegen/`; the generated Fortran is already in the module and the kernel.

```
make                 # the ten programs
make check           # tension and shear tests with ten increments against check/ (a few seconds)
./run_paper.sh       # every run of the paper, outputs in out/   (./run_paper.sh light | maps for a part)
```

`make check` compares strain, stress and hardening variable with the reference files and prints the largest
difference; "agrees" means below 1e-8 of the largest value of the column. The last column of those files is a ratio
of two numbers at round-off level once the Picard iteration has converged; it depends on the compiler and is not
compared.

Tested with gfortran 15.3.0 (Ubuntu) on a 24-core workstation: `./run_paper.sh` took 29 minutes, 19 of them for the
map of AA2090-T3 from the undeformed state; `./run_paper.sh light` takes a few minutes. The programs use every core
(OpenMP), and the results do not depend on the number of threads. For another compiler: `make FC=... FFLAGS=...`.

## Runs and outputs

| Run (see `run_paper.sh`) | Output in `out/` | In the paper |
|---|---|---|
| `pctests tension\|shear <Ns> <Nn>` | `pc_<test>_<Ns>_<Nn>.txt` | Figures 2 and 3; `Nn` = 0 is the converged update |
| `wkmap pc\|pla 200` | `wk_pc.txt`, `wk_pla.txt` | Figures 4 and 7 |
| `isomap pc\|pla virgin 81 0.25 1000` | `pc_virgin.txt`, `pla_virgin.txt` | Table 3, Figures 9, 10(a)-(e) and 11 |
| `basemap pla virgin 81 0.25 1000` | `pla_base.txt` | Figure 10(f) |
| `isomap barlat\|iso point 81 5 1000` | `barlat_point.txt`, `iso_point.txt` | Table 3, Figure 12 |
| `isomap barlat virgin 81 0.1 1000` | `barlat_virgin.txt` | Table 3, Figure 13 |
| `limitscan <law> <map> 720`, `limitscan <law> general 2000` | `limit_*.txt` | Section 8.2 |
| `raydiag pc\|pla` | `raydiag_*.txt` | Section 8.2 |
| `countmap pla\|pc 81 0.25` | `countmap_*.txt` | Section 6.4 |
| `dencheck`, `symcheck` | `dencheck*.txt`, `symcheck.txt` | Section 5.2 |
| `nnsens` | `nnsens.txt` | Section 8.3 |

Each program prints its usage in the header of its source file, and every output file starts with a line that
names its columns.

## Words of the code

`virgin`: increments from the undeformed state. `point`: increments from the balanced biaxial yield point. `iso`:
the isotropic control (Hosford, exponent 8, with the elasticity and hardening of AA2090-T3). `parts`: the equal
parts into which a failed increment is divided (Section 6.5). `info`: the status of an update (0 a stress with the
Picard iteration converged; the other values are listed in the header of the module).
