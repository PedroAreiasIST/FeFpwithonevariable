# FeFpwithonevariable

Sources of the Newton/Picard return mapping for finite-strain plasticity of

P. Areias, T. Rabczuk, N. Silvestre, M.F. Vaz, *A Newton/Picard semi-implicit algorithm for finite-strain
plasticity*, Computational Mechanics.

## `NewtonPicard_revised_code/`: the code of the revised manuscript

Fortran 2003 module and programs that reproduce the material-point results of the revised manuscript: Algorithm 1
(safeguarded scalar Newton iteration for the plastic multiplier at a frozen flow direction, Picard update of the
direction, consistent tangent), the adaptive substepping with the exact tangent of a substepped increment, the
staggered-hardening Newton baseline, and the materials of the examples (polycarbonate, porous PLA, AA2090-T3 with
Yld2004-18p). Requirements: gfortran and GNU make.

```
cd NewtonPicard_revised_code
make            # the ten programs
make check      # a quick check against reference values (a few seconds)
./run_paper.sh  # the runs behind the figures and tables (outputs in out/)
```

Its README lists which program produces which figure, table and number of the paper. The framed beam was computed
with the same module inside the authors' finite element code, which is not part of this repository.

## Files at the top level

Mathematica/AceGen notebooks (`*.nb`) and the Fortran routines generated from them (yield functions, elastic law,
tangent operators), `plasticity1.f90` and `readlatex.py`, uploaded in April 2026.
