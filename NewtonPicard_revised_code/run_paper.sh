#!/usr/bin/env bash
# run_paper.sh [check|light|maps|all] - the runs behind the figures and tables of the paper; outputs in out/.
#   check  the tension and shear tests with ten increments; strain, stress and hardening variable are compared
#          with check/ (a few seconds)
#   light  everything except the plane stress maps (a few minutes on one workstation)
#   maps   the iso-error maps: the five maps of Table 3 and the baseline map (about half an hour on 24 cores; the
#          map of AA2090-T3 from the undeformed state takes most of it), and the six maps of the figures (a few
#          minutes)
#   all    light and maps (default)
# The programs use every core (OpenMP); the results do not depend on the number of threads.
set -e
cd "$(dirname "$0")"
what=${1:-all}
mkdir -p out

tests() {              # Section 8.1.1: tension and shear tests of the polycarbonate
  for t in tension shear; do
    for nn in 1 5 10 0; do ./pctests $t 10 $nn out/pc_${t}_10_${nn}.txt > /dev/null; done
  done
}

light() {
  tests
  for t in tension shear; do ./pctests $t 1000 0 out/pc_${t}_1000_0.txt > /dev/null; done
  ./wkmap pc 200 out/wk_pc.txt                       # Section 8.1: W and kappa as functions of E11, E22
  ./wkmap pla 200 out/wk_pla.txt
  ./countmap pla 81 0.25 > out/countmap_pla.txt      # Section 6.4: iteration counts
  ./countmap pc 81 0.25 > out/countmap_pc.txt
  ./symcheck > out/symcheck.txt                      # Section 5.2: symmetry of the tangent
  ./dencheck > out/dencheck.txt                      # Section 5.2: denominator of the multiplier increment
  ./dencheck 321 > out/dencheck_321.txt
  ./dencheck 641 > out/dencheck_641.txt
  ./raydiag pc out/raydiag_pc.txt                    # Section 8.2: the first ray of a single increment
  ./raydiag pla out/raydiag_pla.txt
  for lm in "pc virgin" "pla virgin" "barlat point" "barlat virgin" "iso point"; do   # Section 8.2: limits
    set -- $lm; ./limitscan $1 $2 720 out/limit_$1_$2.txt
  done
  for l in pc pla barlat iso; do ./limitscan $l general 2000 out/limit_${l}_general.txt; done
  { echo "NNS_PERT=1 NNS_DELTA=1e-7"; NNS_PERT=1 ./nnsens                 # Section 8.3: N_n updates
    echo "NNS_PERT=3 NNS_DELTA=1e-7"; NNS_PERT=3 ./nnsens; } > out/nnsens.txt
}

maps() {               # Section 8.2: iso-error maps
  ./maymap pc virgin 101 0.25 90 out/fig_pc_virgin.txt          # Figures 9 to 13: N equal increments
  ./maymap pla virgin 101 0.25 90 out/fig_pla_virgin.txt
  ./maymap pla virgin 101 0.25 90 out/fig_pla_base.txt base
  ./maymap barlat point 101 5 90 out/fig_barlat_point.txt
  ./maymap iso point 101 5 90 out/fig_iso_point.txt
  ./maymap barlat virgin 101 0.1 90 out/fig_barlat_virgin.txt
  ./isomap pc virgin 81 0.25 1000 out/pc_virgin.txt              # Table 3: single increments
  ./isomap pla virgin 81 0.25 1000 out/pla_virgin.txt
  ./basemap pla virgin 81 0.25 1000 out/pla_base.txt
  ./isomap iso point 81 5 1000 out/iso_point.txt
  ./isomap barlat point 81 5 1000 out/barlat_point.txt
  ./isomap barlat virgin 81 0.1 1000 out/barlat_virgin.txt
}

compare() {            # strain, stress and hardening variable of out/<file> against check/<file>
  awk -v name="$(basename "$1")" '
    function abs(x) { return x < 0 ? -x : x }
    FNR == NR { if ($0 !~ /^#/ && NF >= 3) { n++; for (k = 1; k <= 3; k++) { r[n, k] = $k; if (abs($k) > big[k]) big[k] = abs($k) } }; next }
    $0 !~ /^#/ && NF >= 3 { m++; for (k = 1; k <= 3; k++) { d = abs($k - r[m, k])/big[k]; if (d > worst) worst = d } }
    END { if (m != n) printf "DIFFERENT  %s: %d rows against %d\n", name, m, n
          else if (worst <= 1e-8) printf "agrees     %s: %d rows, largest difference %.1e of the largest value of its column\n", name, n, worst
          else printf "DIFFERENT  %s: largest difference %.1e of the largest value of its column\n", name, worst }' "$1" "out/$(basename "$1")"
}

case $what in
  check) tests
         for f in check/*.txt; do compare "$f"; done ;;
  light) light ;;
  maps)  maps ;;
  all)   light; maps ;;
  *)     echo "usage: $0 [check|light|maps|all]"; exit 1 ;;
esac
