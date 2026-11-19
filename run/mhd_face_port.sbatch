#!/bin/bash
#SBATCH -J cof_port
#SBATCH -p p.shared
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nice=10000
#SBATCH -t 02:00:00
#SBATCH -o /orion/ptmp/uli/orion/enzyme/cofactor_runs/port_%j.out
# Porting check, target mhd_face: cofactor's extraction + harness on the build copies and no-FMA dumps of job 231830
# (enzyme/mhd_audit). Bar (fixed before the run): extraction byte-identical to the old gen/<cfg>, and every results row
# identical to enzyme/mhd_audit/out/results_<cfg>.txt (name, n, max_rel). No GIZMO build or run.
set -u
C=/u/uli/cofactor
A=/orion/ptmp/uli/orion/enzyme/mhd_audit
W=/orion/ptmp/uli/orion/enzyme/cofactor_runs/port_${SLURM_JOB_ID:?sbatch only}
mkdir -p $W
source /etc/profile.d/modules.sh
module purge; module load gcc openmpi gsl hwloc fftw-mpi hdf5-mpi clang/18
export LD_LIBRARY_PATH=$HDF5_HOME/lib:$GSL_HOME/lib:$FFTW_HOME/lib:${LD_LIBRARY_PATH:-}
ENZ=/orion/ptmp/uli/orion/enzyme/src/build18/Enzyme/ClangEnzyme-18.so
MPIF="$(mpicxx -show | cut -d' ' -f2-)"
git -C $C log --oneline -1 2>/dev/null; git -C $C status --short
echo "start $(date)"
comp() { # $1 cfg, $2 extra define
  timeout -k 10 1800 clang++ -std=c++14 -O2 -ffp-contract=off -fplugin=$ENZ -mllvm -enzyme-loose-types=1 -mllvm -enzyme-globals-default-inactive=1 \
      -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unknown-pragmas -Wno-deprecated-declarations $2 \
      -I$W/gen/$1 -I$C/lib -I$A/build/$1/hydro -I$A/build/$1 -I$GSL_HOME/include -I$HDF5_HOME/include -DH5_USE_16_API -DHDF5_DISABLE_VERSION_CHECK \
      $C/targets/mhd_face/harness.cc $A/build/$1/declarations/allvars.cc $MPIF -L$GSL_HOME/lib -lgsl -lgslcblas -lm -o $W/bin_$1 \
      > $W/compile_$1$2.txt 2>&1
  local e=$?; echo "compile $1 $2 exit $e"; grep -m 5 -E " error|Assertion" $W/compile_$1$2.txt; return $e
}
ALL=0
for cfg in mfm_psi_nf mfv_psi_nf mfm_ded_nf mfv_ded_nf mfm_psi_kvp_nf mfm_ded_kvp_nf; do
  echo "=== $cfg"
  timeout 120 python3 $C/targets/mhd_face/spec.py $A/build/$cfg $W/gen/$cfg $W/manifests > $W/extract_$cfg.txt 2>&1; e=$?
  echo "extract exit $e"; grep cofactor $W/extract_$cfg.txt
  same=1; for f in $A/harness/gen/$cfg/ext_*.h; do cmp -s $f $W/gen/$cfg/$(basename $f) || { echo "EXTRACTION DIFFERS: $(basename $f)"; same=0; }; done
  [ $same = 1 ] && echo "extraction byte-identical to the old gen/$cfg" || ALL=1
  comp $cfg "" || comp $cfg "-DFWD_ONLY" || { ALL=1; continue; }
  COFACTOR_BARS=$C/targets/mhd_face/bars.txt timeout -k 10 3000 $W/bin_$cfg $A/runs/$cfg/output/face_dump.bin $cfg > $W/results_$cfg.txt 2>&1
  echo "harness exit $? (0 pass, 1 fail, 2 incomplete)"; grep "^cofactor\|FAIL (bar" $W/results_$cfg.txt
  python3 $C/tools/compare_results.py $A/out/results_$cfg.txt $W/results_$cfg.txt > $W/compare_$cfg.txt; c=$?
  tail -12 $W/compare_$cfg.txt; [ $c = 0 ] || ALL=1
done
echo "porting check overall: $([ $ALL = 0 ] && echo 'PASS (extraction identical, every row identical)' || echo FAIL)"
echo "done $(date)"
