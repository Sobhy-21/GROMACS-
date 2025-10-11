#!/bin/bash
# ================================================
# check_npt.sh — extract and summarize NPT metrics
# Forces use of npt.edr and ignores other .edr files
# Requires: gmx_mpi (GROMACS)
# ================================================

# Force npt.edr
EDR="npt.edr"

if [ ! -f "$EDR" ]; then
  echo "Error: $EDR not found in current directory!"
  exit 1
fi

BASENAME=$(basename "$EDR" .edr)
SUMMARY="${BASENAME}_summary.txt"

echo "Analyzing $EDR ..."
echo "Results will be saved to ${SUMMARY}"
echo "---------------------------------------" > "$SUMMARY"

# Terms to extract
terms=("Temperature" "Pressure" "Potential")

for term in "${terms[@]}"; do
  # Set output filename
  case $term in
    Temperature) out="${BASENAME}_temperature.xvg" ;;
    Pressure)    out="${BASENAME}_pressure.xvg" ;;
    Potential)   out="${BASENAME}_potential.xvg" ;;
  esac

  echo "Extracting $term → $out ..."
  
  # Run gmx_mpi energy non-interactively
  echo "$term" | gmx_mpi energy -f "$EDR" -o "$out" 2>/dev/null

  # Compute average and std deviation
  stats=$(awk '!/^[@#]/ {sum+=$2; sumsq+=$2*$2; n++}
               END {if(n>0){
                 mean=sum/n; sd=sqrt(sumsq/n - mean*mean);
                 printf "mean=%.4f  sd=%.4f  n=%d", mean, sd, n
               }}' "$out")

  echo -e "$term\t$stats" | tee -a "$SUMMARY"
done

echo "---------------------------------------" >> "$SUMMARY"
echo "All .xvg files and $SUMMARY generated."
