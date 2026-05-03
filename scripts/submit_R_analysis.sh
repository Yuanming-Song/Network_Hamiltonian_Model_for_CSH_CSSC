#!/bin/bash
#SBATCH --job-name=NHM_RestFiber_analysis
#SBATCH --account=dtobias_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --output=NHM_RestFiber_analysis.out
#SBATCH --error=NHM_RestFiber_analysis.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yuanmis1@uci.edu
#SBATCH --time=120:00:00
#SBATCH --mem=180G

# Load required modules
module purge
module load R

# Allow the user to define the publication bundle root once.
MAINDIR="${1:-${NHM_MAIN_DIR:-$(cd "$(dirname "$0")/.." && pwd)}}"
DATADIR="${NHM_DATA_DIR:-$MAINDIR/data/Restrained_Fiber_simulation}"
SCRIPT_DIR="$MAINDIR/scripts"

# Set R library path
export R_LIBS_USER="/dfs9/tw/yuanmis1/R_libs/"
export NHM_MAIN_DIR="$MAINDIR"
export NHM_DATA_DIR="$DATADIR"

cd "$SCRIPT_DIR" || exit 1

# Process data first (fix CSX to CSSC naming)
echo "Processing MD data..."
#Rscript process_MD_data.R

# Run atomistic analysis (already completed)
# echo "Starting atomistic analysis..."
# Rscript csh_cssc_1node_estimation.R

# Run coarse-grained analysis (already completed)
# echo "Starting coarse-grained analysis..."
# echo "Detailed output will be written to: Single_Node_Edges/nhm_analysis_results/csh_cssc_1node_fit_CG_with_degcor.log"
# Rscript csh_cssc_1node_estimation_CG.R
# echo "CG analysis complete!"
# echo "Check log file: Single_Node_Edges/nhm_analysis_results/csh_cssc_1node_fit_CG_with_degcor.log"

# Run restrained fiber analysis
echo "Starting restrained fiber analysis..."
echo "Main directory: $NHM_MAIN_DIR"
echo "Data directory: $NHM_DATA_DIR"
echo "Detailed output will be written to: $NHM_MAIN_DIR/results/csh_cssc_1node_fit_RestFiber.log"
echo "Using checkpoint recovery - will skip data processing and baseline model fitting"
Rscript csh_cssc_1node_estimation_RestFiber.R
echo "Restrained fiber analysis complete!"
echo "Check log file: $NHM_MAIN_DIR/results/csh_cssc_1node_fit_RestFiber.log"
