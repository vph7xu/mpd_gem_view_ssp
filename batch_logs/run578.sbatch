#!/bin/bash
#SBATCH --job-name=gem_578
#SBATCH --partition=standard
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/sfs/ceph/standard/sbs/CERN_BEAM_TEST/mpd_gem_view_ssp_cern/mpd_gem_view_ssp/batch_logs/slurm_run578_%j.out
#SBATCH --account=phys_nrf

set -euo pipefail
cd /sfs/ceph/standard/sbs/CERN_BEAM_TEST/mpd_gem_view_ssp_cern/mpd_gem_view_ssp
export OMP_NUM_THREADS=1
/sfs/ceph/standard/sbs/CERN_BEAM_TEST/mpd_gem_view_ssp_cern/mpd_gem_view_ssp/bin/replay -c 0 -t 0 -z 1 -n -1 --tracking off --pedestal /sfs/ceph/standard/sbs/CERN_BEAM_TEST/mpd_gem_view_ssp_cern/mpd_gem_view_ssp/database/gem_ped_472.dat --common_mode /sfs/ceph/standard/sbs/CERN_BEAM_TEST/mpd_gem_view_ssp_cern/mpd_gem_view_ssp/database/CommonModeRange_472.txt /sfs/ceph/standard/sbs/CERN_BEAM_TEST/UVA_DATA_1/data_20Ar80CO2iso/20260726_DUT_GEM_hv_scan_uRWELL370_iso/cern_beamtest_exp1_578.evio.0
python3 /sfs/ceph/standard/sbs/CERN_BEAM_TEST/mpd_gem_view_ssp_cern/mpd_gem_view_ssp/analyze_all_runs.py --repo /sfs/ceph/standard/sbs/CERN_BEAM_TEST/mpd_gem_view_ssp_cern/mpd_gem_view_ssp --collect-output 578 /sfs/ceph/standard/sbs/CERN_BEAM_TEST/mpd_gem_view_ssp_cern/mpd_gem_view_ssp/run_calibrations_Muon_GEM_HVscan_ArCO2Iso.csv
