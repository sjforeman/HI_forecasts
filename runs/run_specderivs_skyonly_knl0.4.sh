#!/bin/bash

run () {
  echo "************"
  echo $1
  echo "************"

  if [ "$COMPUTE_SETUP" = true ]
  then
      mpirun -n $NTASKS python $SCRIPT $1 setup $2
  fi
  if [ "$COMPUTE_REC_DERIVATIVES" = true ]
  then
      mpirun -n $NTASKS python $SCRIPT $1 rec $2
  fi
  if [ "$COMPUTE_FS_DERIVATIVES" = true ]
  then
      mpirun -n $NTASKS python $SCRIPT $1 fs $2
  fi
  if [ "$COMPUTE_LENSING_DERIVATIVES" = true ]
  then
      mpirun -n $NTASKS python $SCRIPT $1 lens $2
  fi
}

NTASKS=6

SCRIPT=../scripts/generate.py

COMPUTE_SETUP=true
COMPUTE_REC_DERIVATIVES=true
COMPUTE_FS_DERIVATIVES=true
COMPUTE_LENSING_DERIVATIVES=true

OUTPUT_NOSUB_DIR=../output/specderivs/skyonly/no_sub_knl0.4/
OUTPUT_SUB_DIR=../output/specderivs/skyonly/sub_knl0.4/

OUTPUT_NOSUB="--output_root_directory $OUTPUT_NOSUB_DIR"
OUTPUT_SUB="--output_root_directory $OUTPUT_SUB_DIR"

PUMA_REDSHIFTS="0.5 6.0 11"
PUMA32K_SURVEY="--HI --fsky 0.5 --Nside 256 --fill_factor 0.5 --t_int 5.0 --dish_diameter 6.0 --hex_pack True --aperture_efficiency 0.7 --sky_coupling 1 --omt_coupling 1 --T_ground 0 --T_ampl 0 --knl_z0 0.4 --dknl_dz 0.0 --dknl_dDinv 0.0 --remove_lowk_delta2_powspec --kmax 3.0"
PUMA5K_SURVEY="--HI --fsky 0.5 --Nside 100 --fill_factor 0.5 --t_int 5.0 --dish_diameter 6.0 --hex_pack True --aperture_efficiency 0.7 --sky_coupling 1 --omt_coupling 1 --T_ground 0 --T_ampl 0 --knl_z0 0.4 --dknl_dz 0.0 --dknl_dDinv 0.0 --remove_lowk_delta2_powspec --kmax 3.0"
PUMA_DERIVS="--deriv_directory ../output/specderivs/chordnoise/no_sub_knl0.4/puma_32k_opt/derivatives/ --deriv_Cl_directory ../output/specderivs/chordnoise/no_sub_knl0.4/puma_32k_opt/derivatives_Cl/ --deriv_recon_directory ../output/specderivs/chordnoise/no_sub_knl0.4/puma_32k_opt/derivatives_recon/"


CHORD_REDSHIFTS="0.2 3.7 7"
CHORD_SURVEY="--HI --fsky 0.5 --Nside 23 --fill_factor 1.0 --t_int 5.0 --dish_diameter 6.0 --hex_pack False --aperture_efficiency 0.7 --sky_coupling 1 --omt_coupling 1 --T_ground 0 --T_ampl 0 --knl_z0 0.4 --dknl_dz 0.0 --dknl_dDinv 0.0 --remove_lowk_delta2_powspec --kmax 3.0"
CHORD_DERIVS="--deriv_directory ../output/specderivs/chordnoise/no_sub_knl0.4/chord_opt/derivatives/ --deriv_Cl_directory ../output/specderivs/chordnoise/no_sub_knl0.4/chord_opt/derivatives_Cl/ --deriv_recon_directory ../output/specderivs/chordnoise/no_sub_knl0.4/chord_opt/derivatives_recon/"


VERBOSE="--verbose"
STOCH="--HI_stoch_file ../input_files/obuljen/TNG_z_meanPerr_stdPerr_kmin0.1_kmax0.3.txt --HI_sampling_file ../input_files/obuljen/TNG_z_Psampling.txt"
BIAS="--b_file ../input_files/obuljen/TNG_z_b1S_fitfunc.txt  --b2_file ../input_files/obuljen/TNG_z_b2S_fitfunc.txt  --bs_file ../input_files/obuljen/TNG_z_bG2S_fitfunc.txt --use_G2_basis"



####### PUMA-32k

NAME="puma_32k_opt"
OUTPUT="$OUTPUT_NOSUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA32K_SURVEY $STOCH $BIAS"
run $NAME "$PARS"

NAME="puma_32k_opt"
OUTPUT="$OUTPUT_SUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA32K_SURVEY $STOCH $BIAS --remove_lowk_delta2_cov $PUMA_DERIVS"
run $NAME "$PARS"

NAME="puma_32k_pess"
OUTPUT="$OUTPUT_NOSUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA32K_SURVEY $STOCH $BIAS  --pessimistic_wedge $PUMA_DERIVS"
run $NAME "$PARS"

NAME="puma_32k_pess"
OUTPUT="$OUTPUT_SUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA32K_SURVEY $STOCH $BIAS  --pessimistic_wedge --remove_lowk_delta2_cov $PUMA_DERIVS"
run $NAME "$PARS"


####### PUMA-5k

NAME="puma_5k_opt"
OUTPUT="$OUTPUT_NOSUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA5K_SURVEY $STOCH $BIAS $PUMA_DERIVS"
run $NAME "$PARS"

NAME="puma_5k_opt"
OUTPUT="$OUTPUT_SUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA5K_SURVEY $STOCH $BIAS --remove_lowk_delta2_cov $PUMA_DERIVS"
run $NAME "$PARS"


NAME="puma_5k_pess"
OUTPUT="$OUTPUT_NOSUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA5K_SURVEY $STOCH $BIAS --pessimistic_wedge $PUMA_DERIVS"
run $NAME "$PARS"


NAME="puma_5k_pess"
OUTPUT="$OUTPUT_SUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA5K_SURVEY $STOCH $BIAS --pessimistic_wedge --remove_lowk_delta2_cov $PUMA_DERIVS"
run $NAME "$PARS"


###### CHORD

NAME="chord_opt"
OUTPUT="$OUTPUT_NOSUB"
PARS="$CHORD_REDSHIFTS $VERBOSE $OUTPUT $CHORD_SURVEY $STOCH $BIAS"
run $NAME "$PARS"

NAME="chord_opt"
OUTPUT="$OUTPUT_SUB"
PARS="$CHORD_REDSHIFTS $VERBOSE $OUTPUT $CHORD_SURVEY $STOCH $BIAS --remove_lowk_delta2_cov $CHORD_DERIVS"
run $NAME "$PARS"


NAME="chord_pess"
OUTPUT="$OUTPUT_NOSUB"
PARS="$CHORD_REDSHIFTS $VERBOSE $OUTPUT $CHORD_SURVEY $STOCH $BIAS --pessimistic_wedge $CHORD_DERIVS"
run $NAME "$PARS"

NAME="chord_pess"
OUTPUT="$OUTPUT_SUB"
PARS="$CHORD_REDSHIFTS $VERBOSE $OUTPUT $CHORD_SURVEY $STOCH $BIAS --pessimistic_wedge --remove_lowk_delta2_cov $CHORD_DERIVS"
run $NAME "$PARS"







