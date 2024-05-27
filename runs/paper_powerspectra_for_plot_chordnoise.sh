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

NTASKS=8

SCRIPT=../scripts/generate.py

COMPUTE_SETUP=true
COMPUTE_REC_DERIVATIVES=false
COMPUTE_FS_DERIVATIVES=false
COMPUTE_LENSING_DERIVATIVES=false

OUTPUT_NOSUB="--output_root_directory ../output/specderivs/pkplot_chordnoise/withnoise/"
OUTPUT_NOTHERMALNOISE="--output_root_directory ../output/specderivs/pkplot_chordnoise/nothermalnoise/"
OUTPUT_NONOISE="--output_root_directory ../output/specderivs/pkplot_chordnoise/nonoise/"

WITHNOISE="--t_int 5.0 --HI_stoch_file ../input_files/obuljen/TNG_z_meanPerr_stdPerr_kmin0.1_kmax0.3.txt"
NOTHERMALNOISE="--t_int 1e40 --HI_stoch_file ../input_files/obuljen/TNG_z_meanPerr_stdPerr_kmin0.1_kmax0.3.txt"
NONOISE="--t_int 1e40 --HI_stoch_file ../input_files/obuljen/TNG_z_meanPerr_stdPerr_kmin0.1_kmax0.3.txt --HI_stoch_multiplier 1e-20"

PUMA_REDSHIFTS="0.25 5.25 10"
PUMA32K_SURVEY="--HI --fsky 0.5 --Nside 256 --fill_factor 0.5 --dish_diameter 6.0 --hex_pack True --aperture_efficiency 0.7 --sky_coupling 1 --omt_coupling 1 --T_ground 0 --T_ampl 30.0 --remove_lowk_delta2_powspec"
PUMA5K_SURVEY="--HI --fsky 0.5 --Nside 100 --fill_factor 0.5 --dish_diameter 6.0 --hex_pack True --aperture_efficiency 0.7 --sky_coupling 1 --omt_coupling 1 --T_ground 0 --T_ampl 30.0 --remove_lowk_delta2_powspec"

CHORD_REDSHIFTS="0.25 3.75 7"
CHORD_SURVEY="--HI --fsky 0.5 --Nside 23 --fill_factor 1.0 --dish_diameter 6.0 --hex_pack False --aperture_efficiency 0.7 --sky_coupling 1 --omt_coupling 1 --T_ground 0 --T_ampl 30.0 --remove_lowk_delta2_powspec"


VERBOSE="--verbose"
STOCH="--HI_sampling_file ../input_files/obuljen/TNG_z_Psampling.txt"
BIAS="--b_file ../input_files/obuljen/TNG_z_b1S_fitfunc.txt  --b2_file ../input_files/obuljen/TNG_z_b2S_fitfunc.txt  --bs_file ../input_files/obuljen/TNG_z_bG2S_fitfunc.txt --use_G2_basis"
KMAX="--kmax 3.0"


####### PUMA-32k

NAME="puma_32k_opt"
OUTPUT="$OUTPUT_NOSUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA32K_SURVEY $STOCH $BIAS $WITHNOISE $KMAX"
run $NAME "$PARS"

NAME="puma_32k_opt"
OUTPUT="$OUTPUT_NOTHERMALNOISE"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA32K_SURVEY $STOCH $BIAS $NOTHERMALNOISE $KMAX"
run $NAME "$PARS"

NAME="puma_32k_opt"
OUTPUT="$OUTPUT_NONOISE"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA32K_SURVEY $STOCH $BIAS $NONOISE $KMAX"
run $NAME "$PARS"


###### PUMA-5k

NAME="puma_5k_opt"
OUTPUT="$OUTPUT_NOSUB"
PARS="$PUMA_REDSHIFTS $VERBOSE $OUTPUT $PUMA5K_SURVEY $STOCH $BIAS $WITHNOISE $KMAX"
run $NAME "$PARS"


##### CHORD

NAME="chord_opt"
OUTPUT="$OUTPUT_NOSUB"
PARS="$CHORD_REDSHIFTS $VERBOSE $OUTPUT $CHORD_SURVEY $STOCH $BIAS $WITHNOISE $KMAX"
run $NAME "$PARS"
