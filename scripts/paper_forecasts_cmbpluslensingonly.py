#!/usr/bin/env python
#
# Script to compute a set of forecasts that only use CMB+lensing information.

import sys, os
import argparse

import numpy as np

from caput import mpiutil

import compute_forecasts as cf

in_root_dir = "../output/specderivs/"
out_root_dir = "../output/forecasts/"


labels = [
    "chordnoise",
]

surveys = [
    "puma_32k_opt",
]

options = {}
for label in labels:
    options[label] = [
        "no_sub_dknldDinv0.3",
    ]
option_out = "cmb_plus_lensing"

verbose = False
verbose_forecasts = True

# Make list of tuples of label, survey, option for each forecast to run
runs = [
    (label, survey, option)
    for label in labels
    for survey in surveys
    for option in options[label]
]

# Partition list of runs based on number of MPI ranks
runs_local = mpiutil.partition_list(runs, mpiutil.rank, mpiutil.size)

# Create root output directory if it doesn't exist already
if mpiutil.rank0:
    if not os.path.exists(out_root_dir):
        os.makedirs(out_root_dir)

# Synchronize all MPI ranks
mpiutil.barrier()

# Print opening statement about how many runs will be processed on each rank
if mpiutil.size > 1:
    print(f"** Process {mpiutil.rank}: Computing {len(runs_local)} runs ***")

# Compute forecasts
for run in runs_local:
    label, survey, option = run
    print(f"* Process {mpiutil.rank}: Computing {label}, {option}, {survey}")

    # Make output directory for this run if it doesn't already exist
    out_dir = os.path.join(out_root_dir, option_out)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Compute cosmological parameter forecasts for LCDM + three 1-parameter extensions
    for par_basis in ["LCDM", "Mnu", "Neff", "OmegaK"]:

        if verbose_forecasts:
            print(f"Process {mpiutil.rank}: Computing {option_out}, {par_basis}, CMB")

        cf.compute_cos_parameter_forecasts(
            out_dir,
            survey,
            os.path.join(in_root_dir, label, option),
            verbose=verbose,
            pars=par_basis,
            Tb_free=True,
            omb_bbn_prior=0,
            lensing=None,
            cmb="Planck",
            exclude_LSS=True,
        )

        cf.compute_cos_parameter_forecasts(
            out_dir,
            survey,
            os.path.join(in_root_dir, label, option),
            verbose=verbose,
            pars=par_basis,
            Tb_free=True,
            omb_bbn_prior=0,
            lensing=None,
            cmb="Planck_SO",
            exclude_LSS=True,
        )

        if verbose_forecasts:
            print(
                f"Process {mpiutil.rank}: Computing {option_out}, {par_basis}, CMB+lensing"
            )

        cf.compute_cos_parameter_forecasts(
            out_dir,
            survey,
            os.path.join(in_root_dir, label, option),
            verbose=verbose,
            pars=par_basis,
            Tb_free=True,
            omb_bbn_prior=0,
            lensing="SO",
            cmb="Planck_SO",
            exclude_LSS=True,
        )
