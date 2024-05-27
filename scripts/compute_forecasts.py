#!/usr/bin/env python
#
# Script to compute forecasts for different cosmological quantities.

import sys, os
from math import ceil
import argparse

from caput import mpiutil

from FishLSS.headers import *
from FishLSS import twoPoint, twoPointNoise, parameter_forecast as pf


def compute_bao_forecasts(
    fc_out_dir, name, out_root_dir="output", verbose=False, overwrite=False
):
    """Compute BAO forecasts for AP parameters.

    Parameters
    ----------
    fc_out_dir : string
        Directory for forecast output file.
    name: string
        Identifier for previously-computed forecasting information.
    out_root_dir: string, optional
        Root output directory that previous forecasting info was saved to.
    verbose : bool, optional
        Print status messages when loading previous forecasting info. Default: False.
    overwrite : bool, optional
        Overwrite output file if it already exists. Default: False.
    """
    # Exit if output file exists and should not be overwritten
    out_file = os.path.join(fc_out_dir, "bao.json")
    if os.path.isfile(out_file) and not overwrite:
        return

    # Load previously-computed forecast information
    fc = pf.load_forecast(
        name, verbose=verbose, out_root_dir=out_root_dir, use_mpi=False
    )

    # Compute BAO forecast
    output = pf.compute_forecast_recon_BAO(fc)

    # Organize forecast output in dict
    data = {
        "onesigma_aperp": output["onesigma_aperp"].tolist(),
        "onesigma_apar": output["onesigma_apar"].tolist(),
        "z_center": fc["Centers of redshift bins"],
        "z_edge": fc["Edges of redshift bins"],
        "F_z": output["F_z"].tolist(),
    }

    # Save outputs to file
    with open(out_file, "w") as f:
        json.dump(data, f, indent=2)


def compute_sigma8_fixedshape_forecasts(
    fc_out_dir,
    name,
    out_root_dir="output",
    verbose=False,
    Tb_free=True,
    overwrite=False,
):
    """Compute sigma_8 forecasts where linear power spectrum shape is held fixed.

    Parameters
    ----------
    fc_out_dir : string
        Directory for forecast output file.
    name: string
        Identifier for previously-computed forecasting information.
    out_root_dir: string, optional
        Root output directory that previous forecasting info was saved to.
    verbose : bool, optional
        Print status messages when loading previous forecasting info. Default: False.
    Tb_free : bool, optional
        Whether to marginalize over mean brightness temparature T_b. Default: True.
    overwrite : bool, optional
        Overwrite output file if it already exists. Default: False.
    """
    # Exit if output file exists and should not be overwritten
    if Tb_free:
        out_file = os.path.join(fc_out_dir, "sigma8_fixedshape_Tbfree.json")
    else:
        out_file = os.path.join(fc_out_dir, "sigma8_fixedshape_Tbfixed.json")
    if os.path.isfile(out_file) and not overwrite:
        return

    # Load previously-computed forecast information
    fc = pf.load_forecast(
        name, verbose=verbose, out_root_dir=out_root_dir, use_mpi=False
    )

    # Compute sigma_8 forecast
    if Tb_free:
        output = pf.compute_forecast_sigma8_fixed_shape(fc, marginalize_over_Tb=True)
    else:
        output = pf.compute_forecast_sigma8_fixed_shape(fc, marginalize_over_Tb=False)

    # Organize forecast output in dict
    data = {
        "onesigma_sigma8": output["onesigma_sigma8"].tolist(),
        "z_center": fc["Centers of redshift bins"],
        "z_edge": fc["Edges of redshift bins"],
        "F_z": output["F_z"].tolist(),
    }

    # Save outputs to file
    with open(out_file, "w") as f:
        json.dump(data, f, indent=2)


def compute_sigma8_fullshape_forecasts(
    fc_out_dir,
    name,
    out_root_dir="output",
    verbose=False,
    Tb_free=True,
    overwrite=False,
):
    """Compute sigma_8 forecasts where all cosmological parameters are free.

    Parameters
    ----------
    fc_out_dir : string
        Directory for forecast output file.
    name: string
        Identifier for previously-computed forecasting information.
    out_root_dir: string, optional
        Root output directory that previous forecasting info was saved to.
    verbose : bool, optional
        Print status messages when loading previous forecasting info. Default: False.
    Tb_free : bool, optional
        Whether to marginalize over mean brightness temparature T_b. Default: True.
    overwrite : bool, optional
        Overwrite output file if it already exists. Default: False.
    """
    # Exit if output file exists and should not be overwritten
    if Tb_free:
        out_file = os.path.join(fc_out_dir, "sigma8_fullshape_Tbfree.json")
    else:
        out_file = os.path.join(fc_out_dir, "sigma8_fullshape_Tbfixed.json")
    if os.path.isfile(out_file) and not overwrite:
        return

    # Load previously-computed forecast information
    fc = pf.load_forecast(
        name, verbose=verbose, out_root_dir=out_root_dir, use_mpi=False
    )

    # Compute sigma_8 forecast
    output = pf.compute_forecast_sigma8_full_cosmology(fc, marginalize_over_Tb=Tb_free)

    # Organize forecast output in dict
    data = {
        "onesigma_sigma8": output["onesigma_sigma8"].tolist(),
        "z_center": fc["Centers of redshift bins"],
        "z_edge": fc["Edges of redshift bins"],
        "F_z": output["F_z"].tolist(),
    }

    # Save outputs to file
    with open(out_file, "w") as f:
        json.dump(data, f, indent=2)


def compute_cos_parameter_forecasts(
    fc_out_dir,
    name,
    out_root_dir="output",
    verbose=False,
    pars="LCDM",
    Tb_free=True,
    omb_bbn_prior=0,
    lensing=None,
    cmb=None,
    overwrite=False,
    exclude_LSS=False,
):
    """Compute cosmological parameter forecasts.

    Parameters
    ----------
    fc_out_dir : string
        Directory for forecast output file.
    name: string
        Identifier for previously-computed forecasting information.
    out_root_dir: string, optional
        Root output directory that previous forecasting info was saved to.
    verbose : bool, optional
        Print status messages when loading previous forecasting info. Default: False.
    pars : string, optional
        Type of parameter set to consider. Must be one of
        ["LCDM", "Mnu", "Neff", "OmegaK"]. Default: "LCDM".
    Tb_free : bool, optional
        Whether to marginalize over mean brightness temparature T_b. Default: True.
    omb_bbn_prior : float, optional
        Gaussian prior to put on omega_b (e.g. from BBN). Default: 0.
    lensing: string, optional
        Identifier for CMB survey for which to include lensing information.
        Default: None.
    cmb: string, optional
        Identifier for CMB survey for which to include primary CMB information.
        Default: None.
    overwrite : bool, optional
        Overwrite output file if it already exists. Default: False.
    exclude_lss : bool, optional
        Whether to exclude LSS from forecast (i.e. perform CMB-only forecast).
        Default: False.
    """
    if pars not in ["LCDM", "Mnu", "Neff", "OmegaK"]:
        raise InputError(f"{pars} not a valid choice for pars argument")

    if pars == "LCDM":
        out_label = "LCDM"
    elif pars == "Mnu":
        out_label = "LCDMplusMnu"
    elif pars == "Neff":
        out_label = "LCDMplusNeff"
    elif pars == "OmegaK":
        out_label = "LCDMplusOmegaK"

    if Tb_free:
        Tb_label = "Tbfree"
    else:
        Tb_label = "Tbfixed"

    if cmb is None:
        cmb_label = "noCMB"
    else:
        cmb_label = f"{cmb}CMB"

    if lensing is None:
        lensing_label = "noLensing"
    else:
        lensing_label = f"{lensing}Lensing"

    if exclude_LSS:
        LSS_label = "_noLSS"
    else:
        LSS_label = ""

    # Exit if output file exists and should not be overwritten
    out_file = os.path.join(
        fc_out_dir,
        f"{out_label}_parameters_{Tb_label}_{cmb_label}_{lensing_label}{LSS_label}.json",
    )
    if os.path.isfile(out_file) and not overwrite:
        return

    # Load previously-computed forecast information
    fc = pf.load_forecast(
        name, verbose=verbose, out_root_dir=out_root_dir, use_mpi=False
    )

    # Compute parameter forecasts
    if pars == "LCDM":
        output = pf.compute_forecast_6parLCDM(
            fc,
            marginalize_over_Tb=Tb_free,
            omb_bbn_prior=omb_bbn_prior,
            cmb=cmb,
            lensing=lensing,
            exclude_LSS=exclude_LSS,
        )
        out_label = "LCDM"
    elif pars == "Mnu":
        output = pf.compute_forecast_Mnu(
            fc,
            marginalize_over_Tb=Tb_free,
            omb_bbn_prior=omb_bbn_prior,
            cmb=cmb,
            lensing=lensing,
            exclude_LSS=exclude_LSS,
        )
        out_label = "LCDMplusMnu"
    elif pars == "Neff":
        output = pf.compute_forecast_Neff(
            fc,
            marginalize_over_Tb=Tb_free,
            omb_bbn_prior=omb_bbn_prior,
            cmb=cmb,
            lensing=lensing,
            exclude_LSS=exclude_LSS,
        )
        out_label = "LCDMplusNeff"
    elif pars == "OmegaK":
        output = pf.compute_forecast_OmegaK(
            fc,
            marginalize_over_Tb=Tb_free,
            omb_bbn_prior=omb_bbn_prior,
            cmb=cmb,
            lensing=lensing,
            exclude_LSS=exclude_LSS,
        )
        out_label = "LCDMplusOmegaK"

    # Organize forecast output in dict
    data = {}
    for key in output.keys():
        if "onesigma" in key:
            data[key] = output[key].tolist()
    data["F_total"] = output["F_total"].tolist()
    data["params_fid"] = fc["forecast"].params_fid

    # Save outputs to file
    with open(out_file, "w") as f:
        json.dump(data, f, indent=2)


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()

    # Define positional command-line arguments
    parser.add_argument(
        "in_root_dir",
        help="Root directory containing pre-computed survey inputs",
    )
    parser.add_argument(
        "out_root_dir",
        help="Root directory for forecast outputs",
    )
    parser.add_argument(
        "survey_name",
        help="Name of survey",
    )
    parser.add_argument(
        "task",
        choices=[
            "bao",
            "sigma8_fixedshape",
            "sigma8_fullshape",
            "LCDM",
            "Mnu",
            "Neff",
            "OmegaK",
        ],
        help="Task: BAO; fixed-shape sigma8; full-cosmology sigma8; LCDM parameters; "
        "LCDM+Mnu parameters; LCDM+Neff parameters; LCDM+OmegaK parameters",
    )

    # Define optional command-line arguments
    parser.add_argument("--verbose", action="store_true", help="Print status messages")
    parser.add_argument(
        "--omb_bbn_prior",
        type=float,
        default=0,
        help=("Prior on omega_b from BBN"),
    )
    parser.add_argument(
        "--cmb",
        help="Label for CMB priors to use",
    )
    parser.add_argument(
        "--lensing",
        help="Label for CMB lensing priors to use",
    )
    parser.add_argument(
        "--Tb_fixed",
        action="store_true",
        help="Whether to fix Tb in the forecasts",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Whether to overwrite output files if they exist already",
    )
    parser.add_argument(
        "--exclude_LSS",
        action="store_true",
        help="Whether to exclude LSS survey from forecasts (i.e. only use CMB and/or lensing)",
    )

    # Parse arguments
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    Tb_free = not args.Tb_fixed

    if args.task == "bao":
        compute_bao_forecasts(
            args.in_root_dir,
            args.survey_name,
            out_root_dir=args.out_root_dir,
            verbose=args.verbose,
            overwrite=args.overwrite,
        )
    elif args.task == "sigma8_fixedshape":
        compute_sigma8_fixedshape_forecasts(
            args.in_root_dir,
            args.survey_name,
            out_root_dir=args.out_root_dir,
            verbose=args.verbose,
            Tb_free=Tb_free,
            overwrite=args.overwrite,
        )
    elif args.task == "sigma8_fullshape":
        compute_sigma8_fullshape_forecasts(
            args.in_root_dir,
            args.survey_name,
            out_root_dir=args.out_root_dir,
            verbose=args.verbose,
            Tb_free=Tb_free,
            overwrite=args.overwrite,
        )
    elif args.task in ["LCDM", "Mnu", "Neff", "OmegaK"]:
        compute_cos_parameter_forecasts(
            args.in_root_dir,
            args.survey_name,
            out_root_dir=args.out_root_dir,
            verbose=args.verbose,
            pars=args.task,
            Tb_free=Tb_free,
            omb_bbn_prior=args.omb_bbn_prior,
            lensing=args.lensing,
            cmb=args.cmb,
            overwrite=args.overwrite,
            exclude_LSS=args.exclude_LSS,
        )
