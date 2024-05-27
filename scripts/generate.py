#!/usr/bin/env python
#
# Script to generate inputs for a survey.
# Based on version by Noah Sailer, but modified by Simon Foreman

import sys
from math import ceil
import argparse

from caput import mpiutil

from FishLSS.headers import *
from FishLSS import twoPoint, twoPointNoise


# Set the default cosmological parameters.
# We take the "Plik best fit" column of Table 1 of 1807.06209.
# Note that with 2 massive neutrinos (N_ncdm=2), N_ur=1.0196 is equivalent to
# N_eff=3.046.
default_cosmo = {
    "A_s": 2.10e-9,
    "n_s": 0.966,
    "alpha_s": 0.0,
    "h": 0.673,
    "N_ur": 1.0196,
    "N_ncdm": 2,
    "m_ncdm": "0.01,0.05",
    "tau_reio": 0.543,
    "omega_b": 0.0224,
    "omega_cdm": 0.120,
    "Omega_k": 0.0,
    "P_k_max_h/Mpc": 2.0,
    "z_pk": "0.0,6",
}


def make_forecast(
    params,
    survey_name,
    task,
    fsky,
    zmin,
    zmax,
    num_zbins,
    verbose=False,
    ell_min=10,
    ell_max=1000,
    kmax=1.,
    out_root_dir="output",
    deriv_dir=None,
    deriv_Cl_dir=None,
    deriv_recon_dir=None,
    b_file=None,
    n_file=None,
    b2_file=None,
    bs_file=None,
    use_G2_basis=False,
    HI=False,
    sigv=None,
    pessimistic_wedge=False,
    HI_stoch_multiplier=1.0,
    HI_stoch_file=None,
    HI_sampling_file=None,
    Nside=256,
    fill_factor=0.5,
    t_int=5.0,
    dish_diameter=6.0,
    hex_pack=True,
    aperture_efficiency=0.7,
    sky_coupling=0.9,
    omt_coupling=0.9,
    T_ground=300.0,
    T_ampl=50.0,
    knl_z0=None,
    dknl_dz=None,
    dknl_dDinv=None,
    remove_lowk_delta2_powspec=False,
    remove_lowk_delta2_cov=False,
):
    """Precompute and store forecasting information."""

    # Store default cosmological parameter values in params dict
    for k in default_cosmo.keys():
        params[k] = default_cosmo[k]

    # Set k_max to input value
    params["P_k_max_h/Mpc"] = kmax

    # Initialize CLASS
    if verbose:
        print("Initializing CLASS")
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()

    if b_file is not None:
        try:
            b_z_in, b_val_in = np.loadtxt(b_file).T
        except:
            raise RuntimeError(f"Error reading from bias file {b_file}")

        b = interp1d(b_z_in, b_val_in, bounds_error=True)
    else:
        b = None

    if b2_file is not None:
        try:
            b2_z_in, b2_val_in = np.loadtxt(b2_file).T
        except:
            raise RuntimeError(f"Error reading from b2 file {b2_file}")

        b2 = interp1d(b2_z_in, b2_val_in, bounds_error=True)
    else:
        b2 = None

    if bs_file is not None:
        try:
            bs_z_in, bs_val_in = np.loadtxt(bs_file).T
        except:
            raise RuntimeError(f"Error reading from bs file {bs_file}")

        bs = interp1d(bs_z_in, bs_val_in, bounds_error=True)
    else:
        bs = None

    if n_file is not None:
        try:
            n_z_in, n_val_in = np.loadtxt(n_file).T
        except:
            raise RuntimeError(f"Error reading from nbar file {n_file}")

        n = interp1d(n_z_in, n_val_in, bounds_error=True)
    else:
        n = None

    if use_G2_basis:
        if b is None or b2 is None or bs is None:
            raise InputError("Must specify all of b, b2, bs to use G2 basis")

        b2_new_val = 2 * b2(b2_z_in) - (4 / 3) * bs(b2_z_in) + (8 / 21) * b(b2_z_in)
        bs_new_val = bs(bs_z_in) - (2 / 7) * b(bs_z_in)

        b2 = interp1d(b2_z_in, b2_new_val, bounds_error=True)
        bs = interp1d(bs_z_in, bs_new_val, bounds_error=True)

    # Set fiducial velocity dispersion responsible for Fingers of God, in km/s
    if sigv is None:
        if HI:
            sigv = 10
        else:
            sigv = 100

    zedges = np.linspace(zmin, zmax, num_zbins + 1, endpoint=True)

    # Set up the experiment object
    if verbose:
        print("Initializing experiment")
    exp = experiment(
        zedges=zedges,
        fsky=fsky,
        b=b,
        n=n,
        b2=b2,
        bs=bs,
        sigv=sigv,
        HI=HI,
        pessimistic=pessimistic_wedge,
        HI_stoch_file=HI_stoch_file,
        HI_stoch_multiplier=HI_stoch_multiplier,
        HI_sampling_file=HI_sampling_file,
        Ndetectors=Nside**2.0,
        fill_factor=fill_factor,
        tint=t_int,
        D=dish_diameter,
        hex_pack=hex_pack,
        aperture_efficiency=aperture_efficiency,
        sky_coupling=sky_coupling,
        omt_coupling=omt_coupling,
        T_ground=T_ground,
        T_ampl=T_ampl,
        knl_z0=knl_z0,
        dknl_dz=dknl_dz,
        dknl_dDinv=dknl_dDinv,
    )

    # Generate the fisherForecast object with directory set by
    # the survey file basename and properties determined by the text file
    setup = task == "setup"
    forecast = fisherForecast(
        experiment=exp,
        cosmo=cosmo,
        name=survey_name,
        setup=setup,
        verbose=verbose,
        kmax=kmax,
        ell=np.arange(ell_min, ell_max, 1),
        remove_lowk_delta2_powspec=remove_lowk_delta2_powspec,
        remove_lowk_delta2_cov=remove_lowk_delta2_cov,
        output_root_directory=out_root_dir,
        deriv_directory=deriv_dir,
        deriv_Cl_directory=deriv_Cl_dir,
        deriv_recon_directory=deriv_recon_dir,
    )

    return forecast


def do_task(
    task,
    survey_name,
    fsky,
    zmin,
    zmax,
    num_zbins,
    verbose=False,
    ell_min=10,
    ell_max=1000,
    kmax=1.,
    b_file=None,
    n_file=None,
    b2_file=None,
    bs_file=None,
    use_G2_basis=False,
    HI=False,
    sigv=None,
    pessimistic_wedge=False,
    HI_stoch_multiplier=1.0,
    HI_stoch_file=None,
    HI_sampling_file=None,
    Nside=256,
    fill_factor=0.5,
    t_int=5.0,
    dish_diameter=6.0,
    hex_pack=True,
    aperture_efficiency=0.7,
    sky_coupling=0.9,
    omt_coupling=0.9,
    T_ground=300.0,
    T_ampl=50.0,
    knl_z0=None,
    dknl_dz=None,
    dknl_dDinv=None,
    remove_lowk_delta2_powspec=False,
    remove_lowk_delta2_cov=False,
    out_root_dir="output",
    deriv_dir=None,
    deriv_Cl_dir=None,
    deriv_recon_dir=None,
):
    """Does the work, performing task 'task' on survey file base name 'sfb'.

    When taking derivatives of P(k,mu) you don't need to get the lensing
    C_ell's from class. So to speed things up we'll use a different CLASS
    object depending on the derivatives being calculated.
    """

    kwargs = {
        "verbose": verbose,
        "ell_min": ell_min,
        "ell_max": ell_max,
        "kmax": kmax,
        "b_file": b_file,
        "n_file": n_file,
        "b2_file": b2_file,
        "bs_file": bs_file,
        "use_G2_basis": use_G2_basis,
        "HI": HI,
        "sigv": sigv,
        "pessimistic_wedge": pessimistic_wedge,
        "HI_stoch_multiplier": HI_stoch_multiplier,
        "HI_stoch_file": HI_stoch_file,
        "HI_sampling_file": HI_sampling_file,
        "Nside": Nside,
        "fill_factor": fill_factor,
        "t_int": t_int,
        "dish_diameter": dish_diameter,
        "hex_pack": hex_pack,
        "aperture_efficiency": aperture_efficiency,
        "sky_coupling": sky_coupling,
        "omt_coupling": omt_coupling,
        "T_ground": T_ground,
        "T_ampl": T_ampl,
        "knl_z0": knl_z0,
        "dknl_dz": dknl_dz,
        "dknl_dDinv": dknl_dDinv,
        "remove_lowk_delta2_powspec": remove_lowk_delta2_powspec,
        "remove_lowk_delta2_cov": remove_lowk_delta2_cov,
        "out_root_dir": out_root_dir,
        "deriv_dir": deriv_dir,
        "deriv_Cl_dir": deriv_Cl_dir,
        "deriv_recon_dir": deriv_recon_dir,
    }

    if task == "setup":
        # Set outputs to unlensed and lensed C_l, along with matter P(k) and
        # all cosmological parameters
        params = {
            "output": "tCl lCl mPk",
            "non linear": "halofit",
            "l_max_scalars": 1000,
            "lensing": "yes",
        }
        forecast = make_forecast(
            params,
            survey_name,
            task,
            fsky,
            zmin,
            zmax,
            num_zbins,
            **kwargs,
        )

    elif task == "rec":
        # Set outputs to matter P(k) and all cosmological parameters
        params = {"output": "mPk"}
        forecast = make_forecast(
            params,
            survey_name,
            task,
            fsky,
            zmin,
            zmax,
            num_zbins,
            **kwargs,
        )

        # Compute derivatives with respect to BAO parameters and linear/quadratic biases,
        # with BAO reconstruction turned on
        basis = np.array(["alpha_perp", "alpha_parallel", "b"])
        forecast.recon = True
        forecast.marg_params = basis
        forecast.compute_derivatives(
            five_point=False,
            verbose=verbose,
            n_partitions=mpiutil.size,
            partition=mpiutil.rank,
        )
        forecast.recon = False

    elif task == "fs":
        # Set outputs to matter P(k) and all cosmological parameters
        params = {"output": "mPk"}
        forecast = make_forecast(
            params,
            survey_name,
            task,
            fsky,
            zmin,
            zmax,
            num_zbins,
            **kwargs,
        )

        # Compute derivatives for full-shape parameters
        basis = np.array(
            [
                "h",
                "log(A_s)",
                "n_s",
                "omega_cdm",
                "omega_b",
                "tau_reio",
                "m_ncdm",
                "f_NL",
                "N_ur",
                "Omega_k",
                "N",
                "alpha0",
                "b",
                "b2",
                "bs",
                "N2",
                "N4",
                "alpha2",
                "alpha4",
                "Tb",
            ]
        )
        forecast.marg_params = basis
        forecast.compute_derivatives(
            five_point=False,
            verbose=verbose,
            n_partitions=mpiutil.size,
            partition=mpiutil.rank,
        )

    elif task == "lens":
        # Set outputs to unlensed and lensed C_l, along with matter P(k) and
        # all cosmological parameters
        params = {
            "output": "tCl lCl mPk",
            "l_max_scalars": 1000,
            "lensing": "yes",
            "non linear": "halofit",
        }

        forecast = make_forecast(
            params,
            survey_name,
            task,
            fsky,
            zmin,
            zmax,
            num_zbins,
            **kwargs,
        )

        # Compute derivatives for parameters entering lensing spectra
        basis = np.array(
            [
                "h",
                "log(A_s)",
                "n_s",
                "omega_cdm",
                "omega_b",
                "tau_reio",
                "m_ncdm",
                # "f_NL",
                "N_ur",
                "Omega_k",
            ]
        )
        forecast.marg_params = basis
        forecast.compute_Cl_derivatives(
            five_point=False,
            verbose=verbose,
            n_partitions=mpiutil.size,
            partition=mpiutil.rank,
            only_kk=True,
        )


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()

    # Define positional command-line arguments
    parser.add_argument(
        "survey_name",
        help="Name of survey, corresponding to existing [survey_name].txt file",
    )
    parser.add_argument(
        "task",
        choices=["setup", "rec", "fs", "lens"],
        help="Task: setup; derivatives for reconstruction; derivatives for full-shape; "
        "derivatives for lensing",
    )
    parser.add_argument("zmin", type=float, help="Minimum redshift of survey")
    parser.add_argument("zmax", type=float, help="Maximum redshift of survey")
    parser.add_argument("num_zbins", type=int, help="Number of redshift bins")

    # Define optional command-line arguments
    parser.add_argument("--verbose", action="store_true", help="Print status messages")
    parser.add_argument(
        "--ellmin", type=int, default=10, help="Minimum ell for computing lensing spectra"
    )
    parser.add_argument(
        "--ellmax", type=int, default=1000, help="Maximum ell for computing lensing spectra"
    )
    parser.add_argument(
        "--kmax", type=float, default=1., help="Maximum k for computing 3d power spectra"
    )
    parser.add_argument(
        "--fsky", type=float, default=0.5, help="Sky fraction covered by survey"
    )
    parser.add_argument(
        "--HI",
        action="store_true",
        help="Whether forecast corresponds to a 21cm survey",
    )
    parser.add_argument("--sigv", type=float, default=None, help="Fiducial FoG velocity dispersion, in km/s")
    parser.add_argument(
        "--pessimistic_wedge",
        action="store_true",
        help="Whether to use Sailer+ prescription for pessimistic 21cm foreground wedge",
    )
    parser.add_argument(
        "--HI_stoch_multiplier",
        type=float,
        default=1.0,
        help="Multiply default HI stochastic noise by this multiplier",
    )
    parser.add_argument(
        "--HI_stoch_file",
        type=str,
        help="HI stochastic noise model file",
    )
    parser.add_argument(
        "--HI_sampling_file",
        type=str,
        help="HI sampling noise (i.e. physical Poisson noise) model file",
    )
    parser.add_argument(
        "--Nside",
        type=int,
        default=256,
        help=(
            "Number of antennas per side of array, "
            "assuming square array with Nside^2 detectors"
        ),
    )
    parser.add_argument(
        "--fill_factor", type=float, default=0.5, help="Array fill factor"
    )
    parser.add_argument(
        "--t_int", type=float, default=5.0, help="Observing time in years"
    )
    parser.add_argument(
        "--dish_diameter", type=float, default=6.0, help="Dish diameter in meters"
    )
    parser.add_argument(
        "--hex_pack",
        type=bool,
        default=True,
        help=("Hex-packed array if True, square-packed array if False"),
    )
    parser.add_argument(
        "--aperture_efficiency",
        type=float,
        default=0.7,
        help="Dish aperture efficiency",
    )
    parser.add_argument(
        "--sky_coupling",
        type=float,
        default=0.9,
        help="Fraction of the primary beam coupled to the sky, vs. the ground",
    )
    parser.add_argument("--omt_coupling", type=float, default=0.9, help="OMT coupling")
    parser.add_argument(
        "--T_ground", type=float, default=300.0, help="Ground temperature in K"
    )
    parser.add_argument(
        "--T_ampl", type=float, default=0.7, help="Amplifier temperature in K"
    )
    parser.add_argument(
        "--knl_z0", type=float, default=None, help="k_nl at z=0, in h/Mpc"
    )
    parser.add_argument("--dknl_dz", type=float, default=None, help="d(k_nl)/dz")
    parser.add_argument(
        "--dknl_dDinv",
        type=float,
        default=None,
        help="d(k_nl)/dDinv, where D is the linear growth factor",
    )
    parser.add_argument(
        "--remove_lowk_delta2_powspec",
        action="store_true",
        help=(
            "Subtract constant contributions to b_2^2, b_s*b_2, and b_s^2 terms in "
            "power spectra"
        ),
    )
    parser.add_argument(
        "--remove_lowk_delta2_cov",
        action="store_true",
        help=(
            "Subtract constant contributions to b_2^2, b_s*b_2, and b_s^2 terms in "
            "computation of covariance"
        ),
    )
    parser.add_argument(
        "--output_root_directory",
        default="output",
        help="Root directory to store forecast output",
    )
    parser.add_argument(
        "--deriv_directory",
        type=str,
        help="Directory from which to read 3d power spectrum derivatives",
    )
    parser.add_argument(
        "--deriv_recon_directory",
        type=str,
        help="Directory from which to read reconstructed 3d power spectrum derivatives",
    )
    parser.add_argument(
        "--deriv_Cl_directory",
        type=str,
        help="Directory from which to read angular power spectrum derivatives",
    )
    parser.add_argument(
        "--b_file",
        default=None,
        help="Text file containing linear bias as a function of z",
    )
    parser.add_argument(
        "--n_file",
        default=None,
        help="Text file containing number density as a function of z",
    )
    parser.add_argument(
        "--b2_file",
        default=None,
        help="Text file containing b2 as a function of z",
    )
    parser.add_argument(
        "--bs_file",
        default=None,
        help="Text file containing bs as a function of z",
    )
    parser.add_argument(
        "--use_G2_basis",
        action="store_true",
        help="Treat input b_s(z) values as corresponding to b_G2(z), and internally "
        "translate b_2, b_s values to correspond to this basis",
    )

    # Parse arguments
    args = parser.parse_args()

    return args


if __name__ == "__main__":

    # Parse command-line arguments
    args = parse_arguments()

    # Check for MPI run
    if mpiutil.size > 1 and args.verbose:
        print(f"Starting MPI process {mpiutil.rank}")

    # Do the work
    do_task(
        args.task,
        args.survey_name,
        args.fsky,
        args.zmin,
        args.zmax,
        args.num_zbins,
        verbose=args.verbose,
        ell_min=args.ellmin,
        ell_max=args.ellmax,
        kmax=args.kmax,
        b_file=args.b_file,
        n_file=args.n_file,
        b2_file=args.b2_file,
        bs_file=args.bs_file,
        use_G2_basis=args.use_G2_basis,
        HI=args.HI,
        sigv=args.sigv,
        pessimistic_wedge=args.pessimistic_wedge,
        HI_stoch_multiplier=args.HI_stoch_multiplier,
        HI_stoch_file=args.HI_stoch_file,
        HI_sampling_file=args.HI_sampling_file,
        Nside=args.Nside,
        fill_factor=args.fill_factor,
        t_int=args.t_int,
        dish_diameter=args.dish_diameter,
        hex_pack=args.hex_pack,
        aperture_efficiency=args.aperture_efficiency,
        sky_coupling=args.sky_coupling,
        omt_coupling=args.omt_coupling,
        T_ground=args.T_ground,
        T_ampl=args.T_ampl,
        knl_z0=args.knl_z0,
        dknl_dz=args.dknl_dz,
        dknl_dDinv=args.dknl_dDinv,
        remove_lowk_delta2_powspec=args.remove_lowk_delta2_powspec,
        remove_lowk_delta2_cov=args.remove_lowk_delta2_cov,
        out_root_dir=args.output_root_directory,
        deriv_dir=args.deriv_directory,
        deriv_Cl_dir=args.deriv_Cl_directory,
        deriv_recon_dir=args.deriv_recon_directory,
    )
