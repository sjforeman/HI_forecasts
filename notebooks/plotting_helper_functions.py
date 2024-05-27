#!/usr/bin/env python
#
# Helper functions for plotting forecasts.

import sys, os
from math import ceil
import argparse

import matplotlib.ticker as mticker

from FishLSS.headers import *
from FishLSS.twoPoint import *
from FishLSS.twoPointNoise import *
# import FishLSS.parameter_forecast as pf


# Get default matplotlib color cycle, so that we can repeat colors easily in plots
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]

# Linestyles for optimistic and pessimistic foregrounds in "shaded region" plots
_LS_OPT = '-'
_LS_PESS = '--'


def fill_stairs(ax, edges, lower, upper, c, alpha, label=None):
    """Analog of fill_between, but for staircase plots."""
    lower_plot = np.concatenate([lower, [lower[-1]]])
    upper_plot = np.concatenate([upper, [upper[-1]]])
    ax.fill_between(edges, lower_plot, upper_plot, step="post", alpha=alpha, color=c)

    if label is not None:
        ax.fill_between([], [], [], color=c, label=label, alpha=alpha)


def plot_AP(
    rec_BAO_output,
    noise,
    p,
    abs_label,
    ratio_ylabel,
    knl_keys,
    knl_titles,
    ylim_abs=[0, 5],
    ylim_ratio=[0.6, 1.01],
    lw=1.5,
    fill_alpha=0.3,
    abs_dtick=None,
    ratio_dtick=None,
    legend_loc="upper left",
    legend2_loc="upper left",
    plot_chord=True,
):
    """Plot forecasts for Alcock-Paczynski parameters, for different kNL assumptions."""
    
    n_knl = len(knl_keys)

    fig, ax = plt.subplots(
        2, n_knl, figsize=(9, 3.7), sharex=True, sharey=False, height_ratios=[1, 0.5]
    )

    for knli, knl_option in enumerate(knl_keys):

        for s, c, label in zip(
            ["puma_32k", "puma_5k", "chord"],
            [colors[0], colors[1], colors[2]],
            ["PUMA-32k", "PUMA-5k", "CHORD"],
        ):

            if s == "chord" and not plot_chord:
                continue

            edges = rec_BAO_output[noise, f"{s}_opt", f"no_sub_{knl_option}"]["z_edge"]
            opt = 100 * np.array(rec_BAO_output[noise, f"{s}_opt", f"no_sub_{knl_option}"][p])
            pess = 100 * np.array(
                rec_BAO_output[noise, f"{s}_pess", f"no_sub_{knl_option}"][p]
            )

            opt_sub = 100 * np.array(rec_BAO_output[noise, f"{s}_opt", f"sub_{knl_option}"][p])
            pess_sub = 100 * np.array(
                rec_BAO_output[noise, f"{s}_pess", f"sub_{knl_option}"][p]
            )

            ax[0, knli].stairs(pess, edges, lw=lw, color=c, ls=_LS_PESS, baseline=None)
            ax[0, knli].stairs(opt, edges, lw=lw, color=c, ls=_LS_OPT, baseline=None)
            fill_stairs(ax[0, knli], edges, opt, pess, c, alpha=fill_alpha, label=label)

            ax[1, knli].stairs(
                pess_sub / pess, edges, lw=lw, color=c, ls=_LS_PESS, baseline=None
            )
            ax[1, knli].stairs(
                opt_sub / opt, edges, lw=lw, color=c, ls=_LS_OPT, baseline=None
            )
            fill_stairs(
                ax[1, knli],
                edges,
                opt_sub / opt,
                pess_sub / pess,
                c,
                alpha=fill_alpha,
                label=label,
            )

        ax[0, knli].set_ylim(ylim_abs)
        if abs_dtick is not None:
            ax[0, knli].set_yticks(
                np.arange(ylim_abs[0], ylim_abs[1] + 1e-10, abs_dtick)
            )
        ax[0, knli].set_xlim(0, 6)
        ax[0, knli].set_ylabel(abs_label)
        ax[0, knli].grid(ls="--")

        ax[1, knli].set_ylim(ylim_ratio)
        if ratio_dtick is not None:
            ax[1, knli].set_yticks(
                np.arange(ylim_ratio[0], ylim_ratio[1] + 1e-10, ratio_dtick)
            )
        ax[1, knli].set_xlim(0, 6)
        ax[1, knli].set_ylabel(ratio_ylabel)
        ax[1, knli].grid(ls="--")
        ax[1, knli].set_xlabel(r"$z$")

    for i, title in enumerate(knl_titles):
        ax[0, i].set_title(title)

    ax[0, 0].legend(loc=legend_loc)
    
    line1, = ax[0, 1].plot([0, 1], [-1e10, -1e10], c='k', ls=_LS_OPT, label="Optimistic foregrounds")
    line2, = ax[0, 1].plot([0, 1], [-1e10, -1e10], c='k', ls=_LS_PESS, label="Pessimistic foregrounds")
    ax[0, 1].legend(handles=[line1, line2], loc=legend2_loc, fontsize=9)

    plt.tight_layout()


def plot_AP_singlekNL(
    rec_BAO_output,
    noise,
    ps,
    abs_labels,
    ratio_ylabels,
    knl_option,
    knl_title,
    ylim_abs=[0, 5],
    ylim_ratio=[0.6, 1.01],
    lw=1.5,
    fill_alpha=0.3,
    abs_dtick=None,
    ratio_dtick=None,
    legend_loc="upper left",
    legend2_loc="upper left",
    plot_chord=True,
):
    """Plot forecasts for Alcock-Paczynski parameters, for a single kNL assumption."""

    n_p = len(ps)

    fig, ax = plt.subplots(
        2, n_p, figsize=(9, 3.7), sharex=True, sharey=False, height_ratios=[1, 0.5]
    )

    # for knli, knl_option in enumerate(knl_keys):

    for pi, (p, abs_label, ratio_ylabel) in enumerate(
        zip(ps, abs_labels, ratio_ylabels)
    ):

        for s, c, label in zip(
            ["puma_32k", "puma_5k", "chord"],
            [colors[0], colors[1], colors[2]],
            ["PUMA-32k", "PUMA-5k", "CHORD"],
        ):

            if s == "chord" and not plot_chord:
                continue

            edges = rec_BAO_output[noise, f"{s}_opt", f"no_sub_{knl_option}"]["z_edge"]
            opt = 100 * np.array(rec_BAO_output[noise, f"{s}_opt", f"no_sub_{knl_option}"][p])
            pess = 100 * np.array(
                rec_BAO_output[noise, f"{s}_pess", f"no_sub_{knl_option}"][p]
            )

            opt_sub = 100 * np.array(rec_BAO_output[noise, f"{s}_opt", f"sub_{knl_option}"][p])
            pess_sub = 100 * np.array(
                rec_BAO_output[noise, f"{s}_pess", f"sub_{knl_option}"][p]
            )

            ax[0, pi].stairs(pess, edges, lw=lw, color=c, ls=_LS_PESS, baseline=None)
            ax[0, pi].stairs(opt, edges, lw=lw, color=c, ls=_LS_OPT, baseline=None)
            fill_stairs(ax[0, pi], edges, opt, pess, c, alpha=fill_alpha, label=label)

            ax[1, pi].stairs(
                pess_sub / pess, edges, lw=lw, color=c, ls=_LS_PESS, baseline=None
            )
            ax[1, pi].stairs(opt_sub / opt, edges, lw=lw, color=c, ls=_LS_OPT, baseline=None)
            fill_stairs(
                ax[1, pi],
                edges,
                opt_sub / opt,
                pess_sub / pess,
                c,
                alpha=fill_alpha,
                label=label,
            )

        ax[0, pi].set_ylim(ylim_abs)
        if abs_dtick is not None:
            ax[0, pi].set_yticks(np.arange(ylim_abs[0], ylim_abs[1] + 1e-10, abs_dtick))
        ax[0, pi].set_xlim(0, 6)
        ax[0, pi].set_ylabel(abs_label)
        ax[0, pi].grid(ls="--")

        ax[1, pi].set_ylim(ylim_ratio)
        if ratio_dtick is not None:
            ax[1, pi].set_yticks(
                np.arange(ylim_ratio[0], ylim_ratio[1] + 1e-10, ratio_dtick)
            )
        ax[1, pi].set_xlim(0, 6)
        ax[1, pi].set_ylabel(ratio_ylabel)
        ax[1, pi].grid(ls="--")
        ax[1, pi].set_xlabel(r"$z$")

    # for i, title in enumerate(knl_titles):
    #     ax[0, i].set_title(title)

    ax[0, 1].legend(loc=legend_loc)
    
    line1, = ax[0, 0].plot([0, 1], [-1e10, -1e10], c='k', ls=_LS_OPT, label="Optimistic foregrounds")
    line2, = ax[0, 0].plot([0, 1], [-1e10, -1e10], c='k', ls=_LS_PESS, label="Pessimistic foregrounds")
    ax[0, 0].legend(handles=[line1, line2], loc=legend2_loc, fontsize=9)

    plt.suptitle(knl_title, y=0.93)

    plt.tight_layout()


def plot_sigma8(
    results_dict,
    noise,
    knl_keys,
    knl_titles,
    ylim_abs=[0, 5],
    ylim_ratio=[0.6, 1.01],
    lw=1.5,
    fill_alpha=0.3,
    abs_dtick=None,
    ratio_dtick=None,
    legend_loc="upper left",
    legend2_loc="upper left",
    figsize=(9, 3.7),
    Tb="freeTb",
    plot_chord=True,
):
    """Plot forecasts for sigma_8."""

    p = "onesigma_sigma8"

    n_knl = len(knl_keys)

    fig, ax = plt.subplots(
        2, n_knl, figsize=figsize, sharex=True, sharey=False, height_ratios=[1, 0.5]
    )

    for knli, knl_option in enumerate(knl_keys):

        for s, c, label in zip(
            ["puma_32k", "puma_5k", "chord"],
            [colors[0], colors[1], colors[2]],
            ["PUMA-32k", "PUMA-5k", "CHORD"],
        ):

            if s == "chord" and not plot_chord:
                continue

            edges = results_dict[noise, f"{s}_opt", f"no_sub_{knl_option}", Tb]["z_edge"]
            opt = 100 * np.array(
                results_dict[noise, f"{s}_opt", f"no_sub_{knl_option}", Tb][p]
            )
            pess = 100 * np.array(
                results_dict[noise, f"{s}_pess", f"no_sub_{knl_option}", Tb][p]
            )

            opt_sub = 100 * np.array(
                results_dict[noise, f"{s}_opt", f"sub_{knl_option}", Tb][p]
            )
            pess_sub = 100 * np.array(
                results_dict[noise, f"{s}_pess", f"sub_{knl_option}", Tb][p]
            )

            ax[0, knli].stairs(pess, edges, lw=lw, color=c, ls=_LS_PESS, baseline=None)
            ax[0, knli].stairs(opt, edges, lw=lw, color=c, ls=_LS_OPT, baseline=None)
            fill_stairs(ax[0, knli], edges, opt, pess, c, alpha=fill_alpha, label=label)

            ax[1, knli].stairs(
                pess_sub / pess, edges, lw=lw, color=c, ls=_LS_PESS, baseline=None
            )
            ax[1, knli].stairs(
                opt_sub / opt, edges, lw=lw, color=c, ls=_LS_OPT, baseline=None
            )
            fill_stairs(
                ax[1, knli],
                edges,
                opt_sub / opt,
                pess_sub / pess,
                c,
                alpha=fill_alpha,
                label=label,
            )

        ax[0, knli].set_ylim(ylim_abs)
        if abs_dtick is not None:
            ax[0, knli].set_yticks(
                np.arange(ylim_abs[0], ylim_abs[1] + 1e-10, abs_dtick)
            )
        ax[0, knli].set_xlim(0, 6)
        ax[0, knli].set_ylabel(r"$\sigma_{\sigma_8} [\%]$")
        ax[0, knli].grid(ls="--")

        ax[1, knli].set_ylim(ylim_ratio)
        if ratio_dtick is not None:
            ax[1, knli].set_yticks(
                np.arange(ylim_ratio[0], ylim_ratio[1] + 1e-10, ratio_dtick)
            )
        ax[1, knli].set_xlim(0, 6)
        ax[1, knli].set_ylabel(r"$\sigma_{\sigma_8}$ ratio")
        ax[1, knli].grid(ls="--")
        ax[1, knli].set_xlabel(r"$z$")

    for i, title in enumerate(knl_titles):
        ax[0, i].set_title(title)

    ax[0, 0].legend(loc=legend_loc)
    
    line1, = ax[0, 1].plot([0, 1], [-1e10, -1e10], c='k', ls=_LS_OPT, label="Optimistic foregrounds")
    line2, = ax[0, 1].plot([0, 1], [-1e10, -1e10], c='k', ls=_LS_PESS, label="Pessimistic foregrounds")
    ax[0, 1].legend(handles=[line1, line2], loc=legend2_loc, fontsize=9)

    plt.tight_layout()


def plot_sigma8_singlepanel(
    results_dict,
    noise, 
    knl_key,
    knl_title,
    ylim_abs=[0, 5],
    ylim_ratio=[0.6, 1.01],
    lw=1.5,
    fill_alpha=0.3,
    abs_dtick=None,
    ratio_dtick=None,
    legend_loc="upper left",
    legend2_loc="upper left",
    figsize=(9, 3.7),
    Tb="freeTb",
    plot_chord=True,
):
    """Plot forecasts for sigma_8, as single panel."""

    p = "onesigma_sigma8"

    fig, ax = plt.subplots(
        2, 1, figsize=figsize, sharex=True, sharey=False, height_ratios=[1, 0.5]
    )

    knl_option = knl_key

    for s, c, label in zip(
        ["puma_32k", "puma_5k", "chord"],
        [colors[0], colors[1], colors[2]],
        ["PUMA-32k", "PUMA-5k", "CHORD"],
    ):

        if s == "chord" and not plot_chord:
            continue

        edges = results_dict[noise, f"{s}_opt", f"no_sub_{knl_option}", Tb]["z_edge"]
        opt = 100 * np.array(results_dict[noise, f"{s}_opt", f"no_sub_{knl_option}", Tb][p])
        pess = 100 * np.array(results_dict[noise, f"{s}_pess", f"no_sub_{knl_option}", Tb][p])

        opt_sub = 100 * np.array(results_dict[noise, f"{s}_opt", f"sub_{knl_option}", Tb][p])
        pess_sub = 100 * np.array(results_dict[noise, f"{s}_pess", f"sub_{knl_option}", Tb][p])

        ax[0].stairs(pess, edges, lw=lw, color=c, ls=_LS_PESS, baseline=None)
        ax[0].stairs(opt, edges, lw=lw, color=c, ls=_LS_OPT, baseline=None)
        fill_stairs(ax[0], edges, opt, pess, c, alpha=fill_alpha, label=label)

        ax[1].stairs(pess_sub / pess, edges, lw=lw, color=c, ls=_LS_PESS, baseline=None)
        ax[1].stairs(opt_sub / opt, edges, lw=lw, color=c, ls=_LS_OPT, baseline=None)
        fill_stairs(
            ax[1],
            edges,
            opt_sub / opt,
            pess_sub / pess,
            c,
            alpha=fill_alpha,
            label=label,
        )

    ax[0].set_ylim(ylim_abs)
    if abs_dtick is not None:
        ax[0].set_yticks(np.arange(ylim_abs[0], ylim_abs[1] + 1e-10, abs_dtick))
    ax[0].set_xlim(0, 6)
    ax[0].set_ylabel(r"$\sigma_{\sigma_8} [\%]$")
    ax[0].grid(ls="--")

    ax[1].set_ylim(ylim_ratio)
    if ratio_dtick is not None:
        ax[1].set_yticks(np.arange(ylim_ratio[0], ylim_ratio[1] + 1e-10, ratio_dtick))
    ax[1].set_xlim(0, 6)
    ax[1].set_ylabel(r"$\sigma_{\sigma_8}$ ratio")
    ax[1].grid(ls="--")
    ax[1].set_xlabel(r"$z$")

    ax[0].set_title(knl_title)

    ax[0].legend(loc=legend_loc)
    
    line1, = ax[1].plot([0, 1], [-1e10, -1e10], c='k', ls=_LS_OPT, label="Optimistic foregrounds")
    line2, = ax[1].plot([0, 1], [-1e10, -1e10], c='k', ls=_LS_PESS, label="Pessimistic foregrounds")
    ax[1].legend(handles=[line1, line2], loc=legend2_loc, fontsize=9)

    plt.tight_layout()


def plot_LCDM(
    parameter_output,
    noise, 
    knl_keys,
    knl_titles,
    ylim_abs=[1e-2, 1e0],
    ylim_ratio=[0.6, 1.01],
    pess_alpha=0.3,
    abs_dtick=None,
    ratio_dtick=None,
    legend_loc="upper left",
    figsize=(9, 3.5),
    Tb_choice="Tbfree",
    plot_chord=True,
    experiment_suffix = ""
):
    """Plot LCDM parameter forecasts."""

    par_list = [
        "h",
        "log(A_s)",
        "n_s",
        "omega_cdm",
        "omega_b",
        # "m_ncdm",
        # "N_ur",
        # "Omega_k"
    ]
    par_list_symbols = [
        r"$h$",
        r"$\log(A_{\rm s})$",
        r"$n_{\rm s}$",
        r"$\omega_{\rm c}$",
        r"$\omega_{\rm b}$",
        # r"$M_\nu$",
        # r"$N_{\rm eff}$",
        # r"$\Omega_k$"
    ]
    w = 0.17
    
    params_fid = parameter_output[
        noise, 
        "puma_32k_opt",
        f"no_sub_{knl_keys[0]}",
        # "no_sub_dknldDinv0.3",
        "LCDM_parameters_Tbfree_Planck_SOCMB_noLensing",
    ]["params_fid"]

    def _onesigma_par(noise, survey, option, priors, p):
        # In percent of fiducial value, except for Omega_K

        if p == "log(A_s)":
            p_fid = np.abs(np.log(params_fid["A_s"]))
        elif p == "Omega_k":
            p_fid = 1
        elif p == "m_ncdm":
            p_fid = 0.06
        else:
            p_fid = params_fid[p]

        if p in ["h", "log(A_s)", "n_s", "omega_cdm", "omega_b"]:
            return (
                100
                * parameter_output[noise, survey, option, f"LCDM_parameters_{priors}"][
                    f"onesigma_{p}"
                ]
                / p_fid
            )
        elif p == "m_ncdm":
            return (
                100
                * parameter_output[noise, survey, option, f"LCDMplusMnu_parameters_{priors}"][
                    f"onesigma_{p}"
                ]
                / p_fid
            )
        elif p == "N_ur":
            return (
                100
                * parameter_output[noise, survey, option, f"LCDMplusNeff_parameters_{priors}"][
                    f"onesigma_{p}"
                ]
                / p_fid
            )
        elif p == "Omega_k":
            return (
                100
                * parameter_output[
                    noise, survey, option, f"LCDMplusOmegaK_parameters_{priors}"
                ][f"onesigma_{p}"]
                / p_fid
            )

    def _noLSS_onesigma_par(priors, p):
        # In percent of fiducial value, except for Omega_K

        if p == "log(A_s)":
            p_fid = np.abs(np.log(params_fid["A_s"]))
        elif p == "Omega_k":
            p_fid = 1
        elif p == "m_ncdm":
            p_fid = 0.06
        else:
            p_fid = params_fid[p]

        if p in ["h", "log(A_s)", "n_s", "omega_cdm", "omega_b"]:
            return (
                100
                * parameter_output[f"LCDM_parameters_{priors}_noLSS"][f"onesigma_{p}"]
                / p_fid
            )
        elif p == "m_ncdm":
            return (
                100
                * parameter_output[f"LCDMplusMnu_parameters_{priors}_noLSS"][
                    f"onesigma_{p}"
                ]
                / p_fid
            )
        elif p == "N_ur":
            return (
                100
                * parameter_output[f"LCDMplusNeff_parameters_{priors}_noLSS"][
                    f"onesigma_{p}"
                ]
                / p_fid
            )
        elif p == "Omega_k":
            return (
                100
                * parameter_output[f"LCDMplusOmegaK_parameters_{priors}_noLSS"][
                    f"onesigma_{p}"
                ]
                / p_fid
            )

    fig = plt.figure(figsize=figsize, layout="constrained")
    subfig = fig.subfigures(3, 1)

    ax = {}
    for si, (s, label) in enumerate(
        zip(
            ["puma_32k", "puma_5k", "chord"], 
            ["PUMA-32k" + experiment_suffix, "PUMA-5k" + experiment_suffix, "CHORD" + experiment_suffix]
        )
    ):
        ax[s] = subfig[si].subplots(
            2, 2, height_ratios=[1, 0.5], sharex=True, sharey="row"
        )

        # subfig[si].suptitle(label, y=1)

        ax[s][0, 0].set_title(" ", fontsize=12)
        
        if s == "chord" and not plot_chord:
            continue
        
        subfig[si].suptitle(label, y=1.01, x=0.57)

        for knli, knl_option in enumerate(knl_keys):

            #         for si, s, label in zip(
            #             [0, 1, 2],
            #             ["puma_32k", "puma_5k", "chord"],
            #             ["PUMA-32k", "PUMA-5k", "CHORD"],
            #         ):

            for pi, p in enumerate(par_list):

                pess = _onesigma_par(
                    noise, f"{s}_pess", f"no_sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
                )
                opt = _onesigma_par(
                    noise, f"{s}_opt", f"no_sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
                )

                pess_sub = _onesigma_par(
                    noise, f"{s}_pess", f"sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
                )
                opt_sub = _onesigma_par(
                    noise, f"{s}_opt", f"sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
                )

                no_LSS = _noLSS_onesigma_par("Tbfree_Planck_SOCMB_SOLensing", p)

                label_noLSS = None
                label_pess_bbn = None
                label_opt_bbn = None
                label_pess_cmb = None
                label_opt_cmb = None

                if pi == len(par_list) - 1 and si == 0 and knli == 0:
                    # if knli == 0:
                    label_noLSS = r"CMB+lensing"
                    label_pess_bbn = r"21cm+BBN, pess. wedge"
                    label_opt_bbn = r"21cm+BBN, opt. wedge"
                    # else:
                    label_pess_cmb = r"21cm+CMB+lensing, pess. wedge"
                    label_opt_cmb = r"21cm+CMB+lensing, opt. wedge"

                ax[s][0, knli].bar(
                    pi - 2.0 * w, no_LSS, w, color="k", alpha=1, label=label_noLSS
                )
                ax[s][0, knli].bar(
                    pi - 1.0 * w,
                    pess,
                    w,
                    color=colors[3],
                    alpha=pess_alpha,
                    label=label_pess_bbn,
                )
                ax[s][0, knli].bar(
                    pi - 0.0 * w, opt, w, color=colors[3], label=label_opt_bbn
                )

                ax[s][1, knli].bar(
                    pi - 1.0 * w,
                    pess_sub / pess,
                    w,
                    color=colors[3],
                    alpha=pess_alpha,
                )
                ax[s][1, knli].bar(
                    pi - 0 * w,
                    opt_sub / opt,
                    w,
                    color=colors[3],
                )

                pess = _onesigma_par(
                    noise, 
                    f"{s}_pess",
                    f"no_sub_{knl_option}",
                    f"{Tb_choice}_Planck_SOCMB_SOLensing",
                    p,
                )
                opt = _onesigma_par(
                    noise, 
                    f"{s}_opt",
                    f"no_sub_{knl_option}",
                    f"{Tb_choice}_Planck_SOCMB_SOLensing",
                    p,
                )

                pess_sub = _onesigma_par(
                    noise, 
                    f"{s}_pess", f"sub_{knl_option}", f"{Tb_choice}_Planck_SOCMB_SOLensing", p
                )
                opt_sub = _onesigma_par(
                    noise, 
                    f"{s}_opt", f"sub_{knl_option}", f"{Tb_choice}_Planck_SOCMB_SOLensing", p
                )

                ax[s][0, knli].bar(
                    pi + 1 * w,
                    pess,
                    w,
                    color=colors[4],
                    alpha=pess_alpha,
                    label=label_pess_cmb,
                )
                ax[s][0, knli].bar(
                    pi + 2 * w, opt, w, color=colors[4], label=label_opt_cmb
                )

                ax[s][1, knli].bar(
                    pi + 1 * w,
                    pess_sub / pess,
                    w,
                    color=colors[4],
                    alpha=pess_alpha,
                )
                ax[s][1, knli].bar(pi + 2 * w, opt_sub / opt, w, color=colors[4])

                # ax[2*si, knli].bar(
                #     pi+1.5*w,
                #     _onesigma_par(((f"{s}_pess", "no_sub", knl_option), Tb_choice, "PlanckSO_prior", "SO_lensing"), p),
                #     w,
                #     color=colors[2],
                #     alpha=pess_alpha,
                # )
                # ax[2*si, knli].bar(
                #     pi+2.5*w,
                #     _onesigma_par(((f"{s}_opt", "no_sub", knl_option), Tb_choice, "PlanckSO_prior", "SO_lensing"), p),
                #     w,
                #     color=colors[2]
                # )

                ax[s][0, knli].set_yscale("log")
                ax[s][0, knli].set_ylim(ylim_abs)
                if abs_dtick is not None:
                    ax[s][0, knli].set_yticks(
                        np.arange(ylim_abs[0], ylim_abs[1] + 1e-10, abs_dtick)
                    )
                ax[s][0, knli].set_xticks(
                    ticks=np.arange(len(par_list)),
                    labels=par_list_symbols[: len(par_list)],
                )
                # ax[2*si, knli].set_xticks([])
                ax[s][0, knli].grid(ls="--", axis="y")

                if knli == 0:
                    ax[s][0, knli].set_ylabel(r"$\sigma_p \;[\%]$")

                ax[s][1, knli].set_ylim(ylim_ratio)
                if ratio_dtick is not None:
                    ax[s][1, knli].set_yticks(
                        np.arange(ylim_ratio[0], ylim_ratio[1] + 1e-10, ratio_dtick)
                    )
                ax[s][1, knli].set_xticks(
                    ticks=np.arange(len(par_list)),
                    labels=par_list_symbols[: len(par_list)],
                )
                ax[s][1, knli].grid(ls="--", axis="y")

                if knli == 0:
                    ax[s][1, knli].set_ylabel(r"$\sigma_p$ ratio")
                    
                ax[s][0, knli].yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
                # ax.ticklabel_format(style='plain', axis='x')

    fig.legend(bbox_to_anchor=(0.5, -0.01), loc="upper center", fontsize=8.2, ncol=2)

    fig.suptitle(
        knl_titles[0] + "    " + knl_titles[1],
        y=1.05,
        x=0.1,
        horizontalalignment="left",
    )


def plot_LCDM_extensions(
    parameter_output,
    noise, 
    knl_keys,
    knl_titles,
    ylim_abs=[1e-2, 1e0],
    ylim_ratio=[0.6, 1.01],
    pess_alpha=0.3,
    abs_dtick=None,
    ratio_dtick=None,
    legend_loc="upper left",
    figsize=(9, 3.5),
    Tb_choice="Tbfree",
    plot_chord=True,
    experiment_suffix = "",
):
    """Plot LCDM-extension parameter forecasts."""

    par_list = [
        # "h",
        # "log(A_s)",
        # "n_s",
        # "omega_cdm",
        # "omega_b",
        "m_ncdm",
        "N_ur",
        "Omega_k",
    ]
    par_list_symbols = [
        # r"$h$",
        # r"$\log(A_{\rm s})$",
        # r"$n_{\rm s}$",
        # r"$\omega_{\rm c}$",
        # r"$\omega_{\rm b}$",
        r"$M_\nu\;[{\rm eV}]$",
        r"$N_{\rm eff}$",
        r"$10^2 \Omega_k$",
    ]
    w = 0.15

    params_fid = parameter_output[
        noise, 
        "puma_32k_opt",
        f"no_sub_{knl_keys[0]}",
        # "no_sub_dknldDinv0.3",
        "LCDM_parameters_Tbfree_Planck_SOCMB_noLensing",
    ]["params_fid"]

    def _onesigma_par(noise, survey, option, priors, p):

        if p in ["h", "log(A_s)", "n_s", "omega_cdm", "omega_b"]:
            return parameter_output[noise, survey, option, f"LCDM_parameters_{priors}"][
                f"onesigma_{p}"
            ]
        elif p == "m_ncdm":
            return parameter_output[noise, survey, option, f"LCDMplusMnu_parameters_{priors}"][
                f"onesigma_{p}"
            ]
        elif p == "N_ur":
            return parameter_output[
                noise, survey, option, f"LCDMplusNeff_parameters_{priors}"
            ][f"onesigma_{p}"]
        elif p == "Omega_k":
            return (
                100
                * parameter_output[
                    noise, survey, option, f"LCDMplusOmegaK_parameters_{priors}"
                ][f"onesigma_{p}"]
            )

    def _noLSS_onesigma_par(priors, p):

        if p in ["h", "log(A_s)", "n_s", "omega_cdm", "omega_b"]:
            return parameter_output[f"LCDM_parameters_{priors}_noLSS"][f"onesigma_{p}"]
        elif p == "m_ncdm":
            return parameter_output[f"LCDMplusMnu_parameters_{priors}_noLSS"][
                f"onesigma_{p}"
            ]
        elif p == "N_ur":
            return parameter_output[f"LCDMplusNeff_parameters_{priors}_noLSS"][
                f"onesigma_{p}"
            ]
        elif p == "Omega_k":
            return (
                100
                * parameter_output[f"LCDMplusOmegaK_parameters_{priors}_noLSS"][
                    f"onesigma_{p}"
                ]
            )

    fig = plt.figure(figsize=figsize, layout="constrained")
    subfig = fig.subfigures(3, 1)

    ax = {}
    for si, (s, label) in enumerate(
        zip(
            ["puma_32k", "puma_5k", "chord"], 
            ["PUMA-32k" + experiment_suffix, "PUMA-5k" + experiment_suffix, "CHORD" + experiment_suffix]
        )
    ):
        ax[s] = subfig[si].subplots(
            2, 2, height_ratios=[1, 0.5], sharex=True, sharey="row"
        )

        ax[s][0, 0].set_title(" ", fontsize=12)

        if s == "chord" and not plot_chord:
            continue
        
        subfig[si].suptitle(label, y=1.01, x=0.57)

        for knli, knl_option in enumerate(knl_keys):

            #         for si, s, label in zip(
            #             [0, 1, 2],
            #             ["puma_32k", "puma_5k", "chord"],
            #             ["PUMA-32k", "PUMA-5k", "CHORD"],
            #         ):

            for pi, p in enumerate(par_list):

                pess = _onesigma_par(
                    noise, f"{s}_pess", f"no_sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
                )
                opt = _onesigma_par(
                    noise, f"{s}_opt", f"no_sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
                )

                pess_sub = _onesigma_par(
                    noise, f"{s}_pess", f"sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
                )
                opt_sub = _onesigma_par(
                    noise, f"{s}_opt", f"sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
                )

                noLSS = _noLSS_onesigma_par("Tbfree_Planck_SOCMB_SOLensing", p)

                label_noLSS = None
                label_pess_bbn = None
                label_opt_bbn = None
                label_pess_cmb = None
                label_opt_cmb = None

                if pi == len(par_list) - 1 and si == 0 and knli == 0:
                    # if knli == 0:
                    label_noLSS = r"CMB+lensing"
                    label_pess_bbn = r"21cm+BBN, pess. wedge"
                    label_opt_bbn = r"21cm+BBN, opt. wedge"
                    # else:
                    label_pess_cmb = r"21cm+CMB+lensing, pess. wedge"
                    label_opt_cmb = r"21cm+CMB+lensing, opt. wedge"

                ax[s][0, knli].bar(
                    pi - 2.0 * w, noLSS, w, color="k", alpha=1, label=label_noLSS
                )
                ax[s][0, knli].bar(
                    pi - 1.0 * w,
                    pess,
                    w,
                    color=colors[3],
                    alpha=pess_alpha,
                    label=label_pess_bbn,
                )
                ax[s][0, knli].bar(
                    pi - 0.0 * w, opt, w, color=colors[3], label=label_opt_bbn
                )

                ax[s][1, knli].bar(
                    pi - 1.0 * w,
                    pess_sub / pess,
                    w,
                    color=colors[3],
                    alpha=pess_alpha,
                )
                ax[s][1, knli].bar(
                    pi - 0.0 * w,
                    opt_sub / opt,
                    w,
                    color=colors[3],
                )

                pess = _onesigma_par(
                    noise, 
                    f"{s}_pess",
                    f"no_sub_{knl_option}",
                    f"{Tb_choice}_Planck_SOCMB_SOLensing",
                    p,
                )
                opt = _onesigma_par(
                    noise, 
                    f"{s}_opt",
                    f"no_sub_{knl_option}",
                    f"{Tb_choice}_Planck_SOCMB_SOLensing",
                    p,
                )

                pess_sub = _onesigma_par(
                    noise, f"{s}_pess", f"sub_{knl_option}", f"{Tb_choice}_Planck_SOCMB_SOLensing", p
                )
                opt_sub = _onesigma_par(
                    noise, f"{s}_opt", f"sub_{knl_option}", f"{Tb_choice}_Planck_SOCMB_SOLensing", p
                )

                ax[s][0, knli].bar(
                    pi + 1 * w,
                    pess,
                    w,
                    color=colors[4],
                    alpha=pess_alpha,
                    label=label_pess_cmb,
                )
                ax[s][0, knli].bar(
                    pi + 2 * w, opt, w, color=colors[4], label=label_opt_cmb
                )

                ax[s][1, knli].bar(
                    pi + 1 * w,
                    pess_sub / pess,
                    w,
                    color=colors[4],
                    alpha=pess_alpha,
                )
                ax[s][1, knli].bar(pi + 2 * w, opt_sub / opt, w, color=colors[4])

                ax[s][0, knli].set_yscale("log")
                ax[s][0, knli].set_ylim(ylim_abs)
                if abs_dtick is not None:
                    ax[s][0, knli].set_yticks(
                        np.arange(ylim_abs[0], ylim_abs[1] + 1e-10, abs_dtick)
                    )
                ax[s][0, knli].set_xticks(
                    ticks=np.arange(len(par_list)),
                    labels=par_list_symbols[: len(par_list)],
                )
                # ax[2*si, knli].set_xticks([])
                ax[s][0, knli].grid(ls="--", axis="y")

                if knli == 0:
                    ax[s][0, knli].set_ylabel(r"$\sigma_p$")

                ax[s][1, knli].set_ylim(ylim_ratio)
                if ratio_dtick is not None:
                    ax[s][1, knli].set_yticks(
                        np.arange(ylim_ratio[0], ylim_ratio[1] + 1e-10, ratio_dtick)
                    )
                ax[s][1, knli].set_xticks(
                    ticks=np.arange(len(par_list)),
                    labels=par_list_symbols[: len(par_list)],
                )
                ax[s][1, knli].grid(ls="--", axis="y")

                if knli == 0:
                    ax[s][1, knli].set_ylabel(r"$\sigma_p$ ratio")
                    
                ax[s][0, knli].yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))

    fig.suptitle(
        knl_titles[0] + "    " + knl_titles[1],
        y=1.05,
        x=0.1,
        horizontalalignment="left",
    )

    fig.legend(bbox_to_anchor=(0.5, -0.01), loc="upper center", fontsize=8.2, ncol=2)


def plot_LCDM_and_extensions(
    parameter_output,
    noise, 
    knl_key,
    knl_title,
    ylim_abs=[1e-2, 1e0],
    ylim_ratio=[0.6, 1.01],
    pess_alpha=0.3,
    abs_dtick=None,
    ratio_dtick=None,
    legend_loc="upper left",
    figsize=(9, 3.5),
    Tb_choice="Tbfree",
):
    """Plot LCDM parameter forecasts and LCDM-extension forecasts in single figure."""

    par_list = [
        "h",
        "log(A_s)",
        "n_s",
        "omega_cdm",
        "omega_b",
        "m_ncdm",
        "N_ur",
        "Omega_k",
    ]
    par_list_symbols = [
        r"$h$",
        r"$\log A_{\rm s} $",
        r"$n_{\rm s}$",
        r"$\omega_{\rm c}$",
        r"$\omega_{\rm b}$",
        r"$M_\nu\;[{\rm eV}]$",
        r"$N_{\rm eff}$",
        r"$10^2 \Omega_k$",
    ]
    w = 0.15

    par_list_A = [
        "h",
        "log(A_s)",
        "n_s",
        "omega_cdm",
        "omega_b",
    ]
    par_list_B = ["m_ncdm", "N_ur", "Omega_k"]

    params_fid = parameter_output[
        noise, 
        "puma_32k_opt",
        f"no_sub_{knl_key}",
        # "no_sub_dknldDinv0.3",
        "LCDM_parameters_Tbfree_Planck_SOCMB_noLensing",
    ]["params_fid"]

    def _onesigma_par(noise, survey, option, priors, p):
        # In percent of fiducial value, except for Mnu, Neff, and OmegaK

        if p == "log(A_s)":
            p_fid = np.abs(np.log(params_fid["A_s"]))
        elif p == "Omega_k":
            p_fid = 1
        elif p == "m_ncdm":
            p_fid = 0.06
        else:
            p_fid = params_fid[p]

        if p in ["h", "log(A_s)", "n_s", "omega_cdm", "omega_b"]:
            return (
                100
                * parameter_output[noise, survey, option, f"LCDM_parameters_{priors}"][
                    f"onesigma_{p}"
                ]
                / p_fid
            )
        elif p == "m_ncdm":
            return parameter_output[noise, survey, option, f"LCDMplusMnu_parameters_{priors}"][
                f"onesigma_{p}"
            ]
        elif p == "N_ur":
            return parameter_output[
                noise, survey, option, f"LCDMplusNeff_parameters_{priors}"
            ][f"onesigma_{p}"]
        elif p == "Omega_k":
            return (
                100
                * parameter_output[
                    noise, survey, option, f"LCDMplusOmegaK_parameters_{priors}"
                ][f"onesigma_{p}"]
            )

    def _noLSS_onesigma_par(priors, p):
        # In percent of fiducial value, except for Mnu, Neff, and OmegaK

        if p == "log(A_s)":
            p_fid = np.abs(np.log(params_fid["A_s"]))
        elif p == "Omega_k":
            p_fid = 1
        elif p == "m_ncdm":
            p_fid = 0.06
        else:
            p_fid = params_fid[p]

        if p in ["h", "log(A_s)", "n_s", "omega_cdm", "omega_b"]:
            return (
                100
                * parameter_output[f"LCDM_parameters_{priors}_noLSS"][f"onesigma_{p}"]
                / p_fid
            )
        elif p == "m_ncdm":
            return parameter_output[f"LCDMplusMnu_parameters_{priors}_noLSS"][
                f"onesigma_{p}"
            ]
        elif p == "N_ur":
            return parameter_output[f"LCDMplusNeff_parameters_{priors}_noLSS"][
                f"onesigma_{p}"
            ]
        elif p == "Omega_k":
            return (
                100
                * parameter_output[f"LCDMplusOmegaK_parameters_{priors}_noLSS"][
                    f"onesigma_{p}"
                ]
            )

    fig = plt.figure(figsize=figsize, layout="constrained")
    subfig = fig.subfigures(3, 1)

    ax = {}
    for si, (s, label) in enumerate(
        zip(["puma_32k", "puma_5k", "chord"], ["PUMA-32k", "PUMA-5k", "CHORD"])
    ):
        ax[s] = subfig[si].subplots(
            2, 2, height_ratios=[1, 0.5], sharex="col"
        )

        # subfig[si].suptitle(label, y=1)

        ax[s][0, 0].set_title(" ", fontsize=12)

        subfig[si].suptitle(label, y=1.01, x=0.57)

        knl_option = knl_key

        for pi, p in enumerate(par_list):

            if p in ["h", "log(A_s)", "n_s", "omega_cdm", "omega_b"]:
                knli = 0
                pii = pi
                # pos = pi-2
            else:
                knli = 1
                pii = pi - 5

            pess = _onesigma_par(
                noise, f"{s}_pess", f"no_sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
            )
            opt = _onesigma_par(
                noise, f"{s}_opt", f"no_sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
            )

            pess_sub = _onesigma_par(
                noise, f"{s}_pess", f"sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
            )
            opt_sub = _onesigma_par(
                noise, f"{s}_opt", f"sub_{knl_option}", f"{Tb_choice}_noCMB_noLensing", p
            )

            no_LSS = _noLSS_onesigma_par("Tbfree_Planck_SOCMB_SOLensing", p)

            label_noLSS = None
            label_pess_bbn = None
            label_opt_bbn = None
            label_pess_cmb = None
            label_opt_cmb = None

            if pi == 0 and si == 0 and knli == 0:
                # if knli == 0:
                label_noLSS = r"CMB+lensing"
                label_pess_bbn = r"21cm+BBN, pess. wedge"
                label_opt_bbn = r"21cm+BBN, opt. wedge"
                # else:
                label_pess_cmb = r"21cm+CMB+lensing, pess. wedge"
                label_opt_cmb = r"21cm+CMB+lensing, opt. wedge"

            ax[s][0, knli].bar(
                pii - 2.0 * w, no_LSS, w, color="k", alpha=1, label=label_noLSS
            )
            ax[s][0, knli].bar(
                pii - 1.0 * w,
                pess,
                w,
                color=colors[3],
                alpha=pess_alpha,
                label=label_pess_bbn,
            )
            ax[s][0, knli].bar(
                pii - 0.0 * w, opt, w, color=colors[3], label=label_opt_bbn
            )

            ax[s][1, knli].bar(
                pii - 1.0 * w,
                pess_sub / pess,
                w,
                color=colors[3],
                alpha=pess_alpha,
            )
            ax[s][1, knli].bar(
                pii - 0 * w,
                opt_sub / opt,
                w,
                color=colors[3],
            )

            pess = _onesigma_par(
                noise, f"{s}_pess", f"no_sub_{knl_option}", f"{Tb_choice}_Planck_SOCMB_SOLensing", p
            )
            opt = _onesigma_par(
                noise, f"{s}_opt", f"no_sub_{knl_option}", f"{Tb_choice}_Planck_SOCMB_SOLensing", p
            )

            pess_sub = _onesigma_par(
                noise, f"{s}_pess", f"sub_{knl_option}", f"{Tb_choice}_Planck_SOCMB_SOLensing", p
            )
            opt_sub = _onesigma_par(
                noise, f"{s}_opt", f"sub_{knl_option}", f"{Tb_choice}_Planck_SOCMB_SOLensing", p
            )

            ax[s][0, knli].bar(
                pii + 1 * w,
                pess,
                w,
                color=colors[4],
                alpha=pess_alpha,
                label=label_pess_cmb,
            )
            ax[s][0, knli].bar(
                pii + 2 * w, opt, w, color=colors[4], label=label_opt_cmb
            )

            ax[s][1, knli].bar(
                pii + 1 * w,
                pess_sub / pess,
                w,
                color=colors[4],
                alpha=pess_alpha,
            )
            ax[s][1, knli].bar(pii + 2 * w, opt_sub / opt, w, color=colors[4])

            ax[s][0, knli].set_yscale("log")
            ax[s][0, knli].set_ylim(ylim_abs)
            if abs_dtick is not None:
                ax[s][0, knli].set_yticks(
                    np.arange(ylim_abs[0], ylim_abs[1] + 1e-10, abs_dtick)
                )

            # ax[2*si, knli].set_xticks([])
            ax[s][0, knli].grid(ls="--", axis="y")
            

            ax[s][1, knli].set_ylim(ylim_ratio)
            if ratio_dtick is not None:
                ax[s][1, knli].set_yticks(
                    np.arange(ylim_ratio[0], ylim_ratio[1] + 1e-10, ratio_dtick)
                )
            ax[s][1, knli].grid(ls="--", axis="y")
                

        ax[s][0, 0].set_xticks(
            ticks=np.arange(len(par_list_A)), labels=par_list_symbols[: len(par_list_A)]# , fontsize=7
        )
        # ax[s][0, 0].set_xticklabels(labels=par_list_symbols[: len(par_list_A)], fontsize=5)
        # ax[s][0, 0].tick_params(axis="x", labelsize=4)
        ax[s][0, 1].set_xticks(
            ticks=np.arange(len(par_list_B)), labels=par_list_symbols[len(par_list_A) :]
        )
        
        ax[s][0, 0].set_ylabel(r"$\sigma_p \;[\%]$")
        ax[s][0, 1].set_ylabel(r"$\sigma_p$")
        
        ax[s][1, 0].set_ylabel(r"$\sigma_p$ ratio")
        ax[s][1, 1].set_ylabel(r"$\sigma_p$ ratio")
        
        ax[s][0, 0].yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
        ax[s][0, 1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))

    fig.legend(bbox_to_anchor=(0.5, -0.01), loc="upper center", fontsize=8.2, ncol=2)

    fig.suptitle(knl_title, y=1.05, x=0.58, horizontalalignment="center")
