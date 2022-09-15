#!/usr/bin/env python3

# Copyright 2022 Asher Preska Steinberg
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
#from mcorr import fit_p2
# from mcorr import fit_p2, read_corr, FitDatas, \
#     write_fitting_results, plot_fit, plot_params, write_fitting_reports, \
#     geom_r1, const_r1
from .fit_data import FitDatas
from .fit import fit_p2, fit_modelopts, geom_r1, const_r1, fit_zerorecombo, zero_r1, \
    solve_zerorecombo, calc_akaike_weights
from .writer import write_fitting_reports, write_fitting_results
from .plot import plot_fit, plot_params, plot_zerorecombo
from .corr_res import read_corr
import numpy as np
import csv
import pandas as pd
import scipy
#from lmfit.printfuncs import report_fit

def main():
    """Fit the data and bootstrap replicates of pair-averaged correlation profiles using coalescent model w/ and w/o recombination"""
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Fit data and bootstrap replicates for pair-averaged correlation profiles using coalescent model w/ and w/o recombination")
    parser.add_argument("corr_file", type = str, help='correlation input file')
    parser.add_argument("output_prefix", type=str, help='output file prefix')
    parser.add_argument('--fit_start', type=int, default=3,
                        help='fitting range starts at')
    parser.add_argument('--fit_end', type=int, default=300,
                        help='fitting range ends at')
    parser.add_argument("--use_geom_frag", action="store_true",
                        help='use geometric distribution for fragment sizes')
    parser.add_argument("--fit_method", type=str, default="least_squares", help="lmfit method (see lmfit documentation)")
    parser.add_argument("--max_nfev", type=int, default=int(1e6),
                        help='max number of function evaluations before lmfit quits')
    parser.add_argument("--genome_length", type=int, default=int(3e4), help="length of viral genome")
    parser.add_argument("--template_switching", type=bool, default=True, help="fit the bootstrap replicates with the template-switching (default) or fragment-incorporation model")
    parser.add_argument('--quiet', action="store_true")
    parser.add_argument("--title", type=str, help="plot title", default="")
    opts = parser.parse_args()
    corr_file = opts.corr_file
    prefix = opts.output_prefix
    fit_start = opts.fit_start
    fit_end = opts.fit_end
    quiet = opts.quiet
    use_geom_frag = opts.use_geom_frag
    title = opts.title
    fit_method = opts.fit_method
    max_nfev = opts.max_nfev
    genome_length = opts.genome_length
    fixed_f = opts.template_switching

    ##for testing fixes
    # dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit/cluster8_cluster221'
    # corr_file = os.path.join(dir, 'cluster8_cluster221_CORE_XMFA_OUT.csv')
    # prefix = 'cluster8_cluster221_CORE_FIT_OUT_0205test'
    # fit_start = 3
    # fit_end = 300
    # quiet = False
    # use_geom_frag = False
    # title=""

    # read correlation results and prepare fitting data
    corr_results = read_corr(corr_file)
    fitdatas = FitDatas(corr_results, fit_start, fit_end)
    ##do fitting
    r1_func = const_r1
    #if you want to use a geometric distribution of fragments
    if use_geom_frag:
        r1_func = geom_r1

    all = fitdatas.get("all")
    x = all.xvalues
    y = all.yvalues
    d_sample = all.d_sample
    # fit with the recombination model w/o fixing fragment size
    fitres = fit_modelopts(x, y, d_sample, r1_func, max_nfev, fit_method, genome_length)
    # fit with the recombination model w/ fixed fragment size
    fixedres = fit_modelopts(x, y, d_sample, r1_func, max_nfev, fit_method, genome_length, fixed_f=True)
    #fit with the zero recombination model
    zdata, zres, zchisq, zred_chisq, zaic, zthetaS, z_dc = solve_zerorecombo(x, y, d_sample)
    ## write a fit report as generated by lmfit (includes chi-squared, uncertainties, etc)
    ## for the recombination model, fixed fragment recombination, and the null recombination model
    params = fitres.params.valuesdict()
    thetaS = fitres.params["thetaS"]
    phiS = fitres.params["phiS"]
    f = fitres.params["f"]
    ##for fixed fragment size model
    fixedparams = fixedres.params.valuesdict()
    fixedthetaS = fixedres.params["thetaS"]
    fixedphiS = fixedres.params["phiS"]
    fixedf = fixedres.params["f"]
    #### also calculate akaike weights ...
    w_v, w_f, w_z = calc_akaike_weights(fitres.aic, fixedres.aic, zaic)
    lmfitfile = prefix + "_comparemodels.csv"
    with open(lmfitfile, "w+") as csvfile:
        lmfit_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        lmfit_writer.writerow(["", "recombo (frag-incorp)", "recombo (template-switch)", "zero recombo"])
        lmfit_writer.writerow(["fit_success", fitres.success, fixedres.success, "n/a"])
        lmfit_writer.writerow(["function_evals", fitres.nfev, fixedres.nfev, "n/a"])
        lmfit_writer.writerow(["data_points", fitres.ndata, fixedres.ndata, zdata])
        lmfit_writer.writerow(["variables", fitres.nvarys, fixedres.nvarys, 1])
        lmfit_writer.writerow(["message", fitres.message, fixedres.message, "n/a"])
        lmfit_writer.writerow(["thetaS (init)", thetaS.init_value, fixedthetaS.init_value, "n/a"])
        lmfit_writer.writerow(["f (init)", f.init_value, fixedf.init_value, 0])
        lmfit_writer.writerow(["phiS (init)", phiS.init_value, fixedphiS.init_value, 0])
        lmfit_writer.writerow([""])
        lmfit_writer.writerow(["recombination", "d_s", "theta_s", "f", "phi_s",
                               "theta_p", "phi_p", "c", "d_theta_p",
                               "d_theta_s", "chisq", "red-chisq", "AIC", "akaike_weight", "evidence_ratio (model/null)"])
        lmfit_writer.writerow(["recombo (vary f)", params["ds"], params["thetaS"], params["f"], params["phiS"],
                               params["thetaP"], params["phiP"], params["c"], params["dp"],
                               params["dc"], fitres.chisqr, fitres.redchi, fitres.aic, w_v])
        lmfit_writer.writerow(["recombo (fixed f)", fixedparams["ds"], fixedparams["thetaS"], fixedparams["f"], fixedparams["phiS"],
                               fixedparams["thetaP"], fixedparams["phiP"], fixedparams["c"], fixedparams["dp"],
                               fixedparams["dc"], fixedres.chisqr, fixedres.redchi, fixedres.aic, w_f])
        lmfit_writer.writerow(["zero_recombo", d_sample, zthetaS, np.NAN, np.NAN,
                               np.NAN, np.NAN, np.NAN, np.NAN,
                               z_dc, zchisq, zred_chisq, zaic, w_z])


    ##save the residuals as a .csv file
    residuals = zres
    resdat = pd.DataFrame(residuals)
    resdat.to_csv(prefix+"_zero-recombo_residuals.csv", header=None)
    #plot the best fit and the residuals
    best_fit_file = prefix + "_zero-recombo_best_fit.svg"
    plot_zerorecombo(all, zres, best_fit_file, title=title)
    if not fixed_f:
        if fitres.success:
            residuals = fitres.residual
            resdat = pd.DataFrame(residuals)
            resdat.to_csv(prefix+"_frag-incorp_residuals.csv", header=None)
            #plot the best fit and the residuals
            best_fit_file = prefix + "_frag-incorp_best_fit.svg"
            plot_fit(all, fitres, best_fit_file, title=title)
        else:
            print("Fitting w/ fragment-incorporation model failed for %s" % corr_file)
    else:
        if fixedres.success:
            residuals = fixedres.residual
            resdat = pd.DataFrame(residuals)
            resdat.to_csv(prefix+"_template-switch_residuals.csv", header=None)
            #plot the best fit and the residuals
            best_fit_file = prefix + "_template-switch_best_fit.svg"
            plot_fit(all, fixedres, best_fit_file, title=title)
        else:
            print("Fitting w/ template-switching model failed for %s" % corr_file)


    ##fit all the bootstraps w/ fragment-incoporation or with template-switching model
    if not fixed_f:
        print("fitting bootstraps with fragment-incorporation model ...")
        fit_results = fit_p2(fitdatas, max_nfev, fit_method, genome_length,
                             r1_func=r1_func, disable_progress_bar=quiet)
        # parameters to report
        model_params = ["group", "d_sample", "theta_pool",
                        "phi_pool", "ratio", "fbar", "c", "d_pool",
                        "d_clonal", 'theta_s', 'phi_s']
        # save the results of fitting the bootstraps into csv file
        csv_file = prefix + "_frag-incorp_fit_results.csv"
        write_fitting_results(fit_results, model_params, csv_file)
        # write fitting report for bootstrapping
        report_file = prefix + "_frag-incorp_fit_report.txt"
        write_fitting_reports(fit_results, model_params[1:7], report_file)
    else:
        print("Fitting bootstraps with template-switching model ...")
        fit_results = fit_p2(fitdatas, max_nfev, fit_method, genome_length,
                             r1_func=r1_func, disable_progress_bar=quiet, fixed_f=True)
        # parameters to report
        model_params = ["group", "d_sample", "theta_pool",
                        "phi_pool", "ratio", "fbar", "c", "d_pool",
                        "d_clonal", 'theta_s', 'phi_s']
        # save the results of fitting the bootstraps into csv file
        csv_file = prefix + "_template-switch_fit_results.csv"
        write_fitting_results(fit_results, model_params, csv_file)
        # write fitting report for bootstrapping
        report_file = prefix + "_template-switch_fit_report.txt"
        write_fitting_reports(fit_results, model_params[1:7], report_file)

if __name__ == "__main__":
    main()





