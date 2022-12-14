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
from distutils.util import strtobool

#from mcorr import fit_p2
# from mcorr import fit_p2, read_corr, FitDatas, \
#     write_fitting_results, plot_fit, plot_params, write_fitting_reports, \
#     geom_r1, const_r1
from .fit_data import FitDatas
from .fit import fit_p2, fit_modelopts, geom_r1, const_r1, fit_zerorecombo, zero_r1, solve_zerorecombo, fit_allpairs
from .writer import write_fitting_reports, write_fitting_results, write_pairfitting_results
from .plot import plot_fit, plot_params, plot_zerorecombo
from .corr_res import read_corr
import numpy as np
import csv
import pandas as pd
import scipy
#from lmfit.printfuncs import report_fit

def main():
    """Fit correlation profiles measured across individual viral sequence pairs"""
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Fit correlation profiles measured across viral sequence pairs to coalescent model w/ and w/o recombination")
    parser.add_argument("corr_file", type = str, help='correlation input file')
    parser.add_argument("output_prefix", type=str, help='output file prefix')
    parser.add_argument('--fit_start', type=int, default=3,
                        help='fitting range starts at')
    parser.add_argument('--fit_end', type=int, default=300,
                        help='fitting range ends at')
    parser.add_argument("--template_switching", type=strtobool,
                        default='True', help="use template-switching (default) or fragment-incorporation model")
    parser.add_argument("--use_geom_frag", action="store_true",
                        help='use geometric distribution for fragment sizes')
    parser.add_argument("--fit_method", type=str, default="least_squares", help="lmfit method (see lmfit documentation)")
    parser.add_argument("--max_nfev", type=int, default=int(1e6),
                        help='max number of function evaluations before lmfit quits')
    parser.add_argument("--genome_length", type=int, default=int(3e4), help="length of viral genome")
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
    # calculate the pair averaged correlation profile and add it to the corr file ....
    # may take this bit out ...
    corrdat = pd.read_csv(corr_file)
    corrdat = corrdat[corrdat["b"] != "all"].copy()
    grouped = corrdat.groupby('l').mean()
    meancorr = grouped.reset_index()
    t = ["Ks"]
    for i in np.arange(0,len(meancorr)-1):
        t.append("P2")
    meancorr["t"] = t
    meancorr["b"] = "all"
    #corrdat = corrdat.append(meancorr)
    corrdat = pd.concat([corrdat, meancorr])
    ##drops rows with nans which screws with read_corr
    corrdat = corrdat.dropna()
    corrdat.to_csv(corr_file, index=False)
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
    ##fit with the recombination model
    if fixed_f:
        fitres = fit_modelopts(x, y, d_sample, r1_func, max_nfev, fit_method, genome_length, fixed_f=True)
    else:
        fitres = fit_modelopts(x, y, d_sample, r1_func, max_nfev, fit_method, genome_length)
    #fitres = fit_modelopts(x, y, d_sample, r1_func, max_nfev, fit_method, genome_length)
    #fit with the zero recombination model
    zdata, zres, zchisq, zred_chisq, zaic, zthetaS, z_dc = solve_zerorecombo(x, y, d_sample)
    ## write a fit report as generated by lmfit (includes chi-squared, uncertainties, etc)
    ## for the recombination model vs the null recombination model
    params = fitres.params.valuesdict()
    thetaS = fitres.params["thetaS"]
    phiS = fitres.params["phiS"]
    f = fitres.params["f"]
    lmfitfile = prefix + "_comparemodels.csv"
    with open(lmfitfile, "w+") as csvfile:
        lmfit_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        lmfit_writer.writerow(["", "recombo", "zero recombo"])
        lmfit_writer.writerow(["fit_success", fitres.success, "n/a"])
        lmfit_writer.writerow(["function_evals", fitres.nfev, "n/a"])
        lmfit_writer.writerow(["data_points", fitres.ndata, zdata])
        lmfit_writer.writerow(["variables", fitres.nvarys, 1])
        lmfit_writer.writerow(["message", fitres.message, "n/a"])
        lmfit_writer.writerow(["thetaS (init)", thetaS.init_value, "n/a"])
        lmfit_writer.writerow(["f (init)", f.init_value, 0])
        lmfit_writer.writerow(["phiS (init)", phiS.init_value, 0])
        lmfit_writer.writerow([""])
        lmfit_writer.writerow(["recombination", "d_s", "theta_s", "f", "phi_s",
                               "theta_p", "phi_p", "c", "d_theta_p",
                               "d_theta_s", "chisq", "red-chisq", "AIC"])
        lmfit_writer.writerow(["recombo", params["ds"], params["thetaS"], params["f"], params["phiS"],
                               params["thetaP"], params["phiP"], params["c"], params["dp"],
                               params["dc"], fitres.chisqr, fitres.redchi, fitres.aic])
        lmfit_writer.writerow(["zero_recombo", d_sample, zthetaS, np.NAN, np.NAN,
                               np.NAN, np.NAN, np.NAN, np.NAN,
                               z_dc, zchisq, zred_chisq, zaic])


    ##save the residuals as a .csv file
    residuals = zres
    resdat = pd.DataFrame(residuals)
    resdat.to_csv(prefix+"_zero-recombo_residuals.csv", header=None)
    #plot the best fit and the residuals
    best_fit_file = prefix + "_zero-recombo_best_fit.svg"
    plot_zerorecombo(all, zres, best_fit_file, title=title)
    if fitres.success:
        residuals = fitres.residual
        resdat = pd.DataFrame(residuals)
        resdat.to_csv(prefix+"_recombo_residuals.csv", header=None)
        #plot the best fit and the residuals
        best_fit_file = prefix + "_recombo_best_fit.svg"
        plot_fit(all, fitres, best_fit_file, title=title)
    else:
        print("Fitting failed for %s" % corr_file)

    ##fit all the pairs ...
    if fixed_f == True:
        print("fitting pairs with template-switching model ...")
        fit_results, all_aic, all_zaic, all_zthetaS = fit_allpairs(fitdatas, max_nfev, fit_method, genome_length,
                                                                   r1_func=r1_func, disable_progress_bar=quiet, fixed_f=True)
    else:
        print("fitting pairs with fragment-incorporation model ...")
        fit_results, all_aic, all_zaic, all_zthetaS = fit_allpairs(fitdatas, max_nfev, fit_method, genome_length,
                                                r1_func=r1_func, disable_progress_bar=quiet)
    # parameters to report
    model_params = ["group", "d_sample", "theta_pool",
                    "phi_pool", "ratio", "fbar", "c", "d_pool",
                    "d_clonal", 'theta_s', 'phi_s']
    # save the results of fitting the bootstraps into csv file
    csv_file = prefix + "_fit_results.csv"
    write_pairfitting_results(fit_results, model_params, all_aic, all_zaic, all_zthetaS, csv_file)
    #write_fitting_results(fit_results, model_params, csv_file)
    # write fitting report for bootstrapping
    # report_file = prefix + "_bootstrapping_report.txt"
    # write_fitting_reports(fit_results, model_params[1:7], report_file)

if __name__ == "__main__":
    main()






