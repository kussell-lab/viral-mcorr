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
from mcorr import fit_p2, read_corr, FitDatas, \
    write_fitting_results, plot_fit, plot_params, write_fitting_reports, \
    geom_r1, const_r1
from mcorr.fit import fit_model, vary_fit
from lmfit import fit_report
import csv
import pandas as pd
#from lmfit.printfuncs import report_fit
from mcorr.lmfitFunctions import perform_lmfit

def main():
    """Run fitting using lmfit, and generate output files and plots"""
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Fit the actual data (not the bootstraps) and return goodness-of fit stats")
    parser.add_argument("corr_file", type = str, help='correlation input file')
    parser.add_argument("output_prefix", type=str, help='output file prefix')
    parser.add_argument('--fit_start', type=int, default=3,
                        help='fitting range starts at')
    parser.add_argument('--fit_end', type=int, default=300,
                        help='fitting range ends at')
    parser.add_argument("--use_geom_frag", action="store_true",
                        help='use geometric distribution for fragment sizes')
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
    fitres = perform_lmfit(x, y, d_sample)
    ## write a fit report as generated by lmfit (includes chi-squared, uncertainties, etc)
    params = fitres.params.valuesdict()
    thetaS = fitres.params["theta_s"]
    phiS = fitres.params["phi_s"]
    f = fitres.params["f"]
    lmfitfile = prefix + "_lmfit_report.csv"
    with open(lmfitfile, "w+") as csvfile:
        lmfit_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        lmfit_writer.writerow(["fit_success", fitres.success])
        lmfit_writer.writerow(["function_evals", fitres.nfev])
        lmfit_writer.writerow(["data_points", fitres.ndata])
        lmfit_writer.writerow(["variables", fitres.nvarys])
        lmfit_writer.writerow(["message", fitres.message])
        lmfit_writer.writerow(["thetaS (init)", thetaS.init_value])
        lmfit_writer.writerow(["f (init)", f.init_value])
        lmfit_writer.writerow(["phiS (init)", phiS.init_value])
        lmfit_writer.writerow([""])
        lmfit_writer.writerow(["d_s", "theta_s", "f", "phi_s",
                               "theta_p", "phi_p", "c", "d_theta_p",
                               "d_theta_s", "chisq", "red-chisq"])
        lmfit_writer.writerow([params["d_s"], params["theta_s"], params["f"], params["phi_s"],
                               params["theta_p"], params["phi_p"], params["c_s"], params["d_theta_p"],
                               params["d_theta_s"], fitres.chisqr, fitres.redchi])
    ##save the residuals as a .csv file
    residuals = fitres.residual
    resdat = pd.DataFrame(residuals)
    resdat.to_csv(prefix+"_residuals.csv", header=None)
    ##plot the best fit and the residuals
    best_fit_file = prefix + "_best_fit.svg"
    plot_fit(all, fitres, best_fit_file, title=title)

if __name__ == "__main__":
    main()






