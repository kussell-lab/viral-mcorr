# Copyright 2022 Asher Preska Steinberg
#
# Initially modified from writer.py of mcorr
# https://github.com/kussell-lab/mcorr

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

from . import FitReport
def write_fitting_results(all_results, model_params, out_file):
    """
    write fitting results into a .csv file.
    """
    # write fitting results.
    sep = ","
    with open(out_file, 'w') as out:
        out.write(sep.join(model_params)+"\n")
        for fit_res in all_results:
            values = fit_res.get_values(model_params)
            out.write(sep.join([str(x) for x in values])+"\n")

def write_pairfitting_results(all_results, model_params, all_aic, all_zaic, all_zthetaS, out_file):
    """
    write fitting results into a .csv file.
    """
    # write fitting results.
    sep = ","
    ##header
    header = list(model_params)
    header.append("aic")
    header.append("z_aic")
    header.append("z_theta_s")
    with open(out_file, 'w') as out:
        out.write(sep.join(header)+"\n")
        i = 0
        for fit_res in all_results:
            values = fit_res.get_values(model_params)
            line = list(values)
            line.append(all_aic[i])
            line.append(all_zaic[i])
            line.append(all_zthetaS[i])
            out.write(sep.join([str(x) for x in line])+"\n")
            i = i + 1

def write_fitting_reports(all_results, model_params, out_file):
    """
    write fitting reports into a .txt file.
    """
    with open(out_file, 'w') as out:
        for param_name in model_params:
            label_name = param_name
            if param_name == "ratio":
                label_name = "gamma/mu"
            report = FitReport(all_results, param_name, label_name)
            out.write(report.report()+"\n")

