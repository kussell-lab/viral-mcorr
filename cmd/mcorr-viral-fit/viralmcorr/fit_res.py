# Copyright 2022 Asher Preska Steinberg
#
# Initially modified from fit_res.py of mcorr
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

class FitRes(object):
    """Fitting results"""
    def __init__(self, group, residual, params, d_sample):
        self.group = group
        self.d_sample = d_sample
        self.residual = residual
        if "thetaP" in params:
            self.theta_pool = params['thetaP']
        if 'phiP' in params:
            self.phi_pool = params['phiP']
        if 'f' in params:
            self.fbar = params['f']
        if 'phiP' in params:
            self.ratio = self.phi_pool / self.theta_pool
            if 'f' in params:
                self.rho = self.phi_pool * self.fbar
        if 'c' in params:
            self.c = params['c']
        if 'dc' in params:
            self.d_clonal = params['dc']
        if 'dp' in params:
            self.d_pool = params['dp']
        if 'phiS' in params:
            self.phi_s = params['phiS']
        if 'thetaS' in params:
            self.theta_s = params['thetaS']

    def get_values(self, attributes):
        """Get attribute values"""
        values = []
        for name in attributes:
            if hasattr(self, name):
                values.append(getattr(self, name))
            else:
                values.append("NA")
        return values


