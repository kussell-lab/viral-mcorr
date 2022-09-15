# Copyright 2022 Asher Preska Steinberg
#
# Initially modified from fit_data.py of mcorr
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

import numpy

class FitData(object):
    """Fitting data"""
    def __init__(self, group, xvalues, yvalues, d_sample):
        self.group = group
        self.xvalues = xvalues
        self.yvalues = yvalues
        self.d_sample = d_sample

class FitDatas(object):
    """Fitting data"""
    def __init__(self, corr_results, fit_start, xmax):
        corr_map = {}
        groups = []
        for row in corr_results:
            rows = corr_map.get(row.group, [])
            if len(rows) == 0:
                groups.append(row.group)
            rows.append(row)
            corr_map[row.group] = rows
        fitdata_map = {}
        for group, items in corr_map.items():
            xvalues, yvalues, d_sample = prepare_fitting_data(
                items, fit_start, xmax)
            fitdata_map[group] = FitData(group, xvalues, yvalues, d_sample)
        self.fitdata_dict = fitdata_map
        self.groups = groups
    def has(self, group):
        """return True if the group is in the data"""
        return group in self.fitdata_dict

    def get(self, group):
        """return fit data"""
        fitdata = self.fitdata_dict.get(group, None)
        return fitdata
    def getall(self):
        """return all"""
        return [self.fitdata_dict[group] for group in self.groups]

def prepare_fitting_data(fitdata, fit_start, xmax):
    """Prepare fitting xvalues and yvalues"""
    xvalues = []
    yvalues = []
    diver = 0
    for row in fitdata:
        if row.corrtype == 'P2' and row.lag >= fit_start and row.lag <= xmax:
            xvalues.append(row.lag)
            yvalues.append(row.value)
        elif row.corrtype == 'Ks':
            diver = row.value
    xvalues = numpy.array(xvalues)
    yvalues = numpy.array(yvalues)
    return (xvalues, yvalues, diver)

class FitGeneDatas(object):
    """Fitting data"""
    def __init__(self, corr_results, fit_start, xmax):
        corr_map = {}
        groups = []
        for row in corr_results:
            rows = corr_map.get(row.group, [])
            if len(rows) == 0:
                groups.append(row.group)
            rows.append(row)
            corr_map[row.group] = rows
        fitdata_map = {}
        for group, items in corr_map.items():
            xvalues, yvalues, d_sample = prepare_fitting_genedata(
                items, fit_start, xmax)
            fitdata_map[group] = FitData(group, xvalues, yvalues, d_sample)
        self.fitdata_dict = fitdata_map
        self.groups = groups
    def has(self, group):
        """return True if the group is in the data"""
        return group in self.fitdata_dict

    def get(self, group):
        """return fit data"""
        fitdata = self.fitdata_dict.get(group, None)
        return fitdata
    def getall(self):
        """return all"""
        return [self.fitdata_dict[group] for group in self.groups]

def prepare_fitting_genedata(fitdata, fit_start, xmax):
    """Prepare fitting xvalues and yvalues"""
    xvalues = []
    yvalues = []
    diver = 0
    for row in fitdata:
        if row.corrtype == 'P2' and row.lag >= fit_start and row.lag <= xmax:
            xvalues.append(row.lag)
            yvalues.append(row.value)
        elif row.lag == 0:
            diver = row.value
    xvalues = numpy.array(xvalues)
    yvalues = numpy.array(yvalues)
    return (xvalues, yvalues, diver)


