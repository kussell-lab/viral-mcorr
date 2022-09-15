# Copyright 2022 Asher Preska Steinberg
#
# Initially modified from corr_res.py of mcorr
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

class CorrRes(object):
    """One correlation result"""
    def __init__(self, terms):
        lag = float(terms[0])
        value = float(terms[1])
        variance = float(terms[2])
        num = float(terms[3])
        corrtype = terms[4]
        group = terms[5]
        self.lag = lag
        self.value = value
        self.variance = variance
        self.num = num
        self.corrtype = corrtype
        self.group = group

def read_corr(csv_file):
    """Read corr results in a csv file"""
    results = []
    with open(csv_file, 'r') as infile:
        for line in infile:
            if line.startswith('#'): continue
            terms = line.rstrip().split(",")
            if terms[0] == 'l': continue
            results.append(CorrRes(terms))
    return results

class GeneCorrRes(object):
    """One correlation result"""
    def __init__(self, terms):
        lag = float(terms[0])
        value = float(terms[1])
        variance = float(terms[3])
        num = float(terms[2])
        corrtype = terms[4]
        group = terms[5]
        self.lag = lag
        self.value = value
        self.variance = variance
        self.num = num
        self.corrtype = corrtype
        self.group = group

def read_genecorr(csv_file):
    """Read corr results in a csv file"""
    results = []
    with open(csv_file, 'r') as infile:
        for line in infile:
            if line.startswith('#'): continue
            terms = line.rstrip().split(",")
            if terms[0] == 'lag': continue
            results.append(CorrRes(terms))
    return results


