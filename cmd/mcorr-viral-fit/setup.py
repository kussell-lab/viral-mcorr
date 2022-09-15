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

from setuptools import setup

# read requirements.
requirements = []
with open("requirements.txt", 'rU') as reader:
    for line in reader:
        requirements.append(line.strip())

setup(name='viralmcorr',
        python_requires='>=3',
        version='20220318',
        description='Inferring recombination rates from correlation profiles for viruses',
        url='https://github.com/apsteinberg/mcorr',
        author='Asher Preska Steinberg',
        author_email='apsteinberg@nyu.edu',
        license='Apache 2.0',
        packages=['viralmcorr'],
        install_requires=requirements,
        entry_points = {
            'console_scripts' : ['mcorr-viral-fit=viralmcorr.FitComparison:main', 'mcorr-pair-fit=viralmcorr.FitPairs:main'],
            },
        zip_safe=False)
