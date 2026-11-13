#!/usr/bin/env python

# Copyright 2024 Informatics Matters Ltd.
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


# Generates MACCS keys fingerprints for molecules, reproducing jaqpotpy 1.0.x.
#
# The MACCS definition is a fixed 167-key set, so there are no value-affecting
# parameters to expose (and none for the predictor to check). 9 of the legacy
# models use MACCS, which `main` had no job for; this fills that gap.


from dm_job_utilities.dm_log import DmLog

import jaqpot_featurizers as jf
from descriptor_generator import AbstractCalculator, get_base_parser


def main():

    # Usage:
    # python src/fp_maccs.py -i data/10.smi -o descriptors-maccs.smi -d tab --write-header

    parser = get_base_parser()

    args = parser.parse_args()
    DmLog.emit_event("descriptor_calc: ", args)

    class MACCSFPCalculator(AbstractCalculator):
        calculator = staticmethod(jf.maccs_bitstring)
        descriptor_names = ["MACCSKeys"]
        manifest_name = jf.MACCS
        binary_output = True

        def calculate(self, *calc_args):
            mol = calc_args[0]
            return (jf.maccs_bitstring(mol),)

    calc = MACCSFPCalculator(args)
    calc.run()


if __name__ == "__main__":
    main()
