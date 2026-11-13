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


# Generates RDKit 2D descriptors for molecules, reproducing jaqpotpy 1.0.x.
#
# LEGACY branch: unlike `main` (which uses MolecularDescriptorCalculator), this
# computes Ipc in its AVERAGED form (jaqpotpy did; the unaveraged value is
# orders of magnitude larger). See jaqpot_featurizers.


from dm_job_utilities.dm_log import DmLog

import jaqpot_featurizers as jf
from descriptor_generator import AbstractCalculator, get_base_parser


def main():
    # Usage:
    # python src/rdkit_descriptors.py -i data/10.smi -o descriptors-rdkit.smi -d tab --write-header

    parser = get_base_parser()
    desc_group = parser.add_argument_group("RDKit descriptor options")
    desc_group.add_argument(
        "--no-fragment-descriptors",
        action="store_true",
        help="Exclude the fr_* fragment-count descriptors (jaqpotpy includes "
        "them, which is the default here)",
    )
    desc_group.add_argument(
        "--no-ipc-average",
        action="store_true",
        help="Compute Ipc unaveraged (jaqpotpy averages it, the default here)",
    )

    args = parser.parse_args()
    DmLog.emit_event("descriptor_calc: ", args)

    use_fragment = not args.no_fragment_descriptors
    ipc_avg = not args.no_ipc_average
    # a missing_val of None keeps an empty cell; the predictor fills it
    missing_val = args.missing_val

    names = jf.rdkit_descriptor_names(use_fragment=use_fragment)

    class DescriptorCalculator(AbstractCalculator):
        calculator = staticmethod(jf.rdkit_descriptor_values)
        descriptor_names = names
        manifest_name = jf.RDKIT_DESCRIPTORS
        binary_output = False

        def calculate(self, *calc_args):
            mol = calc_args[0]
            return jf.rdkit_descriptor_values(
                mol,
                use_fragment=use_fragment,
                ipc_avg=ipc_avg,
                missing_val=missing_val,
            )

        def manifest_params(self):
            return jf.rdkit_descriptor_params(
                use_fragment=use_fragment, ipc_avg=ipc_avg
            )

    calc = DescriptorCalculator(args)
    calc.run()


if __name__ == "__main__":
    main()
