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


# Generates Morgan fingerprints for molecules, reproducing jaqpotpy 1.0.x.
#
# LEGACY branch: this uses the deprecated GetMorganFingerprintAsBitVect (via
# jaqpot_featurizers), NOT the newer rdFingerprintGenerator Morgan generator
# that `main` uses. The two do not always produce the same bits, and the
# legacy models were trained on the former, so this reproduces them exactly.
# The default radius is 2 (jaqpotpy's default), where `main` defaulted to 3.


from dm_job_utilities.dm_log import DmLog

import jaqpot_featurizers as jf
from descriptor_generator import AbstractCalculator, get_base_parser


def main():

    # Usage:
    # python src/fp_morgan.py -i data/10.smi -o descriptors-morganfp.smi -d tab --write-header

    parser = get_base_parser()

    morganfp_group = parser.add_argument_group("Morgan fingerprint options")
    morganfp_group.add_argument(
        "--morganfp-radius",
        default=2,
        type=int,
        help="Morgan radius (jaqpotpy default is 2)",
    )
    morganfp_group.add_argument(
        "--morganfp-size",
        default=2048,
        type=int,
        help="Number of bits in the fingerprint",
    )
    morganfp_group.add_argument(
        "--morganfp-use-chirality",
        action="store_true",
        help="If set, chirality information is added to the fingerprint",
    )
    morganfp_group.add_argument(
        "--morganfp-no-bond-types",
        action="store_true",
        help="If set, bond types are NOT included in the bond invariants "
        "(jaqpotpy includes them, which is the default here)",
    )
    morganfp_group.add_argument(
        "--morganfp-use-features",
        action="store_true",
        help="If set, use feature invariants (FCFP-style) instead of the "
        "default connectivity invariants (ECFP-style)",
    )

    args = parser.parse_args()
    DmLog.emit_event("descriptor_calc: ", args)

    radius = args.morganfp_radius
    n_bits = args.morganfp_size
    chiral = args.morganfp_use_chirality
    bond_types = not args.morganfp_no_bond_types
    features = args.morganfp_use_features

    class MorganFPCalculator(AbstractCalculator):
        calculator = staticmethod(jf.morgan_bitstring)
        descriptor_names = ["MorganFingerprint"]
        manifest_name = jf.MORGAN
        binary_output = True

        def calculate(self, *calc_args):
            mol = calc_args[0]
            return (
                jf.morgan_bitstring(
                    mol,
                    radius=radius,
                    n_bits=n_bits,
                    chiral=chiral,
                    bond_types=bond_types,
                    features=features,
                ),
            )

        def manifest_params(self):
            return jf.morgan_params(
                radius=radius,
                n_bits=n_bits,
                chiral=chiral,
                bond_types=bond_types,
                features=features,
            )

    calc = MorganFPCalculator(args)
    calc.run()


if __name__ == "__main__":
    main()
