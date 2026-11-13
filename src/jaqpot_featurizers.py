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

"""Featurizers that reproduce jaqpotpy 1.0.x's descriptor calculations
bit-for-bit, for the legacy environment (rdkit 2023.9.6, numpy 1.26).

This is the whole reason the ``legacy`` branch exists. The pre-trained
jaqpotpy models we need to serve were trained on features computed by specific,
old rdkit calls; a newer rdkit, or a different call (e.g. the newer
``rdFingerprintGenerator`` Morgan generator in place of the deprecated
``GetMorganFingerprintAsBitVect``), silently changes the numbers. So these
functions deliberately mirror what jaqpotpy did:

- **Morgan**: ``GetMorganFingerprintAsBitVect`` (NOT the generator API),
  radius 2, 2048 bits, ``useBondTypes=True``.
- **MACCS**: ``GetMACCSKeysFingerprint`` (167 keys).
- **RDKit descriptors**: ``Descriptors.descList`` with ``Ipc`` computed in its
  averaged form (``avg=True``) — unaveraged ``Ipc`` is orders of magnitude
  larger and is what ``MolecularDescriptorCalculator`` would give.

Every molecule is normalised the way jaqpotpy did before featurising (strip
explicit Hs, canonicalise to SMILES and re-parse), so SDF input with explicit
hydrogens or 3D coordinates yields the same features as the equivalent SMILES.

These calls match the predictor's proven featurizers, which reproduce the
legacy models' reference predictions exactly (see the predictor repo).
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# the predictor's featurizer names, written into the manifest's featurizer.name
MORGAN = "morgan"
MACCS = "maccs"
RDKIT_DESCRIPTORS = "rdkit_descriptors"


def normalize(mol):
    """Bring a molecule to what parsing its canonical SMILES would give.

    jaqpotpy featurised ``MolFromSmiles(MolToSmiles(RemoveHs(mol)))``, so a
    molecule read from an SDF (often with explicit Hs / 3D coordinates) matches
    the same molecule read from SMILES. Raises if the molecule can't be
    processed, which the caller turns into a flagged (empty) output row.
    """
    normalized = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.RemoveHs(mol)))
    if normalized is None or normalized.GetNumAtoms() == 0:
        raise ValueError("molecule could not be normalised")
    return normalized


# --- Morgan -----------------------------------------------------------------


def morgan_bitstring(
    mol,
    radius=2,
    n_bits=2048,
    chiral=False,
    bond_types=True,
    features=False,
):
    """A single 0/1 bit string for the Morgan fingerprint, jaqpotpy-style."""
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        normalize(mol),
        radius,
        nBits=n_bits,
        useChirality=chiral,
        useBondTypes=bond_types,
        useFeatures=features,
    )
    return fp.ToBitString()


def morgan_params(radius=2, n_bits=2048, chiral=False, bond_types=True, features=False):
    """The value-affecting Morgan knobs, named as the manifest/predictor want
    them (``feature_invariants`` is rdkit's ``useFeatures``)."""
    return {
        "radius": radius,
        "n_bits": n_bits,
        "chiral": chiral,
        "bond_types": bond_types,
        "feature_invariants": features,
    }


# --- MACCS ------------------------------------------------------------------


def maccs_bitstring(mol):
    """A single 167-bit 0/1 bit string for the MACCS keys."""
    return rdMolDescriptors.GetMACCSKeysFingerprint(normalize(mol)).ToBitString()


# --- RDKit descriptors ------------------------------------------------------


def rdkit_descriptor_list(use_fragment=True):
    """(name, function) pairs in the order jaqpotpy produced them, optionally
    dropping the ``fr_*`` fragment-count descriptors."""
    return [
        (name, fn)
        for name, fn in Descriptors.descList
        if use_fragment or not name.startswith("fr_")
    ]


def rdkit_descriptor_names(use_fragment=True):
    """Just the descriptor column names, in order."""
    return [name for name, _ in rdkit_descriptor_list(use_fragment)]


def rdkit_descriptor_values(mol, use_fragment=True, ipc_avg=True, missing_val=None):
    """Descriptor values for one molecule, in ``rdkit_descriptor_names`` order.

    A descriptor that raises yields ``missing_val`` for that cell (the fill the
    predictor records in the manifest), rather than failing the whole row.
    """
    mol = normalize(mol)
    values = []
    for name, fn in rdkit_descriptor_list(use_fragment):
        try:
            if name == "Ipc" and ipc_avg:
                values.append(fn(mol, avg=True))
            else:
                values.append(fn(mol))
        except Exception:  # pylint: disable=broad-except
            values.append(missing_val)
    return values


def rdkit_descriptor_params(use_fragment=True, ipc_avg=True):
    """The value-affecting RDKit-descriptor knobs, as the manifest wants them."""
    return {"use_fragment": use_fragment, "ipc_avg": ipc_avg}
