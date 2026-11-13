"""Reference-value tests for the legacy featurizers.

These freeze the exact output of ``jaqpot_featurizers`` under the branch's
pinned rdkit (2023.9.6), so the build fails loudly if anything drifts from what
jaqpotpy produced — for example an accidental switch to the newer Morgan
generator API, a changed radius, or ``Ipc`` losing its averaging. The golden
on-bit lists and values were generated once with rdkit 2023.9.6, the same rdkit
whose Morgan output reproduces the legacy CYP2D6 model's reference predictions
exactly (50/50) in the predictor repo.
"""

import pytest
from rdkit import Chem

import jaqpot_featurizers as jf

# Golden, frozen under rdkit 2023.9.6.
ETHANOL_MORGAN_ONBITS = [80, 222, 294, 807, 1057, 1410]
ETHANOL_MACCS_ONBITS = [82, 109, 114, 139, 153, 155, 157, 160, 164]


def _onbits(bitstring):
    return [i for i, b in enumerate(bitstring) if b == "1"]


def _mol(smiles):
    return Chem.MolFromSmiles(smiles)


# --- Morgan -----------------------------------------------------------------


def test_morgan_is_2048_bits_radius_2_by_default():
    bits = jf.morgan_bitstring(_mol("CCO"))
    assert len(bits) == 2048
    assert _onbits(bits) == ETHANOL_MORGAN_ONBITS


def test_morgan_radius_changes_the_bits():
    # radius is a value-affecting knob: the whole reason the manifest records it
    aspirin = _mol("CC(=O)Oc1ccccc1C(=O)O")
    assert jf.morgan_bitstring(aspirin, radius=2) != jf.morgan_bitstring(
        aspirin, radius=3
    )


def test_morgan_params_names_match_the_manifest():
    assert jf.morgan_params() == {
        "radius": 2,
        "n_bits": 2048,
        "chiral": False,
        "bond_types": True,
        "feature_invariants": False,
    }


# --- MACCS ------------------------------------------------------------------


def test_maccs_is_167_bits():
    bits = jf.maccs_bitstring(_mol("CCO"))
    assert len(bits) == 167
    assert _onbits(bits) == ETHANOL_MACCS_ONBITS


# --- RDKit descriptors ------------------------------------------------------


def test_rdkit_descriptor_values_for_ethanol():
    names = jf.rdkit_descriptor_names()
    values = dict(zip(names, jf.rdkit_descriptor_values(_mol("CCO"))))
    assert values["MolWt"] == pytest.approx(46.069, abs=1e-3)
    assert values["NumHDonors"] == 1
    assert values["TPSA"] == pytest.approx(20.23, abs=1e-2)


def test_ipc_is_averaged_not_astronomical():
    # unaveraged Ipc for aspirin is ~730; averaged it is ~2. jaqpotpy used the
    # averaged form, and MolecularDescriptorCalculator would NOT — this guards
    # that distinction.
    names = jf.rdkit_descriptor_names()
    values = dict(zip(names, jf.rdkit_descriptor_values(_mol("CC(=O)Oc1ccccc1C(=O)O"))))
    assert values["Ipc"] == pytest.approx(2.0325, abs=1e-3)


def test_use_fragment_drops_the_fr_columns():
    with_fr = jf.rdkit_descriptor_names(use_fragment=True)
    without_fr = jf.rdkit_descriptor_names(use_fragment=False)
    assert any(name.startswith("fr_") for name in with_fr)
    assert not any(name.startswith("fr_") for name in without_fr)
    assert len(without_fr) < len(with_fr)


# --- normalisation ----------------------------------------------------------


def test_normalize_makes_sdf_style_hs_match_smiles():
    # a molecule with explicit Hs must featurise like the plain SMILES
    plain = _mol("CCO")
    explicit = Chem.AddHs(_mol("CCO"))
    assert jf.morgan_bitstring(explicit) == jf.morgan_bitstring(plain)


def test_normalize_raises_on_unusable_molecule():
    with pytest.raises(ValueError):
        jf.normalize(Chem.Mol())  # empty molecule -> empty SMILES -> unusable
