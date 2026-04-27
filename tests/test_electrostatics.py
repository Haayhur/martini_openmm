import math
import tempfile
import unittest
from pathlib import Path

import openmm as mm
import openmm.unit as unit

import martini_openmm as martini


TOPOLOGY = """[ defaults ]
1 1

[ atomtypes ]
; name mass charge ptype c6 c12
D 24.0 0.000 A 0.0 0.0

[ moleculetype ]
; molname nrexcl
TST 1

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1 D 1 TST A 1  1.0 24.0
2 D 1 TST B 2 -1.0 24.0

[ exclusions ]
1 2

[ system ]
test system

[ molecules ]
TST 1
"""


ONE_PARENT_VSITE_TOPOLOGY = """[ defaults ]
1 1

[ atomtypes ]
; name mass charge ptype c6 c12
D 24.0 0.000 A 0.0 0.0

[ moleculetype ]
; molname nrexcl
VST 1

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1 D 1 VST A 1  1.0 24.0
2 D 1 VST V 2 -1.0  0.0

[ virtual_sitesn ]
; site funct constructing atom indices
2 1 1

[ system ]
one-parent vsite system

[ molecules ]
VST 1
"""


def make_topology(tmpdir, epsilon_r=4.0, topology_text=TOPOLOGY):
    path = Path(tmpdir) / "system.top"
    path.write_text(topology_text)
    box = (
        mm.Vec3(2.0, 0.0, 0.0) * unit.nanometer,
        mm.Vec3(0.0, 2.0, 0.0) * unit.nanometer,
        mm.Vec3(0.0, 0.0, 2.0) * unit.nanometer,
    )
    return martini.MartiniTopFile(str(path), periodicBoxVectors=box, epsilon_r=epsilon_r)


def forces_by_type(system):
    return {force.__class__.__name__: force for force in system.getForces()}


class TestElectrostaticsMethods(unittest.TestCase):
    def test_reaction_field_is_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            top = make_topology(tmpdir)
            system = top.create_system()

        force_types = forces_by_type(system)
        self.assertIn("CustomNonbondedForce", force_types)
        self.assertNotIn("NonbondedForce", force_types)
        self.assertIsNone(top.es_force)
        self.assertIsNotNone(top.es_self_excl_force)
        self.assertIsNotNone(top.es_except_force)

    def test_pme_uses_scaled_charges_and_exceptions(self):
        epsilon_r = 4.0
        with tempfile.TemporaryDirectory() as tmpdir:
            top = make_topology(tmpdir, epsilon_r=epsilon_r)
            system = top.create_system(
                electrostatics_method="PME",
                ewald_error_tolerance=1e-5,
            )

        nonbonded_forces = [
            force for force in system.getForces() if isinstance(force, mm.NonbondedForce)
        ]
        self.assertEqual(len(nonbonded_forces), 1)
        pme_force = nonbonded_forces[0]
        self.assertEqual(pme_force.getNonbondedMethod(), mm.NonbondedForce.PME)
        self.assertAlmostEqual(pme_force.getEwaldErrorTolerance(), 1e-5)

        charge0, _, _ = pme_force.getParticleParameters(0)
        charge1, _, _ = pme_force.getParticleParameters(1)
        expected_scale = 1.0 / math.sqrt(epsilon_r)
        self.assertAlmostEqual(
            charge0.value_in_unit(unit.elementary_charge),
            expected_scale,
        )
        self.assertAlmostEqual(
            charge1.value_in_unit(unit.elementary_charge),
            -expected_scale,
        )

        self.assertEqual(pme_force.getNumExceptions(), 1)
        i, j, charge_product, _, _ = pme_force.getExceptionParameters(0)
        self.assertEqual((i, j), (0, 1))
        self.assertAlmostEqual(
            charge_product.value_in_unit(unit.elementary_charge**2),
            0.0,
        )

    def test_electrostatics_aliases(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            top = make_topology(tmpdir)
            top.create_system(electrostatics_method="rf")
            self.assertIsNone(top.es_force)

        with tempfile.TemporaryDirectory() as tmpdir:
            top = make_topology(tmpdir)
            system = top.create_system(electrostatics_method="pme")
            self.assertTrue(
                any(isinstance(force, mm.NonbondedForce) for force in system.getForces())
            )

    def test_one_parent_vsite_exclusion_is_mirrored_to_pme(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            top = make_topology(tmpdir, topology_text=ONE_PARENT_VSITE_TOPOLOGY)
            system = top.create_system(electrostatics_method="pme")

        nonbonded_forces = [
            force for force in system.getForces() if isinstance(force, mm.NonbondedForce)
        ]
        self.assertEqual(len(nonbonded_forces), 1)

        pme_force = nonbonded_forces[0]
        exception_pairs = set()
        charge_products = {}
        for index in range(pme_force.getNumExceptions()):
            i, j, charge_product, _, _ = pme_force.getExceptionParameters(index)
            pair = tuple(sorted((i, j)))
            exception_pairs.add(pair)
            charge_products[pair] = charge_product.value_in_unit(
                unit.elementary_charge**2
            )

        self.assertIn((0, 1), exception_pairs)
        self.assertAlmostEqual(charge_products[(0, 1)], 0.0)


if __name__ == "__main__":
    unittest.main()
