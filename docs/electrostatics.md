# Electrostatics

`MartiniTopFile.create_system()` supports two electrostatics treatments:

- `reaction-field`, the default Martini cutoff-style electrostatics used by the original implementation.
- `pme`, which keeps Martini Lennard-Jones interactions in the existing custom force and moves Coulomb interactions to an OpenMM `NonbondedForce` using PME.

The default remains `reaction-field`, so existing scripts keep the previous behavior unless `electrostatics_method` is explicitly set.

## Usage

```python
system = top.create_system(
    nonbonded_cutoff=1.1 * unit.nanometer,
    electrostatics_method="pme",
    ewald_error_tolerance=5e-4,
)
```

The `electrostatics_method` argument accepts `reaction-field`, `rf`, `reactionfield`, and `pme`, case-insensitively. `ewald_error_tolerance` is only used for PME.

## Implementation

For reaction-field electrostatics, the implementation continues to use the original custom nonbonded expression:

- Lennard-Jones and Coulomb reaction-field terms are evaluated together in a `CustomNonbondedForce`.
- Excluded electrostatic self/reaction-field corrections are handled with custom bond forces.
- Explicit pair exceptions remain custom bond terms.

For PME electrostatics:

- The existing `CustomNonbondedForce` evaluates only the Martini Lennard-Jones potential-shifted interaction.
- A separate OpenMM `NonbondedForce` evaluates Coulomb interactions with PME.
- Charges are scaled by `1 / sqrt(epsilon_r)` before they are added to the PME force. This gives the same Coulomb prefactor as using `q_i q_j / epsilon_r`.
- Topology exclusions and pair exceptions are added to both the Lennard-Jones custom force and the PME force.
- Virtual-site particles are treated as ordinary topology particles for nonbonded parameters. If a virtual-site atom has a charge in the topology, that charge is added to the PME force and its topology exclusions are honored.

## Validation

We compared single-frame rerun energies and forces between GROMACS and OpenMM for several cases.

### Martini PME

Two Martini polarizable systems were tested with PME: pure refPOL water (`refPOL2797`) and refPOL water with polarizable NaCl (`pnapcl_aq_0.10mol`). These systems include off-center charged particles and explicit exclusions.

| System | Ewald tolerance | OpenMM - GROMACS energy (kJ/mol) | Coulomb contribution to mismatch (kJ/mol) | LJ mismatch (kJ/mol) | Force RMS vector delta (kJ/mol/nm) |
|---|---:|---:|---:|---:|---:|
| refPOL2797 | 5e-4 | +2.204 | +2.251 | -0.041 | 0.077 |
| pnapcl_aq_0.10mol | 5e-4 | +1.969 | +2.057 | -0.085 | 0.081 |
| refPOL2797 | 1e-5 | +6.287 | +6.329 | -0.041 | 0.049 |
| pnapcl_aq_0.10mol | 1e-5 | +7.287 | +7.379 | -0.085 | 0.062 |

The residual is dominated by Coulomb/PME terms. Lennard-Jones and angle terms agree closely.

### Stock OpenMM/GROMACS water PME

The same comparison was repeated with stock OpenMM `GromacsTopFile` on a conventional SPC/E water box. This test does not use Martini topology handling.

| System | Ewald tolerance | OpenMM - GROMACS full energy (kJ/mol) | OpenMM - GROMACS Coulomb (kJ/mol) | OpenMM - GROMACS LJ (kJ/mol) | Force RMS vector delta (kJ/mol/nm) |
|---|---:|---:|---:|---:|---:|
| SPC/E water box | 5e-4 | -1.915 | -1.911 | -0.001 | 0.792 |
| SPC/E water box | 1e-5 | +2.135 | +2.137 | -0.001 | 0.492 |

This showed that the observed PME residual is not specific to `martini_openmm`. It is a generic OpenMM/GROMACS PME parity difference under these settings.

### Charged virtual sites and exclusions

A minimal Martini-like topology was constructed with two charged `virtual_sites2` particles and explicit parent-vsite exclusions. Lennard-Jones parameters were set to zero so any residual was electrostatic.

| Test | Ewald tolerance | OpenMM - GROMACS energy (kJ/mol) | Coulomb mismatch (kJ/mol) | LJ mismatch (kJ/mol) | Parent-vsite PME exceptions found | Force RMS vector delta (kJ/mol/nm) |
|---|---:|---:|---:|---:|---:|---:|
| Charged virtual sites | 5e-4 | -0.0053 | -0.0053 | 0.0000 | 4 / 4 | 0.141 |
| Charged virtual sites | 1e-5 | +0.1284 | +0.1284 | 0.0000 | 4 / 4 | 0.427 |

This test did not show evidence of missing virtual-site exclusions. OpenMM created all expected parent-vsite PME exceptions.

### Reaction-field cutoff check

The normal Martini reaction-field path was also checked with `epsilon_r = 15`.

| System | Electrostatics | epsilon_r | OpenMM - GROMACS energy (kJ/mol) | Relative energy difference | Force RMS vector delta (kJ/mol/nm) |
|---|---|---:|---:|---:|---:|
| refPOL2797 | reaction-field | 15 | -0.361 | 1.07e-5 | 0.00106 |
| pnapcl_aq_0.10mol | reaction-field | 15 | -0.440 | 1.23e-5 | 0.00115 |

The reaction-field totals and forces remain close to GROMACS. The PME implementation does not replace or alter the default reaction-field path.
