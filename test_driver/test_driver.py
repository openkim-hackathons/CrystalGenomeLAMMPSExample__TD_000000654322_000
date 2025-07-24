"""
Example OpenKIM Crystal Genome Test Driver
==========================================

This is an example demonstrating usage of the kim-tools package. See
https://kim-tools.readthedocs.io for more information.

"""

import subprocess
import os
from kim_tools import SingleCrystalTestDriver, get_stoich_reduced_list_from_prototype
from tempfile import NamedTemporaryFile


def run_lammps(atoms, species, lmp_cmd, model, atom_style):
    LAMMPS_INPUT = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "in.lammps"
    )
    with NamedTemporaryFile() as f_config:
        atoms.write(
            f_config.name, format="lammps-data", atom_style=atom_style, masses=True
        )
        f_config.seek(0)
        species_str = " ".join(species)
        with NamedTemporaryFile() as f_out:
            subprocess.run(
                f"{lmp_cmd} -i {LAMMPS_INPUT} -var model {model} "
                f'-var datafile {f_config.name} -var species "{species_str}" '
                f"-var outfile {f_out.name}",
                shell=True,
            )
            f_out.seek(0)
            return float(f_out.read())


class TestDriver(SingleCrystalTestDriver):
    def _calculate(self, lmp_cmd="lmp", **kwargs) -> None:
        """
        Recompute binding energy using LAMMPS
        """
        atoms = self._get_atoms()
        species = self._get_nominal_crystal_structure_npt()["stoichiometric-species"][
            "source-value"
        ]
        model = self.kim_model_name

        # Always need to try "atomic" and "charge" atom
        # styles because we don't know what the SM requires
        try:
            energy = run_lammps(
                atoms=atoms,
                species=species,
                lmp_cmd=lmp_cmd,
                model=model,
                atom_style="atomic",
            )
        except Exception:
            energy = run_lammps(
                atoms=atoms,
                species=species,
                lmp_cmd=lmp_cmd,
                model=model,
                atom_style="charge",
            )

        # Subtract out isolated energy
        energy_per_atom = energy / len(atoms)
        prototype_label = self._get_nominal_crystal_structure_npt()["prototype-label"][
            "source-value"
        ]
        stoichiometry = get_stoich_reduced_list_from_prototype(prototype_label)
        num_atoms_in_formula = sum(stoichiometry)
        energy_per_formula = energy_per_atom * num_atoms_in_formula
        isolated_energy_per_formula = 0
        for species, num in zip(species, stoichiometry):
            spec_isol_energ = self.get_isolated_energy_per_atom(species)
            isolated_energy_per_formula += spec_isol_energ * num

        binding_energy_per_formula = energy_per_formula - isolated_energy_per_formula
        binding_energy_per_atom = binding_energy_per_formula / num_atoms_in_formula

        # WRITE OUT PROPERTY
        self._add_property_instance_and_common_crystal_genome_keys(
            property_name="binding-energy-crystal",
            write_stress=False,
            write_temp=False,
        )
        self._add_key_to_current_property_instance(
            "binding-potential-energy-per-atom", binding_energy_per_atom, "eV"
        )
        self._add_key_to_current_property_instance(
            "binding-potential-energy-per-formula", binding_energy_per_formula, "eV"
        )
