"""
Example OpenKIM Crystal Genome Test Driver
==========================================

This is an example demonstrating usage of the kim-tools package. See
https://kim-tools.readthedocs.io for more information.

"""
from kim_tools import (
    SingleCrystalTestDriver
)
from tempfile import NamedTemporaryFile
from lammps import lammps
from mpi4py import MPI


class TestDriver(SingleCrystalTestDriver):
    def _calculate(
        self, lmp_cmd="lmp", **kwargs
    ) -> None:
        """
        Recompute binding energy using LAMMPS
        """
        with NamedTemporaryFile() as f_config:
            self._get_atoms().write(
                f_config.name,
                format="lammps-data",
                atom_style="charge",
                masses=True
            )
            f_config.seek(0)
            species = self._get_nominal_crystal_structure_npt()["stoichiometric-species"]["source-value"]
            with NamedTemporaryFile() as f_out:
                lmp_in = f"""
                boundary        p p p
                kim             init {self.kim_model_name} metal unit_conversion_mode
                read_data       {f_config.name}
                change_box      all triclinic 
                kim             interactions {" ".join(species)}

                variable        E_eV equal pe/${{_u_energy}}
                thermo_modify   norm yes

                run 0

                print           "${{E_eV}}" file {f_out.name}
                """
                l = lammps()
                l.commands_string(lmp_in)
                f_out.seek(0)
                energy = float(f_out.read())

        print(f"ENERGY!!! {energy}")
            