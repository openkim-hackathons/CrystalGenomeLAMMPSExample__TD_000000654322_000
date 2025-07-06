"""
Example OpenKIM Crystal Genome Test Driver
==========================================

This is an example demonstrating usage of the kim-tools package. See
https://kim-tools.readthedocs.io for more information.

"""

import shutil
from os import path

from kim_tools import (
    SingleCrystalTestDriver,
    get_stoich_reduced_list_from_prototype,
    minimize_wrapper,
)


class TestDriver(SingleCrystalTestDriver):
    def _calculate(
        self, max_volume_scale: float = 1e-2, num_steps: int = 10, **kwargs
    ) -> None:
        """
        Computes the energy vs. volume curve for isotropic expansion and compression.

        Args:
            max_volume_scale:
                Maximum fractional change in volume to investigate
            num_steps:
                Number of steps to take in each direction
        """
        # The base class provides ``self._get_atoms()``, which provides a starting
        # point for your calculations. This Atoms object is a primitive cell in the
        # setting defined in https://doi.org/10.1016/j.commatsci.2017.01.017.
        #
        # Perform your calculations using this copy. If using ASE, the calculator is
        # already attached. Note that if you invoke the Test Driver using an Atoms
        # object, this is not the same object. It is a re-generated object produced
        # from the symmetry analysis of the input Atoms object. Additionally, any
        # changes made to this Atoms object will not be automatically reflected in the
        # output of this Test Driver. If your Test Driver changes the nominal crystal
        # structure (through equilibration or relaxation), see how to report that in
        # the comments preceding the invocation of
        # ``self._update_nominal_parameter_values()`` below.
        original_atoms = self._get_atoms()
        original_cell = original_atoms.get_cell()
        num_atoms = len(original_atoms)

        # The base class also provides an instance of the OpenKIM
        # ``crystal-structure-npt`` property. This is a symmetry-reduced description of
        # the nominal crystal structure. See the definition of this property at
        # https://openkim.org/properties/show/2023-02-21/staff@noreply.openkim.org/crystal-structure-npt
        # For a general explanation of the structure of a KIM Property Instance, see
        # https://openkim.org/doc/schema/properties-framework/
        # To access a copy of the property instance as a Python dictionary, use
        # ``self._get_nominal_crystal_structure_npt()```

        # Here we will use the prototype label to infer the number of atoms per
        # stoichiometric formula
        prototype_label = self._get_nominal_crystal_structure_npt()["prototype-label"][
            "source-value"
        ]

        num_atoms_in_formula = sum(
            get_stoich_reduced_list_from_prototype(prototype_label)
        )

        # Because temperature and stress are such commonly used parameters, we provide
        # utility functions to access them. This example Test Driver does not use them,
        # this is for demonstration. The units requested here are the defaults, they
        # will be what you get if you omit the ``unit`` argument.
        print(f'Temperature (K): {self._get_temperature(unit="K")}')
        print(
            "Stress (eV/angstrom^3): "
            + str(self._get_cell_cauchy_stress(unit="eV/angstrom^3"))
        )

        binding_potential_energy_per_atom = []
        binding_potential_energy_per_formula = []
        volume_per_atom = []
        volume_per_formula = []
        step_size = max_volume_scale / num_steps
        disclaimer = None

        print("\nPerforming energy scan...\n")

        for i in range(-num_steps, num_steps + 1):
            volume_scale = 1 + step_size * i
            linear_scale = volume_scale ** (1 / 3)
            # Get atoms again
            atoms = self._get_atoms()
            atoms.set_cell(original_cell * linear_scale, scale_atoms=True)
            volume = atoms.get_volume()
            current_volume_per_atom = volume / num_atoms
            problem_occurred = False
            try:
                # The Atoms object returned by ``self._get_atoms()`` comes
                # pre-initialized with an ASE calculator. If you need to access the
                # name of the KIM model (for example, if you are exporting the atomic
                # configuration to run an external simulator like LAMMPS), it can
                # be accessed at ``self.kim_model_name``
                minimize_wrapper(atoms, variable_cell=False)

                # ``self._verify_unchanged_symmetry()`` is used to check that the
                # material has not undergone a phase transition.
                # It is OpenKIM convention that a run of SingleCrystalTestDriver should
                # not involve phase transitions, so here any points where the crystal
                # symmetry has changed are not recorded in the curve.
                if self._verify_unchanged_symmetry(atoms):
                    potential_energy = atoms.get_potential_energy()
                    stress = atoms.get_stress()
                    print(f"Volume: {volume} Energy: {potential_energy}")
                    # If your Test Driver changes the nominal crystal structure (e.g.
                    # through relaxation or MD) AND you intend to write the changed
                    # structure as a property instance, you must update the nominal
                    # crystal structure before writing properties. This is done by
                    # providing an Atoms object to
                    # ``self._update_nominal_crystal_structure()``.
                    # SingleCrystalTestDriver expects the crystal to maintain the same
                    # space group and occupied Wyckoff positions, meaning only the free
                    # parameters of the crystal unconstrained by symmetry are allowed
                    # to change. You must provide an Atoms object that is a primitive
                    # cell in the same setting. Translations, permutations, and rigid
                    # rotations of the unit cell are permissible, but the identity of
                    # the lattice vectors may not change during the property computation
                    # (e.g. lattice vectors a and c may not exchange places, even if
                    # allowed by symmetry). ``_update_nominal_parameter_values``
                    # raises an exception if a symmetry change is detected. Here, it
                    # will be caught and the Test Driver will move on to the next
                    # structure evaluation. If your Test Driver only deals with one
                    # structure and it changes symmetry, you should allow it to fail.
                    self._update_nominal_parameter_values(atoms)
                    # Normally, you would use either ``_verify_unchanged_symmetry``
                    # (if you don't intend to write the changed structure), or
                    # ``_update_nominal_parameter_values`` (if you do), not both, as
                    # done here for demonstration.
                else:
                    print("Atoms underwent symmetry change")
                    problem_occurred = True
            except Exception:
                print(
                    "Failed to get energy, stress, or parameter values "
                    f"at volume {volume}"
                )
                problem_occurred = True

            if problem_occurred:
                disclaimer = (
                    "At least one of the requested deformations of the unit cell "
                    "underwent a phase transformation or otherwise failed to evaluate."
                )
            else:
                current_binding_potential_energy_per_atom = potential_energy / num_atoms
                volume_per_atom.append(current_volume_per_atom)
                volume_per_formula.append(
                    current_volume_per_atom * num_atoms_in_formula
                )
                binding_potential_energy_per_atom.append(
                    current_binding_potential_energy_per_atom
                )
                binding_potential_energy_per_formula.append(
                    current_binding_potential_energy_per_atom * num_atoms_in_formula
                )

                # Because we have evaluated a structure at some combination of
                # temperature and/or stress (only stress in this case), we might as
                # well write an instance of the "crystal-structure-npt" property at
                # every loop iteration. The newly defined
                # "energy-vs-volume-isotropic-crystal" property will be written
                # at the end of the loop.

                # This method initializes the Property Instance and adds the keys
                # common to all Crystal Genome properties. `property_name`` can be the
                # full "property-id" field in the Property Definition, or the
                # "Property Name", which is just the short name after the slash, as
                # used here. The options `write_stress` and `write_temp` can be set to
                # `False` (default), `True`, in which case the internally stored nominal
                # values will be written, or set to a value in case you wish to write
                # a different value. Here we tell the Driver to write a provided stress
                # value and write the nominal temperature (zero).
                # If you specify the value for stress, you must also specify
                # `stress_unit` for the value you have provided.
                # For temperature, `temp_unit` defaults to K. See the API documentation
                # for the method for more info.
                self._add_property_instance_and_common_crystal_genome_keys(
                    property_name="crystal-structure-npt",
                    write_stress=stress,
                    write_temp=True,
                    stress_unit="eV/angstrom^3",
                )
                # Because "crystal-structure-npt" has no additional custom fields, we
                # are done with this property

        # Now it is time to write the output in the format you created in your Property
        # Definition. Here we omit stress and temperature because we did not include
        # them in our Definition, and optionally add a disclaimer. First, since we have
        # been changing the nominal crystal structure, we must reset it to the original
        # structure.
        self._update_nominal_parameter_values(original_atoms)
        self._add_property_instance_and_common_crystal_genome_keys(
            property_name="energy-vs-volume-isotropic-crystal",
            write_stress=False,
            write_temp=False,
            disclaimer=disclaimer,
        )
        # This method adds additional fields to your property instance by specifying the
        # key names you defined in your property definition and providing units if
        # necessary.
        self._add_key_to_current_property_instance(
            "volume-per-atom", volume_per_atom, unit="angstrom^3"
        )
        self._add_key_to_current_property_instance(
            "volume-per-formula", volume_per_formula, unit="angstrom^3"
        )

        # You may also provide a dictionary supplying uncertainty information. It is
        # optional, and normally would not be reported for a deterministic calculation
        # like this, only one involving some statistics, such as molecular dynamics or
        # fitting. There are many possible ways to report uncertainty information,
        # detailed at https://openkim.org/doc/schema/properties-framework/
        uncertainty_info = {
            "source-std-uncert-value": [0] * len(binding_potential_energy_per_atom)
        }
        self._add_key_to_current_property_instance(
            "binding-potential-energy-per-atom",
            binding_potential_energy_per_atom,
            unit="eV",
            uncertainty_info=uncertainty_info,
        )
        self._add_key_to_current_property_instance(
            "binding-potential-energy-per-formula",
            binding_potential_energy_per_formula,
            unit="eV",
            uncertainty_info=uncertainty_info,
        )
        # This adds a "file" type key to the current Property Instance. It will
        # automatically be numbered according to the current Property Instance and
        # moved to the output directory if it is not already there
        # (e.g. cat.jpg will be automatically moved to the path output/cat-1.jpg
        # and reported in the property accordingly.) For this example, we just copy
        # the picture packaged with this Test Driver.
        cat_path = path.join(path.dirname(path.realpath(__file__)), "data", "cat.jpg")
        save_path = "cat.jpg"
        shutil.copy(cat_path, save_path)
        self._add_file_to_current_property_instance("cat-picture", save_path)

        # If your Test Driver reports multiple Property Instances, repeat the process
        # above for each one.
