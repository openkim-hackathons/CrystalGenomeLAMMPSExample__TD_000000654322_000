"""
Example OpenKIM Crystal Genome Test Driver
==========================================

This is an example demonstrating usage of the kim-tools package. See https://kim-tools.readthedocs.io for more information.

"""

from kim_tools import SingleCrystalTestDriver
from kim_tools import get_stoich_reduced_list_from_prototype

class TestDriver(SingleCrystalTestDriver):
    def _calculate(self, max_volume_scale: float = 1e-2, num_steps: int = 10, **kwargs):
        """
        Computes the energy vs. volume curve for isotropic expansion and compression. 

        Args:
            max_volume_scale:
                Maximum fractional change in volume to investigate
            num_steps:
                Number of steps to take in each direction
        """
        # The base class provides self._get_nominal_atoms(), which returns a copy of the internal Atoms object (not accessible directly).
        # This Atoms object is a primitive cell in the setting defined in https://doi.org/10.1016/j.commatsci.2017.01.017
        # Perform your calculations using this copy. If using ASE, the calculator is already attached.
        atoms = self._get_nominal_atoms()
        original_cell = atoms.get_cell()
        num_atoms = len(atoms)

        # The base class also provides an instance of the OpenKIM `crystal-structure-npt` property (not accessible directly).
        # This is a symmetry-reduced description of the crystal structure that is kept in sync with the internal Atoms object.
        # See the definition of this property at 
        # https://openkim.org/properties/show/2023-02-21/staff@noreply.openkim.org/crystal-structure-npt
        # For a general explanation of the structure of a KIM Property Instance, see
        # https://openkim.org/doc/schema/properties-framework/
        # To access a copy of the property instance as a Python dictionary, use `self._get_nominal_crystal_structure_npt()`
        
        # Here we will use the prototype label to infer the number of atoms per stoichiometric formula
        prototype_label = self._get_nominal_crystal_structure_npt()['prototype-label']['source-value']
        
        num_atoms_in_formula = sum(get_stoich_reduced_list_from_prototype(prototype_label))

        binding_potential_energy_per_atom = []
        binding_potential_energy_per_formula = []
        volume_per_atom = []
        volume_per_formula = []
        step_size = max_volume_scale/num_steps
        disclaimer = None

        print('\nPerforming energy scan...\n')

        for i in range(-num_steps,num_steps+1):
            volume_scale = 1 + step_size*i
            linear_scale = volume_scale ** (1/3)
            atoms.set_cell(original_cell*linear_scale,scale_atoms=True)
            volume = atoms.get_volume()
            current_volume_per_atom = volume/num_atoms

            problem_occurred = False
            try:
                # The self.atoms object comes pre-initialized with the calculator set to the interatomic model
                # the Test Driver was called with. If you need to access the name of the KIM model (for example,
                # if you are exporting the atomic configuration to run an external simulator like LAMMPS), it can
                # be accessed at self.kim_model_name
                potential_energy = atoms.get_potential_energy()                
                print('Volume: %5.5f Energy: %5.5f'%(volume,potential_energy))
            except:
                print('Failed to get energy at volume %f'%volume)
                problem_occurred = True
                disclaimer = "At least one of the requested deformations of the unit cell failed to compute a potential energy."                
            
            if not problem_occurred:
                current_binding_potential_energy_per_atom = potential_energy/num_atoms
                volume_per_atom.append(current_volume_per_atom)
                volume_per_formula.append(current_volume_per_atom*num_atoms_in_formula)
                binding_potential_energy_per_atom.append(current_binding_potential_energy_per_atom)
                binding_potential_energy_per_formula.append(current_binding_potential_energy_per_atom*num_atoms_in_formula)

        # Now it is time to write the output in the format you created in your Property Definition. The base class provides utility methods
        # to facilitate this process.
        
        # If your Test Driver changes the nominal crystal structure (e.g. through relaxation or MD), you must update the nominal
        # crystal structure before writing properties. This is done by providing an Atoms object to self._update_nominal_crystal_structure().
        # SingleCrystalTestDriver expects the crystal to maintain the same space group and occupied Wyckoff positions, meaning only
        # the free parameters of the crystal unconstrained by symmetry are allowed to change. You must provide an Atoms object that
        # is a primitive cell in the same setting. Translations, permutations, and rigid rotations of the unit cell are permissible,
        # but the identity of the lattice vectors may not change during the property computation (e.g. lattice vectors a and c may
        # not exchange places, even if allowed by symmetry)
        
        # This Test Driver does not actually need to do this, because the nominal structure does not change -- the energy-vs-volume
        # curve is defined w.r.t. the original undeformed structure. For demonstration, we return the atoms object to its original 
        # cell and update
        atoms.set_cell(original_cell,scale_atoms=True)
        self._update_nominal_crystal_structure
        
        # This method initializes the Property Instance and adds the keys common to all Crystal Genome properties.
        # property_name can be the full "property-id" field in your Property Definition, or the "Property Name",
        # which is just the short name after the slash, as used here. You can also specify whether your property
        # includes stress and temperature (no by default), and have the option to specify a disclaimer.
        self._add_property_instance_and_common_crystal_genome_keys("energy-vs-volume-isotropic-crystal",
                                                                   write_stress=False, write_temp=False, disclaimer=disclaimer)

        # This method adds additional fields to your property instance by specifying the key names you defined
        # in your property definition and providing units if necessary.
        self._add_key_to_current_property_instance("volume-per-atom",volume_per_atom,
                                                   units="angstrom^3")
        self._add_key_to_current_property_instance("volume-per-formula",volume_per_formula,
                                                   units="angstrom^3")

        # You may also provide a dictionary supplying uncertainty information. It is optional, and normally
        # would not be reported for a deterministic calculation like this, only one involving some statistics,
        # such as molecular dynamics or fitting. There are many possible ways to report uncertainty information,
        # detailed at https://openkim.org/doc/schema/properties-framework/
        uncertainty_info = {
            "source-std-uncert-value": [0]*len(binding_potential_energy_per_atom)
        }
        self._add_key_to_current_property_instance("binding-potential-energy-per-atom",binding_potential_energy_per_atom,
                                                   units="eV",uncertainty_info=uncertainty_info)
        self._add_key_to_current_property_instance("binding-potential-energy-per-formula",binding_potential_energy_per_formula,
                                                   units="eV",uncertainty_info=uncertainty_info)
        
        # If your Test Driver reports multiple Property Instances, repeat the process above for each one.
