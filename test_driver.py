"""
Example OpenKIM Crystal Genome Test Driver
==========================================

This is an example demonstrating usage of the kim-test-utils package. See https://kim-test-utils.readthedocs.io for more information.

"""

from kim_test_utils import CrystalGenomeTestDriver
from kim_test_utils import get_stoich_reduced_list_from_prototype

class TestDriver(CrystalGenomeTestDriver):
    def _calculate(self, max_volume_scale: float = 1e-2, num_steps: int = 10, **kwargs):
        """
        Computes the energy vs. volume curve for isotropic expansion and compression. 

        Args:
            max_volume_scale:
                Maximum fractional change in volume to investigate
            num_steps:
                Number of steps to take in each direction
        """

        # The base class provides self.atoms, an ase.Atoms object representing the initial configuration of the crystal.
        # Use this configuration to evaluate the material property of interest
        original_cell = self.atoms.get_cell()
        num_atoms = len(self.atoms)

        # Besides temperature, stress, and atoms, you may wish to access other attributes of the base class for information about 
        # the material, such as its symmetry-reduced AFLOW prototype label. Here we use it to get information about the stoichiometry of the crystal.
        # See the documentation below this method definition and the API documentation for CrystalGenomeTestDriver for more information.
        num_atoms_in_formula = sum(get_stoich_reduced_list_from_prototype(self.prototype_label))

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
            self.atoms.set_cell(original_cell*linear_scale,scale_atoms=True)
            volume = self.atoms.get_volume()
            current_volume_per_atom = volume/num_atoms

            problem_occurred = False
            try:
                # The self.atoms object comes pre-initialized with the calculator set to the interatomic model
                # the Test Driver was called with. If you need to access the name of the KIM model (for example,
                # if you are exporting the atomic configuration to run an external simulator like LAMMPS), it can
                # be accessed at self.kim_model_name
                potential_energy = self.atoms.get_potential_energy()                
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
