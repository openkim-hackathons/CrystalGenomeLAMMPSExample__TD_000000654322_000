#!/usr/bin/python

"""
test_driver.py
==============

Example usage of the :mod:`kim_test_utils` package to make a Crystal Genome Test Driver using ASE. Use this file as a tutorial and template to write your own Crystal Genome Test Driver.

.. note::

    The comments in this file are written to be rendered using `Sphinx-Gallery <https://sphinx-gallery.github.io>`_. Unless you wish to document your Test Driver in the same way, feel free to use a simpler format to comment your code!

You must create a class named ``TestDriver`` inheriting from  :class:`~kim_test_utils.CrystalGenomeTestDriver`.
In your ``TestDriver`` class, you must overload the function :func:`~kim_test_utils.KIMTestDriver._calculate`. 
Besides ``self``, the function must also accept ``**kwargs``. Before ``**kwargs``, you may add any additional arguments that 
you wish users or the OpenKIM Pipeline to be able to vary. 

.. note::

  Temperature and stress are commonly used, so they are built-in attributes and do not need to be added as additional arguments: 

    *  :attr:`~kim_test_utils.CrystalGenomeTestDriver.temperature_K`
    *  :attr:`~kim_test_utils.CrystalGenomeTestDriver.cell_cauchy_stress_eV_angstrom3`

The ``_calculate`` function implemented below computes the energy vs. volume curve for isotropic expansion and compression of a crystal
at zero temperature. You can use it as a starting point for your Test Driver. See the comments for explanations. Documentation regarding 
more complex usage of :mod:`kim_test_utils` can be found below the function definition.

``_calculate`` method
=====================

"""

from kim_test_utils import CrystalGenomeTestDriver
from kim_test_utils import get_stoich_reduced_list_from_prototype

class TestDriver(CrystalGenomeTestDriver):
    def _calculate(self, max_volume_scale: float = 1e-2, num_steps: int = 10, **kwargs):
        """
        Computes the energy vs. volume curve for isotropic expansion and compression. Example Test Driver demonstrating
        usage of the kim-test-utils package. 

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

# %%
# See the code and comments for usage examples of the various functions of the :class:`~kim_test_utils.CrystalGenomeTestDriver`.
# Click on the link to see full API documentation of the class and its members.


from kim_test_utils import query_crystal_genome_structures
from ase.build import bulk

if __name__ == "__main__":        
    ####################################################
    # if called directly, do some debugging examples
    ####################################################
    kim_model_name = "MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002"

    # For initialization, only pass a KIM model name or an ASE calculator
    test_driver = TestDriver(kim_model_name)

    # To do a calculation, you can pass an ASE.Atoms object or a Crystal Genome prototype designation.
    # Atoms object example:
    atoms = bulk('Fe','bcc',a=2.863,cubic=True)
    test_driver(atoms,example_arg="my example argument")

    # You can get a list of dictionaries of the results like this:
    print(test_driver.get_property_instances())

    # Or write it to a file (by default `output/results.edn`) like this:
    test_driver.write_property_instances_to_file()

    # Alternatively, you can pass a Crystal Genome designation. You can automatically query for all equilibrium structures for a given 
    # species and prototype label like this:
    cg_des_list = query_crystal_genome_structures("MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002",['Fe','P'],'AB_oP8_62_c_c')

    # IMPORTANT: cg_des is a LIST. Pass only one element of it to the test, as keywords (i.e. using **):
    for cg_des in cg_des_list:
        test_driver(**cg_des,example_arg="my example argument")

    # Now both results are in the property instances:
    print(test_driver.get_property_instances())

    # Here are some other crystal prototypes supported by the current model you can try:
    # ["Fe", "P"], "A2B_hP9_189_fg_ad"
    # ["Fe", "P"], "A3B_tI32_82_3g_g"
    # ["Fe", "P"], "AB_oP8_62_c_c"
    # ["Fe", "P"], "AB2_oP6_58_a_g"
    # ["Fe", "P"], "AB4_mC40_15_ae_4f"
    # ["Fe", "P"], "AB4_mP30_14_ae_6e"
    # ["Fe", "P"], "AB4_oC20_20_a_2c"
    # ["Fe"], "A_cF4_225_a"
    # ["Fe"], "A_cI2_229_a"
    # ["Fe"], "A_hP2_194_c"
    # ["Fe"], "A_tP28_136_f2ij"
    # ["P"], "A_aP24_2_12i"
    # ["P"], "A_cP1_221_a"
    # ["P"], "A_mC16_12_2ij"
    # ["P"], "A_oC8_64_f"
    # ["P"], "A_tI4_139_e"

# %%
