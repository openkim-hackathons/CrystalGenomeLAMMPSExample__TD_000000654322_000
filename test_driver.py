"""
CrystalGenomeASEExample
=======================

Example usage of the kim-test-utils package to make a Crystal Genome Test Driver using ASE. Use this file as a tutorial and template to write your own Crystal Genome Test Driver.

.. note::
    The comments in this file are written to be rendered using `Sphinx-Gallery <https://sphinx-gallery.github.io>`_. Unless you wish to document your Test Driver in the same way, feel free to use a simpler format to comment your code!

"""

#!/usr/bin/python

from kim_test_utils import CrystalGenomeTestDriver
from kim_python_utils.ase import get_isolated_energy_per_atom
from kim_test_utils import get_stoich_reduced_list_from_prototype

# %%
# You must create a class named ``TestDriver`` inheriting from  :class:`~kim_test_utils.CrystalGenomeTestDriver`.
# In your ``TestDriver`` class, you must overload the function :func:`~kim_test_utils.CrystalGenomeTestDriver._calculate`. 
# Besides `self`, the function must also accept ``**kwargs``. Before ``**kwargs``, you may add any additional arguments that 
# you wish users or the OpenKIM Pipeline to be able to vary. 
# .. note::
#   Temperature and stress are commonly used, so they are built-in attributes and do not need to be added as additional arguments: 
#     *  :attr:`~kim_test_utils.CrystalGenomeTestDriver.temperature_K`
#     *  :attr:`~kim_test_utils.CrystalGenomeTestDriver.cell_cauchy_stress_eV_angstrom3`

class TestDriver(CrystalGenomeTestDriver):
    def _calculate(self, nu, **kwargs):
        """
        Example calculate method. Doesn't actually do any calculations, just demonstrates example functionality.

        You may add arbitrary arguments, which will be passed to this method when the test driver is invoked.

        You must include **kwargs in the argument list, but you don't have to do anything with it

        Args:
            example_arg:
                An example argument. This would be an additional variable needed to define your simulation besides the crystal
        """

        ####################################################
        # ACTUAL CALCULATION BEGINS 
        ####################################################
        # calculate potential energy and do the required stuff to figure out per-formula and per-atom, and subtract isolated energy
        potential_energy = self.atoms.get_potential_energy()
        potential_energy_per_atom = potential_energy/len(self.atoms)
        reduced_stoichiometry = get_stoich_reduced_list_from_prototype(self.prototype_label) # i.e. "AB3\_...." -> [1,3]        
        binding_energy_per_formula = potential_energy_per_atom * sum(reduced_stoichiometry)
        for num_in_formula,species in zip(reduced_stoichiometry,self.stoichiometric_species):
            binding_energy_per_formula -= num_in_formula*get_isolated_energy_per_atom(self.kim_model_name,species)
        binding_energy_per_atom = binding_energy_per_formula/sum(reduced_stoichiometry)
        print("I was passed the following string argument as an example:\n\n%s\n\n"%example_arg)
        ####################################################
        # ACTUAL CALCULATION ENDS 
        ####################################################

        ####################################################
        # SOME USAGE EXAMPLES NOT NECESSARY FOR THE PRESENT TEST 
        ####################################################
        # If your self.atoms object has changed, this is how you update the Crystal Genome designation in your class instance:
        self._update_crystal_genome_designation_from_atoms()

        # If you just need to check that the symmetry hasn't changed without re-writing the Crystal Genome designation, do this:
        # Triclinic and monoclinic crystals may not be able to be matched perfectly, there is an option for loose matching (False by default)
        self._get_crystal_genome_designation_from_atoms_and_verify_unchanged_symmetry(loose_triclinic_and_monoclinic=True)

        # You can also check the symmetry of a passed atoms object
        atoms = self.atoms
        # Let's pretend we did something, like MD with LAMMPS
        # This is only meant to work with unit cells or small supercells, 
        # if you pass this function a large supercell, it will be intractably slow!
        self._get_crystal_genome_designation_from_atoms_and_verify_unchanged_symmetry(atoms)
        ####################################################
        # USAGE EXAMPLES END
        ####################################################

        ####################################################
        # PROPERTY WRITING
        ####################################################
        # add property instance and automatically pre-fill it with common Crystal Genome keys
        self._add_property_instance_and_common_crystal_genome_keys("binding-energy-crystal",write_stress=False, write_temp=False)

        # add the fields unique to this property
        self._add_key_to_current_property_instance("binding-potential-energy-per-atom",binding_energy_per_atom,"eV")
        self._add_key_to_current_property_instance("binding-potential-energy-per-formula",binding_energy_per_formula,"eV")


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
