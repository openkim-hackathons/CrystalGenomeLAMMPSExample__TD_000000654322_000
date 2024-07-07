#!/usr/bin/python

"""
==============
debug.py
==============

.. note::

    The comments in this file are written to be rendered using `Sphinx-Gallery <https://sphinx-gallery.github.io>`_. Unless you wish to document your Test Driver in the same way, feel free to use a simpler format to comment your code!

"""

from kim_test_utils import query_crystal_genome_structures
from ase.build import bulk
from test_driver import TestDriver

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
