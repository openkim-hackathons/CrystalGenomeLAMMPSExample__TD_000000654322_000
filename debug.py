#!/usr/bin/python

"""
=====================================================
Testing and Debugging Your Crystal Genome Test Driver
=====================================================

.. note::

    This Python file has comments written to be rendered using `Sphinx-Gallery <https://sphinx-gallery.github.io>`_ as part
    of the `documentation for the kim-test-utils package <https://kim-test-utils.readthedocs.io>`_. 
    A rendered version of this file should be available as part of that documentation.
    
When your Test Driver is integrated into the OpenKIM Processing Pipeline, the auxiliary files supplied with your
Test Driver will be used to pass inputs and parse outputs from the calculation you have implemented in ``test_driver.py``.
This functionality is also available for testing in the KIM Developer Platform. More information about this is avalable
here: :doc:`../../tutorial_pipeline`, however it is likely that you will first wish to test your Driver by invoking it directly.

Examples for how to do that are found in this file, ``debug.py`` from |example_url|. 
Everything you need to run this file is containerized in the :ref:`doc.KDP`,
or alternatively will be installed if you follow the :ref:`doc.standalone_installation`.

We will use a model from OpenKIM.org to test our Driver. The KIM Model needs to first be installed using 

.. code-block:: bash

    kim-api-collections-management install user SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_003

or (KDP only)

.. code-block:: bash

    kimitems install SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_003

more models can be found at https://openkim.org/browse/models/by-species, just replace the model name in the above
commands and it will be automatically downloaded and installed.

First, import your Test Driver and instantiate it by passing it a KIM Model Name:
"""
from test_driver import TestDriver
kim_model_name = "SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_003"
test_driver = TestDriver(kim_model_name)

###############################################################################
# Testing Using an :class:`~ase.Atoms` Object
# ===========================================
#
# You can test your Driver by directly passing it an :class:`ase.Atoms` object. The base :class:`~kim_test_utils.CrystalGenomeTestDriver`
# (actually, its base :class:`~kim_test_utils.KIMTestDriver`) provides the option to ``optimize`` (``False`` by default) the atomic positions
# and lattice vectors before running the simulation. The base class will automatically perform a symmetry analysis on the structure and populate
# the common Crystal Genome fields. Let's build a bulk wurtzite structure and run our Driver on it, setting the ``_calculate()`` argument ``max_volume_scale`` 
# to 0.1 and leaving the other argument as default. We are also demonstrating how to pass temperature and stress, even if our Test Driver doesn't use it.
from ase.build import bulk
atoms = bulk('ZnS','wurtzite',a=3.8)

print ('\nRUNNING TEST DRIVER ON WURTZITE ATOMS OBJECT\n')
test_driver(atoms,optimize=True,max_volume_scale=0.1,temperature_K=0,cell_cauchy_stress_eV_angstrom3=[0,0,0,0,0,0])
###############################################################################
# You can now access the results of the calculation in the format defined by your Property Definition 
# and the `KIM Properties Framework <https://openkim.org/doc/schema/properties-framework/>`_. 
# It can be accessed as a list of dictionaries, each corresponding to a Property Instance. 
# Each caclulation may produce multiple Property Instances, and any additional calls to 
# ``test_driver`` will append to this list. 
#
# This is how you access the value of the key ``a`` (corresponding to the lattice constant), 
# found in every Crystal Genome property, from the first Property Instance. You can access the 
# keys you defined in your Property Definition in the same way.
print('\n--------------------------------------')
print('Lattice constant a: %f %s'%(
    test_driver.property_instances[0]['a']['source-value'],
    test_driver.property_instances[0]['a']['source-unit']
))
print('--------------------------------------\n')
###############################################################################
# You can also dump the instances to a file, by default ``output/results.edn``.
test_driver.write_property_instances_to_file()

###############################################################################
# Testing Using a Prototype Label
# =============================== 
#
# In the KIM Processing Pipeline, your ``TestDriver`` will automatically run
# on of thousands of different crystal structures under the Crystal Genome 
# testing framework. These are precomputed relaxations 
# of each structure with each compatible interatomic potential in OpenKIM. 
# 
# Instead of passing your ``TestDriver`` an :class:`ase.Atoms` object,
# the pipeline will pass a set of keyword arguments containing the symmetry-reduced
# AFLOW prototype description of the crystal. The base class will then
# construct ``self.atoms``. You can replicate this functionality using
# the utility method :func:`kim_test_utils.query_crystal_genome_structures` 
# to query for relaxed structures:

from kim_test_utils import query_crystal_genome_structures
list_of_queried_structures = query_crystal_genome_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=['Zn','S'],
    prototype_label='AB_hP4_186_b_b')

###############################################################################
# ``AB_hP4_186_b_b`` is the AFLOW prototype label describing the symmetry of the 
# wurtzite structure. To fully specify the crystal structure, this label must
# be combined with one or more free parameters (unit cell and internal atom 
# degrees of freedom). The equilibrium values of these free parameters will depend
# on the potential being used, which is why the query function takes the KIM model
# name as an argument.
#
# .. todo::
#
#   Currently all structures queried for are zero temperature, zero stress.
#   We will add support for finite temperature and stress queries in the future.
# 
# A detailed description of AFLOW prototype labels can be found in
# Part IIB of https://arxiv.org/abs/2401.06875. A practical guide for testing your
# Driver on a variety of crystal structures is found lower on this page.
# 
# A single set of the arguments given above (model, species, and AFLOW prototype label)
# may or may not correspond to multiple local minima (i.e. multiple sets of free parameters),
# so :func:`kim_test_utils.query_crystal_genome_structures` returns a list of 
# dictionaries. You can then run your ``TestDriver`` by passing one of these
# dictionaries as keyword arguments instead of an :class:`ase.Atoms` object. 
# Do not use the ``optimize`` option.

for queried_structure in list_of_queried_structures:
    print ('\nRUNNING TEST DRIVER ON QUERIED STRUCTURE\n')
    test_driver(**queried_structure,max_volume_scale=0.1)

###############################################################################
# Remember that this will append the results to ``test_driver.property_instances``.
#
# Curated Set of Test Cases
# -------------------------
# 
# Because your Test Driver will run on a wide variety of crystals and interatomic potentials, 
# it is important to test it on a diverse selection of both. Below is a curated table
# of test cases that you can query for. These are arranged in increasing symmetry order
# (from triclinic to cubic).
# 
# You can explore more prototypes at 
# http://aflow.org/prototype-encyclopedia/, but it is not guaranteed that OpenKIM
# will have results or a compatible interatomic potential 
# (https://openkim.org/browse/models/by-species).
#
# Every time you use a new model, you will need to install the model and re-instantiate
# your ``TestDriver`` class.
#
# .. csv-table:: 
#    :header-rows: 1
#    :file: ../../structure_table.csv
#    :widths: 10, 10, 40, 40
#    :delim: tab