#!/usr/bin/python

"""
Invoking a Crystal Genome Test Driver Directly
==============================================

.. note::

    This Python file has comments written to be rendered using `Sphinx-Gallery <https://sphinx-gallery.github.io>`_ as part
    of the `documentation for the kim-tools package <https://kim-tools.readthedocs.io>`_. 
    A rendered version of this file should be available as part of that documentation. The rendered documentation
    is contained entirely in the comments, so this file can be run just like any other Python file from your shell, e.g. 

    .. code-block:: bash

        python CrystalGenomeASEExample__TD_000000654321_000/run.py
    
We will use a model from OpenKIM.org to run our Driver. For this, `kimpy <https://github.com/openkim/kimpy>`_
must be installed (see Note below regarding using non-KIM ASE calculators). The KIM Model needs to first be installed 
using the following command (in your shell, not in Python):

.. code-block:: bash

    kim-api-collections-management install user SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_003

or (`KIM Developer Platform <https://openkim.org/doc/evaluation/kim-developer-platform/>`_ only)

.. code-block:: bash

    kimitems install SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_003

more models can be found at https://openkim.org/browse/models/by-species, just replace the model name in the above
commands and it will be automatically downloaded and installed.

First, import your Test Driver and instantiate it by passing it a KIM Model Name:
"""
from test_driver.test_driver import TestDriver
kim_model_name = "SW_ZhouWardMartin_2013_CdTeZnSeHgS__MO_503261197030_003"
test_driver = TestDriver(kim_model_name)

###############################################################################
# .. note::
#
#   If the Test Driver uses only ASE for computations, you can use any ASE calculator,
#   not just KIM models. Simply instantiate the ``TestDriver`` class by passing it
#   the :class:`~ase.calculators.calculator.Calculator` instance you wish to use instead
#   of the string identifying the KIM model.
#
# Running Using an :class:`~ase.Atoms` Object
# -------------------------------------------
#
# You can run your Driver by directly passing it an :class:`ase.Atoms` object. The base :class:`~kim_tools.CrystalGenomeTestDriver`
# (actually, its base :class:`~kim_tools.KIMTestDriver`) provides the option to ``optimize`` (``False`` by default) the atomic positions
# and lattice vectors before running the simulation. The base class will automatically perform a symmetry analysis on the structure and populate
# the common Crystal Genome fields. Let's build a bulk wurtzite structure and run our Driver on it, setting the ``_calculate()`` argument ``max_volume_scale`` 
# to 0.1 and leaving the other argument as default. We are also demonstrating how to pass temperature and stress, even if our Test Driver doesn't use it.
from ase.build import bulk
atoms = bulk('ZnS','wurtzite',a=3.8)

print ('\nRUNNING TEST DRIVER ON WURTZITE ATOMS OBJECT\n')
test_driver(atoms,optimize=True,max_volume_scale=0.1,temperature_K=0,cell_cauchy_stress_eV_angstrom3=[0,0,0,0,0,0])
###############################################################################
# You can now access the results of the calculation in the format defined by the Property Definitions
# that the Driver uses and the `KIM Properties Framework <https://openkim.org/doc/schema/properties-framework/>`_. 
# It can be accessed as a list of dictionaries, each corresponding to a Property Instance. 
# Each caclulation may produce multiple Property Instances, and any additional calls to 
# ``test_driver`` will append to this list. 
#
# This is how you access the value of the key ``a`` (corresponding to the lattice constant), 
# found in every Crystal Genome property, from the first Property Instance. You can access the 
# other keys in the same way.
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
# In the KIM Processing Pipeline, Test Drivers automatically run
# on of thousands of different crystal structures under the Crystal Genome 
# testing framework. These are precomputed relaxations 
# of each structure with each compatible interatomic potential in OpenKIM. 
# 
# Instead of passing your ``TestDriver`` an :class:`ase.Atoms` object,
# the pipeline will pass a set of keyword arguments containing the symmetry-reduced
# AFLOW prototype description of the crystal. The base class will then
# construct ``self.atoms``. You can replicate this functionality using
# the utility method :func:`kim_tools.query_crystal_genome_structures` 
# to query for relaxed structures:

from kim_tools import query_crystal_genome_structures
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
# Part IIB of https://arxiv.org/abs/2401.06875. 
# 
# A single set of the arguments given above (model, species, and AFLOW prototype label)
# may or may not correspond to multiple local minima (i.e. multiple sets of free parameters),
# so :func:`kim_tools.query_crystal_genome_structures` returns a list of 
# dictionaries. You can then run your ``TestDriver`` by passing one of these
# dictionaries as keyword arguments instead of an :class:`ase.Atoms` object. 
# Do not use the ``optimize`` option.

for queried_structure in list_of_queried_structures:
    print ('\nRUNNING TEST DRIVER ON QUERIED STRUCTURE\n')
    test_driver(**queried_structure,max_volume_scale=0.1)

###############################################################################
# Remember that this will append the results to ``test_driver.property_instances``.