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

This file is intended to demonstrate how to directly run Test Drivers that derive from :class:`~kim_tools.test_driver.core.SingleCrystalTestDriver`
for debugging. When developing your Test Driver, copy this file into its top-level directory and modify it (probably removing
all of these complicated comments as well).

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
# You can run your Driver by directly passing it an :class:`ase.Atoms` object. The base class will automatically
# perform a symmetry analysis on the structure and store a symmetry-reduced description of it. Note that the Atoms
# object you pass will not itself be passed to the ``_calculate()`` method, the crystal will be re-created from the
# symmetry-reduced description. Let's build a bulk zincblende structure and run our Driver on it, setting the
# ``_calculate()`` argument ``max_volume_scale`` to 0.1 and leaving the other argument as default. When testing a different
# Test Driver, this is where you would instead pass the specific arguments your ``_calculate()`` method uses instead.
# We are also demonstrating how to pass temperature and stress, even if our Test Driver doesn't use it.

from ase.build import bulk

atoms = bulk("ZnS", "zincblende", a=5.4093)

print("\nRUNNING TEST DRIVER ON ZINCBLENDE ATOMS OBJECT\n")
computed_property_instances = test_driver(
    atoms,
    max_volume_scale=0.1,
    num_steps=1,
    temperature_K=0,
    cell_cauchy_stress_eV_angstrom3=[0, 0, 0, 0, 0, 0],
)
###############################################################################
# The results of the calculation is returned in the format defined by the Property Definitions
# that the Driver uses and the `KIM Properties Framework <https://openkim.org/doc/schema/properties-framework/>`_.
# It can be accessed as a list of dictionaries, each corresponding to a Property Instance.
# Each caclulation may produce multiple Property Instances.
#
# For example, this is how you access the value of the key ``a`` (corresponding to the lattice constant),
# found in every Crystal Genome property, from the first Property Instance.

print("\n--------------------------------------")
print(
    "Lattice constant a: %f %s"
    % (
        computed_property_instances[0]["a"]["source-value"],
        computed_property_instances[0]["a"]["source-unit"],
    )
)
print("--------------------------------------\n")
###############################################################################
# You probably received a warning that the configuration you provided had a
# non-negligible force or stress. This is because the atoms object was not minimized
# with the potential we are using. One option to minimize the configuration
# is to use the ``EquilibriumCrystalStructure`` Test Driver, available in the
# `kimvv <https://github.com/openkim/kimvv>`_ package, which can simply be installed
# with
#
# .. code-block:: bash
#
#     pip install kimvv
#
# Once your Test Driver is published on openkim.org, it will be incorporated into
# ``kimvv`` and available to users just as easily!
#
# Instead an :class:`ase.Atoms` object, your Test Driver can use a KIM Property
# Instance containing a symmetry-reduced description of the crystal.
# ``EquilibriumCrystalStructure`` returns
# 3 different properties for each crystal (structure, energy, and density),
# but all of them contain the required fields to build a crystal, so we can
# pass any of them to our Test Driver to use as a relaxed structure.

from kimvv import EquilibriumCrystalStructure

ecs = EquilibriumCrystalStructure(kim_model_name)
print("\nMINIMIZING STRUCTURE\n")
relaxed_structure = ecs(atoms)[0]

print("\nRUNNING TEST DRIVER ON EquilibriumCrystalStructure OUTPUT\n")
computed_property_instances = test_driver(
    relaxed_structure,
    max_volume_scale=0.1,
    num_steps=1,
)

###############################################################################
# Testing Using a Prototype Label
# ===============================
#
# In the KIM Processing Pipeline, Test Drivers automatically run
# on of thousands of different crystal structures under the Crystal Genome
# testing framework. These are precomputed relaxations
# of each structure with each compatible interatomic potential in OpenKIM
# expressed as KIM Property Instances.
#
# You can replicate this functionality using
# the utility method :func:`kim_tools.test_driver.core.query_crystal_structures`
# to query for relaxed structures:

from kim_tools import query_crystal_structures

list_of_queried_structures = query_crystal_structures(
    kim_model_name=kim_model_name,
    stoichiometric_species=["Zn", "S"],
    prototype_label="AB_hP4_186_b_b",
)
# do something with computed_property_instances if you want

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
#   While the option to query for finite temperature and pressure exists, no such structures
#   are in the database yet, plus the query function needs to have tolerances implemented
#
# .. todo::
#
#   query_crystal_structures should be moved out of kim-tools and become a query endpoint on query.openkim.org
#
# A detailed description of AFLOW prototype labels can be found in
# Part IIB of https://arxiv.org/abs/2401.06875.
#
# A single set of the arguments given above (model, species, and AFLOW prototype label)
# may or may not correspond to multiple local minima (i.e. multiple sets of free parameters),
# so :func:`kim_tools.test_driver.core.query_crystal_structures` returns a list of
# dictionaries. You can then run your ``TestDriver`` by passing one of these
# dictionaries instead of an :class:`ase.Atoms` object.

for queried_structure in list_of_queried_structures:
    print("\nRUNNING TEST DRIVER ON QUERIED STRUCTURE\n")
    computed_property_instances = test_driver(
        queried_structure, max_volume_scale=0.1, num_steps=1
    )
    # do something with computed_property_instances if you want

###############################################################################
# You can also omit the ``kim_model_name`` argument, in which case
# instead of querying for potential-specific minima, reference data will be queried instead
# (typically DFT-relaxed). This is useful if you are using a model that is not on OpenKIM.org.
# In this case you should minimize the structure first. Just like any other Test Driver,
# :class:`kimvv.EquilibriumCrystalStructure` can take the dictionaries returned by
# :func:`kim_tools.test_driver.core.query_crystal_structures`.
# Here we are also demonstrating the ability to query by the crystal's name instead of
# or in addition to the AFLOW Prototype Label. This will be searched as a case-insensitive
# regex, so partial matches will work. Please note that, like any human-curated set of names,
# this is inevitably limited and incomplete. For example, "Face-Centered Cubic" is recognized,
# but not "FCC" or "Face Centered Cubic".
#
# Note that when querying for Reference Data, it is quite likely that the query will return duplicate structures
# for this, the :func:`kim_tools.test_driver.core.detect_unique_crystal_structures` utility is provided.
#
# .. todo::
#
#   Expand the database of short names

list_of_queried_structures = query_crystal_structures(
    stoichiometric_species=["Zn", "S"], short_name="Wurtzite"
)

from kim_tools import detect_unique_crystal_structures

unique_structure_indices = detect_unique_crystal_structures(list_of_queried_structures)

print(
    f"\n{len(unique_structure_indices)} of {len(list_of_queried_structures)} queried structures were found to be unique.\n"
)

for i in unique_structure_indices:
    # Minimize the structure with the EquilibriumCrystalStructure we previously instantiated
    print("\nMINIMIZING STRUCTURE\n")
    relaxed_structure = ecs(list_of_queried_structures[i])[0]
    print("\nRUNNING TEST DRIVER ON MINIMIZED QUERIED STRUCTURE\n")
    computed_property_instances = test_driver(
        relaxed_structure,
        max_volume_scale=0.1,
        num_steps=1,
    )
    # do something with computed_property_instances if you want

###############################################################################
# In addition to returning the computed property instances on each run,
# Test Drivers accumulate all the property instances computed over all calls,
# which can be accessed from :attr:`~kim_tools.test_driver.core.KIMTestDriver.property_instances` as
# a list of Python dictionaries, or written to a file (default ``output/results.edn``)
#

print(f"\nI've accumulated {len(test_driver.property_instances)} Property Instances\n")
test_driver.write_property_instances_to_file()
