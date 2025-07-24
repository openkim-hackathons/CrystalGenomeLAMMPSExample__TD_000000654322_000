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
from ase.build import bulk

kim_model_name = "Sim_LAMMPS_ReaxFF_RaymandVanDuinBaudin_2008_ZnOH__SM_449472104549_001"
test_driver = TestDriver(kim_model_name)
atoms = bulk("Zn")
test_driver(atoms)