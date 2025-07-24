#!/usr/bin/python

"""
Ultra-minimal example of running a LAMMPS Test Driver
"""
from test_driver.test_driver import TestDriver
from ase.build import bulk

kim_model_name = "Sim_LAMMPS_ReaxFF_RaymandVanDuinBaudin_2008_ZnOH__SM_449472104549_001"
test_driver = TestDriver(kim_model_name)
atoms = bulk("Zn")
test_driver(atoms, lmp_cmd="mpirun -np 2 lmp")

test_driver.write_property_instances_to_file()
