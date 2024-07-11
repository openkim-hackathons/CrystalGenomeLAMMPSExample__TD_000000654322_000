#!/usr/bin/env python3

from test_driver.test_driver import TestDriver
from ast import literal_eval

test_level_optional_inputs = {}
model_name = input('Model name?\n')
print(model_name)
stoichiometric_species = literal_eval(input('Stoichiometric species (literal list of strings)?\n'))
print(stoichiometric_species)
prototype_label = input('Prototype label?\n')
print(prototype_label)
temperature_K = input('Temperature (K)?\n')
if temperature_K != '':
    print(temperature_K)
    test_level_optional_inputs['temperature_K'] = int(temperature_K)
else:
    print('No temperature given')
cell_cauchy_stress_eV_angstrom3 = input('Cauchy stress (literal list of floats, Voigt order xx,yy,zz,yz,xz,xy, eV/A^3)?\n')
if cell_cauchy_stress_eV_angstrom3 != '':
    print(cell_cauchy_stress_eV_angstrom3)
    test_level_optional_inputs['cell_cauchy_stress_eV_angstrom3'] = literal_eval(cell_cauchy_stress_eV_angstrom3)
else:
    print('No stress given')
query_result = literal_eval(input('Query result (literal list of dicts)?\n'))
print(query_result)
runtime_args = literal_eval(input('Runtime arguments (literal dictonary)?\n'))
print(runtime_args)
test = TestDriver(model_name)
# TODO: generalize querying to finite temperature and pressure. Right now these just let the
# driver know what to evaluate the property at, the query is always for 0K 0bar
for query_element in query_result:
    structure_level_optional_inputs = test_level_optional_inputs
    parameter_values_angstrom = [query_element['a.si-value']*1e10]
    if 'parameter-values.source-value' in query_element:
        parameter_values_angstrom += query_element['parameter-values.source-value']
    if 'parameter-names.source-value' in query_element:
        structure_level_optional_inputs['parameter_names'] = query_element['parameter-names.source-value']
    if 'library-prototype-label.source-value' in query_element:
        structure_level_optional_inputs['library_prototype_label'] = query_element['library-prototype-label.source-value']
    if 'short-name.source-value' in query_element:
        short_name = query_element['short-name.source-value']
        if not isinstance(short_name,list):
            short_name = [short_name]
        structure_level_optional_inputs['short_name'] = short_name
    if 'crystal-genome-source-structure-id.source-value' in query_element:
        crystal_genome_source_structure_id = query_element['crystal-genome-source-structure-id.source-value']
    else:
        crystal_genome_source_structure_id = [query_element['meta.uuid']+':'+str(query_element['instance-id'])]
    test(
        stoichiometric_species = stoichiometric_species,
        prototype_label = prototype_label,
        parameter_values_angstrom = parameter_values_angstrom,
        crystal_genome_source_structure_id = crystal_genome_source_structure_id,
        **structure_level_optional_inputs,
        **runtime_args
        )
test.write_property_instances_to_file()