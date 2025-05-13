import importlib
import inspect
from pathlib import Path

FORBIDDEN_NAMES = [
    'As'
]  # List of forbidden names to avoid conflicts with vba keywords
MODULES = [
    'structuralcodes.codes.ec2_2023._annexB_time_dependent',
    'structuralcodes.codes.ec2_2023._section5_materials',
    'structuralcodes.codes.ec2_2023._section9_sls',
    'structuralcodes.codes.ec2_2004._concrete_creep_and_shrinkage',
    'structuralcodes.codes.ec2_2004._concrete_material_properties',
    'structuralcodes.codes.ec2_2004._reinforcement_material_properties',
    'structuralcodes.codes.ec2_2004._section_7_3_crack_control',
    'structuralcodes.codes.ec2_2004.shear',
]
OUTPUT = Path('xlfunc.py')


def generate_xlwings_udfs(modules, output_path):
    """Generates an xlwings UDF wrapper file for the given modules."""
    with open(output_path, 'w', encoding='utf-8') as f:
        # Write header imports
        f.write('import xlwings as xw\n')
        f.write('import typing as t\n')
        f.write('from typing import Tuple,Dict\n')
        for mod in modules:
            f.write(f'import {mod}\n')
        f.write('\n')

        # Write the example main function
        f.write('def main():\n')
        f.write('    """Example function to demonstrate xlwings UDFs."""\n')
        f.write('    wb = xw.Book.caller()\n')
        f.write('    sheet = wb.sheets[0]\n')
        f.write('    # Toggle the cell value between two states\n')
        f.write('    if sheet["A1"].value == "It works!":\n')
        f.write('        sheet["A1"].value = "Still works!"\n')
        f.write('    else:\n')
        f.write('        sheet["A1"].value = "It works!"\n')
        f.write('\n')

        # Dynamically generate UDF wrappers
        for mod_name in modules:
            module = importlib.import_module(mod_name)
            prefix = ''.join(
                mod_name.rsplit('.', 2)[-2:]
            )  # Get the last two parts of the module name

            for func_name, func_obj in inspect.getmembers(
                module, inspect.isfunction
            ):
                sig = inspect.signature(func_obj)
                params_with_types = []
                for name, param in sig.parameters.items():
                    if name in FORBIDDEN_NAMES:
                        name += '_'
                    if param.annotation is inspect._empty:
                        params_with_types.append(name)
                    else:
                        ann_str = inspect.formatannotation(param.annotation)
                        if ann_str.startswith('Literal'):
                            ann_str = f't.{ann_str}'
                        # params_with_types.append(f'{name}: {ann_str}') TODO
                        params_with_types.append(f'{name}')

                # return
                if sig.return_annotation is inspect._empty:
                    ret_ann = 't.Any'
                else:
                    ret_ann = inspect.formatannotation(sig.return_annotation)

                # docstring
                doc = inspect.getdoc(func_obj) or ''

                # write wrapper method
                f.write('@xw.func\n')
                f.write(
                    f"def SC_{prefix}_{func_name}"
                    f"({', '.join(params_with_types)}):\n"
                    # f" -> {ret_ann}:\n" TODO
                )
                f.write(f'    """{doc}"""\n')
                f.write(
                    f"    return {mod_name}.{func_name}"
                    f"({', '.join(p.split(':')[0] for p in params_with_types)})\n\n"
                )

        #  __main__
        f.write('if __name__ == "__main__":\n')
        f.write('    xw.Book("pruebaxlwings.xlsm").set_mock_caller()\n')
        f.write('    main()\n')


if __name__ == '__main__':
    generate_xlwings_udfs(MODULES, OUTPUT)
