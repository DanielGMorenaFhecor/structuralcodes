import inspect
import os
import sys

import structuralcodes.codes.ec2_2023

print(structuralcodes.codes.ec2_2023)

# remove previous file
if os.path.exists('xlfunc.py'):
    os.remove('xlfunc.py')

# imports
imports = ['structuralcodes.codes.ec2_2023']

file_xlwings = open('xlfunc.py', 'x')
file_xlwings.write('import xlwings as xw' + '\n')
file_xlwings.write('import typing as t' + '\n')
for module in imports:
    file_xlwings.write('import ' + module + '\n')

# main xlwings
file_xlwings.write('def main():' + '\n')
file_xlwings.write('    wb = xw.Book.caller()' + '\n')
file_xlwings.write('    sheet = wb.sheets[0]' + '\n')
file_xlwings.write('    if sheet["A1"].value == "Hello xlwings!":' + '\n')
file_xlwings.write('        sheet["A1"].value = "Bye xlwings!"' + '\n')
file_xlwings.write('    else:' + '\n')
file_xlwings.write('        sheet["A1"].value = "Hello xlwings!"' + '\n')

# search functions
for module in imports:
    prefix = module.split('.')[-1] + '_'
    directory = os.getcwd() + '\\' + module.replace('.', '\\')
    for subdir, dirs, files in os.walk(directory):
        for file in files:
            file_dir = subdir + '\\' + os.path.splitext(file)[0]
            file_ext = os.path.splitext(file)[1]
            file = file_dir + file_ext
            if file_ext == '.py':
                with open(file, 'r') as f:
                    lines = f.readlines()
                    for i, row in enumerate(lines):
                        if row.find('def ') != -1:
                            def_name = row.split('(')[0]
                            def_name = def_name.split('def ')[1]

                            # def_param_types -> list of parameters with type
                            def_param_types = row.split('(', 1)[1]
                            j = 0
                            while def_param_types.find(')') == -1:
                                def_param_types = (
                                    def_param_types + lines[i + j + 1]
                                )
                                j = j + 1

                            def_param_types = def_param_types.split(')')[0]
                            def_param_types = def_param_types.split(',')
                            def_param_types = [
                                s.strip() for s in def_param_types
                            ]  # remove blanks and  \n
                            def_param_types = [
                                s.lstrip('_') for s in def_param_types
                            ]  # remove initial "_" in parameters
                            def_param_types = [
                                item.replace('As', 'As_')
                                for item in def_param_types
                            ]  # remove protected names in vba

                            # def_param_types -> list of params without type
                            def_param = []
                            for param in def_param_types:
                                def_param.append(param.split(':')[0])

                            # get docstrings
                            my_def = getattr(sys.modules[module], def_name)
                            docstring = inspect.getdoc(my_def)

                            # Write @xw function
                            file_xlwings.write('\n@xw.func' + '\n')
                            file_xlwings.write(
                                'def '
                                + prefix
                                + def_name
                                + '('
                                + ','.join(def_param_types)
                                + '):\n'
                            )
                            file_xlwings.write('    """' + docstring + '"""\n')
                            file_xlwings.write(
                                '    return '
                                + module
                                + '.'
                                + def_name
                                + '('
                                + ','.join(def_param)
                                + ')\n'
                            )


file_xlwings.write('\nif __name__ == "__main__":' + '\n')
file_xlwings.write(
    '    xw.Book("pruebaxlwings.xlsm").set_mock_caller()' + '\n'
)
file_xlwings.write('    main()' + '\n')

file_xlwings.close
