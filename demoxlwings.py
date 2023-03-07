import xlwings as xw

import structuralcodes.codes.mc2010 as mc2010


def main():
    wb = xw.Book.caller()
    sheet = wb.sheets[0]
    if sheet["A1"].value == "Hello xlwings!":
        sheet["A1"].value = "Bye xlwings!"
    else:
        sheet["A1"].value = "Hello xlwings!"


@xw.func
def fcm(fck: float, delta_f: float = 8.0):
    return mc2010.fcm(fck, delta_f)


if __name__ == "__main__":
    xw.Book("pruebaxlwings.xlsm").set_mock_caller()
    main()


# import structuralcodes.material.concrete as conc
