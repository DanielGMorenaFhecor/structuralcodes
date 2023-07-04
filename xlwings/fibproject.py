import xlwings as xw


def main():
    wb = xw.Book.caller()
    sheet = wb.sheets[0]
    if sheet["A1"].value == "it works":
        sheet["A1"].value = "it works!!!!!"
    else:
        sheet["A1"].value = "it works"


# if __name__ == "__main__":
#    xw.Book("pruebaxlwings.xlsm").set_mock_caller()
#    main()
