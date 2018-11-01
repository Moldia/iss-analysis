# find oligos from the oligo database file

import openpyxl

wb = openpyxl.load_workbook(r'C:\Users\Xiaoyan\Downloads\Stockholm ordered oligos.xlsx')

allOligos = wb.get_sheet_by_name('Blad1')

nameOligos = allOligos.columns[1]
typeOligos = allOligos.columns[3]
hitOligos = []
for c, i in enumerate(nameOligos):
    try:
        if 'AP' in i.value and 'padlock' not in typeOligos[c].value.lower():
            hitOligos.append(c)
    except:     # name empty
        pass

with open(r'C:\Users\Xiaoyan\Downloads\oligosFound.csv', 'w') as f:
    for i in hitOligos:
        for j in allOligos.rows[i][:15]:
            try:
                f.write('%s,' % j.value)
            except:
                f.write("UnicodeEncodeError,")
        f.write('\n')

