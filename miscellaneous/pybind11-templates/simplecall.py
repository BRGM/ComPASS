import simplecall as sc

for i, j in [(2, 2), (3, 3)]:
    print('Sum ', i, ' + ', j)
    print('C++ call from python')
    res = sc.add(i, j)
    print('result: ', res)
    print('Fortran call from python')
    res = sc.fadd(i, j)
    print('result: ', res)
    print('Fortran call from python - by reference')
    res = sc.fadd_byref(i, j)
    print('result: ', res)
