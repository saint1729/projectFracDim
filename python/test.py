from numpy import zeros

a = int(raw_input())
b = int(raw_input())

c = [[0 for col in range(a)] for row in range(b)]

c[0][0] = 1

print c[0][0]

"""from array import array

a = array('l')
b = array('c', 'hello world')
c = array('u', u'hello \u2641')
d = array('l', [1, 2, 3, 4, 5])
e = array('d', [1.0, 2.0, 3.14])

print d[0]"""
