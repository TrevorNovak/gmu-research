import numpy as np

l1 = [1, 2, 3, 4, 5]
l2 = [-1,-2,-3,-4,-5]
l3 = [10, 11, 12, 13, 14, 15]

ll = np.matrix([l1, l2, l3])

outfile = "matrix.txt"
ll.tofile(outfile, sep=" ", format="%s")
