# Questão 1 b)

import numpy as np

# e1 = n
# e2 = e

n = np.array([0.1666666666,0.6666666666,0.1666666666])

e = np.array([0.1666666666,0.1666666666,0.6666666666])

coord = np.array([
    [0,0],
    [6,-3],
    [6,3],
    [2.5981,-1.5],
    [6,0],
    [2.5981,1.5]
])


w = np.array([0.1666666666,0.1666666666,0.1666666666])

K = np.zeros((6,6))

F_s = np.zeros((6,1))

one = np.ones((6,1))

for i in range(3):

    N6 = np.array([
    [2*(n[i])**2 - n[i]],
    [2*(e[i])**2 - e[i]],
    [2*(n[i])**2 + 2*(e[i])**2 +4*n[i]*e[i] - 3*n[i] - 3*e[i] + 1],
    [4*n[i]*e[i]],
    [4*e[i] - 4*n[i]*e[i] - 4*(e[i])**2],
    [4*n[i] - 4*n[i]*e[i] - 4*(n[i])**2]
])

    GN6 = np.array([
        [4*n[i] - 1, 0, 4*n[i] + 4*e[i] - 3, 4*e[i], -4*e[i], 4 - 8*n[i] - 4*e[i]],
        [0, 4*e[i] - 1, 4*e[i] + 4*n[i] - 3, 4*n[i], 4 - 4*n[i] - 8*e[i], -4*n[i]]
    ])

    j = np.matmul(GN6,coord)

    detj = np.linalg.det(j)

    j_inv = np.linalg.inv(j)

    B = np.matmul(j_inv,GN6)

    K = K + (w[i]*detj*np.matmul(B.transpose(),B))

    F_s = F_s + (w[i]*detj*N6)

K = 5*K

F_s = 6*F_s

print("\nJ\n",j)
print("det J\n",detj)
print("\nJ inver\n",j_inv)
print("\nB\n",B)
print("\nK\n",K)
print("\nF_s\n",F_s)