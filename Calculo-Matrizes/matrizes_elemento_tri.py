import numpy as np

j = np.array([
    [0,1.5],
    [-0.5,1.5]
])

deltah = np.array([
    [-1,1,0],
    [-1,0,1]
])

detj = j[0][0]*j[1][1] - j[0][1]*j[1][0]

j_inv = (1/detj)*np.array([[j[1][1],-j[0][1]],[-j[1][0],j[0][0]]])

B = np.matmul(j_inv,deltah)
K = 0.5*np.matmul(B.transpose(),B)*detj

print("det J\n",detj)
print("\nJ inver\n",j_inv)
print("\nB\n",B)
print("\nK\n",K)
