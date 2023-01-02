import numpy as np

n = np.array([-0.57735,0.57735])

e = np.array([-0.57735,0.57735])

coord = np.array([
    [0,1],
    [0,0],
    [2,0.5],
    [2,1]
])

w = np.array([1,1])

K = np.zeros((4,4))

for i in range(2):
    for k in range(2):
        GN4 = (1/4)*np.array([
            [(n[i]-1),(1-n[i]),(1+n[i]),(-n[i]-1)],
            [(e[i]-1),(-e[i]-1),(e[i]+1),(1-e[i])]
        ])

        j = np.matmul(GN4,coord)

        detj = np.linalg.det(j)

        j_inv = np.linalg.inv(j)

        B = np.matmul(j_inv,GN4)

        K = K + (w[i]*w[k]*detj*np.matmul(B.transpose(),B)) 

K = 5*K

print("\nJ\n",j)
print("det J\n",detj)
print("\nJ inver\n",j_inv)
print("\nB\n",B)
print("\nK\n",K)
