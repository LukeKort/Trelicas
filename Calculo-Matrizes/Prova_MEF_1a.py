#Questão 1 a)

#Calculo de matrizes para método dos elementos finitos
#Matriz B, K, F_forte, N
import numpy as np

#Dados do problema
k = 5
s = 6

#Dados do elemento 1
coord = np.array([
    [0,0],
    [6,-3],
    [6,3]
    ])

x1 = coord[0][0]
x2 = coord[1][0]
x3 = coord[2][0]
y1 = coord[0][1]
y2 = coord[1][1]
y3 = coord[2][1]

A = ((x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1))

B = (1/A)*(np.array([[(y2-y3),(y3-y1),(y1-y2)],[(x3-x2),(x1-x3),(x2-x1)]]))

K = k*(A/2)*(np.matmul(B.transpose(),B))

f_source = (s*(A/2)/3)*np.array([[1],[1],[1]])

# n1 = (1/A)*(x2*y3 - x3*y2 + (y2-y3)*x + (x3-x2)*y)
# n2 = (1/A)*(x3*y1 - x1*y3 + (y3-y1)*x + (x1-x3)*y)
# n3 = (1/A)*(x1*y2 - x2*y1 + (y1-y2)*x + (x2-x1)*y)

n1 = (1/A)*np.array([(x2*y3 - x3*y2),(y2-y3),(x3-x2)])
n2 = (1/A)*np.array([(x3*y1 - x1*y3),(y3-y1),(x1-x3)])
n3 = (1/A)*np.array([(x1*y2 - x2*y1),(y1-y2),(x2-x1)])

N_1 = np.array([[n1],[n2],[n3]])

print("\n###Elemento 1 ###\n")
print("\nÁrea do elemento\n",A/2)
print("\nMatriz B\n",B)
print("\nMatriz K\n",K)
print("\nMatriz termo fonte\n",f_source)

#----------------------------------------------------------------------------
#Segunda parte - tem que compeltar com dados obtidos na primeira parte

f_fluxo = np.array([
    [0],
    [0],
    [0],
])


print("\nF fluxo Global\n", f_fluxo)

#----------------------------------------------------------------------------
#Terceira parte - tem que compeltar com dados obtidos na segunda parte

D = np.array([
    [14.4],
    [0],
    [0]
])

q = -k*np.matmul(B,D)

print("\nMatriz q\n", q)