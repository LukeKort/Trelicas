#Calculo de matrizes para método dos elementos finitos
#Matriz B, K, F_forte, N
import numpy as np

#Dados do problema
k = 5
s = 6

#Dados do elemento 1
coord = np.array([
    [0,0],
    [2,0.5],
    [0,1]
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

f_s_1 = f_source
B1 = B
K1 = K

print("\n###Elemento 1 ###\n")
print("Área do elemento\n",A/2)
print("Matriz B\n",B)
print("Matriz K\n",K)

#Dados do elemento 2

coord = np.array([
    [2,0.5],
    [2,1],
    [0,1]
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

n1 = (1/A)*np.array([(x2*y3 - x3*y2),(y2-y3),(x3-x2)])
n2 = (1/A)*np.array([(x3*y1 - x1*y3),(y3-y1),(x1-x3)])
n3 = (1/A)*np.array([(x1*y2 - x2*y1),(y1-y2),(x2-x1)])

N_2 = np.array([[n1],[n2],[n3]])

f_s_2 = f_source
B2 = B
K2 = K

print("\n###Elemento 2 ###\n")
print("Área do elemento\n",A/2)
print("Matriz B\n",B)
print("Matriz K\n",K)

L1 = np.array([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,1,0]
])

L2 = np.array([
    [0,1,0,0],
    [0,0,0,1],
    [0,0,1,0]
])

KG = (
    np.matmul(np.matmul(L1.transpose(),K1),L1) +
    np.matmul(np.matmul(L2.transpose(),K2),L2)
)

F_FonteG = (
    np.matmul(L1.transpose(),f_s_1) +
    np.matmul(L2.transpose(),f_s_2)
)

print("\nMatriz Global\n", KG)
print("\nF fonte Global\n", F_FonteG)
print("\nMatriz N1\n"," [0][X][Y]\n", N_1)
print("\nMatriz N2\n"," [0][X][Y]\n", N_2)

#----------------------------------------------------------------------------
#Segunda parte - tem que compeltar com dados obtidos na primeira parte

f_fluxo_1 = np.array([
    [0],
    [0],
    [0],
])

f_fluxo_2 = np.array([
    [0],
    [-20],
    [-20],
])

F_FluxoG = (
    np.matmul(L1.transpose(),f_fluxo_1) +
    np.matmul(L2.transpose(),f_fluxo_2)
)

print("\nF fluxo Global\n", F_FluxoG)

#----------------------------------------------------------------------------
#Terceira parte - tem que compeltar com dados obtidos na segunda parte

D = np.array([
    [0],
    [0],
    [0],
    [-1.788]
])

D1 = np.matmul(L1,D)
D2 = np.matmul(L2,D)

q1 = -k*np.matmul(B1,D1)
q2 = -k*np.matmul(B2,D2)

print("\nMatriz D\n", D)
print("\nMatriz D1\n", D1)
print("\nMatriz D2\n", D2)
print("\nMatriz q1\n", q1)
print("\nMatriz q2\n", q2)