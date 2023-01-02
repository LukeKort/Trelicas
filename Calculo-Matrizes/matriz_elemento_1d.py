# Calculo de matrizes para problema 1D

import numpy as np

X = np.array([4-(2/np.sqrt(3)),4+(2/np.sqrt(3))])

k11 = 0
k12 = 0
k13 = 0
k22 = 0
k23 = 0
k33 = 0

for i in range(2):
    x = X[i]
    k11 = k11 + x*(x-5)**2

    k12 = k12 + x*(x-5)*(8-2*x)

    k13 = k13 + x*(x-5)*(x-3)

    k22 = k22 + x*(8-2*x)**2

    k23 = k23 + x*(8-2*x)*(x-3)

    k33 = k33 + x*(x-3)**2

K = 2*np.array([
    [k11, k12, k13],
    [k12, k22, k23],
    [k13, k23, k33]
])

fb1 = 0
fb2 = 0
fb3 = 0

for j in range(2):
    x = X[j]
    fb1 = fb1 + 0.125*(x-4)*(x-6)
    fb2 = fb2 - 0.25*(x-2)*(x-6)
    fb3 = fb3 + 0.125*(x-2)*(x-4)


# O 8 multiplicando veio do b e o 2 do método de integração numérica
F_body = 8*2*np.array([[fb1],[fb2],[fb3]])

p = 5 #ponto de aplicação do forçamento P

fp1 = 0.125*(p-4)*(p-6)
fp2 = -0.25*(p-2)*(p-6)
fp3 = 0.125*(p-2)*(p-4)

#O 24 multiplicando vei do P
F_ext = 24*np.array([[fp1],[fp2],[fp3]])

F = F_body + F_ext

print('\nMatriz de rigidez\n', K)
print('\nMatriz de força\n', F)