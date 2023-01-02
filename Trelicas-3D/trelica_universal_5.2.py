# Code by Kort, nov. 8, 2022
# Available at github/lukekort

import numpy as np
import matplotlib.pyplot as plt

def matriz_ID(nos_presos,n_nos):
    # Essa função monda a matriz os dos restritos (1) e irrestritos (0)
    # Entradas: nós restritos, número de nós
    
    ID = np.zeros((n_nos,n_gl)) # pre-alocagem com zeros
    for i in range(n_nos):
        if n_gl == 2:
            if (i+1) in nos_presos:
                ID[i,0] = 1
                ID[i,1] = 1
        if n_gl == 3:
            if (i+1) in nos_presos:
                ID[i,0] = 1
                ID[i,1] = 1
                ID[i,2] = 1

    # Controle de qualidade
    # print('ID\n',ID)
    return(ID)

def matriz_propriedades(E,A,n_el):
    # Essa função monda a matriz de propriedades dos elementos
    # Entradas: Modo de Young, Área da seção transversão, número de elementos

    prop = np.zeros((2,n_el)) # pre-alocagem com zeros
    for i in range(n_el):
        prop[0][i] = E
        prop[1][i] = A

    # Controle de qualidade
    # print('Prop\n',prop)
    return(prop)

def matriz_IDM(ID):
    IDM = np.zeros((size_ID,n_gl),int) #pre-alocagem com zeros
    neq = 0
    for i in range(size_ID):
        for j in range(n_gl):
            aux = ID[i][j]
            if aux == 0:
                neq = neq + 1
                IDM[i][j] = neq #Nó inrestrito, recebe um número de equação
    
    # Controle de qualidade
    # print('IDM\n',IDM)
    return(IDM,neq)

def matriz_elementar(E,A,inci1,inci2,coordxi,coordxj,coordyi,coordyj,coordzi,coordzj):
    # Essa função monta a matriz elementar
    # Entradas:Módulo de Young,
    #     Área da seção transveral
    #     primeiro nó do elemento
    #     segundo nó do elemento
    #     coordenada x do nó inicial
    #     coordenada x do nó final
    #     coordenada y do nó inicial
    #     coordenada y do nó final
    #     coordenada z do nó inicial
    #     coordenada z do nó final

    if n_gl == 2:
        l = np.sqrt((coordxj-coordxi)**2 + (coordyj-coordyi)**2)
        cos_theta = (coordxj-coordxi)/l
        sin_theta = (coordyj-coordyi)/l
        R = np.array(
            [[cos_theta**2, cos_theta*sin_theta, -cos_theta**2, -cos_theta*sin_theta],
            [cos_theta*sin_theta, sin_theta**2, -cos_theta*sin_theta, -sin_theta**2],
            [-cos_theta**2, -cos_theta*sin_theta, cos_theta**2, cos_theta*sin_theta],
            [-cos_theta*sin_theta, -sin_theta**2, cos_theta*sin_theta, sin_theta**2]]
            )
    else:
        l = np.sqrt((coordxj-coordxi)**2 + (coordyj-coordyi)**2 + (coordzj-coordzi)**2)
        c_x = (coordxj-coordxi)/l
        c_y = (coordyj-coordyi)/l
        c_z = (coordzj-coordzi)/l
        R = np.array(
            [[c_x**2, c_x*c_y, c_x*c_z, - c_x**2, -c_x*c_y, -c_x*c_z],
            [c_x*c_y, c_y**2, c_y*c_z, -c_x*c_y, -c_y**2, -c_y*c_z],
            [c_x*c_z, c_y*c_z, c_z**2, -c_x*c_z, -c_y*c_z, -c_z**2],
            [-c_x**2, -c_x*c_y, -c_y*c_z, c_x**2, c_x*c_y, c_x*c_z],
            [-c_x*c_y, -c_y**2, -c_y*c_z, c_x*c_y, c_y**2, c_y*c_z],
            [-c_x*c_z, -c_y*c_z, -c_z**2, c_x*c_z, c_y*c_z, c_z**2]])
 
    K = ((E*A)/l)*R

    # Controle de qualidade
    # print(inci1,inci2,coordxi,coordxj,coordyi,coordyj)
    # print('l',l)
    # print('cos:',cos_theta)
    # print('sin:',sin_theta)
    # print('K\n',K)
    return(K)

# Entrada de dados-----------------------------------------------------------------
# Coordenadas dos nós
coord = np.transpose(np.array(
    [[0,4,0],
    [0,2,5],
    [-3,0,0],
    [3,0,0],
    [0,0,5]]))

# Mapa dos nós em cada elemento
inci = np.transpose(np.array(
    [[1,2],
    [1,3],
    [1,4],
    [2,4],
    [2,3],
    [2,5]]))     

# Módulo de young
E = 250e6
# Área da seção transversal
A = 1e-2

# Nós presos (restritos)
no_preso = [3,4,5]

# Aplicação de força

force = np.transpose(np.array(
    [[0,0,-6e3],
    [0,0,0],
    [0,0,0],
    [0,0,0],
    [0,0,0]]))

# Controle do gráfico

# Plotar gráfico [y/n] - (defalt='y')

plt_q = 'y'

# Fator de aumento do dano (para efeitos de visualização)

fator = 5 

# Não tem interação do usuário a partir desse ponto--------------------------------

# Montagem das matrizes
n_gl = np.size(coord,0) #Número de graus de liberdade
size_ID = np.size(coord,1)    # número de nós(colunas de coord)
size_inci = np.size(inci,1)  # número de elementos (colunas de inci)
ID = matriz_ID(no_preso,size_ID) #nós com restrição          
prop = matriz_propriedades(E,A,size_inci) #matriz de propriedades
IDM,neq = matriz_IDM(ID)

# Controle de qualidade
# elemento = 6
# e = elemento - 1
# K = matriz_elementar(
#     prop[0][e],
#     prop[1][e],
#     inci[0][e],
#     inci[1][e],
#     coord[0][inci[0][e]-1],
#     coord[0][inci[1][e]-1],
#     coord[1][inci[0][e]-1],
#     coord[1][inci[1][e]-1])
# print('Ke\n',K)
# print('gl',n_gl)
# print('size_ID',size_ID)
# print('size_inci', size_inci)
# print('ID\n', ID)
# print('IMD\n', IDM)
# print('neq', neq)
# print('coord\n', coord)
# print('inci\n', inci)
# print('força', force)


# Vetores de força

F = np.zeros((neq,1))    #pre-alocagem (-1 pq começa na posição zero)
for i in range(size_ID):
    for j in range(n_gl):
        eq = IDM[i][j]
        if eq != 0:
            F[eq-1] = force[j][i]   #-1 pq começa na posição zero

# Controle de qualidade
# print('F',F)

# Montagem do sistema (forma eficiente)-------------------------------------------

K = np.zeros((neq,neq))     # pre-alocagem com zeros
LM = np.zeros((n_gl,2),int)      # pre-alocagem com zeros

# Controle de qualidade
# print('neq\n',neq)

for e in range(size_inci):
    for j in range(2):      # Número de nós por elemento
        no = inci[j][e]
        
        # Controle de qualidade
        # print(no)
        
        for i in range(n_gl):  # Range até o número de gl por nó
            eq = IDM[no-1][i]
            LM[i][j] = eq

            # Controle de qualidade
            # print(eq)
            #print(LM)

    for j in range(1,3):
        ii = 2*(j-1)
        for i in range(1,n_gl+1):
            l = LM[i-1][j-1]

            # Controle de qualidade
            # print(l)

            ll = ii + i
            for n in range(1,3):
                if l != 0:
                    jj = 2*(n-1)
                    for k in range(1,n_gl+1):
                        m = LM[k-1][n-1]
                        ml = jj + k
                        if n_gl == 2:
                            coordzi = 1
                            coordzj = 1
                        else:
                            coordzi = coord[2][inci[0][e]-1]
                            coordzj = coord[2][inci[1][e]-1]
                        Ke = matriz_elementar(
                        prop[0][e],
                        prop[1][e],
                        inci[0][e],
                        inci[1][e],
                        coord[0][inci[0][e]-1],
                        coord[0][inci[1][e]-1],
                        coord[1][inci[0][e]-1],
                        coord[1][inci[1][e]-1],
                        coordzi,
                        coordzj)

                        # Controle de qualidade
                        # print(l,m,ll,ml)

                        if m != 0:
                            K[l-1][m-1] = K[l-1][m-1] + Ke[ll-1][ml-1]

# Controle de qualidade
# print('K\n',K)

# Solução do problema ------------------------------------------------------------

U = np.matmul(np.linalg.inv(K),F)

# Controle de qualidade
# print('U\n',U)

V = np.zeros((size_ID,n_gl))
for i in range(size_ID):
    for j in range(n_gl): #graus de liberdade
        eq = IDM[i][j]
        if eq != 0:
            V[i][j] = U[eq-1]



if n_gl == 2:
    print('\n# Deslocamentos por nó #\n [x] [y]\n',V)
else:
    print('\n# Deslocamentos por nó #\n [x] [y] [z]\n',V)

# Calculo do sistema deformado----------------------------------------------------

coord_f = coord + fator*V.transpose()

# Controle de qualidade
# print(coord_f)

# Carregamento interno------------------------------------------------------------

Fe = np.zeros((size_inci,2))
for i in range(size_inci):
    noI = inci[0][i] - 1
    noJ = inci[1][i] - 1
    xI = coord[0][noI]
    xJ = coord[0][noJ]
    yI = coord[1][noI]
    yJ = coord[1][noJ]
    if n_gl == 2:
        l = np.sqrt((xJ-xI)**2 + (yJ-yI)**2)
        cos_theta = (xJ-xI)/l
        sin_theta = (yJ-yI)/l
        sigma_e = -cos_theta*V[noI][0] - sin_theta*V[noI][1] + cos_theta*V[noJ][0] + sin_theta*V[noJ][1]
    else:
        zI = coord[2][noI]
        zJ = coord[2][noJ]
        l = np.sqrt((xJ-xI)**2 + (yJ-yI)**2 + (zJ-zI)**2)
        c_x = (xJ-xI)/l
        c_y = (yJ-yI)/l
        c_z = (zJ-zI)/l
        sigma_e = -c_x*V[noI][0] -c_y*V[noI][1] -c_z*V[noI][2] + c_x*V[noJ][0] + c_y*V[noJ][1] + c_z*V[noJ][2]

    Fe[i][0] = i + 1
    Fe[i][1] = ((E/l)*sigma_e)/1000

# Controle de qualidade
# print(Fe)

if plt_q.lower() == 'y':
    plt.title('Tensão em cada elemento')
    plt.xlabel('N° do elemento')
    plt.ylabel('Tensão [1e3]')
    plt.bar_label(plt.bar(Fe[:,0],Fe[:,1],1))
    plt.show()