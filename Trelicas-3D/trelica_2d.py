import numpy as np
import matplotlib.pyplot as plt

def matriz_ID(nos_presos,n_nos):
    # Essa função monda a matriz os dos restritos (1) e irrestritos (0)
    # Entradas: nós restritos, número de nós
    
    ID = np.zeros((n_nos,2)) # pre-alocagem com zeros
    for i in range(n_nos):
        if (i+1) in nos_presos:
            ID[i,0] = 1
            ID[i,1] = 1

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
    IDM = np.zeros((size_ID,2),int) #pre-alocagem com zeros
    neq = 0
    for i in range(size_ID):
        for j in range(2):
            aux = ID[i][j]
            if aux == 0:
                neq = neq + 1
                IDM[i][j] = neq #Nó inrestrito, recebe um número de equação
    
    # Controle de qualidade
    # print('IDM\n',IDM)
    return(IDM,neq)

def matriz_elementar(E,A,inci1,inci2,coordxi,coordxj,coordyi,coordyj):
    # Essa função monda a matriz elementar
    # Entradas:Módulo de Young,
    #     Área da seção transveral
    #     primeiro nó do elemento
    #     segundo nó do elemento
    #     coordenada x do nó inicial
    #     coordenada x do nó final
    #     coordenada y do nó inicial
    #     coordenada y do nó final

    l = np.sqrt((coordxj-coordxi)**2 + (coordyj-coordyi)**2)
    cos_theta = (coordxj-coordxi)/l
    sin_theta = (coordyj-coordyi)/l
    R = np.array([[cos_theta**2, cos_theta*sin_theta, -cos_theta**2, -cos_theta*sin_theta],
           [cos_theta*sin_theta, sin_theta**2, -cos_theta*sin_theta, -sin_theta**2],
           [-cos_theta**2, -cos_theta*sin_theta, cos_theta**2, cos_theta*sin_theta],
           [-cos_theta*sin_theta, -sin_theta**2, cos_theta*sin_theta, sin_theta**2]])
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
coord = np.transpose(np.array([[0,0],[3,0],[0,3],[3,3],[6,3]]))

# Mapa dos nós em cada elemento
inci = np.transpose(np.array([[1,2],[1,4],[3,4],[4,5],[2,4],[2,5]]))     

# Módulo de young
E = 1e6
# Área da seção transversal
A = 1e-2

# Nós presos (restritos)
no_preso = [1,3]

# Aplicação de força

force = np.transpose(np.array([[0,0],[0,0],[0,0],[0,0],[0,-40]]))

# Controle do gráfico

# Plotar gráfico (y/n)

plt_q = 'y'

# Fator de aumento do dano (para efeitos de visualização)

fator = 5 

# Não tem interação do usuário a partir desse ponto--------------------------------

# Montagem das matrizes
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
    for j in range(2):
        eq = IDM[i][j]
        if eq != 0:
            F[eq-1] = force[j][i]   #-1 pq começa na posição zero

# Controle de qualidade
# print('F',F)

# Montagem do sistema (forma eficiente)-------------------------------------------

neq = 2*(size_ID-np.size(no_preso,0)) #número de equações
K = np.zeros((neq,neq))     # pre-alocagem com zeros
LM = np.zeros((2,2),int)      # pre-alocagem com zeros

# Controle de qualidade
# print('neq\n',neq)

for e in range(size_inci):
    for j in range(2):      # Número de nós por elemento
        no = inci[j][e]
        
        # Controle de qualidade
        # print(no)
        
        for i in range(2):  # Range até o número de gl por nó
            eq = IDM[no-1][i]
            LM[i][j] = eq

            # Controle de qualidade
            # print(eq)
            #print(LM)

    for j in range(1,3):
        ii = 2*(j-1)
        for i in range(1,3):
            l = LM[i-1][j-1]

            # Controle de qualidade
            # print(l)

            ll = ii + i
            for n in range(1,3):
                if l != 0:
                    jj = 2*(n-1)
                    for k in range(1,3):
                        m = LM[k-1][n-1]
                        ml = jj + k
                        Ke = matriz_elementar(
                        prop[0][e],
                        prop[1][e],
                        inci[0][e],
                        inci[1][e],
                        coord[0][inci[0][e]-1],
                        coord[0][inci[1][e]-1],
                        coord[1][inci[0][e]-1],
                        coord[1][inci[1][e]-1])

                        # Controle de qualidade
                        # print(l,m,ll,ml)
                        
                        if m != 0:
                            K[l-1][m-1] = K[l-1][m-1] + Ke[ll-1][ml-1]

# Controle de qualidade
# print('K\n',K)

# Solução do problema ------------------------------------------------------------

U = np.matmul(np.linalg.inv(K),F)

V = np.zeros((size_ID,2))
for i in range(size_ID):
    for j in range(2): #graus de liberdade
        eq = IDM[i][j]
        if eq != 0:
            V[i][j] = U[eq-1]

# Controle de qualidade
# print('U\n',U)
# print('V\n',V)

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
    l = np.sqrt((xJ-xI)**2 + (yJ-yI)**2)
    cos_theta = (xJ-xI)/l
    sin_theta = (yJ-yI)/l
    sigma_e = -cos_theta*V[noI][0] - sin_theta*V[noI][1] + cos_theta*V[noJ][0] + sin_theta*V[noJ][1]
    Fe[i][0] = i + 1
    Fe[i][1] = (E*A/l)*sigma_e

# Controle de qualidade
# print(Fe)

plt.title('Força em cada elemento')
plt.xlabel('N° do elemento')
plt.ylabel('Força em N')
plt.bar_label(plt.bar(Fe[:,0],Fe[:,1]))
plt.show()

# Plotagem de gráficos (Não implementado)-----------------------------------------

# xy1 = np.zeros((size_ID+1,2))
# xy2 = np.zeros((size_ID+1,2))

# for i in range(size_ID+1):
#     a = inci[0][i]-1
#     b = inci[1][i]-1
#     xy1[i][0] = coord[0][a]
#     xy1[i][1] = coord[1][a]
#     xy2[i][0] = coord[0][b]
#     xy2[i][1] = coord[1][b]
    
# # print('xy1\n',xy1,'\nxy2\n',xy2)
# plt.plot(xy1[:,0],xy1[:,1])
# plt.show()
