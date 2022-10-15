#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do parâmetros de DH do manipulador pioneer7DOF e algumas funções auxiliares
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

from  math import pi
import numpy as np
from funcoes import matriz_homogenea, S, Esfera, deteccao_de_colisao
from random import uniform

def getDH_paramaters(q):
    elos = getLengthElos()
    base = elos[0]
    d = [elos[1] + base,
        0, elos[2] + elos[3], 0, elos[4] + elos[5], 0, elos[6] + elos[7],
        0, 0]
    a = [0,
        0, 0, 0, 0, 0, 0,
       elos[8],0]

    alpha = [pi/2,
        -pi/2, pi/2, -pi/2, pi/2, -pi/2, pi/2,
        pi/2, pi/2]

    theta = [pi/2 + q[0],q[1],q[2],
        q[3], q[4], q[5], q[6],
        pi/2 + q[7], pi/2 + q[8]]

    L = 0.075 #distância da ultima junta a extremidade do efetuador
    p_n = np.array([[0,0,L,1]]).T #ponto de atuação do manipulador no sistema de coordenadas on xn yn zn
    return d,a,alpha,theta,p_n.copy()

def getRaio():
    return 0.025

def getLimits():
    qlim = [2.6179,
        1.6144, 2.6179, 1.6144, 2.6179, 1.6144, 2.6179,
        1.8413, 1.7889]
    return qlim

def getNumberJoints():
    n = 9 #número de juntas
    return n

def getAlphasNegatives():
    return [False,
            True,False,True,False,True,False,
            False,False] 

#Calcula a posição das juntas a partir da configuração q
def Cinematica_Direta(q,orientacao = False):
    #Pontos de interesse
    elos = getLengthElos()
    p = np.array([[0,0,0,1]]).T #Base
    p1_1 = np.array([[0,-elos[0],0,1]]).T #junta1

    p2_2 = p #junta2
    p3_3 = np.array([[0,-elos[2],0,1]]).T #junta3
    p4_4 = p #junta4
    p5_5 = np.array([[0,-elos[4],0,1]]).T #junta5

    p6_6 = p #junta6
    p7_7 = np.array([[0,-elos[6],0,1]]).T #junta7

    p8_8 = np.array([[-elos[7],0,0,1]]).T #junta8
    p9_9 = p #junta9

    d,a,alpha,theta,pn_9 = getDH_paramaters(q)

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])
    A8 = matriz_homogenea(d[7],a[7],alpha[7],theta[7])
    A9 = matriz_homogenea(d[8],a[8],alpha[8],theta[8])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    T8 = T7@A8
    T9 = T8@A9

    p1_0 = T1@p1_1
    p2_0 = T2@p2_2
    p3_0 = T3@p3_3
    p4_0 = T4@p4_4
    p5_0 = T5@p5_5
    p6_0 = T6@p6_6
    p7_0 = T7@p7_7
    p8_0 = T8@p8_8
    p9_0 = T9@p9_9
    pn_0 = T9@pn_9

    pontos = np.array([p1_0[0:3,0],p2_0[0:3,0],p3_0[0:3,0],p4_0[0:3,0]\
                    ,p5_0[0:3,0],p6_0[0:3,0],p7_0[0:3,0],p8_0[0:3,0]
                    ,p9_0[0:3,0],pn_0[0:3,0]]).T

    if(orientacao):
        return pontos,T9[0:3,0:3]

    return pontos

#Calcula as posições das juntas e seus eixos de atuação
def Cinematica_Direta2(q):
    #Pontos de interesse
    elos = getLengthElos()
    p = np.array([[0,0,0,1]]).T #Base
    p1_1 = np.array([[0,-elos[0],0,1]]).T #junta1

    p2_2 = p #junta2
    p3_3 = np.array([[0,-elos[2],0,1]]).T #junta3
    p4_4 = p #junta4
    p5_5 = np.array([[0,-elos[4],0,1]]).T #junta5

    p6_6 = p #junta6
    p7_7 = np.array([[0,-elos[6],0,1]]).T #junta7

    p8_8 = np.array([[-elos[7],0,0,1]]).T #junta8
    p9_9 = p #junta9

    d,a,alpha,theta,pn_9 = getDH_paramaters(q)

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])
    A8 = matriz_homogenea(d[7],a[7],alpha[7],theta[7])
    A9 = matriz_homogenea(d[8],a[8],alpha[8],theta[8])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    T8 = T7@A8
    T9 = T8@A9

    p1_0 = T1@p1_1
    p2_0 = T2@p2_2
    p3_0 = T3@p3_3
    p4_0 = T4@p4_4
    p5_0 = T5@p5_5
    p6_0 = T6@p6_6
    p7_0 = T7@p7_7
    p8_0 = T8@p8_8
    p9_0 = T9@p9_9
    pn_0 = T9@pn_9

    pontos = np.array([p1_0[0:3,0],p2_0[0:3,0],p3_0[0:3,0],p4_0[0:3,0]\
                    ,p5_0[0:3,0],p6_0[0:3,0],p7_0[0:3,0],p8_0[0:3,0]
                    ,p9_0[0:3,0],pn_0[0:3,0]]).T

    vetores = np.array([T1[0:3,1],T2[0:3,1],T3[0:3,1],T4[0:3,1],T5[0:3,1],T6[0:3,1]\
                        ,T7[0:3,1],T8[0:3,1],T9[0:3,1]]).T
    return [pontos,vetores]

#Calcula a posição e a matriz de rotação do efetuador final
def Cinematica_Direta3(q):
    #Pontos de interesse
    d,a,alpha,theta,p10_9 = getDH_paramaters(q)

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])
    A8 = matriz_homogenea(d[7],a[7],alpha[7],theta[7])
    A9 = matriz_homogenea(d[8],a[8],alpha[8],theta[8])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    T8 = T7@A8
    T9 = T8@A9

    p10_0 = T9@p10_9
    return [p10_0[0:3] ,T9[0:3,0:3]]

#Gera uma pose alcançável posição + matriz de rotação
def random_pose(esferas = []): 
    #valor maximo que a junta pode assumir
    qlim = getLimits() 
    n = getNumberJoints()
    #angulos de juntas iniciais
    q = np.zeros([n,1])
    for a in range(np.size(q)):
        q[a] = uniform(-qlim[a],qlim[a])

    colidiu = True
    while(colidiu):
        colidiu = False
        
        for a in range(np.size(q)):
            q[a] = uniform(-qlim[a],qlim[a])

        pontos = Cinematica_Direta(q)
        end = np.shape(pontos)[1] -1

        for esfera in esferas:
            for i in range(end):
                r = esfera.get_raio() + getRaio()
                if(deteccao_de_colisao(pontos[:,i],pontos[:,i+1],esfera.get_centro(),r)):
                    colidiu = True
                    break

    return Cinematica_Direta3(q)

#Calcula o jacobiano considerando apenas a posição
def jacobianoGeometrico(q):

    #vetores colunas do sistema de coordenadas global
    z = np.array([[0,0,1,1]]).T
    o = np.array([[0,0,0,1]]).T #origem
    n = 9
    #Parâmetros de DH e ponto de atuação
    d,a,alpha,theta,p10_9 = getDH_paramaters(q)

    #Matrizes homogêneas
    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])
    A8 = matriz_homogenea(d[7],a[7],alpha[7],theta[7])
    A9 = matriz_homogenea(d[8],a[8],alpha[8],theta[8])

    #Definindo os pontos de interesse em seus sistemas locais
    o1_1 = o #origem do SC 1
    o2_2 = o #origem do SC 2
    o3_3 = o #origem do SC 3
    o4_4 = o #origem do SC 4
    o5_5 = o #origem do SC 5
    o6_6 = o #origem do SC 6
    o7_7 = o #origem do SC 7
    o8_8 = o #origem do SC 8
    o9_9 = o #origem do SC 9

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    T8 = T7@A8
    T9 = T8@A9

    o1_0 = T1@o1_1
    o2_0 = T2@o2_2
    o3_0 = T3@o3_3
    o4_0 = T4@o4_4
    o5_0 = T5@o5_5
    o6_0 = T6@o6_6
    o7_0 = T7@o7_7
    o8_0 = T8@o8_8
    o9_0 = T9@o9_9

    p10_0 = T9@p10_9

    #os vetores z serao transformados em vetores  no R^3
    z0_0 = z[0:3]
    z1_0 = (T1@z)[0:3]
    z2_0 = (T2@z)[0:3]
    z3_0 = (T3@z)[0:3]
    z4_0 = (T4@z)[0:3]
    z5_0 = (T5@z)[0:3]
    z6_0 = (T6@z)[0:3]
    z7_0 = (T7@z)[0:3]
    z8_0 = (T8@z)[0:3]
    #z9_0 = (T9@z)[0:3] nao eh usado.

    #cálculo do Jacobiano geométrico
    J = np.zeros([3,n])

    #produto vetorial de Z0_0 por (o7_0 - o) 
    J[:,0] = S(z0_0)@(o9_0[0:3] - o[0:3])[:,0]
    J[:,1] = S(z1_0)@(o9_0[0:3] - o1_0[0:3])[:,0]
    J[:,2] = S(z2_0)@(o9_0[0:3] - o2_0[0:3])[:,0]
    J[:,3] = S(z3_0)@(o9_0[0:3] - o3_0[0:3])[:,0]
    J[:,4] = S(z4_0)@(o9_0[0:3] - o4_0[0:3])[:,0]
    J[:,5] = S(z5_0)@(o9_0[0:3] - o5_0[0:3])[:,0]
    J[:,6] = S(z6_0)@(o9_0[0:3] - o6_0[0:3])[:,0]
    J[:,7] = S(z7_0)@(o9_0[0:3] - o7_0[0:3])[:,0]
    J[:,8] = S(z8_0)@(o9_0[0:3] - o8_0[0:3])[:,0]

    return J,p10_0

#Calcula o jacobiano considerando apenas a posição e aorientação
def jacobianoGeometrico2(q):

    #vetores colunas do sistema de coordenadas global
    z = np.array([[0,0,1,1]]).T
    o = np.array([[0,0,0,1]]).T #origem

    #Parâmetros de DH e ponto de atuação
    d,a,alpha,theta,p10_9 = getDH_paramaters(q)

    #Matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])
    A8 = matriz_homogenea(d[7],a[7],alpha[7],theta[7])
    A9 = matriz_homogenea(d[8],a[8],alpha[8],theta[8])

    #Definindo os pontos de interesse em seus sistemas locais
    o1_1 = o #origem do SC 1
    o2_2 = o #origem do SC 2
    o3_3 = o #origem do SC 3
    o4_4 = o #origem do SC 4
    o5_5 = o #origem do SC 5
    o6_6 = o #origem do SC 6
    o7_7 = o #origem do SC 7
    o8_8 = o #origem do SC 8
    o9_9 = o #origem do SC 9
    
    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    T8 = T7@A8
    T9 = T8@A9

    o1_0 = T1@o1_1
    o2_0 = T2@o2_2
    o3_0 = T3@o3_3
    o4_0 = T4@o4_4
    o5_0 = T5@o5_5
    o6_0 = T6@o6_6
    o7_0 = T7@o7_7
    o8_0 = T8@o8_8
    o9_0 = T9@o9_9
    p10_0 = T9@p10_9  

    #os vetores z serao transformados em vetores  no R^3
    z0_0 = z[0:3]
    z1_0 = (T1@z)[0:3]
    z2_0 = (T2@z)[0:3]
    z3_0 = (T3@z)[0:3]
    z4_0 = (T4@z)[0:3]
    z5_0 = (T5@z)[0:3]
    z6_0 = (T6@z)[0:3]
    z7_0 = (T7@z)[0:3]
    z8_0 = (T8@z)[0:3]
    #z9_0 = (T9@z)[0:3] nao eh usado.

    #cálculo do Jacobiano geométrico
    J = np.zeros([6,13])

    #produto vetorial de Z0_0 por (o7_0 - o) 
    J[0:3,0] = S(z0_0)@(o9_0[0:3] - o[0:3])[:,0]
    J[3:6,0] = z0_0[:,0]
    J[0:3,1] = S(z1_0)@(o9_0[0:3] - o1_0[0:3])[:,0]
    J[3:6,1] = z1_0[:,0]
    J[0:3,2] = S(z2_0)@(o9_0[0:3] - o2_0[0:3])[:,0]
    J[3:6,2] = z2_0[:,0]
    J[0:3,3] = S(z3_0)@(o9_0[0:3] - o3_0[0:3])[:,0]
    J[3:6,3] = z3_0[:,0]
    J[0:3,4] = S(z4_0)@(o9_0[0:3] - o4_0[0:3])[:,0]
    J[3:6,4] = z4_0[:,0]
    J[0:3,5] = S(z5_0)@(o9_0[0:3] - o5_0[0:3])[:,0]
    J[3:6,5] = z5_0[:,0]
    J[0:3,6] = S(z6_0)@(o9_0[0:3] - o6_0[0:3])[:,0]
    J[3:6,6] = z6_0[:,0]
    J[0:3,7] = S(z7_0)@(o9_0[0:3] - o7_0[0:3])[:,0]
    J[3:6,7] = z7_0[:,0]
    J[0:3,8] = S(z8_0)@(o9_0[0:3] - o8_0[0:3])[:,0]
    J[3:6,8] = z8_0[:,0]

    return J,p10_0,T9

def getLengthElos():
    k = 1
    return  k*np.array([0.05,
            0.075,0.075,0.075,0.075,0.075,
            0.075,0.075,0.075,0.075]) 

def getTypeJoints():
    return ['pivot',
            'hinge','pivot','hinge','pivot',
            'hinge','pivot','hinge','hinge']