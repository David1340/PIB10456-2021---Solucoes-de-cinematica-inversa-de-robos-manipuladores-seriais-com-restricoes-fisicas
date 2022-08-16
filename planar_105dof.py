#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do parâmetros de DH do manipulador pioneer7DOF e algumas funções auxiliares
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

from  math import pi
import numpy as np
from funcoes import matriz_homogenea
from random import uniform

def getDH_paramaters(q):

    d = 7*[0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0]
    a = 7*[0.075, 0.075, 0.075, 0.075, 0.075,
        0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075, 0.075,
        0.075, 0.075]

    alpha = 7*[0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0]

    theta = [-pi/2 + q[0]]
    for i in range(104):
        theta = theta + [q[i+1]]

    L = 0 #distância da ultima junta a extremidade do efetuador
    p_n = np.array([[0,0,L,1]]).T #ponto de atuação do manipulador no sistema de coordenadas on xn yn zn
    return d,a,alpha,theta,p_n.copy()

def getLimits():
    qlim = 7*[2.6179,1.6144,2.6179,1.6144,2.6179,
        1.6144, 2.6179, 1.6144, 2.6179, 1.6144, 2.6179, 1.6144, 2.6179,
        1.8413, 1.7889]
    return qlim

def getNumberJoints():
    n = 105 #15 #número de juntas
    return n

#Calcula a posição e a matriz de rotação do efetuador final
def Cinematica_Direta3(q):
    #Pontos de interesse
    d,a,alpha,theta,p61_60 = getDH_paramaters(q)
    Tn = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    #Calculando as matrizes homogêneas
    for i in range(getNumberJoints() -1):
        Tn = Tn@matriz_homogenea(d[i+1],a[i+1],alpha[i+1],theta[i+1])

    p61_0 = Tn@p61_60
    return [p61_0[0:3] ,Tn[0:3,0:3]]

#Gera uma pose alcançável posição + matriz de rotação
def random_pose(): 
    #valor maximo que a junta pode assumir
    qlim = getLimits() 
    n = getNumberJoints()
    #angulos de juntas iniciais
    q = np.zeros([n,1])
    for a in range(np.size(q)):
        q[a] = uniform(-qlim[a],qlim[a])

    return Cinematica_Direta3(q)

