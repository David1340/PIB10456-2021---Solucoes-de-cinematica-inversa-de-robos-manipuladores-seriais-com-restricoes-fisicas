#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação de funções auxiliares à implementação de técnicas de Cinemática Inversa
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022)

import random
import numpy as np
from math import pi, sin, cos, sqrt, atan2

#Calcula a distancia Euclidiana entre dois pontos no R^n
def distancia(a,b,n):
    d = 0
    for i in range(n):
        d = d + (a[i] - b[i])**2      
    return sqrt(d)

#Retorna a Matriz de transformacao Homogeneada usando como entrada os parametros de DH
def matriz_homogenea(d,a,alfa,theta):
    L1 = np.array([cos(theta), -sin(theta)*cos(alfa),\
                sin(theta)*sin(alfa),a*cos(theta)])
    L2 = np.array([sin(theta), cos(theta)*cos(alfa),\
                -cos(theta)*sin(alfa),a*sin(theta)])
    L3 = np.array([0.0,sin(alfa), cos(alfa), d])
    L4 = np.array([0.0,0.0,0.0,1.0])
    A = np.array([L1,L2,L3,L4])
    return A

#Calcula os angulos de RPY a partir de uma matriz de Rotacao
def orientacao(A):
    #calcular os ângulos de orientação na conversão Z -> Y -> X
    R = atan2(A[1,0],A[0,0]) #Roll
    P = atan2(-A[2,0],sqrt((A[2,1]**2)+(A[2,2]**2))) #Pitch
    Y = atan2(A[2,1],A[2,2]) #Yaw
    result = np.array([[R,P,Y]]).T
    return result

#matriz_antissimetrica
def S(a):
    #A = [0,-az,ay ; az,0,-ax ; -ay,ax,0]
    #Uso para calcular produto vetorial entre vetores u x v = S(u) * v
    A = np.zeros((3,3))
    A[0,1] = -a[2]
    A[0,2] = a[1]
    A[1,2] = -a[0]
    A[1,0] = - A[0,1]
    A[2,0] = - A[0,2]
    A[2,1] = - A[1,2]
    return A

#Retorna a Matriz de rotação associada ao ângulos Roll, Pitch Yaw
def matriz_RPY(v):
    #v = (R,P,Y)
    A = np.zeros([3,3])
    A[0,0] = cos(v[0])*cos(v[1])
    A[0,1] = -sin(v[0])*cos(v[2]) + cos(v[0])*sin(v[1])*sin(v[2])
    A[0,2] = sin(v[0])*sin(v[2]) + cos(v[0])*sin(v[1])*cos(v[2])
    A[1,0] = sin(v[0])*cos(v[1])
    A[1,1] = cos(v[0])*cos(v[2]) + sin(v[0])*sin(v[1])*sin(v[2])
    A[1,2] = -cos(v[0])*sin(v[2]) + sin(v[0])*sin(v[1])*cos(v[2])
    A[2,0] = -sin(v[1])
    A[2,1] = cos(v[1])*sin(v[2])
    A[2,2] = cos(v[1])*cos(v[2]) 
    return A

#Gera uma pose alcançável posição + matriz de rotação
def random_pose(): 
    #valor maximo que a junta pode assumir
    qlim = [2.6179,1.5358,2.6179,1.6144,2.6179,1.8413,1.7889]   
    #angulos de juntas iniciais
    q = np.zeros([7,1])
    for a in range(np.size(q)):
        q[a] = random.uniform(-qlim[a],qlim[a])

    #Parâmetros Físicos do manipulador [m]
    base = 0.05 #5 cm
    L = 0.075 #distância da ultima junta a extremidade do efetuador

    #parametros de DH constantes
    d = [0.075 + base,0,0.15,0,0.145,0,0]
    a = [0,0,0,0,0,0.075,0]
    alpha = [pi/2,-pi/2,pi/2,-pi/2,pi/2,pi/2,pi/2]
    # parametros de DH variáveis
    theta = [pi/2 + q[0],q[1],q[2],q[3],q[4],pi/2 + q[5],pi/2 + q[6]]
    #Matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])
    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    p_7 = np.array([[0,0,L,1]]).T
    p_0 = T7@p_7
    return p_0[0:3] , T7[0:3,0:3]

#Calcula a posição das juntas a partir da configuração q
def Cinematica_Direta(q):
    #Pontos de interesse
    L = 0.075 #distância da ultima junta a extremidade do efetuador
    p = np.array([[0,0,0,1]]).T #Base
    p1_1 = np.array([[0,-0.05,0,1]]).T #junta1
    p2_2 = p #junta2
    p3_3 = np.array([[0,-0.075,0,1]]).T #junta3
    p4_4 = p #junta4
    p5_5 = np.array([[0,-0.0725,0,1]]).T #junta5
    p6_6 = np.array([[-0.075,0,0,1]]).T #junta6
    p7_7 = p #junta7
    p8_7 = np.array([[0,0,L,1]]).T #contato de atuação do efetuador

    #Parâmetros Físicos do manipulador [m]
    base = 0.05 #5 cm
    #parametros de DH constantes
    d = [0.075 + base,0,0.15,0,0.145,0,0]
    a = [0,0,0,0,0,0.075,0]
    alpha = [pi/2,-pi/2,pi/2,-pi/2,pi/2,pi/2,pi/2]
    theta = [pi/2 + q[0],q[1],q[2],q[3],q[4],pi/2 + q[5],pi/2 + q[6]]

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    p1_0 = T1@p1_1
    p2_0 = T2@p2_2
    p3_0 = T3@p3_3
    p4_0 = T4@p4_4
    p5_0 = T5@p5_5
    p6_0 = T6@p6_6
    p7_0 = T7@p7_7
    p8_0 = T7@p8_7
    pontos = np.array([p1_0[0:3,0],p2_0[0:3,0],p3_0[0:3,0],p4_0[0:3,0]\
                    ,p5_0[0:3,0],p6_0[0:3,0],p7_0[0:3,0],p8_0[0:3,0]]).T
    return pontos

def Cinematica_Direta2(q):
    #Pontos de interesse
    L = 0.075 #distância da ultima junta a extremidade do efetuador
    p = np.array([[0,0,0,1]]).T #Base
    p1_1 = np.array([[0,-0.05,0,1]]).T #junta1
    p2_2 = p #junta2
    p3_3 = np.array([[0,-0.075,0,1]]).T #junta3
    p4_4 = p #junta4
    p5_5 = np.array([[0,-0.0725,0,1]]).T #junta5
    p6_6 = np.array([[-0.075,0,0,1]]).T #junta6
    p7_7 = p #junta7
    p8_7 = np.array([[0,0,L,1]]).T #contato de atuação do efetuador

    #Parâmetros Físicos do manipulador [m]
    base = 0.05 #5 cm
    #parametros de DH constantes
    d = [0.075 + base,0,0.15,0,0.145,0,0]
    a = [0,0,0,0,0,0.075,0]
    alpha = [pi/2,-pi/2,pi/2,-pi/2,pi/2,pi/2,pi/2]
    theta = [pi/2 + q[0],q[1],q[2],q[3],q[4],pi/2 + q[5],pi/2 + q[6]]

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    p1_0 = T1@p1_1
    p2_0 = T2@p2_2
    p3_0 = T3@p3_3
    p4_0 = T4@p4_4
    p5_0 = T5@p5_5
    p6_0 = T6@p6_6
    p7_0 = T7@p7_7
    p8_0 = T7@p8_7
    pontos = np.array([p1_0[0:3,0],p2_0[0:3,0],p3_0[0:3,0],p4_0[0:3,0]\
                    ,p5_0[0:3,0],p6_0[0:3,0],p7_0[0:3,0],p8_0[0:3,0]]).T
    vetores = np.array([T1[0:3,1],T2[0:3,1],T3[0:3,1],T4[0:3,1],T5[0:3,1],T6[0:3,1]\
                        ,T7[0:3,1]]).T
    return [pontos,vetores]


def Cinematica_Direta3(q):
    #Pontos de interesse
    L = 0.075 #distância da ultima junta a extremidade do efetuador
    p = np.array([[0,0,0,1]]).T #Base
    p8_7 = np.array([[0,0,L,1]]).T #contato de atuação do efetuador

    #Parâmetros Físicos do manipulador [m]
    base = 0.05 #5 cm
    #parametros de DH constantes
    d = [0.075 + base,0,0.15,0,0.145,0,0]
    a = [0,0,0,0,0,0.075,0]
    alpha = [pi/2,-pi/2,pi/2,-pi/2,pi/2,pi/2,pi/2]
    theta = [pi/2 + q[0],q[1],q[2],q[3],q[4],pi/2 + q[5],pi/2 + q[6]]


    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    p8_0 = T7@p8_7
    return [p8_0[0:3,0] ,T7[0:3,0:3]]

#Projeta um ponto em um plano
def projecao_ponto_plano(normal,p0,p):
    #normal -> vetor normal ao plano
    #p0 -> ponto pertecente ao plano
    #p -> ponto a ser projetado
    #constante d da equação do plano
    d = -normal.T@p0 #produto escala 
    #distancia do ponto ao plano
    alpha = (-d - normal.T@p)/(normal[0,0]**2 +  normal[1,0]**2 + normal[2,0]**2)
    ponto_projetado = p + alpha*normal
    return ponto_projetado

#realiza a operação de multiplicação de quaternios
def multiplicacao_quaternios(q,q2):  
    resultado = np.zeros([4,1])
    resultado[0,0] = q[0,0]*q2[0,0] -q[1,0]*q2[1,0] -q[2,0]*q2[2,0] -q[3,0]*q2[3,0] 
    resultado[1,0] = q[0,0]*q2[1,0] +q[1,0]*q2[0,0] +q[2,0]*q2[3,0] -q[3,0]*q2[2,0] 
    resultado[2,0] = q[0,0]*q2[2,0] -q[1,0]*q2[3,0] +q[2,0]*q2[0,0] +q[3,0]*q2[1,0] 
    resultado[3,0] = q[0,0]*q2[3,0] +q[1,0]*q2[2,0] -q[2,0]*q2[1,0] +q[3,0]*q2[0,0] 
    return resultado

#gira p em torno de v em th rad
def rotationar_vetor(p,v,th):
    a = cos(th/2)
    if(norm(v) > 0.0001): v = v/norm(v)
    b = v[0,0]*sin(th/2)
    c = v[1,0]*sin(th/2)
    d = v[2,0]*sin(th/2)
    p_aumentado = np.zeros([4,1])
    p_aumentado[1:4,0] = p[:,0]
    h = np.array([[a,b,c,d]]).T
    hx = np.array([[a,-b,-c,-d]]).T
    p_r = multiplicacao_quaternios(h,p_aumentado)
    q_r = multiplicacao_quaternios(p_r,hx)
    return q_r[1:4]

#Calcula a norma de um vetor 
def norm(v):
    return sqrt(v[[0]]**2 + v[[1]]**2 + v[[2]]**2)

def matriz_homogenea_final(d,a,alpha,theta,metodo = None):
    Tn = np.eye(4)
    n = len(a)
    if(metodo == 'CCD'):
        Vy = []
        for i in range(n):
            Tn = Tn@matriz_homogenea(d[i],a[i],alpha[i],theta[i])
            Vy.append(Tn[0:3,1])
        return Tn,np.array(Vy).T
    else:
        for i in range(n):
            Tn = Tn@matriz_homogenea(d[i],a[i],alpha[i],theta[i])
        
        return Tn


