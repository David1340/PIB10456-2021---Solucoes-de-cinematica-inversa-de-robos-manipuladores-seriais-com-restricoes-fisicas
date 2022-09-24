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

def deteccao_de_colisao(pr1,pr2,p,raio):
    v = pr2 - pr1
    u = p - pr1
    #Seja w a projeção de u em v então w = k*v
    #em que w = dot(u,v)/norm(v)^2
    k = np.dot(u,v)/np.sum(np.square(v))
    if(k > 1 or k < 0): #Se a projeção de p não está no segmento de reta pr1->pr2
        return 0
    #d = |u x v|/|v|^2 distância entre ponto e reta
    d = np.sqrt(np.sum(np.square(np.cross(u,v))))/np.sum(np.square(v))
    if(d <= raio):
        return 1
    else:
        return 0

class Esfera:
    def __init__(self,x,y,z,r):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
    def get_centro(self):
        return np.array([self.x,self.y,self.z])
    def get_raio(self):
        return self.r



