#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do Fully Resampled Particle Swarm Optimizarion
#para encontrar encontrar uma configuração q
#dada uma posição (x,y,z) e uma orientacao 
#no espaço para o Pioneer 7DOF.
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import do modulo funcoes.py
import sys
sys.path.append('C:\PIBIC 2022 - Python')
from funcoes import matriz_homogenea, distancia, orientacao


#Import das bibliotecas python
from random import uniform
from math import exp
import numpy as np
from pioneer_7dof import getDH_paramaters, getLimits

class particle:
    def __init__(self,position,dimension):
        self.p = position #posição atual da particula/configuração do robô
        self.v = np.zeros(dimension) #velocidade atual da particula
        self.bp = position.copy() #melhor posição que a particula ja esteve
        self.n = dimension #dimensão da particula
        self.d = 0 #Diferença em módulo da distância atual para a desejada
        self.o = np.array([0,0,0]) #Diferença em módulo da orientacao atual para a desejada
        self.f = np.Inf #Função de custo/fitnees atual da particula
        self.bf = self.f #Melhor valor de função de custo da obtida pela particula
                   
    def update_fuction(self,o,o2): #Calcula a função de custo/fitness da particula
        #(posição,orientacao) da pose desejada

        #Parâmetros de DH e ponto de atuação
        d,a,alpha,theta,p = getDH_paramaters(self.p)

        #Matrizes homogêneas
        A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
        A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
        A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
        A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
        A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
        A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
        A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])
        A = A1@(A2@(A3@(A4@(A5@(A6@A7)))))

        #posição do efetuador em relação ao sistema de coordenadas global
        p = A@p

        #calculo do erro em módulo da orientacao desejada e da particula
        self.o = distancia(orientacao(A),o2,3)

        #calculo da distancia euclidiana da posição do efetuador em relação ao objetivo
        self.d = distancia(p,o,3)

        #Calculo da função de custo       
        k1 = 0 #orientacao
        k2 = 1 #posição
        self.f = (k1*self.o) + (k2*self.d)
        if(self.f < self.bf):
            self.bf = self.f
            self.bp = self.p.copy()

def FRPSO2(o,o2,number,n,L,erro_min,Kmax):
    #numero limite de interações
    k = Kmax     
    q = []
    Nbests = 5
    tau = 20
    #criando as particulas de dimensão n e calculando o valor de sua função de custo
    for i in range(number):
        p = np.array([uniform(-L[0],L[0]),uniform(-L[1],L[1]),uniform(-L[2],L[2]),uniform(-L[3],L[3])\
                      ,uniform(-L[4],L[4]),uniform(-L[5],L[5]),uniform(-L[6],L[6])])
        q.append(particle(p,n))
        q[i].update_fuction(o,o2)

    #Criando as configurações qbests e sua funções de custos
    qbests = []
    qvalues = []
    f = np.inf
    for i in range(Nbests):      
        qbests.append(q[i].p.copy())
        qvalues.append(q[i].f)
        if(f <  q[i].f):
            f = q[i].f

    for i in range(number):
        if(max(qvalues) > q[i].f):
            for j in range(Nbests):
                if(qvalues[j] == max(qvalues)):
                    qvalues[j] = q[i].f
                    qbests[j] = q[i].p.copy()
                    break
            f = min(qvalues)
            
    #Executando PSO
    for j in range(k):
        q = []

        sig = np.sqrt(1 - 0.9999 * exp(-f/tau))
        for N in range(Nbests):
            for i in range(int(number/(Nbests +1))):
                p = sig*np.random.randn(n)
                for i2 in range(n):
                    p[i2] = qbests[N][i2] + p[i2]
                q.append(particle(p,n))

        for i in range(number - len(q)):
            p = np.array([uniform(-L[0],L[0]),uniform(-L[1],L[1]),uniform(-L[2],L[2]),uniform(-L[3],L[3])\
                ,uniform(-L[4],L[4]),uniform(-L[5],L[5]),uniform(-L[6],L[6])])
            q.append(particle(p,n))

        for i in range(number):          
            q[i].update_fuction(o,o2)

            if(max(qvalues) > q[i].f):
                for i2 in range(Nbests):
                    if(qvalues[i2] == max(qvalues)):
                        qvalues[i2] = q[i].f
                        qbests[i2] = q[i].p.copy()
                        break
                f = min(qvalues)   
        #Atualiza a configuração do robô no Rviz
        #Critério de parada
        if(f <= erro_min):
            break;   

    return [f,j+1]



def FRPSO(posicaod,orientacaod,erro_min,Kmax):

    orientacaod = orientacao(orientacaod)
    numero_particulas = 200
    dimensao = 7 #dimensão do robô
    #restrições de cada ângulo
    L = [2.6179,1.5358,2.6179,1.6144,2.6179,1.8413,1.7889]
    return FRPSO2(posicaod,orientacaod,numero_particulas,dimensao,L,erro_min,Kmax)