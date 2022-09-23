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

#Import das bibliotecas python
from random import uniform
from math import exp
import numpy as np

#Import das minhas funções
from funcoes import distancia, orientacao
#Import das funções associadas ao manipulador usado
from pioneer_7dof import *

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
        limits = getLimits()
        for (qi,li) in zip(self.p, limits):
            if(np.abs(qi) > np.abs(li)):
                self.f = np.Inf
                return 

        p,orient = Cinematica_Direta3(self.p)

        #calculo do erro em módulo da orientacao desejada e da particula
        self.o = distancia(orientacao(orient),o2,3)

        #calculo da distancia euclidiana da posição do efetuador em relação ao objetivo
        self.d = distancia(p,o,3)

        #Calculo da função de custo       
        k1 = 0.1 #orientacao
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
        p = []
        for i2 in range(n):
            p.append(uniform(-L[i2],L[i2]))

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

        sig = np.sqrt(1 - exp(-f/tau))
        #sig = f/tau
        for N in range(Nbests):
            for i in range(int(number/(Nbests +1))):
                p = sig*np.random.randn(n)
                #p = -2*sig*np.random.random(n) + sig
                for i2 in range(n):
                    p[i2] = qbests[N][i2] + p[i2]
                q.append(particle(p,n))

        for i in range(number - len(q)):
            p = []
            for i2 in range(n):
                p.append(uniform(-L[i2],L[i2]))
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
                qBest = qbests[np.argmin(qvalues)]   
        #Atualiza a configuração do robô no Rviz
        #Critério de parada
        if(f <= erro_min):
            break;   

    return [f,j+1,qBest]

def FRPSO(posicaod,orientacaod,erro_min,Kmax):

    orientacaod = orientacao(orientacaod)
    numero_particulas = 200
    dimensao = getNumberJoints() #dimensão do robô
    #restrições de cada ângulo
    L = getLimits()

    f,k,qBest = FRPSO2(posicaod,orientacaod,numero_particulas,dimensao,L,erro_min,Kmax)
    #plot(qBest,posicaod)

    return [f,k]