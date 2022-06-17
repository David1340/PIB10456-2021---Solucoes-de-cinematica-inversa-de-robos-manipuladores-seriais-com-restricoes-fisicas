#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do Particle Swarm Optimizarion - P
#para encontrar encontrar uma configuração q
#dada uma posição (x,y,z) e uma orientacao 
#no espaço para o Pioneer 7DOF.
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import das bibliotecas python
from random import random,uniform
import numpy as np

#Import do modulo funcoes.py
from funcoes import distancia, orientacao
from manipulador_15dof import *

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

    def update_position(self,qbest,L): #Atualiza a posição da particula/configuração do robô
        c1 = 1 #grupo
        c2 = 1 #individual
        for i in range(self.n):
            w = 0.5 #+ random()/2
            vmax = np.inf
            #w = random()
            self.v[i] = w*self.v[i] + c1*random()*(qbest[i] - self.p[i])+ c2*random()*(self.bp[i] - self.p[i])
            if(self.v[i] > vmax):
                self.v[i] = vmax
            elif(self.v[i] < -vmax):
                self.v[i] = -vmax
            self.p[i] = self.p[i] + self.v[i]
            if(self.p[i] > L[i]):
                self.p[i] = L[i]
            elif(self.p[i] < -L[i]):
                self.p[i] = -L[i]
       
            
    def update_fuction(self,o,o2): #Calcula a função de custo/fitness da particula
        #(posição,orientacao) da pose desejada
        p,orient = Cinematica_Direta3(self.p)

        #calculo do erro em módulo da orientacao desejada e da particula
        self.o = distancia(orientacao(orient),o2,3)

        #calculo da distancia euclidiana da posição do efetuador em relação ao objetivo
        self.d = distancia(p,o,3)

        #Calculo da função de custo       
        k1 = 0 #orientacao
        k2 = 1 #posição
        self.f = (k1*self.o) + (k2*self.d)
        if(self.f < self.bf):
            self.bf = self.f
            self.bp = self.p.copy()

def PSO2(o,o2,number,n,L,erro_min,Kmax):
    #numero limite de interações
    k = Kmax     
    q = []
    
    #criando as particulas de dimensão n e calculando o valor de sua função de custo
    for i in range(number):
        p = []
        for i2 in range(n):
            p.append(uniform(-L[i2],L[i2]))

        q.append(particle(p,n))
        q[i].update_fuction(o,o2)

    #Criando a configuração qbest e sua função de custo 
    qbest = q[0].p.copy()
    f = q[0].f
    for i in range(number):
        if(f > q[i].f):
            qbest = q[i].p.copy()
            f = q[i].f
            
    #Executando PSO
    for j in range(k):
        for i in range(number):          
            q[i].update_position(qbest,L)
            q[i].update_fuction(o,o2)
            
            #Se alguma particula possui função de custo menor do que qbest, ela se torna a qbest
            if(f > q[i].f): 
                qbest = q[i].p.copy()
                f = q[i].f
        #Atualiza a configuração do robô no Rviz
        #Critério de parada
        if(f <= erro_min):
            break;   

    return [f,j+1]



def PSO(posicaod,orientacaod,erro_min,Kmax):

    orientacaod = orientacao(orientacaod)
    numero_particulas = 200
    dimensao = getNumberJoints() #dimensão do robô
    #restrições de cada ângulo
    L = getLimits()
    return PSO2(posicaod,orientacaod,numero_particulas,dimensao,L,erro_min,Kmax)