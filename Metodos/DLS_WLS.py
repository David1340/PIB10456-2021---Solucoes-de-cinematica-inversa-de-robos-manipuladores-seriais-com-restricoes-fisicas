#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do Damped Last Square usando a solução de norma mínima ponderada
#para encontrar encontrar uma configuração q
#dada uma posição (x,y,z) 
#no espaço para o Pioneer 7DOF.
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import das bibliotecas python
import numpy as np

#Import das minhas funções
from funcoes import distancia
#Import das funções associadas a manipulador usado
from pioneer_7dof import *

def DLS_WLS(posicaod,q,erro_min,Kmax):
    #variaveis do método
    n = getNumberJoints()
    #Constante alpha para melhorar  a aproximação
    alfa = 1
    #Constante de amortecimento
    lbd = 0.005
    I = np.eye(3)
    #passo maximo entre atualizacoes das juntas
    qmax = 0.1 
    K = Kmax #número máximo de iterações

    #valor maximo que a junta pode assumir
    qlim = getLimits()
    thmax = np.array(qlim)
    thmin = -thmax

    for k in range(K):
        
        J,pn_0 = jacobianoGeometrico(q)  
        
        #Calcula a distancia entre o efetuador a o objetivo(Posição)
        erro = distancia(pn_0,posicaod,3)

        #Condição de parada
        if(erro < erro_min):
            break  
        #erro
        f = posicaod - pn_0[0:3]

        #Matriz de pesos
        W = np.zeros([n,n])

        for i in range(n):
            num = ((thmax[i] - thmin[i])**2)*(2*q[i] - thmax[i]-thmin[i]) #numerador
            den = 4*((thmax[i]  - q[i])**2)*((q[i]-thmin[i])**2) #denominador
            W[i,i] = 1 + np.abs(num/den)

        #Equação do DLS com WLS
        Wi = np.linalg.inv(W)
        dq = alfa*Wi@J.T@np.linalg.inv(J@Wi@J.T + lbd*I)@f

        #limitando o delta q
        for i in range(np.size(dq)):
            if(dq[i] > qmax):
                dq[i] = qmax
            elif(dq[i] < -qmax):
                dq[i] = -qmax 

        #Atualizando a configuração do robô
        q = q + dq

        #limitando os limites das juntas, (nos meus testes nunca precisou)
        for i in range(np.size(q)):
            if(q[i] > qlim[i]):
                q[i] = qlim[i] - 0.001
            elif(q[i] < -qlim[i]):
                q[i] = -qlim[i] + 0.001 
    return [erro,k+1]

