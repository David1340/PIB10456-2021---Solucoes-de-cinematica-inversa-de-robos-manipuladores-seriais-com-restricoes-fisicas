#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do Damped Last Square 
#para encontrar encontrar uma configuração q
#dada uma posição (x,y,z) .
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import das bibliotecas
import numpy as np

#Import das minhas funções
from funcoes import distancia
#Import das funções associadas ao manipulador usado
from pioneer_7dof import *

def DLS(posicaod,q,erro_min,Kmax):
    #variaveis do método   
    alfa = 1 #Constante alpha para melhorar  a aproximação  
    lbd = 0.005 #Constante de amortecimento
    I = np.eye(3)   
    qmax = 0.1 #passo maximo entre atualizacoes das juntas
    K = Kmax #número máximo de iterações

    #valor maximo que a junta pode assumir
    qlim = getLimits()

    #vetores colunas do sistema de coordenadas global
    z = np.array([[0,0,1,1]]).T
    o = np.array([[0,0,0,1]]).T #origem

    for k in range(K):
        
        J,pn_0 = jacobianoGeometrico(q)
        
        #Calcula a distancia entre o efetuador a o objetivo(Posição)
        erro = distancia(pn_0,posicaod,3)

        #Condição de parada
        if(erro < erro_min):
            break  

        #erro
        f = posicaod - pn_0[0:3] 
        #Equação do DLS
        dq =  alfa*((J.T@np.linalg.inv(J@J.T + lbd*I))@f)
              
        #limitando o delta q
        for i in range(np.size(dq)):
            if(dq[i] > qmax):
                dq[i] = qmax
            elif(dq[i] < -qmax):
                dq[i] = -qmax 

        #Atualizando a configuração do robô
        q = q + dq

        for i in range(np.size(q)):
            if(q[i] > qlim[i]):
                q[i] = qlim[i]
            elif(q[i] < -qlim[i]):
                q[i] = -qlim[i]       
    return [erro,k+1]