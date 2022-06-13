#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do Método do Gradiente Descendente/Transposto
#para encontrar encontrar uma configuração q
#dada uma posição (x,y,z) com orientação definida
#no espaço para o Pioneer 7DOF.
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import das bibliotecas python
import numpy as np
from pioneer_7dof import *

#Import do modulo funcoes.py
from funcoes import  distancia, S, orientacao

def GradDesc_Completo(posicaod,orientacaod,q,erro_min,Kmax):

    #variaveis do método
    #passo maximo entre atualizacoes das juntas
    qmax = 0.1 
    K = Kmax #número máximo de iterações
    #valor maximo que a junta pode assumir
    qlim = getLimits()

    rpyd = orientacao(orientacaod)

    for k in range(K):

        J,pn_0,Tn = jacobianoGeometrico2(q)     

        #Condição de parada   
        errop =  distancia(pn_0,posicaod,3) #erro de posiçao
        rpy = orientacao(Tn[0:3,0:3]) #angulos Roll, Pitch Yall
        erroa = distancia(rpy,rpyd,3) #erro angular
        c1 = 1
        c2 = 0.1
        erro = c1*errop + c2*erroa #composição de erro 
        
        if(erro < erro_min):
            break      

        #Calculo da erro
        f = np.zeros([6,1])
        f[0:3] = posicaod - pn_0[0:3]
        #a parte angular peguei do artigo A closed-loop inverse kinematic scheme
        #online joint based robot control
        f[3:6] = 0.5*(S(Tn[0:3,0:1])@orientacaod[:,0:1]+ S(Tn[0:3,1:2])@orientacaod[:,1:2]  + \
            S(Tn[0:3,2:3])@orientacaod[:,2:3])
        
        #cálculo da constante de passo c
        e = np.array([f[0,0],f[1,0],f[2,0],f[3,0],f[4,0],f[5,0]])
        c = e.dot(J@J.T@e)/((J@J.T@e).dot(J@J.T@e))#peguei essa equacao do
        #artigo Inverse Kinematics Techniques in Computer Graphics: A Survey

        #Equação do Jacobiano Transposto
        dq = c*(J.T@f)

        #limitando o delta q
        for i in range(np.size(dq)):
            if(dq[i] > qmax):
                dq[i] = qmax
            elif(dq[i] < -qmax):
                dq[i] = -qmax 

        #Atualizando a cofiguração
        q = q + dq

        #Limitando os valores das juntas
        for i in range(np.size(q)):
            if(q[i] > qlim[i]):
                q[i] = qlim[i]
            elif(q[i] < -qlim[i]):
                q[i] = -qlim[i]

    return [erro,k+1]

