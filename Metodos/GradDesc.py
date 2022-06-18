#Autor David Oliveira
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS
#Implementação do Método Gradiente Descendente/Jacobiano Transposto 
#para encontrar uma configuação q 
#dada uma posição (x,y,z) no espaço para o Pioneer 7DOF
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).


#Import das bibliotecas python
import numpy as np

#Import das minhas funções
from funcoes import distancia
#Import das funções associadas ao manipulador usado
from pioneer_7dof import *

def GradDesc(posicaod,q,erro_min,Kmax):

    #variaveis do método
    qmax = 0.1 #passo maximo entre atualizacoes das juntas

    #valor maximo que a junta pode assumir
    qlim = getLimits()

    K = Kmax #número máximo de iterações

    for k in range(K):

        J,pn_0 = jacobianoGeometrico(q)
        
        #Calcula a distancia entre o efetuador a o objetivo(Posição)
        erro = distancia(pn_0,posicaod,3)

        #Condição de parada
        if(erro < erro_min):
            break  

        #erro
        f =  posicaod - pn_0[0:3] 
        
        #cálculo da constante de passo c
        e = np.array([f[0,0],f[1,0],f[2,0]])
        c = e.dot(J@J.T@e)/((J@J.T@e).dot(J@J.T@e))#peguei essa equacao do
        #artigo Inverse Kinematics Techniques in Computer Graphics: A Survey

        #Equação do Jacobiano Transposto
        dq =  c*(J.T@f)
                
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

    return [erro, k+1]

