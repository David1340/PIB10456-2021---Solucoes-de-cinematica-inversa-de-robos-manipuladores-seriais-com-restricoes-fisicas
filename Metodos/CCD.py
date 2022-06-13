#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do Cyclic Coordinate Descent para encontrar uma configuação q 
#dada uma posição (x,y,z) no espaço para o Pioneer 7DOF.
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import das bibliotecas python
from math import sqrt, pi, acos
import numpy as np

#Import das minhaa funções criadas 
from funcoes import distancia, matriz_homogenea_final
from funcoes import projecao_ponto_plano, rotationar_vetor
#from pioneer_7dof import *
from manipulador_15dof import *

#Calcula o acos de um angulo arrendondado em 10 casas decimais
def acosr(x):
    return acos(np.round(x,10))

#Calcula a norma de um vetor 
def norm(v):
    return sqrt(v[[0]]**2 + v[[1]]**2 + v[[2]]**2)

#criar um vetor coluna a partir de uma lista
def vetor(v):  
    return np.array([[v[0],v[1],v[2]]]).T


def CCD(posicaod,q,erro_min,Kmax,Solucao = False):
 
    qlim = getLimits() #valor maximo que a junta pode assumir
    n = getNumberJoints() #número de juntas
    alphasNegatives = getAlphasNegatives()

    for k in range(Kmax):
        # parametros de DH e ponto de atuação
        d,a,alpha,theta,p_n = getDH_paramaters(q)

        #Vy são os vetores de atuação das juntas
        Tn,Vy = matriz_homogenea_final(d,a,alpha,theta,'CCD')
        #pontos são as coordenadas das juntas 
        pontos = Cinematica_Direta(q) 
        erro = distancia(pontos[:,n],posicaod,3)
        for i in range(n-1,-1,-1):
            pontos = Cinematica_Direta(q) 
            proj = projecao_ponto_plano(vetor(Vy[:,i]),pontos[:,i],posicaod[:])
            va = proj - vetor(pontos[:,i])
            va = va/norm(va)
            proj = projecao_ponto_plano(vetor(Vy[:,i]),pontos[:,i],vetor(pontos[:,n]))
            vb = proj - vetor(pontos[:,i])
            vb = vb/norm(vb)
            th = acosr(va.T@vb)
            if(alphasNegatives[i]):#Se for a junta 2 ou 4, aparentemente as juntas que alpha < 0
                v = rotationar_vetor(va,vetor(Vy[:,i]),pi/2)
            else:
                v = rotationar_vetor(va,vetor(Vy[:,i]),-pi/2) 

            if(vb.T@v < 0): th = -th
            q[i,0] = q[i,0] + th
            
            if(q[i] > qlim[i]):
                q[i] = qlim[i]
            elif(q[i] < -qlim[i]):
                q[i] = -qlim[i]

            pontos = Cinematica_Direta(q)
            erro = distancia(pontos[:,n],posicaod,3) 

        pontos = Cinematica_Direta(q)
        erro = distancia(pontos[:,n],posicaod,3) 

        if(erro < erro_min):
            break

    if(Solucao):
        return [erro,k+1,pontos[:,n]]        
    return [erro,k+1]

