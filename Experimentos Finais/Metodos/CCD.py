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

from funcoes import matriz_homogenea, distancia, Cinematica_Direta
from funcoes import projecao_ponto_plano, rotationar_vetor
from pioneer_7dof import getDH_paramaters, getLimits




#Calcula o acos de um angulo arrendondado em 10 casas decimais
def acosr(x):
    return acos(np.round(x,10))

#Calcula a norma de um vetor 
def norm(v):
    return sqrt(v[[0]]**2 + v[[1]]**2 + v[[2]]**2)

#criar um vetor coluna a partir de uma lista
def vetor(v):  
    return np.array([[v[0],v[1],v[2]]]).T


def CCD(posicaod,q,erro_min,Kmax):

    #valor maximo que a junta pode assumir
    qlim = getLimits()  
    n = 7 #número de juntas

    for k in range(Kmax):
        # parametros de DH e ponto de atuação
        d,a,alpha,theta,p_7 = getDH_paramaters(q)

        #Matrizes homogêneas
        A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
        A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
        A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
        A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
        A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
        A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
        A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

        #Calculando os pontos de interesse e vetores no sistema Global
        T1 = A1
        T2 = T1@A2
        T3 = T2@A3
        T4 = T3@A4
        T5 = T4@A5
        T6 = T5@A6
        T7 = T6@A7

        #vetores de atuação das juntas
        Vy = np.array([T1[0:3,1],T2[0:3,1],T3[0:3,1],T4[0:3,1],T5[0:3,1],T6[0:3,1],T7[0:3,1]]).T

        pontos = Cinematica_Direta(q) 
        erro = distancia(pontos[:,7],posicaod,3) 

        for i in range(n-1,-1,-1):
            pontos = Cinematica_Direta(q) 
            proj = projecao_ponto_plano(vetor(Vy[:,i]),pontos[:,i],posicaod[:])
            va = proj - vetor(pontos[:,i])
            va = va/norm(va)
            proj = projecao_ponto_plano(vetor(Vy[:,i]),pontos[:,i],vetor(pontos[:,7]))
            vb = proj - vetor(pontos[:,i])
            vb = vb/norm(vb)
            th = acosr(va.T@vb)
            j = i + 1
            if(j == 4 or j == 2):#Se for a junta 2 ou 4
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
            erro_anterior = erro
            erro = distancia(pontos[:,7],posicaod,3) 

        pontos = Cinematica_Direta(q)
        erro = distancia(pontos[:,7],posicaod,3) 

        if(erro < erro_min):
            break

    return [erro,k+1]

