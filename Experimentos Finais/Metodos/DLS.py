#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do Damped Last Square 
#para encontrar encontrar uma configuração q
#dada uma posição (x,y,z) .
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import do modulo funcoes.py
import sys
sys.path.append('C:\PIBIC 2022 - Python')
from funcoes import matriz_homogenea, distancia, S

#Import das bibliotecas
import numpy as np
from pioneer_7dof import getDH_paramaters, getLimits

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
        
        #Parâmetros de DH e ponto de atuação
        d,a,alpha,theta,p_7 = getDH_paramaters(q)

        #Matrizes homogêneas
        A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
        A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
        A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
        A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
        A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
        A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
        A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

        #Definindo os pontos de interesse em seus sistemas locais
        o1_1 = o #origem do SC 1
        o2_2 = o #origem do SC 2
        o3_3 = o #origem do SC 3
        o4_4 = o #origem do SC 4
        o5_5 = o #origem do SC 5
        o6_6 = o #origem do SC 6
        o7_7 = o #origem do SC 7
        

        #Calculando os pontos de interesse no sistema Global
        T1 = A1
        T2 = T1@A2
        T3 = T2@A3
        T4 = T3@A4
        T5 = T4@A5
        T6 = T5@A6
        T7 = T6@A7
        o1_0 = T1@o1_1
        o2_0 = T2@o2_2
        o3_0 = T3@o3_3
        o4_0 = T4@o4_4
        o5_0 = T5@o5_5
        o6_0 = T6@o6_6
        o7_0 = T7@o7_7
        p_0 = T7@p_7      
        
        #Calcula a distancia entre o efetuador a o objetivo(Posição)
        erro = distancia(p_0,posicaod,3)

        #Condição de parada
        if(erro < erro_min):
            break  

        #os vetores z serao transformados em vetores  no R^3
        z0_0 = z[0:3]
        z1_0 = (T1@z)[0:3]
        z2_0 = (T2@z)[0:3]
        z3_0 = (T3@z)[0:3]
        z4_0 = (T4@z)[0:3]
        z5_0 = (T5@z)[0:3]
        z6_0 = (T6@z)[0:3]
        #z7_0 = (T7@z)[0:3] nao eh usado.

        #cálculo do Jacobiano geométrico
        J = np.zeros([3,7])
        #produto vetorial de Z0_0 por (o7_0 - o) 
        J[:,0] = S(z0_0)@(o7_0[0:3] - o[0:3])[:,0]
        J[:,1] = S(z1_0)@(o7_0[0:3] - o1_0[0:3])[:,0]
        J[:,2] = S(z2_0)@(o7_0[0:3] - o2_0[0:3])[:,0]
        J[:,3] = S(z3_0)@(o7_0[0:3] - o3_0[0:3])[:,0]
        J[:,4] = S(z4_0)@(o7_0[0:3] - o4_0[0:3])[:,0]
        J[:,5] = S(z5_0)@(o7_0[0:3] - o5_0[0:3])[:,0]
        J[:,6] = S(z6_0)@(o7_0[0:3] - o6_0[0:3])[:,0]

        #erro
        f = posicaod - p_0[0:3]
        
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