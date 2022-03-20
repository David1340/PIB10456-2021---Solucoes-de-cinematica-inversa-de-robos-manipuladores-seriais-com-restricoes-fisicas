#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do Forward And Backward Reaching Inverse Kinematic 
#para encontrar uma configuação q dada uma posição (x,y,z) no espaço para o Pioneer 7DOF
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import do modulo funcoes.py
import sys
sys.path.append('C:\PIBIC 2022 - Python')
from funcoes import Cinematica_Direta2, norm, distancia, random_pose, S \
                    , rotationar_vetor, projecao_ponto_plano

from math import pi,acos
import numpy as np

#criar um vetor coluna a partir de uma lista
def vetor(v):   
    return np.array([[v[0],v[1],v[2]]]).T

#Interação básica do Fabrik restringida no plano
def iteracao_Fabrik(p,p0,d,normal):
    #normal -> vetor normal ao plano
    #p0 -> ponto pertecente ao plano
    #p -> ponto a ser projetado
    #d -> distância que o ponto que será projetado terá de p0
    normal = vetor(normal)
    p0 = vetor(p0)
    p = vetor(p)
    proj = projecao_ponto_plano(normal,p0,p)
    r = distancia(p0,proj,3)
    delta = d/r
    pnew = (1- delta)*p0 + delta*proj
    return pnew

#Calcula o acos de um angulo arrendondado em 10 casas decimais
def acosr(x):
    return acos(np.round(x,10))

#Calcula um novo p1 que respeite o limite ângular da junta
#Serve apenas para junta do tipo Hinge
def restringe_junta(p1,p2,p3,i):

    qlim = [2.6179,1.5358,2.6179,1.6144,2.6179,1.8413,1.7889] 
    v1 = vetor(p1 - p2)
    v1 = v1/norm(v1)
    v2 = vetor(p3 - p2)
    v2 = v2/norm(v2)

    th = acosr(v1.T@v2)  

    d = distancia(p1,p2,3)

    v3 = S(v1)@v2 #vetor ortogonal a v1 e v2
    #Se V3 não for o vetor nulo, normalize V3
    if(norm(v3) > 0.0001): v3 = v3/norm(v3)
    #Se o ângulo é maior que o permitido sature ele em qlim[i] - 0.1
    if(th < pi - qlim[i]):
        v1 = rotationar_vetor(v1,v3,th-(pi - qlim[i] + 0.1))
        p1 = vetor(p2) + d*v1
    else:
        p1 = vetor(p1)
    return p1[:,0]

#Calcula um novo v1 que respeite o limite ângular da junta
#Serve apenas para junta do tipo Pivot
def restringe_junta2(v1,v2,i): 
    qlim = [2.6179,1.5358,2.6179,1.6144,2.6179,1.8413,1.7889] 
    th = acosr(v1.T@v2)
    v3 = S(v2)@v1 #vetor ortogonal a v1 e v2
    #Se o ângulo é maior que o permitido sature ele em qlim[i] - 0.01
    if(th > qlim[i]):
        v1 = rotationar_vetor(v2,v3,qlim[i]-0.01)
    return v1[:,0]/norm(v1)

def FABRIK(posicaod,q,erro_min,Kmax):
    n = 7 #numero de juntas
    x  = vetor([1,0,0])
    y  = vetor([0,1,0])
    z  = vetor([0,0,1])


    destino = posicaod

    #tamanho dos elos
    b = np.array([0.05,0.075,0.075,0.0725,0.0725,0.075,0.075])  

    #Calcula os pontos e os vetores de atuação do robô atuais
    p,D = Cinematica_Direta2(q)

    #No Backward a junta 1 é sempre iniciada no mesmo lugar 
    pcte = p[:,0].copy()

    #pl e Dl são p_linha e D_linha respectivamente
    pl = p.copy()
    Dl = D.copy()

    #Calculo o erro inicial (distância euclidiana)
    erro = distancia(p[:,n],destino,3)
    K = Kmax #número máximo de iterações
    k = 0 #iteração inicial
    erromin = erro_min #erro minimo usado como um dos critérios de parada

    while(erro > erromin and k < K):

        #Forward
        for i in range(n-1,0,-1):
            if(i == 6): #Se for a junta 7    
                pl[:,i+1] = destino[:,0]#coloca o efetuador no destino          
                v1 = vetor(pl[:,i+1] - p[:,i])#p6 -> p7' (efetuador)
                v1 = v1/norm(v1)
                v2 = vetor(p[:,i-1] - p[:,i])#p6 -> p5
                v2 = v2/norm(v2)
                naux = (S(v1)@v2)#produto vetorial            
                if(norm(naux) > 0.00001): #Se p7',p6 e p5 não forem colineares
                    Dl[:,i] = naux[:,0]/norm(naux)
                else: #Caso não seja mantém o vetor de direção da iteração anterior
                    Dl[:,i] = D[:,i].copy()

                pl[:,i] = iteracao_Fabrik(p[:,i],pl[:,i+1],b[i],Dl[:,i])[:,0] 
                

            elif(i == 1 or i == 3 or i == 5): #Se atual for junta Hinge (2,4 e 6)

                if(i == 1 or i == 3):#Se a junta prev for pivot (2 e 4)
                    pl[:,i] = pl[:,i+1] - Dl[:,i+1]*b[i]
                    #paux é pl é o ponto da próxima iteração
                    #eu calculo ele aqui porque como proposto na abordagem FABRIK-R
                    #Dl(:,i) é cálculo de forma que pl(:,i+1) seja alcancável
                    paux = iteracao_Fabrik(p[:,i-1],pl[:,i],b[i],Dl[:,i+1])[:,0]
                    paux = restringe_junta(paux,pl[0:3,i],pl[0:3,i+1],i+1)
                    v1 = vetor(paux - pl[:,i])                
                    v1 = v1/norm(v1) #vetor de refência
                    v2 = vetor(Dl[:,i+2])
                    v2 = v2/norm(v2)
                    th = np.real(acosr(v1.T@v2))

                    #v3 é um vetor ortogonal ao vetor de referência (v1)
                    v3 = rotationar_vetor(v1,vetor(Dl[:,i+1]),pi/2)[:,0]
                    if(v3.T@v2 < 0):
                        th = -th

                    Dl[:,i] = rotationar_vetor(v2,vetor(Dl[:,i+1]),(pi/2) - th)[:,0] 
                    Dl[:,i] = Dl[:,i]/norm(D[:,i])

                else: #Se a junta prev for hinge (6)
                    pl[:,i] = iteracao_Fabrik(p[:,i-1],pl[:,i+1],b[i],Dl[:,i+1])[:,0]
                    pl[:,i] = restringe_junta(pl[0:3,i],pl[0:3,i+1],pl[0:3,i+2],i+1)
                    v1 = vetor(pl[:,i] - pl[:,i+1])
                    v1 = v1/norm(v1) 
                    Dl[:,i] = rotationar_vetor(vetor(Dl[:,i+1]),v1,pi/2)[:,0] 
                    Dl[:,i] = Dl[:,i]/norm(D[:,i])
                    
            elif(i == 2 or i == 4): #Se a junta for pivot (3 e 5)
                pl[:,i] = iteracao_Fabrik(p[:,i-1],pl[:,i+1],b[i],Dl[:,i+1])[:,0]
                pl[:,i] = restringe_junta(pl[0:3,i],pl[0:3,i+1],pl[0:3,i+2],i+1)
                v1 = pl[:,i] - pl[:,i+1]
                v1 = v1/norm(v1)
                Dl[:,i] = -v1
                
    
        Dl[:,0] = [0,0,1] #O vetor de atuação da junta 1 no forward não muda
        pl[:,0] = pl[:,1] - b[0]*Dl[:,0] 
        #Atualiza os valores de p e D
        D = Dl.copy()
        p = pl.copy()

        #Backward
        for i in range(0,n):
            if(i == 1 or i == 3 or i == 5 or i == 6): #Se for junta Hinge
                if(i == 1 or i == 3 or i == 5):#Se a junta prev for pivot (2,4 e 6)
                    pl[:,i] = pl[:,i-1] + Dl[:,i-1]*b[i-1]
                    paux = iteracao_Fabrik(p[:,i+1],pl[:,i],b[i],Dl[:,i-1])[:,0]
                    v1 = vetor(paux - pl[:,i])
                    v1 = v1/norm(v1)
                    if(i != 1): 
                        v2 = vetor(Dl[:,i-2])
                    #Como a junta 2 é primeira pivot da cadeia
                    #eu uso como direção iniciao o vetor -i
                    else: v2 = np.array([[1,0,0]]).T 
                    th = np.real(acosr(v1.T@v2))

                    #v3 é um vetor ortogonal ao vetor de referência (v1)
                    v3 = rotationar_vetor(v1,vetor(Dl[:,i-1]),pi/2)[:,0]
                    if(v3.T@v2 < 0):
                        th = -th

                    Dl[:,i] = rotationar_vetor(v2,vetor(Dl[:,i-1]),(pi/2) - th)[:,0] 
                    Dl[:,i] = Dl[:,i]/norm(D[:,i])
                    if(i != 1): Dl[:,i]  = restringe_junta2(vetor(Dl[:,i]),vetor(Dl[:,i-2]),i-1) #1-3-5
                    else: Dl[:,i]  = restringe_junta2(vetor(Dl[:,i]),vetor([1,0,0]),i-1)
                else: #Se a junta prev for Hinge (7)
                    pl[:,i] = iteracao_Fabrik(p[:,i+1],pl[:,i-1],b[i-1],Dl[:,i-1])[:,0]
                    pl[:,i] = restringe_junta(pl[:,i],pl[0:3,i-1],pl[0:3,i-2],i-1)
                    paux = iteracao_Fabrik(p[:,i+1],pl[:,i],b[i],Dl[:,i-1])[:,0]
                    v1 = vetor(paux - pl[:,i])
                    v1 = v1/norm(v1)
                    Dl[:,i] = rotationar_vetor(vetor(Dl[:,i-1]),v1,pi/2)[:,0]
                    Dl[:,i] = Dl[:,i]/norm(D[:,i])
                    #efetuador
                    pl[:,7]  = iteracao_Fabrik(p[:,7],pl[:,6],b[6],Dl[:,6])[:,0] 
                    pl[:,7] = restringe_junta(pl[:,7],pl[0:3,6],pl[0:3,5],6) 

            elif(i == 2 or i == 4): #Se a junta for pivot (3 e 5)
                pl[:,i] = iteracao_Fabrik(p[:,i+1],pl[:,i-1],b[i-1],Dl[:,i-1])[:,0]
                pl[:,i] = restringe_junta(pl[0:3,i],pl[0:3,i-1],pl[0:3,i-2],i-1)
                v1 = pl[:,i] - pl[:,i-1]
                v1 = v1/norm(v1)
                Dl[:,i] = v1
            elif(i == 0): #Primeira junta eixo de atuação e posição são fixos
                pl[:,0] = pcte #[0,0,0.075]
            #END-Backward
        
        #Atualiza os valores de p e D
        p = pl.copy()
        D = Dl.copy()
        
        #cálculo da distância euclidiana
        erro = distancia(p[:,n],destino,3)
        k = k +1

    ####Conversão da solução gráfica em um vetor de ângulo####
    #Primeiro normalizo os vetores, porque as vezes por causa de erro numérico eles chegam
    #com norma maior do que 1
    for i in range(7):
        D[:,i] = D[:,i]/norm(D[:,i])
        
    #Calcula um vetor v que seja ortogonal ao vetor de refência para o cálculo dos ângulos 
    for i in range(7):
        if(i == 0):
            vref = x
        elif(i == 6):
            vref = vetor(p[:,6] - p[:,5])
            vref = vref/norm(vref)
        else:
            vref = vetor(D[:,i-1])

        #v é o vetor ortogonal ao vetor de refência para o cálculo dos ângulos 
        v = rotationar_vetor(vref,vetor(D[:,i]),pi/2) 
        
        #cálculo o ângulo
        if(i == 5):
            vaux = vetor(p[:,6] - p[:,5])
            vaux = vaux/norm(vaux)
            q[5] = acosr(vaux.T@vref)
            if(vaux.T@v < 0): q[5] = - q[5]
        elif(i == 6):
            vaux = vetor(p[:,7] - p[:,6])
            vaux = vaux/norm(vaux)
            q[6] = acosr(vaux.T@vref)
            if(vaux.T@v < 0): q[6] = - q[6]
        else:    
            q[i] = acosr(D[:,i+1].T@vref)
            if(D[:,i+1].T@v < 0): q[i] = - q[i]

    return [erro,k]

 

