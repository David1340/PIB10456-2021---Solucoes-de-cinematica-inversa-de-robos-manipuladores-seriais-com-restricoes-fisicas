#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Implementação do parâmetros de DH do manipulador pioneer7DOF e algumas funções auxiliares
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

from  math import pi
import numpy as np
from funcoes import matriz_homogenea,S,Esfera
from random import uniform
import matplotlib.pyplot as plt

def getDH_paramaters(q):
    base = 0.05
    d = [0.075 + base,0,0.15,0,0.145,0,0]
    a = [0,0,0,0,0,0.075,0]
    alpha = [pi/2,-pi/2,pi/2,-pi/2,pi/2,pi/2,pi/2]
    theta = [pi/2 + q[0],q[1],q[2],q[3],q[4],pi/2 + q[5],pi/2 + q[6]]
    #distância da ultima junta a extremidade do efetuador
    L = 0.075 
    #ponto de atuação do manipulador no sistema de coordenadas on xn yn zn
    p_n = np.array([[0,0,L,1]]).T 
    return d,a,alpha,theta,p_n.copy()

def getRaio():
    return 0.025

def getLimits():
    qlim = [2.6179,1.5358,2.6179,1.6144,2.6179,1.8413,1.7889]
    return qlim

def getNumberJoints():
    n = 7 #número de juntas
    return n

def getAlphasNegatives():
    return [False,True,False,True,False,False,False] 

#Calcula a posição das juntas a partir da configuração q
def Cinematica_Direta(q,orientacao = False):
    #Pontos de interesse
    p = np.array([[0,0,0,1]]).T #Base
    p1_1 = np.array([[0,-0.05,0,1]]).T #junta1
    p2_2 = p #junta2
    p3_3 = np.array([[0,-0.075,0,1]]).T #junta3
    p4_4 = p #junta4
    p5_5 = np.array([[0,-0.0725,0,1]]).T #junta5
    p6_6 = np.array([[-0.075,0,0,1]]).T #junta6
    p7_7 = p #junta7

    d,a,alpha,theta,p8_7 = getDH_paramaters(q)

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    p1_0 = T1@p1_1
    p2_0 = T2@p2_2
    p3_0 = T3@p3_3
    p4_0 = T4@p4_4
    p5_0 = T5@p5_5
    p6_0 = T6@p6_6
    p7_0 = T7@p7_7
    p8_0 = T7@p8_7
    pontos = np.array([p1_0[0:3,0],p2_0[0:3,0],p3_0[0:3,0],p4_0[0:3,0]\
                    ,p5_0[0:3,0],p6_0[0:3,0],p7_0[0:3,0],p8_0[0:3,0]]).T
    
    if(orientacao):
        return pontos,T7[0:3,0:3]
    
    return pontos

#Calcula as posições das juntas e seus eixos de atuação
def Cinematica_Direta2(q):
    #Pontos de interesse
    p = np.array([[0,0,0,1]]).T #Base
    p1_1 = np.array([[0,-0.05,0,1]]).T #junta1
    p2_2 = p #junta2
    p3_3 = np.array([[0,-0.075,0,1]]).T #junta3
    p4_4 = p #junta4
    p5_5 = np.array([[0,-0.0725,0,1]]).T #junta5
    p6_6 = np.array([[-0.075,0,0,1]]).T #junta6
    p7_7 = p #junta7

    d,a,alpha,theta,p8_7 = getDH_paramaters(q)

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    p1_0 = T1@p1_1
    p2_0 = T2@p2_2
    p3_0 = T3@p3_3
    p4_0 = T4@p4_4
    p5_0 = T5@p5_5
    p6_0 = T6@p6_6
    p7_0 = T7@p7_7
    p8_0 = T7@p8_7
    pontos = np.array([p1_0[0:3,0],p2_0[0:3,0],p3_0[0:3,0],p4_0[0:3,0]\
                    ,p5_0[0:3,0],p6_0[0:3,0],p7_0[0:3,0],p8_0[0:3,0]]).T
    vetores = np.array([T1[0:3,1],T2[0:3,1],T3[0:3,1],T4[0:3,1],T5[0:3,1],T6[0:3,1]\
                        ,T7[0:3,1]]).T
    return [pontos,vetores]

#Calcula a posição e a matriz de rotação do efetuador final
def Cinematica_Direta3(q):
    #Pontos de interesse
    d,a,alpha,theta,p8_7 = getDH_paramaters(q)

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    p8_0 = T7@p8_7
    return [p8_0[0:3] ,T7[0:3,0:3]]

#Gera uma pose alcançável posição + matriz de rotação
def random_pose(): 
    #valor maximo que a junta pode assumir
    qlim = getLimits() 
    n = getNumberJoints()
    #angulos de juntas iniciais
    q = np.zeros([n,1])
    for a in range(np.size(q)):
        q[a] = uniform(-qlim[a],qlim[a])

    return Cinematica_Direta3(q)

#Calcula o jacobiano considerando apenas a posição
def jacobianoGeometrico(q):

    #vetores colunas do sistema de coordenadas global
    z = np.array([[0,0,1,1]]).T
    o = np.array([[0,0,0,1]]).T #origem

    #Parâmetros de DH e ponto de atuação
    d,a,alpha,theta,p8_7 = getDH_paramaters(q)

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
    p8_0 = T7@p8_7   

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

    return J,p8_0

#Calcula o jacobiano considerando apenas a posição e aorientação
def jacobianoGeometrico2(q):

    #vetores colunas do sistema de coordenadas global
    z = np.array([[0,0,1,1]]).T
    o = np.array([[0,0,0,1]]).T #origem

    #Parâmetros de DH e ponto de atuação
    d,a,alpha,theta,p8_7 = getDH_paramaters(q)

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
    p8_0 = T7@p8_7   

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
    J = np.zeros([6,7])
    #produto vetorial de Z0_0 por (o7_0 - o) 
    J[0:3,0] = S(z0_0)@(o7_0[0:3] - o[0:3])[:,0]
    J[3:6,0] = z0_0[:,0]
    J[0:3,1] = S(z1_0)@(o7_0[0:3] - o1_0[0:3])[:,0]
    J[3:6,1] = z1_0[:,0]
    J[0:3,2] = S(z2_0)@(o7_0[0:3] - o2_0[0:3])[:,0]
    J[3:6,2] = z2_0[:,0]
    J[0:3,3] = S(z3_0)@(o7_0[0:3] - o3_0[0:3])[:,0]
    J[3:6,3] = z3_0[:,0]
    J[0:3,4] = S(z4_0)@(o7_0[0:3] - o4_0[0:3])[:,0]
    J[3:6,4] = z4_0[:,0]
    J[0:3,5] = S(z5_0)@(o7_0[0:3] - o5_0[0:3])[:,0]
    J[3:6,5] = z5_0[:,0]
    J[0:3,6] = S(z6_0)@(o7_0[0:3] - o6_0[0:3])[:,0]
    J[3:6,6] = z6_0[:,0]

    return J,p8_0,T7

def getLengthElos():
    return  np.array([0.05,0.075,0.075,0.0725,0.0725,0.075,0.075]) 

def getTypeJoints():
    return ['pivot','hinge','pivot','hinge','pivot','hinge','hinge']

def plot_junta_revolucao(A,p,c,ax,h,r,cor = 'blue', offset = 0):
    #A matriz de Rotação, p origem da junta no seu sistema de coordenadas, c eixo do

    theta = np.arange(0,2*pi + 0.12,0.8) #theta de 0 a 2pi com passos de 0.1
    if(c == 'z'):  
        z = np.linspace(0,h,np.size(theta)) + offset#z de 0.1 a 0.1 com o numero de elementos iguais ao de theta
        z, theta = np.meshgrid(z, theta) #transforma em vetores 2D por causa do plot de superficie
        #[x,y,z].T = A@([x,y,z].T + p)
        x = (r*np.cos(theta) + p[0,0])*A[0,0] + (r*np.sin(theta) + p[1,0])*A[0,1] + \
         (z + p[2,0])*A[0,2] + np.ones_like(z)*A[0,3]
        y = (r*np.cos(theta) + p[0,0])*A[1,0] + (r*np.sin(theta) + p[1,0])*A[1,1] + \
         (z + p[2,0])*A[1,2] + np.ones_like(z)*A[1,3]
        z = (r*np.cos(theta) + p[0,0])*A[2,0] + (r*np.sin(theta) + p[1,0])*A[2,1] + \
         (z + p[2,0])*A[2,2] + np.ones_like(z)*A[2,3]
        ax.plot_surface(x, y, z,color = cor, alpha = 1)
        
    elif(c == 'y'):
        y = np.linspace(-h,h,np.size(theta)) + offset#z de 0.1 a 0.1 com o numero de elementos iguais ao de theta
        y, theta = np.meshgrid(y, theta) #transforma em vetores 2D por causa do plot de superficie
        #[x,y,z].T = A@([x,y,z].T + p)
        x = (r*np.cos(theta) + p[0,0])*A[0,0] + (r*np.sin(theta) + p[2,0])*A[0,2] + \
         (y + p[1,0])*A[0,1] + np.ones_like(y)*A[0,3]
        z = (r*np.cos(theta) + p[0,0])*A[2,0] + (r*np.sin(theta) + p[2,0])*A[2,2] + \
         (y + p[1,0])*A[2,1] + np.ones_like(y)*A[2,3]
        y = (r*np.cos(theta) + p[0,0])*A[1,0] + (r*np.sin(theta) + p[2,0])*A[1,2] + \
         (y + p[1,0])*A[1,1] + np.ones_like(y)*A[1,3]
        ax.plot_surface(x, y, z,color = cor, alpha = 1)
        
    elif(c == 'x'):
        x = np.linspace(-h,h,np.size(theta)) + offset#z de 0.1 a 0.1 com o numero de elementos iguais ao de theta
        x, theta = np.meshgrid(x, theta) #transforma em vetores 2D por causa do plot de superficie
        #[x,y,z].T = A@([x,y,z].T + p)
        y = (r*np.cos(theta) + p[2,0])*A[1,2] + (r*np.sin(theta) + p[1,0])*A[1,1] + \
         (x + p[0,0])*A[1,0] + np.ones_like(x)*A[1,3]
        z = (r*np.cos(theta) + p[2,0])*A[2,2] + (r*np.sin(theta) + p[1,0])*A[2,1] + \
         (x + p[0,0])*A[2,0] + np.ones_like(x)*A[2,3]
        x = (r*np.cos(theta) + p[2,0])*A[0,2] + (r*np.sin(theta) + p[1,0])*A[0,1] + \
         (x + p[0,0])*A[0,0] + np.ones_like(x)*A[0,3]
        ax.plot_surface(x, y, z,color = cor, alpha = 1)

#Função que plota o manipulador, quando recebe os pontos de interesse
def plot(q,t,esferas):
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    #Pontos de interesse
    p = np.array([[0,0,0,1]]).T #Base
    p1_1 = np.array([[0,-0.05,0,1]]).T #junta1
    p2_2 = p #junta2
    p3_3 = np.array([[0,-0.075,0,1]]).T #junta3
    p4_4 = p #junta4
    p5_5 = np.array([[0,-0.0725,0,1]]).T #junta5
    p6_6 = np.array([[-0.075,0,0,1]]).T #junta6
    p7_7 = p #junta7

    w = 0.05 #largura do efetuador
    h = 0.025 #altura do efetuador
    L = 0.075 #comprimento do elo que liga a junta 7 ao efetuador

    #O efetuador será um u, para plotarmos um precisamos de 4 pontos
    #topo esquerdo do u
    
    e1_7 = np.array([[0,-w/2,L + h,1]]).T
    #base esquerda do u
    e2_7 = np.array([[0,-w/2,L,1]]).T
    #base direita do u
    e3_7 = np.array([[0,w/2,L,1]]).T
    #topo direito do u
    e4_7 = np.array([[0,w/2,L+ h,1]]).T 

    d,a,alpha,theta,p8_7 = getDH_paramaters(q)

    #Calculando as matrizes homogêneas
    A1 = matriz_homogenea(d[0],a[0],alpha[0],theta[0])
    A2 = matriz_homogenea(d[1],a[1],alpha[1],theta[1])
    A3 = matriz_homogenea(d[2],a[2],alpha[2],theta[2])
    A4 = matriz_homogenea(d[3],a[3],alpha[3],theta[3])
    A5 = matriz_homogenea(d[4],a[4],alpha[4],theta[4])
    A6 = matriz_homogenea(d[5],a[5],alpha[5],theta[5])
    A7 = matriz_homogenea(d[6],a[6],alpha[6],theta[6])

    #Calculando os pontos de interesse no sistema Global
    T1 = A1
    T2 = T1@A2
    T3 = T2@A3
    T4 = T3@A4
    T5 = T4@A5
    T6 = T5@A6
    T7 = T6@A7
    e1_0 = T7@e1_7
    e2_0 = T7@e2_7
    e3_0 = T7@e3_7
    e4_0 = T7@e4_7   
   
    #Plotando Elos 
    ax.clear()   
   
    #Plotando efetuador
    plt.plot([e1_0[0,0],e2_0[0,0],e3_0[0,0],e4_0[0,0],e3_0[0,0]]\
             ,[e1_0[1,0],e2_0[1,0],e3_0[1,0],e4_0[1,0],e3_0[1,0]]\
             ,[e1_0[2,0],e2_0[2,0],e3_0[2,0],e4_0[2,0],e3_0[2,0]],'green')
    
    #Plotando objetivo    
    ax.scatter(t[0,0],t[1,0],t[2,0],color = 'red')
    
    #Plotando Juntas e base   
    r = getRaio()
    plot_junta_revolucao(np.eye(4),p,'z',ax,0.025,r,'green')
    plot_junta_revolucao(T1,p1_1,'y',ax,0.05,r,'green',-0.025*0.5)
    plot_junta_revolucao(T2,p2_2,'y',ax,0.025,r,'green')
    plot_junta_revolucao(T3,p3_3,'y',ax,0.05,r,'green')
    plot_junta_revolucao(T4,p4_4,'y',ax,0.025,r,'green')
    plot_junta_revolucao(T5,p5_5,'y',ax,0.05,r,'green')
    plot_junta_revolucao(T6,p6_6,'y',ax,0.025,r,'green')
    plot_junta_revolucao(T6,p6_6,'x',ax,0.5*0.025,r,'green',1.5*0.025)
    plot_junta_revolucao(T7,p7_7,'y',ax,0.025,r,'green')
    plot_junta_revolucao(T7,p7_7,'z',ax,0.04,r,'green',0.025)

    for esfera in esferas:
        plot_esfera(esfera,ax)

    ax.set_xlabel('Eixo x(m)')
    ax.set_ylabel('Eixo y(m)')
    ax.set_zlabel('Eixo z(m)')

    #titulo
    plt.title('Pioneer 7DOF')
    ax.set_xlim3d(-0.5,0.5)
    ax.set_ylim3d(-0.5,0.5)
    ax.set_zlim3d(0,1)
    fig.canvas.draw() #mostra o plot
    plt.pause(15)

def plot_esfera(esfera: Esfera,ax):
    x0 = esfera.x
    y0 = esfera.y
    z0 = esfera.z
    r = esfera.r

    u = np.linspace(0,2*np.pi,10)
    v = np.linspace(0,2*np.pi,10)
    u,v = np.meshgrid(u,v)
    x = r*np.cos(u)*np.cos(v) + x0
    y = r*np.cos(u)*np.sin(v) + y0
    z = r*np.sin(u) + z0
    ax.plot_surface(x,y,z,color = 'blue', alpha = 1)
