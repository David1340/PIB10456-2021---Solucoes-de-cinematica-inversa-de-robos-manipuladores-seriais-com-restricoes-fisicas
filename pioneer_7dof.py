from  math import pi
import numpy as np

def getDH_paramaters(q):
    base = 0.05
    d = [0.075 + base,0,0.15,0,0.145,0,0]
    a = [0,0,0,0,0,0.075,0]
    alpha = [pi/2,-pi/2,pi/2,-pi/2,pi/2,pi/2,pi/2]
    theta = [pi/2 + q[0],q[1],q[2],q[3],q[4],pi/2 + q[5],pi/2 + q[6]]
    L = 0.075 #distância da ultima junta a extremidade do efetuador
    p_7 = np.array([[0,0,L,1]]).T #ponto de atuação do manipulador no sistema de coordenadas o7x7y7z7
    return d,a,alpha,theta,p_7.copy()

def getLimits():
    qlim = [2.6179,1.5358,2.6179,1.6144,2.6179,1.8413,1.7889]
    return qlim