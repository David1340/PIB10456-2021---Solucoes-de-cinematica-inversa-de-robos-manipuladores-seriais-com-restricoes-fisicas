#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Experimentos apenas a posição .
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import das bibliotecas
import numpy as np
from random import uniform
import sys
import os
import platform

#Import diretorio
diretorio_atual = os.getcwd()
sys.path.append(diretorio_atual)
if(platform.system == 'Windows'):
    sys.path.append(diretorio_atual + '\Metodos')
else:
    sys.path.append(diretorio_atual + '/Metodos')

#Import dos métodos
from DLS import DLS
from DLS_WLS import DLS_WLS
from GradDesc import GradDesc
from PSO import PSO
from CCD import CCD
from FABRIK import FABRIK
from FRPSO import FRPSO
from pioneer_7dof import *

#Configurações do experimento
Kmax = 1000
erro_min = 0.001
repeticoes = 1000

#parâmetros do manipulador
qlim = getLimits() 
n = getNumberJoints()
q = np.zeros([n,1])

kDLS = []
kDLS_WLS = []
kGradDesc = []
kPSO = []
kFRPSO = []
kCCD = []
kFABRIK = []
tc = [0,0,0,0,0,0,0]
mi = tc.copy()

for i in range(repeticoes):
    print('i:',i)
    
    #Gerando a configuração inicial
    for i2 in range(np.size(q)):
        q[i2] = uniform(-qlim[i2],qlim[i2])

    [posicaod,orientacaod] = random_pose()

    [erro,k] = DLS(posicaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kDLS.append(k)

    [erro,k] = DLS_WLS(posicaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kDLS_WLS.append(k)

    [erro,k] = GradDesc(posicaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kGradDesc.append(k)

    [erro,k] = PSO(posicaod,orientacaod,erro_min,Kmax)
    if(erro < erro_min):
        kPSO.append(k)

    [erro,k] = FRPSO(posicaod,orientacaod,erro_min,Kmax)
    if(erro < erro_min):
        kFRPSO.append(k)

    [erro,k] = CCD(posicaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kCCD.append(k)

    [erro,k] = FABRIK(posicaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kFABRIK.append(k)

tc[0] = (len(kDLS)/repeticoes) * 100
tc[1] = (len(kDLS_WLS)/repeticoes) * 100
tc[2] = (len(kGradDesc)/repeticoes) * 100
tc[3] = (len(kPSO)/repeticoes) * 100
tc[4] = (len(kFRPSO)/repeticoes) * 100
tc[5] = (len(kCCD)/repeticoes) * 100
tc[6] = (len(kFABRIK)/repeticoes) * 100

mi[0] = np.mean(kDLS)
mi[1] = np.mean(kDLS_WLS)
mi[2] = np.mean(kGradDesc)
mi[3] = np.mean(kPSO)
mi[4] = np.mean(kFRPSO)
mi[5] = np.mean(kCCD)
mi[6] = np.mean(kFABRIK)

print(tc)
print(np.round(mi,2))
metodos = ["DLS","DLS_WLS","GradDesc","PSO-P","FRPSO","CCD","FABRIK"]
arquivo = open("ExperimentoA.txt", "w")
arquivo.write("Metodos: " + str(metodos) + "\n")
arquivo.write("tc: " + str(np.round(tc,2)) + "\n")
arquivo.write("mi: " + str(np.round(mi,2))+ "\n")
arquivo.write("Kmax: " + str(Kmax))
arquivo.close()
