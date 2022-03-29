#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Experimentos envolvedo a pose completa.
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import das bibliotecas
import numpy as np
import random 
import sys
import os

#Import diretorio e funções
diretorio_atual = os.getcwd()
sys.path.append(diretorio_atual)
sys.path.append(diretorio_atual + '\Metodos')

from funcoes import random_pose
from pioneer_7dof import getLimits
from DLS_Completo import DLS_Completo
from DLS_WLS_Completo import DLS_WLS_Completo
from GradDesc_Completo import GradDesc_Completo
from PSO import PSO
from FRPSO import FRPSO
from FABRIK_Completo import FABRIK_Completo

#valor maximo que a junta pode assumir
qlim = getLimits()

Kmax = 1000
erro_min = 0.001
repeticoes = 1000

kDLS = []
kDLS_WLS = []
kGradDesc = []
kPSO = []
kFRPSO = []
kFABRIK = []
tc = [0,0,0,0,0,0]
mi = tc.copy()
q = np.zeros([7,1])

for i in range(repeticoes):
    print('i:',i)
    #Gerando a cofiguração inicial
    for i2 in range(np.size(q)):
        q[i2] = random.uniform(-qlim[i2],qlim[i2])
    
    [posicaod,orientacaod] = random_pose()

    [erro,k] = DLS_Completo(posicaod,orientacaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kDLS.append(k)

    [erro,k] = DLS_WLS_Completo(posicaod,orientacaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kDLS_WLS.append(k)

    [erro,k] = GradDesc_Completo(posicaod,orientacaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kGradDesc.append(k)

    [erro,k] = PSO(posicaod,orientacaod,erro_min,Kmax)
    if(erro < erro_min):
        kPSO.append(k)

    [erro,k] = FRPSO(posicaod,orientacaod,erro_min,Kmax)
    if(erro < erro_min):
        kFRPSO.append(k)

    [erro,k] = FABRIK_Completo(posicaod,orientacaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kFABRIK.append(k)

tc[0] = (len(kDLS)/repeticoes) * 100
tc[1] = (len(kDLS_WLS)/repeticoes) * 100
tc[2] = (len(kGradDesc)/repeticoes) * 100
tc[3] = (len(kPSO)/repeticoes) * 100
tc[4] = (len(kFRPSO)/repeticoes) * 100
tc[5] = (len(kFABRIK)/repeticoes) * 100

mi[0] = np.mean(kDLS)
mi[1] = np.mean(kDLS_WLS)
mi[2] = np.mean(kGradDesc)
mi[3] = np.mean(kPSO)
mi[4] = np.mean(kFRPSO)
mi[5] = np.mean(kFABRIK)

print(tc)
print(np.round(mi,2))
metodos = ["DLS","DLS_WLS","GradDesc","PSO-P","FRPSO","FABRIK"]
arquivo = open("ExperimentoB2.txt", "w")
arquivo.write("Metodos: " + str(metodos) + "\n")
arquivo.write("tc: " + str(tc) + "\n")
arquivo.write("mi: " + str(np.round(mi,2))+ "\n")
arquivo.write("Kmax: " + str(Kmax))