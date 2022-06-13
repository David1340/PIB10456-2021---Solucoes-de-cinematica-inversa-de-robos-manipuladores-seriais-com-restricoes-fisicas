#Autor David Oliveira.
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS.
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS.
#Experimentos apenas a posição .
#Implementações feitas durante durante a iniciação científica intitulada:
#PIB10456-2021 - Soluções de cinemática inversa de robôs manipuladores seriais com restrições físicas
#Durante o período: PIBIC 2021/2022 (01/09/2021 a 31/08/2022).

#Import das bibliotecas
import numpy as np
import random 
import sys
import os
import platform

#Import diretorio e funções
diretorio_atual = os.getcwd()
sys.path.append(diretorio_atual)
if(platform.system == 'Windows'):
    sys.path.append(diretorio_atual + '\Metodos')
else:
    sys.path.append(diretorio_atual + '/Metodos')


from DLS import DLS
from CCD import CCD
from DLS_WLS import DLS_WLS
from DLS_Completo import DLS_Completo
from DLS_WLS_Completo import DLS_WLS_Completo
from GradDesc import GradDesc
from GradDesc_Completo import GradDesc_Completo
from pioneer_7dof import *

#Configurações do experimento
Kmax = 1000
erro_min = 0.001
repeticoes = 1000

#valor maximo que a junta pode assumir
qlim = getLimits() 
n = getNumberJoints()
q = np.zeros([n,1])
kk =[]
for i in range(repeticoes):
    print('i:',i)
    #Gerando a cnfiguração inicial
    for i2 in range(np.size(q)):
        q[i2] = random.uniform(-qlim[i2],qlim[i2])

    [posicaod,orientacaod] = random_pose()

    [erro,k] = GradDesc_Completo(posicaod,orientacaod,q.copy(),erro_min,Kmax)
    if(erro < erro_min):
        kk.append(k)

mi = np.mean(kk)
tc = len(kk)
print('mi: ', mi)
print('tc', tc/repeticoes *100)
