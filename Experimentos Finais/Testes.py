#Autor David Oliveira
#Estudante de Engenharia Eletrônica da Universidade Federal de Sergipe-UFS
#Membro do Grupo de Pesquisa em Robotica da UFS-GPRUFS
#Experimentos apenas a posição 

#Import das bibliotecas
import numpy as np
import random 
import sys

#Import diretorio e funções
sys.path.append('C:\PIBIC 2022 - Python')
sys.path.append('C:\PIBIC 2022 - Python\Experimentos Finais\Metodos')

from funcoes import random_pose
from FABRIK_Completo import FABRIK_Completo
#Configurações do experimento
Kmax = 50
erro_min = 0.001
repeticoes = 100
tentativas = 5
#valor maximo que a junta pode assumir
qlim = [2.6179,1.5358,2.6179,1.6144,2.6179,1.8413,1.7889] 
q = np.zeros([7,1])
ks = []
tc = []
mi = []

for i in range(repeticoes):
    print('i:',i)

    k2 = 0  
    [posicaod,orientacaod] = random_pose()
    for i2  in range(tentativas):
        #Gerando a cnfiguração inicial
        for i3 in range(np.size(q)):
            q[i3] = random.uniform(-qlim[i3],qlim[i3])
        [erro,k] = FABRIK_Completo(posicaod,orientacaod,q.copy(),erro_min,Kmax)
        k2 = k + k2
        if(erro < erro_min):
            ks.append(k2)
            break
    

tc.append( (len(ks)/repeticoes) * 100)
mi.append( np.mean(ks))

print(tc)
print(np.round(mi,2))
metodos = ["Undefined,"]
arquivo = open("ExperimentoC.txt", "w")
arquivo.write("Metodos: " + str(metodos) + "\n")
arquivo.write("tc: " + str(tc) + "\n")
arquivo.write("mi: " + str(np.round(mi,2)))