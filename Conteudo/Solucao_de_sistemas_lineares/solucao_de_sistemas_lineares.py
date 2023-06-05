import random

import numpy as np
import pandas as pd


def gerarMatriz(dimensao: int, densidade_matriz: float):
    # Encontra a quantidade de valores não-nulos na matriz
    num_nao_nulo = int((dimensao**2) * densidade_matriz)

    # Gerar indices aleatorios
    linha_indices = [random.randint(0, dimensao-1) for _ in range(num_nao_nulo)]

    coluna_indices = [random.randint(0, dimensao-1) for _ in range(num_nao_nulo)]

    # Gera valores aleatorios para os valores não nulos
    valores = [random.randint(1, 10) for _ in range(num_nao_nulo)]

    # Gera a matriz esparsa
    esparsa_matriz = [[0 for _ in range(dimensao)] for _ in range(dimensao)]

    for i in range(num_nao_nulo):
        linha = linha_indices[i]
        coluna = coluna_indices[i]
        val = valores[i]
        esparsa_matriz[linha][coluna] = val

    return esparsa_matriz

def solucaoSuperior(data):
    coluna = len(data[0])-1

    resposta = []
    for linha in range(len(data)-1, -1, -1):
        if (linha == len(data)-1):

            aux = data[linha][coluna]/data[linha][linha]
            resposta.append(aux)

        else:

            auxVaAnterior = 0
            tempsoma = 0
            for soma in range(coluna-1, linha, -1):

                tempsoma += data[linha][soma]*resposta[auxVaAnterior]
                auxVaAnterior += 1

            aux = (data[linha][coluna] - tempsoma)/data[linha][linha]
            resposta.append(aux)
    return resposta

def solucaoInferior(data):
    coluna = len(data[0])-1

    resposta = []
    for linha in range(len(data)):
        if (linha == 0):
            aux = data[linha][coluna]/data[linha][linha]
            resposta.append(aux)

        else:

            auxVaAnterior = 0
            tempsoma = 0
            for soma in range(linha):

                tempsoma += data[linha][soma]*resposta[auxVaAnterior]
                auxVaAnterior += 1

            aux = (data[linha][coluna] - tempsoma)/data[linha][linha]
            resposta.append(aux)
    return resposta

def trianguloSuperior(dimensao, data):
    aux = 0
    aux2 = 0
    for i in range(dimensao):
        aux += i

    for linha in range(dimensao):
        for coluna in range(dimensao):
            if (linha == coluna):
                break
            elif (data[linha][coluna] == 0):
                aux2 += 1

    return True if aux == aux2 else False

def maiorPivo(data, coluna):
    maior = [0, 0]
    maior[0] = data[coluna][coluna]
    maior[1] = coluna

    for linha in range(coluna, len(data)):
        if data[linha][coluna] > maior[0]:
            maior[0] = data[linha][coluna]
            maior[1] = linha
    return maior[1]

def pivoteamento(iteracao, data, fatores=[]):

    auxpivo = iteracao

    if (data[iteracao][iteracao] != 0):
        pivo = data[iteracao][iteracao]
        aux = 1

        for linha in range(iteracao+1, len(data)):
            fator = data[linha][iteracao]/pivo
            fatores.append(fator)
            for coluna in range(len(data[0])):
                data[linha][coluna] = round((data[linha][coluna] - fator * data[linha-aux][coluna]), 4)

            aux += 1
        return 1
    else:
        while (data[iteracao][iteracao] == 0 and auxpivo != len(data)):
            auxpivo += 1
            temp = data[iteracao].copy()
            data[iteracao] = data[auxpivo].copy()
            data[auxpivo] = temp

        return 0

def pivoteamentoParcial(iteracao, data, fatores=[]):

    maiorlinha = maiorPivo(data, iteracao)

    if (data[iteracao][iteracao] != 0 and iteracao == maiorlinha):

        pivo = data[iteracao][iteracao]
        aux = 1

        for linha in range(iteracao+1, len(data)):
            fator = data[linha][iteracao]/pivo
            fatores.append(fator)
            for coluna in range(len(data[0])):
                data[linha][coluna] = round((data[linha][coluna] - fator * data[linha-aux][coluna]), 4)

            aux += 1
        return 1

    else:

        temp = data[iteracao].copy()
        data[iteracao] = data[maiorlinha].copy()
        data[maiorlinha] = temp

        return 0

def eliminacaoGaussianaParcial(data):
    iteracao = 0

    while (iteracao < len(data)):
        if trianguloSuperior(len(data), data) is False:
            iteracao += pivoteamentoParcial(iteracao, data)

        else:
            break

    # 3. Substituindo e resolvendo a equação
    resposta = solucaoSuperior(data)

    dataf = pd.DataFrame(data)
    return dataf, iteracao, resposta

def eliminacaoGaussiana(data):
    iteracao = 0

    while (iteracao < len(data)):
        if trianguloSuperior(len(data), data) is False:
            iteracao += pivoteamento(iteracao, data)

        else:
            break
    # 3. Substituindo e resolvendo a equação
    resposta = solucaoSuperior(data)
    dataf = pd.DataFrame(data)
    return dataf, iteracao, resposta

def menoresPrincipais(data):

    for iteracao in range(len(data)):
        auxM = []
        for linha in range(iteracao+1):
            auxL = []
            for coluna in range(iteracao+1):
                auxL.append(data[linha][coluna])
            # print("Lista:",auxL)
            auxM.append(auxL)

        if np.linalg.det(auxM) == 0:
            return False

    return True

def identidade(dimensao):
    auxm = []
    for linha in range(dimensao):
        auxl = []
        for coluna in range(dimensao):
            if (linha == coluna):
                auxl.append(1)
            else:
                auxl.append(0)
        auxm.append(auxl)
    return auxm

def metodoLU(data):
    if menoresPrincipais(data) is True:
        U = data.copy()
        linha = len(data)
        L = identidade(linha)
        for j in range(linha-1):
            for i in range(j+1, linha):
                L[i][j] = U[i][j]/U[j][j]
                for k in range(j+1, linha):
                    U[i][k] = U[i][k] - L[i][j]*U[j][k]
                U[i][j] = 0
        return L, U

    else:
        print("A matriz não tem decomposição LU.")

def solucaoLU(data, b):

    L, U = metodoLU(data)

    # 1. solucao Ly = b
    auxy = L.copy()
    for i in range(len(data)):
        auxy[i].append(b[i])

    y = solucaoInferior(auxy)

    # 1. Ux = y:
    auxx = U.copy()
    for i in range(len(data)):
        auxx[i].append(y[i])

    x = solucaoSuperior(auxx)

    return y, x, pd.DataFrame(L), pd.DataFrame(U)

def metodoJacobiano(data, chute, erro = 0.05):
    d = []
    B = []
    #1. escrevemos as matrizes B e d do processo iterativo:
    for linha in range(len(data)):
        pivo = data[linha][linha]
        auxl = []
        for coluna in range(len(data)+1):


            if(coluna == linha):
                data[linha][coluna] = 0
                auxl.append(0)
            elif(coluna == len(data[0])-1):
                data[linha][coluna] /= pivo
                d.append(data[linha][coluna] )
            else:
                data[linha][coluna] = (data[linha][coluna]/pivo)*(-1)
                auxl.append(data[linha][coluna])
        B.append(auxl)
                
    
    x = chute
    x0 = chute
    #2. daremos um chute inicial
    while True:
        x = np.dot(B,x) + d
        print(x)
        erroR = max(abs(x-x0))/max(abs(x))
        print("Erro relativo:",erroR)
        if(erroR > erro): break
        x0 = x.copy()
        
    return data

def metodoGaussSeidel(data, chute, erro = 0.05, iteracao = 5):
    d = []
    B = []
    #1. escrevemos as matrizes B e d do processo iterativo:
    for linha in range(len(data)):
        pivo = data[linha][linha]
        auxl = []
        for coluna in range(len(data)+1):


            if(coluna == linha):
                data[linha][coluna] = 0
                auxl.append(0)
            elif(coluna == len(data[0])-1):
                data[linha][coluna] /= pivo
                d.append(data[linha][coluna] )
            else:
                data[linha][coluna] = (data[linha][coluna]/pivo)*(-1)
                auxl.append(data[linha][coluna])
        B.append(auxl)
                

    x = chute
    x0 = chute.copy()
    
    while True:
        for i in range(len(B[0])):
            x[i] = np.dot(B[i],x) + d[i]

        print(x)
        erroR = max(abs(np.array(x)-np.array(x0))) / max(abs(np.array(x)))
        print("Erro relativo:",erroR)
        if(erroR > erro): break
        x0 = x.copy()
    return data