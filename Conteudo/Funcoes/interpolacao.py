def matriz_inter_quadratica(data):
    matriz = []
    for i in range(len(data.x)):
        linha=[]
        for j in range(len(data.x)):
            if (j==0):
                linha.append(1)
            else:   
                linha.append(data.x[i]**j)
        matriz.append(linha)
    return matriz


P1 = lambda x, x0, x1, y1, y0 : y0 + ((y1 - y0)/(x1 - x0))*(x - x0)