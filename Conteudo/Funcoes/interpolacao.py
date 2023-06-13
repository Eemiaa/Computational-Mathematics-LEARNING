import pandas as pd
import operator
import numpy as np

def matriz_inter_quadratica(data):
    matriz = []
    for i in range(len(data.x)):
        linha=[]
        for j in range(len(data.x)):
            if (j==0):
                linha.append(1)
            else:   
                linha.append(data.x[i]**j)
        linha.append(data.y[i])
        matriz.append(linha)
    return matriz

def Pquadratica(coeficientes, x):
    resultado = 0
    for i in range(len(coeficientes)):
        resultado += coeficientes[i]*x**i
    return resultado

Plinear = lambda x, x0, x1, y1, y0 : y0 + ((y1 - y0)/(x1 - x0))*(x - x0)

def PolinomioLagrange(x, y, z):
    S = 0
    m = len(x)
    for i in range(1, m):
        c=1
        d=1
        for j in range(1, m):
            if i!=j:
                c = c*(z-x[j])
                d = d*(x[i]-x[j]) 
        S = S + y[i] * c/d
    return S

def tabela_diferencas_divididas(x, y):
    tab = pd.DataFrame({
    }, index=[i for i in x], columns=["Ordem {}".format(i) for i in range(len(x)-1)])

    #tabela:
    for j in range(len(x)-1):
        auxcolumn = []
        for i in range(1, len(x)):
            if(i < len(x)-j):
                if(j == 0):
                    auxcolumn.append((y[i] - y[i-1])/(x[i+j] - x[i-1]))
                else:
                    auxcolumn.append((tab["Ordem {}".format(j-1)][x[i]] - tab["Ordem {}".format(j-1)][x[i-1]])/(x[i+j] - x[i-1]))
            else:
                auxcolumn.append(None)
        auxcolumn.append(None)
        tab["Ordem {}".format(j)] = auxcolumn

    
    #b's:
    bs = operator.add([y[0]], tab.loc[x[0]].tolist())
    
    return tab, bs

def polinomioNewton(bs, xinterpolar, x):
    resposta = 0
    for i in range(len(bs)):
        multiplicando = 1
        #print(bs[i])
        for j in range(i):
            multiplicando *= (xinterpolar - x[j])
            #print("({} - {})".format(xinterpolar , x[j]))
        resposta += bs[i]*multiplicando
    
    return resposta

def neville( x, y, X, verbose=False ):
  """
  Calcula a matriz de Neville em X
  Recebe:
    x,y     => Vetores com dados
    X       => Posição onde interpolar
    verbose => Mostra passo-a-passo
  Retorna:
    Matriz de Neville
  """
  N = len(x)
  
  assert len(y)==N, "Vetores de dados não têm o mesmo tamanho!"
  assert x[0] <= X <= x[-1] , "Ponto fora do domínio"
  
  Q = np.zeros( (N, N) )
  for i in range(N):
    Q[i,0] = y[i]
  
  for i in range(1,N):
    for j in range(1,i+1):
      Q[i,j] = ( (X-x[i-j])*Q[i,j-1] - (X-x[i])*Q[i-1,j-1] )/float(x[i]-x[i-j])
      
      if verbose:
        print("Q%d%d = (x-x%d)Q%d%d - (x-x%d)Q%d%d" % (i,j,i-j,i,j-1,i,i-1,j-1))
        print("      ---------------------"                                    )
        print("            (x%d - x%d)" % (i,i-j)                              )

  return Q

def solve_spline( x, y, derivatives=[], verbose=False ):
  """
  Resolve o sistema linear das splines para encontrar os coeficientes.

  Recebe:
    x,y         => Vetores com dados
    derivatives => Se recebe alguma coisas, calcula spline fixa, caso contrário, livre.
    verbose     => Mostra passo-a-passo
  Retorna:
    coefficients => Matriz com coeficientes da spline interpolada
                    Linha: intervalos
                    Coluna: potências
                    
                    [ [ a0, b0, c0, d0 ],
                      [ a1, b1, c1, d1 ],
                      ... ]
  """
  n = len(x)  # Número de pontos
  assert len(y) == n, "Vetores de dados com tamanhos diferentes!"

  a = y.copy()
  h = np.zeros(n-1)
  for i in range(n-1):
    h[i] = x[i+1] - x[i] 

  # Matrizes para resolver
  A = np.zeros( (n,n) )
  B = np.zeros( n )

  # Parte comum à Spline Fixa e Natural
  for i in range(1,n-1):
    A[i][i-1] = h[i-1]
    A[i][i]   = 2*( h[i-1] + h[i] )
    A[i][i+1] = h[i]
  
    B[i] = 3.0*(a[i+1] - a[i])/float(h[i]) - 3.0*(a[i] - a[i-1])/float(h[i-1])
  
  # Parte que depende do tipo da Spline
  # Tenta usar as derivadas, se não der, faz a spline natural
  try:
    df0 = derivatives[0]
    df1 = derivatives[1]

    A[0][0]   = 2*h[0]
    A[0][1]   = h[0]
    A[-1][-2] = h[-2]
    A[-1][-1] = 2*h[-2]
  
    B[0] = 3.0*(a[1] - a[0])/float(h[0]) - 3*df0
    B[-1] = 3*df1 - 3.0*(a[-1] - a[-2])/float(h[-2])
  
  # Spline natural
  except:
    A[0][0]   = 1.0
    A[-1][-1] = 1.0

  # Encontra os coeficientes
  c = np.linalg.solve( A, B )
  b = np.zeros(n-1)
  d = np.zeros(n-1)
  for i in range(n-1):
    b[i] = (a[i+1]-a[i])/float(h[i]) - h[i]*( 2*c[i] + c[i+1] )/3.0
    d[i] = (c[i+1]-c[i])/( 3.0 * h[i] )

  if verbose:
    
    print("Matriz de splines:")
    print()
    for i in range(len(A)):
      for j in range(len(A[i])):
        print(A[i][j],end=" ")
      print("| ",B[i])
    
    print("Vetores:")
    print(" i  |      x      |      y      |       h     |       a     |      b      |      c      |      d      |")
    fmt="%2d  " + "|  %9.6g  "*7
    for i in range(n-1):
      print(fmt%(i,x[i],y[i],h[i],a[i],b[i],c[i],d[i]))

  coef = np.zeros( (n-1,4) )
  for i in range(n-1):
    coef[i][0] = a[i]
    coef[i][1] = b[i]
    coef[i][2] = c[i]
    coef[i][3] = d[i]


  return coef

def calc_spline( x, X, coef ):
  """
  Calcula a spline no ponto X
  Os coeficientes estão na matriz coef
  
  Recebe:
    x => Vetor com dados
    X => Ponto (escalar ou vetor) onde calcular a spline
    coef => Coeficientes da spline, conforme calculados pela função solve_spline
  
  Retorna:
    valor ou vetor da spline calculada em X
  """
  
  y = 0
  
  try:
    n = len(X)
    
    y = np.zeros(n)
    for i in range(n):
      y[i] = calc_spline( x, X[i], coef )
  
  except:
    
    # Encontra a posição de X dentro do vetor x
    k = x.searchsorted( X )
  
    if k>0      : k -= 1  # Use a função S do ponto da esquerda
    if k==len(x): k -= 1  # X[i] > x[n]
    
    H = X - x[k]
    ak = coef[k][0]
    bk = coef[k][1]
    ck = coef[k][2]
    dk = coef[k][3]
    y = ak + H*( bk + H*( ck + H*dk ) )
  
  
  return y