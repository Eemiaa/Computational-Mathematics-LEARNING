import numpy as np
import math 
from numpy.linalg import solve

class MQ:
  def __init__(self):
    self.alfas = []

  def fit_exp(self, x, y):
    self.alfas=[]
    lnY = np.log(y)
    self.fit(x, lnY, [lambda x:1, lambda x:x])
    self.alfas[0] = math.e**self.alfas[0]
    self.alfas[1] =- self.alfas[1]

  def fit(self,x,y,G):
    self.G=G
    A=[]
    B=[]
    j=0
    for g_lin in G:
      b=0
      for i in range(0,len(x)):
        b+=g_lin(x[i])*y[i]
      B.append(b)
      A.append([])
      for g_col in G:
        a=0
        for i in range(0,len(x)):
          a+=g_lin(x[i])*g_col(x[i])
        A[j].append(a)
      j+=1

    mat = np.append(A, np.array([B]).T,axis=1)
    #print('A: ', A)
    #print(mat)
    #print('B: ', B)
    #print('mat: ', mat)
    self.alfas = solve(A, B)
    #print(f"Alfas: {self.alfas}")

  def calc(self, x):
    s = 0
    #print(self.alfas)
    for i in range(0,len(self.G)):
      s+=self.alfas[i]*self.G[i](x)
    return s
  
  def calc_exp(self, x):
    return self.alfas[0]*(math.e**(-self.alfas[1]*x))
  

  