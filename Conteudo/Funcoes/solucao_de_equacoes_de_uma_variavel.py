import numpy as np
from numpy import cos, sin, exp
import sympy as sp

def funcEmIntervalo(func, a = -3, b = 3):
  aux = {"x":[],"y":[]}
  for i in range(a, b+1):
    aux["y"].append(func(i))
    aux["x"].append(i)

  return aux

def minLocal(f= lambda x : x**2, a = 1, b = 3):
  aux = f(a)
  for i in range(abs(a-b)+1):
    if(a<b):
      if aux >= f(a+i): aux = f(a+i)
    else:
      if aux >= f(a+i): aux = f(b+i)
  return aux

def maxLocal(f= lambda x : x**2, a = 1, b = 3):
  aux = f(a)
  for i in range(abs(a-b)+1):
    if(a<b):
      if aux <= f(a+i): aux = f(a+i)
    else:
      if aux <= f(a+i): aux = f(b+i)
  return aux

def f(x):
  return np.exp(x)-x-2

def fl(x): 
  return np.exp(x)-1

def bissecao(fx, a=-2, b=0, aproximacao=5):
  f = lambda x: eval(fx)

  cont=0
  while(True):
    x = (a + b)/2
    if(round(f(x), aproximacao) == 0):
      print("\n\nResultado:{}; f(x): {}".format(x, round(f(x), aproximacao)))
      break

    if f(a)*f(x) < 0 : b = x
    else: a = x
    print("\n{}:\n- x: {}\na: {} \nb: {}".format(cont, x, a, b))   

    cont+=1
  
def pontoFixo(fx, gx, a=0.05, x=1.7, aproximacao=5):
  gx = gx.replace('f(x)', "({})".format(fx))
  g = lambda x, a: eval(gx)
  f = lambda x: eval(fx)
  cont=0

  while(True):
    if round(f(x), aproximacao) == 0: 
      print("\n\nResultado:{}; f(x): {}".format(x, round(f(x), aproximacao)))
      break
    print("\n{}: {}".format(cont, x))
    x = g(x, a)

    cont+=1

def newtonRaphson(fx , aproximacao=5, x1=1):
  f = lambda x : eval(fx)
  d1 = lambda x : eval(str(sp.diff(fx, sp.Symbol('x'))))
  d2 = lambda x : eval(str(sp.diff(fx, sp.Symbol('x'), 2)))
  g = lambda x : x - f(x)/d1(x)
  x0 = 0
  cont = 0
  
  while(True):
    x1 = g(x1)
    if round(x0, aproximacao)==round(x1, aproximacao):
      print("\n{}: {}".format(cont, x1))
      print("\n\nResultado:{}; x0:{}".format(x1,x0))
      break
    print("\n{}: {}".format(cont, x1))
    x0 = x1
    cont+=1

def secante(fx, x1=1, x2=2, aproximacao=5):
  f = lambda x : eval(fx)
  g = lambda x2 : x2 - f(x2)*( (x2-x1)  / (f(x2)-f(x1)) )
  x0 = 0 
  cont=0

  while(True):
    x2 = g(x2)
    if round(x0, aproximacao)==round(x2, aproximacao):
      print("\n{}: {}".format(cont, x2))
      print("\n\nResultado:{}\n x0:{}\n".format(x2, x0))
      break
    print("\n{}: {}".format(cont, x2))
    x0 = x2
    cont+=1


