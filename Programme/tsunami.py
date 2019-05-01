# -*- coding: utf-8 -*-
#
# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# Canevas de départ
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 

import numpy as np
print('coucou depuis tsunami.py')
# -------------------------------------------------------------------------

def readMesh(fileName) :
  print('coucou depuis readMesh')
  with open(fileName,"r") as f :
    nNode = int(f.readline().split()[3])
    xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(nNode)))
    nElem = int(f.readline().split()[3])
    elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(nElem)))
  X = xyz[:,0]
  Y = xyz[:,1]
  H = xyz[:,2] 
  return [nNode,X,Y,H,nElem,elem]

# -------------------------------------------------------------------------
  
def readResult(fileBaseName,iter,nElem) :
  fileName = fileBaseName % iter
  with open(fileName,"r") as f :
    nSize = int(f.readline().split()[3])
    if (nElem != nSize) :
      print(" ==== Error : incoherent sizes : %d != %d" % (nElem,nSize))     
    E = np.array(list(list(float(w) for w in f.readline().split()[2:5]) for i in range(nElem)))
    #print(" === iteration %6d : reading %s ===" % (iter,fileName))
  return E

# -------------------------------------------------------------------------

def writeResult(fileBaseName,iter,E) :
  fileName = fileBaseName % iter
  nElem = E.shape[0]  
  with open(fileName,"w") as f :
    f.write("Number of elements %d\n" % nElem)
    for i in range(nElem):
      f.write("%6d : %14.7e %14.7e %14.7e\n" % (i,*E[i,:]))
    print(" === iteration %6d : writing %s ===" % (iter,fileName))
 
# -------------------------------------------------------------------------

def initialConditionOkada(x,y) :
  R = 6371220;
  x3d = 4*R*R*x / (4*R*R + x*x + y*y);
  y3d = 4*R*R*y / (4*R*R + x*x + y*y);
  z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
  lat = np.arcsin(z3d/R)*180/np.pi;
  lon = np.arctan2(y3d,x3d)*180/np.pi;
  lonMin = 142;
  lonMax = 143.75;
  latMin = 35.9;
  latMax = 39.5;
  olon = (lonMin+lonMax)/2;
  olat = (latMin+latMax)/2;
  angle = -12.95*np.pi/180; 
  lon2 = olon + (lon-olon)*np.cos(angle) + (lat-olat)*np.sin(angle);
  lat2 = olat - (lon-olon)*np.sin(angle) + (lat-olat)*np.cos(angle);
  return np.all([lon2 <= lonMax,lon2 >= lonMin,lat2 >= latMin,lat2 <= latMax],axis=0).astype(int)
 
  # -------------------------------------------------------------------------
# -------------------------------------------------------------------------
  # -------------------------------------------------------------------------
# -------------------------------------------------------------------------
  # -------------------------------------------------------------------------
# -------------------------------------------------------------------------
  # -------------------------------------------------------------------------
# -------------------------------------------------------------------------
  # -------------------------------------------------------------------------
# -------------------------------------------------------------------------
  
#def IntegGL(X1, X2, f):
#    eta = 0.577350269189626
#    x1 = -1 * eta * ((X2-X1)/2) + ((X2+X1)/2)
#    x2 = eta * ((X2-X1)/2) + ((X2+X1)/2)
#    
#    return (X2-X1)/2 * (f(x1)+f(x2))
#
#def hammerIntegrate(x,y,f) :
#    A = np.array([[0.666666666666667,0.166666666666667,0.166666666666667],
#                  [0.166666666666667,0.166666666666667,0.666666666666667],
#                  [0.166666666666667,0.666666666666667,0.166666666666667]])
#    
#    xInt = np.dot(A,x)
#    yInt = np.dot(A,y)
#    
#    Jacobien = abs((x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]))
#    I = (Jacobien/6) * ( f(xInt[0], yInt[0]) + f(xInt[1], yInt[1]) + f(xInt[2], yInt[2]) )
#        
#    return I


def addIntegralsTriangles(X, Y, H, U, V, nElem, elem):
    print('coucou depuis addIntegralsTriangles')
    
    B = np.zeros(3*nElem)
    R = 6371220
    
    #C = theProblem.C
    print('===1')
    Xsi    = np.array([0.166666666666667,0.666666666666667,0.166666666666667])
    Eta    = np.array([0.166666666666667,0.166666666666667,0.666666666666667])
    weight = np.array([0.166666666666667,0.166666666666667,0.166666666666667])
    print('===2')
    dphidxsi = np.array([ 1.0, 0.0,-1.0])
    dphideta = np.array([ 0.0, 1.0,-1.0])
    print('===3')
    for iElem in range(nElem):
        mapElem = [3*iElem+j for j in range(3)]
        nodes = elem[iElem]
        
        x = X[nodes]
        y = Y[nodes]
        
        dxdxsi = x @ dphidxsi
        dxdeta = x @ dphideta
        dydxsi = y @ dphidxsi
        dydeta = y @ dphideta
                
        jac = abs(dxdxsi*dydeta - dxdeta*dydxsi)
        
        dphidx = (dphidxsi * dydeta - dphideta * dydxsi) / jac
        dphidy = (dphideta * dxdxsi - dphidxsi * dxdeta) / jac
        print('===4')
        for j in range(3):
            xsi = Xsi[j]
            eta = Eta[j]   
            print('===5')
            u = U[mapElem[0]]*(1-xsi-eta)+U[mapElem[1]]*(xsi)+U[mapElem[2]]*(eta)
            v = V[mapElem[0]]*(1-xsi-eta)+V[mapElem[1]]*(xsi)+V[mapElem[2]]*(eta)
            h = 5
            p = (4*R*R + x[j]**2 + y[j]**2) / (4*R*R)
            print('===6')
            for k in range(3):
                
                a = mapElem[k]
                w = weight[j]
                dxk = dphidx[k]
                dyk = dphidy[k]
                truc = jac * p * w * h
                print('===7')
                print(truc)
                machin = (u * dxk + v * dyk)
                print('===7')
                print(machin)
                print('===7')
                B[a] = truc @ machin
                print('===8')
    print('===9')            
    return B


def compute(theMeshFile,theResultFiles,U,V,E,dt,nIter,nSave):
    print('coucou depuis compute')
    
    [nNode,X,Y,H,nElem,elem] = readMesh(theMeshFile)
    
    B = addIntegralsTriangles(X, Y, H, U, V, nElem, elem)
    print(B)
    
    return [U,V,E]