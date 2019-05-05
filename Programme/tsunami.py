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
from timeit import default_timer as timer
#print('coucou depuis tsunami.py')
# -------------------------------------------------------------------------
def tic(message = ''):
  global startTime
  startTime = timer()

def toc(message = ''):
  global startTime
  stopTime = timer()
  if message:
    message = ' (' + message + ')' ;
  print("Elapsed time is %.6f seconds %s" % ((stopTime - startTime),message) )
  elapsedTime = stopTime - startTime;
  startTime = timer()
  return elapsedTime
# -------------------------------------------------------------------------

def readMesh(fileName) :
#  print('coucou depuis readMesh')
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
  z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y)
  lat = np.arcsin(z3d/R)*180/np.pi
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
  
def ComputeEdges(mesh): #vient de fem.py
    [nNode,X,Y,H,nElem,elem] = readMesh(mesh)
    nEdges = nElem * 3
    nBoundary = 0
    edges = [[0 for i in range(4)] for i in range(nEdges)]
    for i in range (nElem) :
      for j in range(3) :
        id = i*3 + j
        edges[id][0] = elem[i][j]
        edges[id][1] = elem[i][(j+1)%3]
        edges[id][2] = i
        edges[id][3] = -1
    edges.sort(key = lambda item : -(min(item[0:2])*nEdges)-max(item[0:2]))
    index = 0
    for i in range(nEdges) :
      if (edges[i][0:2] != edges[i-1][1::-1]) :
         edges[index] = edges[i]
         index += 1
      else :
         edges[index-1][3] = edges[i][2]
    del edges[index:]
    edges.sort(key = lambda item : item[3])
    nBoundary = 2*index - nEdges
    nEdges = index
    
    return nEdges, nBoundary, edges

def mapEdge(theMeshFile, elem): #vient du prof
    [nEdges, nBoundary, edges] = ComputeEdges(theMeshFile)
    
    mapEdgeLeft = np.zeros((nEdges,2),dtype=np.int)
    mapEdgeRight = np.zeros((nEdges,2),dtype=np.int)
    for iEdge in range(nBoundary,nEdges):
        myEdge = edges[iEdge]
        elementLeft  = myEdge[2]
        elementRight = myEdge[3]
        nodesLeft    = elem[elementLeft]
        nodesRight   = elem[elementRight]
        mapEdgeLeft[iEdge,:]  = [np.nonzero(nodesLeft  == myEdge[j])[0][0] for j in range(2)]
        mapEdgeRight[iEdge,:] = [np.nonzero(nodesRight == myEdge[j])[0][0] for j in range(2)]
    return[mapEdgeLeft, mapEdgeRight]

def addIntegralsTriangles(X, Y, H, Ei, Ui, Vi, E, U, V, nElem, elem):
#    print('coucou depuis addIntegralsTriangles')
    R = 6371220
    gamma = 10**(-7)
    g = 9.81
    
    Xsi    = np.array([0.166666666666667,0.666666666666667,0.166666666666667])
    Eta    = np.array([0.166666666666667,0.166666666666667,0.666666666666667])
    weight = 0.166666666666667
    
    dphidxsi = np.array([ 1.0, 0.0,-1.0])
    dphideta = np.array([ 0.0, 1.0,-1.0])
    
    Phi = np.array([1-Xsi-Eta, Xsi, Eta])
    
    for iElem in range(nElem):
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
        
        u = U[iElem] @ Phi
        v = V[iElem] @ Phi
        e = E[iElem] @ Phi
        h =  H[nodes] @ Phi
        xi = x @ Phi
        yi = y @ Phi
        z3d = R*(4*R*R - xi*xi - yi*yi) / (4*R*R + xi*xi + yi*yi)
        f = (4*np.pi*z3d) / (R*86400)   
        p1 = (4*R*R + xi**2 + yi**2) / (4*R*R)
        p2 = (h*(xi*u+yi*v)) / (R**2)
        p3 = (g*xi*e) / (2*(R**2))
        for j in range(3):
            for k in range(3):
                coefIntegr = jac*weight
                Ei[iElem, k] += coefIntegr * (h[j]*p1[j]*(u[j] * dphidx[k] + v[j] * dphidy[k]) + Phi[j,k]*p2[j])
                Ui[iElem, k] += coefIntegr * (Phi[j,k]*(f[j]*v[j]-u[j]*gamma) + dphidx[k]*g*e[j]*p1[j] + Phi[j,k]*p3[j])
                Vi[iElem, k] += coefIntegr * (Phi[j,k]*(-f[j]*v[j]-u[j]*gamma) + dphidy[k]*g*e[j]*p1[j] + Phi[j,k]*p3[j])


def addIntegralsEdges(X, Y, H, Ei, Ui, Vi, E, U, V, nEdges, nBoundary, edges, mapEdgeLeft, mapEdgeRight):
#    print('coucou depuis addIntegralsEdges')
    R = 6371220
    g = 9.81
    Xsi = np.array([-0.5773502691896257, 0.5773502691896257])
    Phi = np.array([(1-Xsi)/2, (1+Xsi/2)])
    
    for iEdge in range(nBoundary):
        elemL = edges[iEdge][2]
        nodes = edges[iEdge][0:2]
        nodeLeft =  mapEdgeLeft[iEdge]
        x = X[nodes]   
        y = Y[nodes] 
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        jac = np.sqrt(dx*dx+dy*dy)
        nx = dy / jac
        ny = -dx / jac
        jac = jac/2
        eL = E[elemL, nodeLeft] @ Phi
        uL = U[elemL, nodeLeft] @ Phi
        vL = V[elemL, nodeLeft] @ Phi
        h = H[nodes] @ Phi
        xi = x @ Phi
        yi = y @ Phi
        p1 = (4*R*R + xi**2 + yi**2) / (4*R*R)
        unL = nx*uL + ny*vL
        unR = -unL
        etaEt = eL + np.sqrt(h/g) * (-unL)
        for j in range(2):
            for k in range(2):
                Ui[elemL, k] -= jac * (Phi[j,k]*nx*g*etaEt[j]*p1[j])
                Vi[elemL, k] -= jac * (Phi[j,k]*ny*g*etaEt[j]*p1[j])
        
    for iEdge in range(nBoundary, nEdges):
        elemL = edges[iEdge][2]
        elemR = edges[iEdge][3]
        nodeLeft  =  mapEdgeLeft[iEdge]
        nodeRight =  mapEdgeRight[iEdge]
        nodes = edges[iEdge][0:2] #a modif ?
        x = X[nodes]   
        y = Y[nodes]   
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        jac = np.sqrt(dx*dx+dy*dy)
        nx = dy / jac
        ny = -dx / jac
        jac = jac/2
        eL = E[elemL, nodeLeft] @ Phi
        uL = U[elemL, nodeLeft] @ Phi
        vL = V[elemL, nodeLeft] @ Phi
        eR = E[elemR, nodeRight] @ Phi
        uR = U[elemR, nodeRight] @ Phi
        vR = V[elemR, nodeRight] @ Phi
        h = H[nodes] @ Phi
        xi = x @ Phi
        yi = y @ Phi
        p1 = (4*R*R + xi**2 + yi**2) / (4*R*R)
        unL = uL*nx + vL*ny
        unR = uR*nx + vR*ny
        etaEt = (eL+eR)/2 + np.sqrt(h/g) * (unR-unL)/2
        unEt  = (unL+unR)/2 + np.sqrt(g/h) * (eR-eL)/2 #peut etre switch entre eL et eR (de base)
        
        for j in range(2):
            for k in range(2):
                Ei[elemL, k] -= jac * (Phi[j,k]*h[j]*unEt[j]*p1[j])
                Ei[elemR, k] += jac * (Phi[j,k]*h[j]*unEt[j]*p1[j])
                Ui[elemL, k] -= jac * (Phi[j,k]*nx*g*etaEt[j]*p1[j])
                Ui[elemR, k] += jac * (Phi[j,k]*nx*g*etaEt[j]*p1[j])
                Vi[elemL, k] -= jac * (Phi[j,k]*ny*g*etaEt[j]*p1[j])
                Vi[elemR, k] += jac * (Phi[j,k]*ny*g*etaEt[j]*p1[j])
                
                

def multiplyInverseMatrix(elem, nElem, Ei, Ui, Vi, X, Y):
    Ainverse = np.array([[18.0,-6.0,-6.0],[-6.0,18.0,-6.0],[-6.0,-6.0,18.0]])
    
    for iElem in range(nElem) :
        nodes = elem[iElem]
        x = X[nodes]
        y = Y[nodes]
        jac = abs((x[0]-x[1]) * (y[0]-y[2]) - (x[0]-x[2]) * (y[0]-y[1]))
        Ei[iElem] = Ainverse @ Ei[iElem] / jac
        Ui[iElem] = Ainverse @ Ui[iElem] / jac
        Vi[iElem] = Ainverse @ Vi[iElem] / jac
    

def compute(theMeshFile,theResultFiles,U,V,E,dt,nIter,nSave):
    
    for s in range(nIter):
        tic()
        
        nNode,X,Y,H,nElem,elem = readMesh(theMeshFile)
        nEdges, nBoundary, edges = ComputeEdges(theMeshFile)
        mapEdgeLeft, mapEdgeRight = mapEdge(theMeshFile, elem)
        
        Ei = np.zeros(E.shape)
        Ui = np.zeros(U.shape)
        Vi = np.zeros(V.shape)
        addIntegralsTriangles(X, Y, H, Ei, Ui, Vi, E, U, V, nElem, elem)
        addIntegralsEdges(X, Y, H, Ei, Ui, Vi, E, U, V, nEdges, nBoundary, edges, mapEdgeLeft, mapEdgeRight)
        multiplyInverseMatrix(elem, nElem, Ei, Ui, Vi, X, Y)
        
        E += dt*Ei
        U += dt*Ui
        V += dt*Vi
        if s % nSave == 0:  
            writeResult(theResultFiles,s,E)
            
        toc()
    
    return [U,V,E]