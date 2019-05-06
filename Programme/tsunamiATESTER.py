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
    
    
    for iElem in range(nElem):
        nodes = elem[iElem]
        
        x = X[nodes]
        y = Y[nodes]
        
        jac = abs((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]))
        dphidx = [y[1]-y[2], y[2]-y[0], y[0]-y[1]] / jac
        dphidy = [x[2]-x[1], x[0]-x[2], x[1]-x[0]] / jac
        
        for j in range(3):
            xsi = Xsi[j]
            eta = Eta[j]   
            phi = [(1-xsi-eta), xsi, eta]
            u = U[iElem,0]*(1-xsi-eta)+U[iElem,1]*(xsi)+U[iElem,2]*(eta)
            v = V[iElem,0]*(1-xsi-eta)+V[iElem,1]*(xsi)+V[iElem,2]*(eta)
            e = E[iElem,0]*(1-xsi-eta)+E[iElem,1]*(xsi)+E[iElem,2]*(eta)
            h = H[nodes[0]]*(1-xsi-eta)+H[nodes[1]]*(xsi)+H[nodes[2]]*(eta)
            xi = x[0]*(1-xsi-eta)+x[1]*(xsi)+x[2]*(eta)
            yi = y[0]*(1-xsi-eta)+y[1]*(xsi)+y[2]*(eta)
            z3d = R*(4*R*R - xi*xi - yi*yi) / (4*R*R + xi*xi + yi*yi)
#            lat = np.arcsin(z3d/R)*180/np.pi
#            f = 2 * (np.pi/43200) * np.sin(lat)
            f = (4*np.pi*z3d) / (R*86400)   
            p1 = (4*R*R + xi**2 + yi**2) / (4*R*R)
            p2 = (h*(xi*u+yi*v)) / (R**2)
            p3 = (g*xi*e) / (2*(R**2))
            p4 = (g*yi*e) / (2*(R**2))
            
            for k in range(3):
                coefIntegr = jac*weight
                Ei[iElem, k] += coefIntegr * ((h*u*dphidx[k] + h*v*dphidy[k])*p1) + coefIntegr * (phi[k]*p2)
                Ui[iElem, k] += coefIntegr * (phi[k]*(f*v-u*gamma) + dphidx[k]*g*e*p1 + phi[k]*p3)
                Vi[iElem, k] += coefIntegr * (phi[k]*((-1)*f*u-v*gamma) + dphidy[k]*g*e*p1 + phi[k]*p4)
                


def addIntegralsEdges(X, Y, H, Ei, Ui, Vi, E, U, V, nEdges, nBoundary, edges, mapEdgeLeft, mapEdgeRight):
#    print('coucou depuis addIntegralsEdges')
    R = 6371220
    g = 9.81
    Xsi = np.array([-0.5773502691896257, 0.5773502691896257])
    
    for iEdge in range(nBoundary):
        elemL = edges[iEdge][2]
        nodes = edges[iEdge][0:2]
        nodeLeft =  mapEdgeLeft[iEdge]
#        nodeLeft =  [np.nonzero(nodesl == edges[iedg][j])[0][0] for j in range(2)]
        x = X[nodes]   
        y = Y[nodes] 
        dx = x[1] - x[0]
        dy = y[1] - y[0]
        jac = np.sqrt(dx*dx+dy*dy)
        nx = dy / jac
        ny = -dx / jac
        jac = jac/2
        
        for j in range(2):
            xsi = Xsi[j]
            phi = np.array([(1-xsi)/2,(1+xsi)/2])
            eL = E[elemL, nodeLeft[0]]*phi[0]+E[elemL, nodeLeft[1]]*phi[1]
            uL = U[elemL, nodeLeft[0]]*phi[0]+U[elemL, nodeLeft[1]]*phi[1]
            vL = V[elemL, nodeLeft[0]]*phi[0]+V[elemL, nodeLeft[1]]*phi[1]
            xi = x[0]*phi[0]+x[1]*phi[1]
            yi = y[0]*phi[0]+y[1]*phi[1]
            h = H[nodes[0]]*phi[0] + H[nodes[1]]*phi[1]
            p1 = (4*R*R + xi**2 + yi**2) / (4*R*R)
            unL = nx*uL + ny*vL
#            unR = -unL
            etaEt = eL + np.sqrt(h/g) * (unL)                 #peut etre g/h #pas -unL
            
            for k in range(2):
                Ui[elemL, nodeLeft[k]] -= jac * (phi[k]*nx*g*etaEt*p1) #nl[k] et pas k
                Vi[elemL, nodeLeft[k]] -= jac * (phi[k]*ny*g*etaEt*p1)
        
    for iEdge in range(nBoundary, nEdges):
        elemL = edges[iEdge][2]
        elemR = edges[iEdge][3]
        nodeLeft =  mapEdgeLeft[iEdge]
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
        
        
        for j in range(2):
            xsi = Xsi[j]
            phi = [(1-xsi)/2,(1+xsi)/2]
            eL = E[elemL, nodeLeft[0]] *phi[0]+E[elemL, nodeLeft[1]] *phi[1]
            eR = E[elemR, nodeRight[0]]*phi[0]+E[elemR, nodeRight[1]]*phi[1]
            uL = U[elemL, nodeLeft[0]] *phi[0]+U[elemL, nodeLeft[1]] *phi[1]
            uR = U[elemR, nodeRight[0]]*phi[0]+U[elemR, nodeRight[1]]*phi[1]
            vL = V[elemL, nodeLeft[0]] *phi[0]+V[elemL, nodeLeft[1]] *phi[1]
            vR = V[elemR, nodeRight[0]]*phi[0]+V[elemR, nodeRight[1]]*phi[1]
            xi = x[0]*phi[0]+x[1]*phi[1]
            yi = y[0]*phi[0]+y[1]*phi[1]
            h = H[nodes[0]]*phi[0] + H[nodes[1]]*phi[1]
            p1 = (4*R*R + xi**2 + yi**2) / (4*R*R)
            unL = uL*nx + vL*ny
            unR = uR*nx + vR*ny        
            
            etaEt = (eL+eR)/2 + np.sqrt(h/g) * (unL-unR)/2 #je change R-L en L-R
            unEt  = (unL+unR)/2 + np.sqrt(g/h) * (eL-eR)/2 #peut etre switch entre eL et eR (de base)
            
            for k in range(2):
                Ei[elemL, nodeLeft[k]]  -= jac * (phi[k]*h*unEt*p1)
                Ei[elemR, nodeRight[k]] += jac * (phi[k]*h*unEt*p1)
                Ui[elemL, nodeLeft[k]]  -= jac * (phi[k]*nx*g*etaEt*p1)
                Ui[elemR, nodeRight[k]] += jac * (phi[k]*nx*g*etaEt*p1)
                Vi[elemL, nodeLeft[k]]  -= jac * (phi[k]*ny*g*etaEt*p1)
                Vi[elemR, nodeRight[k]] += jac * (phi[k]*ny*g*etaEt*p1)
                
                

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
    
    nNode,X,Y,H,nElem,elem = readMesh(theMeshFile)
    nEdges, nBoundary, edges = ComputeEdges(theMeshFile)
    mapEdgeLeft, mapEdgeRight = mapEdge(theMeshFile, elem)
    
    for s in range(nIter):
        tic()
        
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
            writeResult(theResultFiles ,s,E)
        
        toc()
    
    return [U,V,E]