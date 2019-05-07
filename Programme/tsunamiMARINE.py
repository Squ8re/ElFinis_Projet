# -------------------------------------------------------------------------
#
# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# -------------------------------------------------------------------------
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
#
# Sylvie Van Schendel (44841700)
# Marine Van Renterghem (31621700)
# FSA1BA
#
# -------------------------------------------------------------------------

import numpy as np
from timeit import default_timer as timer

# -------------------------------------------------------------------------
def tic(message = ''):
  global startTime
  startTime = timer()

def toc(Es = None,IsEs = False,message = ''):
  global startTime
  stopTime = timer()
  if message:
    message = ' (' + message + ')' ;
#  print("Elapsed time is %.6f seconds %s" % ((stopTime - startTime),message) )
  if(IsEs == True):
      for iElem in [27] :
          print(" == Elevations for element %d : %14.7e %14.7e %14.7e " % (iElem,*Es[iElem][:]) )
  elapsedTime = stopTime - startTime;
  startTime = timer()
  return elapsedTime
#---------------------------------------------------------------------------
 
def readMesh(fileName) :
  with open(fileName,"r") as f :
    nNode = int(f.readline().split()[3])
    xyz   = np.array(list(list(float(w) for w in f.readline().split()[2:]) for i in range(nNode)))
    nElem = int(f.readline().split()[3])
    elem  = np.array(list(list(int(w)   for w in f.readline().split()[2:]) for i in range(nElem)))
  X = xyz[:,0]
  Y = xyz[:,1]
  H = xyz[:,2]
  return [nNode,X,Y,H,nElem,elem]
#--------------------------------------------------------------------------
def computeEdge(nElem,elem):
    nEdges=nElem*3
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
    return [nEdges,nBoundary,edges]
# -------------------------------------------------------------------------
 
def readResult(fileBaseName,iter,nElem) :
  fileName = fileBaseName % iter
  with open(fileName,"r") as f :
    nSize = int(f.readline().split()[3])
    if (nElem != nSize) :
      print(" ==== Error : incoherent sizes : %d != %d" % (nElem,nSize))    
    E = np.array(list(list(float(w) for w in f.readline().split()[2:5]) for i in range(nElem)))
    print(" === iteration %6d : reading %s ===" % (iter,fileName))
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

 
#def mapEdge(elem,edges,iEdge) : 
#    
#  myEdge = edges[iEdge]
#  elementLeft  = myEdge[2]
#  nodesLeft    = elem[elementLeft]
#  mapEdgeLeft  = [np.nonzero(nodesLeft == myEdge[j])[0][0]  for j in range(2)]
#  elementRight = myEdge[3]
#  if (elementRight != -1) :
#      nodesRight   = elem[elementRight]
#      mapEdgeRight = [np.nonzero(nodesRight == myEdge[j])[0][0] for j in range(2)]
#  else :
#    mapEdgeRight = []                  
#  return [mapEdgeLeft,mapEdgeRight]
def mapEdge(elem,nEdges,edges): 
    mapEdgeLeft = np.zeros((nEdges,2),dtype=int)
    mapEdgeRight = np.zeros((nEdges,2),dtype=int)
    for iEdge in range(nEdges):
        myEdge = edges[iEdge]
        elementLeft  = myEdge[2]
        nodesLeft    = elem[elementLeft]
        mapEdgeLeft[iEdge][:]  = [np.nonzero(nodesLeft  == myEdge[j])[0][0] for j in range(2)]
        
        elementRight = myEdge[3]
        if (elementRight != -1) :
            nodesRight   = elem[elementRight]
            mapEdgeRight[iEdge][:] = [np.nonzero(nodesRight == myEdge[j])[0][0] for j in range(2)]
    return[mapEdgeLeft, mapEdgeRight]
 
# -------------------------------------------------------------------------
 
def compute(theMeshFile,theResultFiles,U,V,E,dt,nIter,nSave):
    [nNode,X,Y,H,nElem,elem]=readMesh(theMeshFile)
    [nEdges,nBoundary,edges]=computeEdge(nElem,elem)
    [mapLeft,mapRight]=mapEdge(elem,nEdges,edges)
   
    R = 6371220
    #save=nSave
    xsiT=np.array([0.166666666666667,0.666666666666667,0.166666666666667])
    etaT=np.array([0.166666666666667,0.166666666666667,0.666666666666667])
    #wT=np.array([0.166666666666667,0.166666666666667,0.166666666666667])
    wT=0.166666666666667
    xsiE=np.array([-0.5773502691896257, 0.5773502691896257],dtype=float)
    #wE=1.0
    #wE=np.array([1.0,1.0],dtype=float)
    
    funformT=np.array([1-xsiT-etaT, xsiT, etaT])
    funformE=np.array([(1-xsiE),(1+xsiE)])/2
    gamma=1e-7
    g=9.81
   
    for iter in range (nIter) :
        tic()
        Be=np.zeros((nElem,3))
        Bu=np.zeros((nElem,3))
        Bv=np.zeros((nElem,3))
       
        #Integrales de triangles
        for iElem in range (nElem) :
              dphidxsi = np.array([ 1.0, 0.0,-1.0])
              dphideta = np.array([ 0.0, 1.0,-1.0])
              nodes = elem[iElem]
              x = X[nodes]
              y = Y[nodes]
              h = H[nodes]
              dxdxsi = x @ dphidxsi
              dxdeta = x @ dphideta
              dydxsi = y @ dphidxsi
              dydeta = y @ dphideta
              jac = abs(dxdxsi*dydeta - dxdeta*dydxsi)
              dphidx = (dphidxsi * dydeta - dphideta * dydxsi) / jac
              dphidy = (dphideta * dxdxsi - dphidxsi * dxdeta) / jac
              u=funformT@U[iElem]
              
              v=funformT@V[iElem]
              e=funformT@E[iElem]
              h=funformT@h
              x=funformT@x
              y=funformT@y
              Term=(4*R*R+x*x+y*y)/(4*R*R)
              z3d= R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
              f=4*np.pi*z3d/(86400*R)
              
#              Be[iElem]+=sum((np.outer(u*h*(4*R**2+x**2+y**2)/(4*R**2),dphidx) + np.outer(v*h*(4*R**2+x**2+y**2)/(4*R**2),dphidy))*wT*jac) + sum(np.multiply((np.outer(h*x*u,np.ones(3)) + np.outer(v*h*y,np.ones(3)))*wT*jac/R**2,funformT))
#              
#              Bu[iElem]+=sum(np.multiply((np.outer(v*f,np.ones(3))-np.outer(gamma*u,np.ones(3))+np.outer(e*x,g/(2*R**2)*np.ones(3)))*wT*jac,funformT))+ sum(np.outer(e*g*(4*R**2+x**2+y**2)/(4*R**2),dphidx)*wT*jac) 
#              
#              Bv[iElem]+=sum(np.multiply(np.outer(g/(2*R**2)*e*y,np.ones(3))+np.outer(u*(-f),np.ones(3))-np.outer(gamma*v,np.ones(3)),funformT)*wT*jac)+ sum(np.outer(g*e*(4*R**2+x**2+y**2)/(4*R**2),dphidy)*jac*wT) 
              
              Be[iElem]+=sum(np.outer(u*h*Term,dphidx) + np.outer(v*h*Term,dphidy))*wT*jac + sum(np.multiply(np.outer(h*(x*u+v*y),np.ones(3)),funformT))*wT*jac/(R*R)
              Bu[iElem]+=sum(np.multiply(np.outer(v*f-gamma*u+(g*x*e)/(2*R*R),np.ones(3)),funformT))*wT*jac+ sum(np.outer(e*g*Term,dphidx))*wT*jac 
              Bv[iElem]+=sum(np.multiply(np.outer(-u*f-gamma*v+(g*e*y)/(2*R*R),np.ones(3)),funformT))*wT*jac+ sum(np.outer(g*e*Term,dphidy))*jac*wT 
    
#              for i in range(3):
#                  for j in range(3):
#                      
#                      Be[iElem,i]+=(dphidx[i]*h[j]*u[j]+dphidy[i]*h[j]*v[j])*(4*R**2+x[j]**2+y[j]**2)/(4*R**2)*jac*wT[j]
#                      Be[iElem,i]+=funformT[i][j]*(h[j]*(x[j]*u[j]+y[j]*v[j]))/R**2*jac*wT[j]
#                      Bu[iElem,i]+=(funformT[i][j]*(f[j]*v[j]-gamma*u[j])+dphidx[i]*g*e[j]*(4*R**2+x[j]**2+y[j]**2)/(4*R**2))*jac*wT[j]
#                      Bu[iElem,i]+=funformT[i][j]*(g*x[j]*e[j])/(2*R**2)*jac*wT[j]
#                      Bv[iElem,i]+=(funformT[i][j]*(-f[j]*u[j]-gamma*v[j])+dphidy[i]*g*e[j]*(4*R**2+x[j]**2+y[j]**2)/(4*R**2))*jac*wT[j]
#                      Bv[iElem,i]+=funformT[i][j]*(g*y[j]*e[j])/(2*R**2)*jac*wT[j]
              
        #print(Bu[36])            
        #Integrales sur les lignes frontieres
        for iEdge in range(nBoundary) :
            nodes = edges[iEdge][0:2]
            x = X[nodes] 
            y = Y[nodes]
            h = H[nodes]
            dx = x[1] - x[0]
            dy = y[1] - y[0]
            jac = np.sqrt(dx*dx+dy*dy)
            nx = dy / jac
            ny = -dx / jac
            jac=jac/2
            mapEdgeLeft=mapLeft[iEdge]
            
            x=funformE@x
            y=funformE@y
            h=funformE@h
            

            uL=funformE@U[edges[iEdge][2],mapEdgeLeft]
            vL=funformE@V[edges[iEdge][2],mapEdgeLeft]
            unL=uL*nx+vL*ny
            EL=funformE@E[edges[iEdge][2],mapEdgeLeft]
            Estar= EL + np.sqrt(h/g) * unL
            Term=(4*R*R+x*x+y*y) / (4*R*R)
            
            Bu[edges[iEdge][2],mapEdgeLeft] -= (funformE)@(nx*g*Estar*Term) *jac
            Bv[edges[iEdge][2],mapEdgeLeft] -=  (funformE)@(ny*g*Estar* Term)*jac 
#            for i in range(2):
#                for j in range(2):
##                    Be[mapEdgeLeft[i]]+= funformE[i][j]*h[i]*unstar[j] * ( (4*R**2+x[i]**2+y[i]**2) / (4*R**2) ) * jac * wE[j]
##                    Be[mapEdgeRight[i]] -= funformE[i][j]*h[i]*unstar[j] * ( (4*R**2+x[i]**2+y[i]**2) / (4*R**2) ) *jac*wE[j]
#                    Bu[edges[iEdge][2],mapEdgeLeft[i]] -= funformE[i][j]*nx*g*Estar[j] * ( (4*R**2+x[j]**2+y[j]**2) / (4*R**2) ) *jac*wE[j]
#                    #Bu[edges[iEdge][3],mapEdgeRight[i]] -= funformE[i][j]*nx*g*Estar[j] * ( (4*R**2+x[i]**2+y[i]**2) / (4*R**2) ) *jac*wE[j]
#                    Bv[edges[iEdge][2],mapEdgeLeft[i]] -= funformE[i][j]*ny*g*Estar[j] * ( (4*R**2+x[j]**2+y[j]**2) / (4*R**2) ) *jac*wE[j]
#                    #Bv[edges[iEdge][3],mapEdgeRight[i]] -= funformE[i][j]*ny*g*Estar[j] * ( (4*R**2+x[i]**2+y[i]**2) / (4*R**2) ) *jac*wE[j]
#               
        
        #Integrales sur les lignes interieures
        for iEdge in range(nBoundary,nEdges):
            nodes = edges[iEdge][0:2]
            x = X[nodes] 
            y = Y[nodes]
            h = H[nodes]
            dx = x[1] - x[0]
            dy = y[1] - y[0]
            jac = np.sqrt(dx*dx+dy*dy)
            nx = dy / jac
            ny =  -dx / jac
            jac=jac/2
            mapEdgeLeft=mapLeft[iEdge]
            mapEdgeRight=mapRight[iEdge]
            uL=funformE@U[edges[iEdge][2],mapEdgeLeft]
            uR=funformE@U[edges[iEdge][3],mapEdgeRight]
            vL=funformE@V[edges[iEdge][2],mapEdgeLeft]
            vR=funformE@V[edges[iEdge][3],mapEdgeRight]
            unL=uL*nx+vL*ny
            unR=uR*nx+vR*ny
            EL=funformE@E[edges[iEdge][2],mapEdgeLeft]
            ER=funformE@E[edges[iEdge][3],mapEdgeRight]
            
            x=funformE@x
            y=funformE@y
            h=funformE@h
            
            Estar= (EL+ER)/2 + np.sqrt(h/g) * (unL-unR)/2
            unstar= (unL+unR)/2 + np.sqrt(g/h) * (EL-ER)/2
            Term=(4*R*R+x*x+y*y) / (4*R*R)
           
            Be[edges[iEdge][2],mapEdgeLeft]-=  (funformE)@(h*unstar*Term)* jac 
            Be[edges[iEdge][3],mapEdgeRight] += (funformE)@(h*unstar*Term)* jac 
            Bu[edges[iEdge][2],mapEdgeLeft] -= (funformE)@(nx*g*Term*Estar) *jac
            Bu[edges[iEdge][3],mapEdgeRight] += (funformE)@(nx*g*Term*Estar) *jac 
            Bv[edges[iEdge][2],mapEdgeLeft] -= (funformE)@(ny*g*Estar*Term )*jac
            Bv[edges[iEdge][3],mapEdgeRight] += (funformE)@(ny*g*Estar*Term)*jac       
            
#            for i in range(2):
#                for j in range(2):
#                    Be[edges[iEdge][2],mapEdgeLeft[i]]-= funformE[i][j]*h[j]*unstar[j] * ( (4*R**2+x[j]**2+y[j]**2) / (4*R**2) ) * jac * wE[j]
#                    Be[edges[iEdge][3],mapEdgeRight[i]] += funformE[i][j]*h[j]*unstar[j] * ( (4*R**2+x[j]**2+y[j]**2) / (4*R**2) ) *jac*wE[j]
#                    Bu[edges[iEdge][2],mapEdgeLeft[i]] -= funformE[i][j]*nx*g*Estar[j] * ( (4*R**2+x[j]**2+y[j]**2) / (4*R**2) ) *jac*wE[j]
#                    Bu[edges[iEdge][3],mapEdgeRight[i]] += funformE[i][j]*nx*g*Estar[j] * ( (4*R**2+x[j]**2+y[j]**2) / (4*R**2) ) *jac*wE[j]
#                    Bv[edges[iEdge][2],mapEdgeLeft[i]] -= funformE[i][j]*ny*g*Estar[j] * ( (4*R**2+x[j]**2+y[j]**2) / (4*R**2) ) *jac*wE[j]
#                    Bv[edges[iEdge][3],mapEdgeRight[i]] += funformE[i][j]*ny*g*Estar[j] * ( (4*R**2+x[j]**2+y[j]**2) / (4*R**2) ) *jac*wE[j]
        
        #print(Bu[36])
        Ainverse = np.array([[18.0,-6.0,-6.0],[-6.0,18.0,-6.0],[-6.0,-6.0,18.0]])
        for iElem in range(nElem) :
                nodes = elem[iElem]
                x = X[nodes]
                y = Y[nodes]
                jac = abs((x[0]-x[1]) * (y[0]-y[2]) - (x[0]-x[2]) * (y[0]-y[1])) 
                Bu[iElem] = Ainverse @ Bu[iElem] / jac
                Bv[iElem] = Ainverse @ Bv[iElem] / jac
                Be[iElem] = Ainverse @ Be[iElem]/ jac
        #print(Bu[36])
        U+=dt*Bu
        V+=dt*Bv
        E+=dt*Be
           
        #Sauver le fichier
        if (iter+1)==nIter:
            writeResult(theResultFiles,iter+1,E)
        elif (iter+1)%nSave==0:
            writeResult(theResultFiles,iter+1,E)
#            save=nSave
        print("Iteration = %d" % (iter+1))
        print(" == Elevations for element %d : %14.7e %14.7e %14.7e " % (27,*E[27][:]) )
        print("eta_max={0:.15f}".format(np.amax(E)))
        toc(Es = E, IsEs = True)
# 
 
 
    return [U,V,E]