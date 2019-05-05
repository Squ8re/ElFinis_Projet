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
# -------------------------------------------------------------------------
def tic(message = ''):
    """Cette fonction lance un chronomètre. L'argument n'est pas utile.
    """
    global startTime
    startTime = timer()

def toc(message = ''):
    """Cette fonction arrête le chronomètre lancé avec tic() et imprime le temps
       écoulé sur la sortie standart.
       @pre: startTime doit être initialisé.
    """
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
    """Cette fonction lit le fichier et sépare les différents éléments.
       @return: La fonction retourne un tableau contenant, dans l'ordre
                   - Le nombre de noeuds;
                   - les coordonnées des points en x;
                   - les coordonnées des points en y;
                   - les coordonnées des points en z (la profondeur);
                   - le nombre d'éléments;
                   - les éléments.
    """
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
    """Cette fonction lit un fichier contenant les résultats d'une itération et 
       renvoie les résultats sous forme de tableau.
       @return: tableau numpy contenant les infos sur chaque élément.
    """
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
    """Cette fonction écrit un fichier contenant les résultats d'une itération.
    """
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
    """Cette fonction retrie le maillage et supprime les doublons. Une fois
       appliquée, cette fonction permet de repérer les éléments frontières et
       de dissocier ces derniers des autres.
       @return: [nombre de sommets, nombre de frontières, sommets]
    """
    [nNode,X,Y,H,nElem,elem] = readMesh(mesh)                                   # Récupération du maillage.
    nEdges = nElem * 3                                                          # Calcul du nombre de sommets.
    nBoundary = 0                                                               # Initialisation du nombre de frontières.
    edges = [[0 for i in range(4)] for i in range(nEdges)]                      # Initialisation des sommets à 0.
    for i in range (nElem) :
        for j in range(3) :
            id = i*3 + j 
            edges[id][0] = elem[i][j]
            edges[id][1] = elem[i][(j+1)%3]
            edges[id][2] = i
            edges[id][3] = -1
    edges.sort(key = lambda item : -(min(item[0:2])*nEdges)-max(item[0:2]))     # Tri des éléments (voir cours).
    index = 0                                                                   # Initialisation de l'index.
    for i in range(nEdges) :
        if (edges[i][0:2] != edges[i-1][1::-1]) :
            edges[index] = edges[i]
            index += 1
        else :
            edges[index-1][3] = edges[i][2]                                                # On place l'élément juste au dessus
    del edges[index:]                                                           # On supprime tous les éléments en dessous de l'index.
    edges.sort(key = lambda item : item[3])                                     # On trie à nouveau les éléments
    nBoundary = 2*index - nEdges                                                # On trouve le nombre de frontières
    nEdges = index                                                              # On met à jour le nombre de sommets.
    
    return nEdges, nBoundary, edges                                             # On retourne le nombre de sommet, le nombre de frontières
                                                                                # et les sommets.

def mapEdge(theMeshFile, elem): #vient du prof
    """Cette fonction crée une carte reprenant la liste des noeuds.
    """
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

def addIntegralsTriangles(X, Y, H, Ei, Ui, Vi, nElem, elem):
    """
    Cette fonction calcul les intégrales sur les triangles et les ajoute aux termes Ei, Ui et Vi
    @pre: Les contiennent les informations, dans l'ordre
            - X contient un vecteur contenant les abscisses de chaque élément;
            - Y contient un vecteur contenant les ordonnées de chaque élément;
            - H contient un vecteur contenant la hauteur de chaque élément;
    """
    R = 6371220                                                                 # Rayon de la terre.
    gamma = 10**(-7)                                                            # Facteur dissipatif proportionnel à la vitesse.
    g = 9.81                                                                    # Gravité.
    #TODO: Vérifier que les Xsi et les Eta sont bons.
    Xsi    = np.array([0.166666666666667,0.666666666666667,0.166666666666667])  # Xsi prédéfini.
    Eta    = np.array([0.166666666666667,0.166666666666667,0.666666666666667])  # Eta prédéfini.
    weight = 0.166666666666667                                                  # Poids prédéfini.
    
    dphidxsi = np.array([ 1.0, 0.0,-1.0])
    dphideta = np.array([ 0.0, 1.0,-1.0])
    
    for iElem in range(nElem):                                                  # Pour tous les éléments...
        nodes = elem[iElem]                                                     # Prendre tous les noeuds d'un élément
        
        x = X[nodes]                                                            # Récupérer les abscisses associées
        y = Y[nodes]                                                            # Récupérer les ordonnées associées
        
        
        dxdxsi = x @ dphidxsi                                                   # On effectue le calcul intégral 
                                                                                #   $ x\frac{\partial phi}{\partial xsi} $ 
        dxdeta = x @ dphideta                                                   # On effectue le calcul intégral 
                                                                                #   $ x\frac{\partial phi}{\partial eta} $ 
        dydxsi = y @ dphidxsi                                                   # On effectue le calcul intégral 
                                                                                #   $ y\frac{\partial phi}{\partial xsi} $ 
        dydeta = y @ dphideta                                                   # On effectue le calcul intégral 
                                                                                #   $ y\frac{\partial phi}{\partial eta} $
        jac = abs(dxdxsi*dydeta - dxdeta*dydxsi)                                # Calcul du Jacobien
        
        dphidx = (dphidxsi * dydeta - dphideta * dydxsi) / jac                  # Calcul de $ \frac{\partial phi}{\partial x} $
        dphidy = (dphideta * dxdxsi - dphidxsi * dxdeta) / jac                  # Calcul de $ \frac{\partial phi}{\partial y} $
        
        for j in range(3):
            xsi = Xsi[j]                                                        # Prendre le Xsi
            eta = Eta[j]                                                        # Prendre le Eta
            phi = [(1-xsi-eta), xsi, eta]                                       # Créer le phi
            u = Ui[iElem,0]*(1-xsi-eta)+Ui[iElem,1]*(xsi)+Ui[iElem,2]*(eta)     # Trouver u
            v = Vi[iElem,0]*(1-xsi-eta)+Vi[iElem,1]*(xsi)+Vi[iElem,2]*(eta)     # Trouver v
            e = Ei[iElem,0]*(1-xsi-eta)+Ei[iElem,1]*(xsi)+Ei[iElem,2]*(eta)     # TODO: Où est e?
            h = H[nodes[0]]*(1-xsi-eta)+H[nodes[1]]*(xsi)+H[nodes[2]]*(eta)     # Trouver h
            xi = x[0]*(1-xsi-eta)+x[1]*(xsi)+x[2]*(eta)                         # Calcul de x_i
            yi = y[0]*(1-xsi-eta)+y[1]*(xsi)+y[2]*(eta)                         # Calcul de y_i
            z3d = R*(4*R*R - xi*xi - yi*yi) / (4*R*R + xi*xi + yi*yi)           # Calcul de z (la hauteur en 3D)
            #TODO: Vérifier que f est correct (n'est pas comme la formule donnée dans l'énoncé?)
            f = (4*np.pi*z3d) / (R*86400)                                       # Calcul du terme de Coriolis
            p1 = (4*R*R + xi*xi + yi*yi) / (4*R*R)                              # TODO: trouver le sens de p1, p2, p3
            p2 = (h*(xi*u+yi*v)) / (R**2)
            p3 = (g*xi*e) / (2*(R**2))
            
            for k in range(3):
                coefIntegr = jac*weight                                         # Calcul de coéfficient d'intégration
                Ei[iElem, k] += coefIntegr* (h*p1*(u * dphidx[k] + v * dphidy[k]) + phi[k]*p2)
                Ui[iElem, k] += coefIntegr * (phi[k]*(f*v-u*gamma) + dphidx[k]*g*e*p1 + phi[k]*p3)
                Vi[iElem, k] += coefIntegr * (phi[k]*(-f*v-u*gamma) + dphidy[k]*g*e*p1 + phi[k]*p3)
                
    return Ei, Ui, Vi
                


def addIntegralsEdges(X, Y, H, E, U, V, nEdges, nBoundary, edges, mapEdgeLeft, mapEdgeRight):
    R = 6371220
    g = 9.81
    #TODO: Vérifier que les Xsi et les Eta sont bons.
    Xsi = np.array([-0.5773502691896257, 0.5773502691896257])
    
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
        
        for j in range(2):
            xsi = Xsi[j]
            phi = [(1-xsi)/2,(1+xsi)/2]
            eL = E[elemL, nodeLeft[0]]*(1-xsi)/2+E[elemL, nodeLeft[1]]*(1+xsi)/2
            uL = U[elemL, nodeLeft[0]]*(1-xsi)/2+U[elemL, nodeLeft[1]]*(1+xsi)/2
            vL = V[elemL, nodeLeft[0]]*(1-xsi)/2+V[elemL, nodeLeft[1]]*(1+xsi)/2
            xi = x[0]*phi[0]+x[1]*phi[1]
            yi = y[0]*phi[0]+y[1]*phi[1]
            h = H[nodes[0]]*phi[0] + H[nodes[1]]*phi[1]
            p1 = (4*R*R + xi**2 + yi**2) / (4*R*R)
            unL = nx*uL + ny*vL
            unR = -unL
            etaEt = eL + np.sqrt(h/g) * (-unL)
            
            for k in range(2):
                U[elemL, k] -= jac * (phi[k]*nx*g*etaEt*p1)
                V[elemL, k] -= jac * (phi[k]*ny*g*etaEt*p1)
        
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
            eL = E[elemL, nodeLeft[0]] *(1-xsi)/2+E[elemL, nodeLeft[1]] *(1+xsi)/2
            eR = E[elemR, nodeRight[0]]*(1-xsi)/2+E[elemR, nodeRight[1]]*(1+xsi)/2
            uL = U[elemL, nodeLeft[0]] *(1-xsi)/2+U[elemL, nodeLeft[1]] *(1+xsi)/2
            uR = U[elemR, nodeRight[0]]*(1-xsi)/2+U[elemR, nodeRight[1]]*(1+xsi)/2
            vL = V[elemL, nodeLeft[0]] *(1-xsi)/2+V[elemL, nodeLeft[1]] *(1+xsi)/2
            vR = V[elemR, nodeRight[0]]*(1-xsi)/2+V[elemR, nodeRight[1]]*(1+xsi)/2
            xi = x[0]*phi[0]+x[1]*phi[1]
            yi = y[0]*phi[0]+y[1]*phi[1]
            h = H[nodes[0]]*phi[0] + H[nodes[1]]*phi[1]
            p1 = (4*R*R + xi**2 + yi**2) / (4*R*R)
            unL = uL*nx + vL*ny
            unR = uR*nx + vR*ny        
            
            etaEt = (eL+eR)/2 + np.sqrt(h/g) * (unR-unL)/2
            unEt  = (unL+unR)/2 + np.sqrt(g/h) * (eR-eL)/2 #peut etre switch entre eL et eR (de base)
            
            for k in range(2):
                E[elemL, k] -= jac * (phi[k]*h*unEt*p1)
                E[elemR, k] += jac * (phi[k]*h*unEt*p1)
                U[elemL, k] -= jac * (phi[k]*nx*g*etaEt*p1)
                U[elemR, k] += jac * (phi[k]*nx*g*etaEt*p1)
                V[elemL, k] -= jac * (phi[k]*ny*g*etaEt*p1)
                V[elemR, k] += jac * (phi[k]*ny*g*etaEt*p1)
                

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
    save = 0
    
    for i in range(nIter):
        tic()
        #TODO: Imposer les conditions initiales.
        nNode,X,Y,H,nElem,elem = readMesh(theMeshFile)
        nEdges, nBoundary, edges = ComputeEdges(theMeshFile)
        mapEdgeLeft, mapEdgeRight = mapEdge(theMeshFile, elem)
        print("This is X")
        print(H)
        Ei = np.zeros(E.shape)
        Ui = np.zeros(U.shape)
        Vi = np.zeros(V.shape)
        print("Ei before")
        print(Ei)
        Ei, Ui, Vi = addIntegralsTriangles(X, Y, H, Ei, Ui, Vi, nElem, elem)
        print("Ei after")
        print(Ei)
        addIntegralsEdges(X, Y, H, Ei, Ui, Vi, nEdges, nBoundary, edges, mapEdgeLeft, mapEdgeRight)
        multiplyInverseMatrix(elem, nElem, Ei, Ui, Vi, X, Y)
        
        E += dt*Ei
        U += dt*Ui
        V += dt*Vi
        if save==nSave:  
            save=0
            writeResult(theResultFiles ,i,E)
        
        save += 1
        toc()
    
    return [U,V,E]