#
# PYTHON for FEM DUMMIES 18-19
# Projet "tsunami"
#
# Visualisation des résultats : OpenGL old-fashion
# Une implémentation avec des shaders sera présentée en S11
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 


from OpenGL.GL   import *
from OpenGL.GLU  import *
from OpenGL.GLUT import *

import numpy as np
import tsunami
import sys

# -------------------------------------------------------------------------
print('coucou depuis tsunamiAnimate.py')

def draw():  
  global E,TransSide,TransVert,theFlagBathymetry,theMouseSide,theMouseVertical,theRatio, EdgeFLAG

  glClearColor( 0.9, 0.9, 0.8, 0.0 );
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(65.0,theRatio,1.0,100.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0.0,1.0,0.0,0.0,20.0,0.0,0.0,0.0,1.0);  
  glTranslatef(0.0,14.0,0.0);
  glTranslatef(0.3*TransSide,0.0,0.0);
  glTranslatef(0.0,0.0,0.3*TransVert);
  glRotatef(0.3*theMouseSide,0.0,0.0,1.0);
  glRotatef(0.3*theMouseVertical, 1.0, 0.0, 0.0);
  
  quadratic = gluNewQuadric();         
  gluQuadricNormals(quadratic, GLU_SMOOTH); 
  glColor3f(0.2,1.0,0.3);
  gluSphere(quadratic,5.95,400,200);
 
  n = 9*nElem
  index  = np.ravel(elem)
  colors = np.zeros(n) 
  coord  = np.zeros(n)

  if (theFlagBathymetry == False) :
    value = 10* np.ravel(E)
  else :
    value = H[index]/BathMax
    
  value = np.clip(value,0,1)        
  colors[0:n:3] = 3.5*(value)*(value)
  colors[1:n:3] = (1-value)*(value)*3.5
  colors[2:n:3] = (1-value)*(1-value)
  
  x = X[index]
  y = Y[index]
  factor = (4*R*R + x*x + y*y)*(R/6)
  coord[0:n:3] = 4*R*R * x / factor
  coord[1:n:3] = 4*R*R * y / factor
  coord[2:n:3] = (4*R*R - x*x - y*y)*R / factor   
 
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(3,GL_FLOAT,0,coord);
  glNormalPointer(GL_FLOAT,0,coord);
  glColorPointer(3,GL_FLOAT,0,colors);
  glDrawArrays(GL_TRIANGLES,0,nElem*3);
  glDisableClientState(GL_NORMAL_ARRAY);    
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);      
    
  coord *= 1.001;

  glColor3f(0.0, 0.0, 0.0);
  glEnableClientState(GL_VERTEX_ARRAY);    
  glVertexPointer(3,GL_FLOAT,0,coord);
  if(EdgeFLAG):
    for i in range(nElem):
      glDrawArrays(GL_LINE_LOOP,3*i,3)
  glDisableClientState(GL_VERTEX_ARRAY) 
  glutSwapBuffers() 
 
# -------------------------------------------------------------------------
    
def reshape(width, height):   
  global theRatio
  
  glViewport(0,0,width,height)
  theRatio = width / float(height)
  glutPostRedisplay()

# -------------------------------------------------------------------------
 
def keyboard(key,x,y):
  global theFlagBathymetry, iter, PauseFlag, EdgeFLAG, TransSide, TransVert
  
  key = key.decode()
  if ord(key) == 27: # Escape
    sys.exit(0)
  elif key == 'b':
    theFlagBathymetry = True
  elif key == 'e':
    theFlagBathymetry = False
  elif key == 'r':
    iter = 0
  elif key == 'p':
    PauseFlag = not PauseFlag
  elif key == 's':
    EdgeFLAG = not EdgeFLAG
  elif key == 'm':
    TransSide += 0.5
  elif key == 'k':
    TransSide -= 0.5
  elif key == 'o':
    TransVert += 0.5
  elif key == 'l':
    TransVert -= 0.5
  else:
    return
  glutPostRedisplay()

# -------------------------------------------------------------------------

def special(symbol,x,y):
  global theMouseSide, theMouseVertical
  
  if symbol == GLUT_KEY_UP :        
    theMouseVertical -= 5
  elif symbol == GLUT_KEY_DOWN :
    theMouseVertical += 5
  elif symbol == GLUT_KEY_LEFT :
    theMouseSide -= 5
  elif symbol == GLUT_KEY_RIGHT :
    theMouseSide += 5
  else:
    return
  glutPostRedisplay()
  
# -------------------------------------------------------------------------

def idle():
  global iter,delta,E,theResultFiles,PauseFlag
  
  if(not PauseFlag):
    iter += delta 
    try :
      E = tsunami.readResult(theResultFiles ,iter,nElem)
      glutPostRedisplay()
    except FileNotFoundError: 
      pass
  elif(PauseFlag):
      pass
  else:
      pass

# -------------------------------------------------------------------------
  
iter = 0; delta = 1;
R = 6371220;
BathMax = 9368;
theMeshFile = "PacificTriangleTiny.txt"
#theResultFiles = "results/eta-%06d.txt" #===========================
theResultFiles = "nosResult/eta-%06d.txt"
theFlagBathymetry = False
PauseFlag = False
EdgeFLAG = True
theMouseSide = 389
theMouseVertical = 0
theRatio = 1.0
TransSide = 1.0
TransVert = 1.0

glutInit(sys.argv)
glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH)
glutInitWindowPosition(0, 0)
glutInitWindowSize(1000, 1000)                                # Initialement 500, 500
glutCreateWindow("MECA1120 : the 2019 project :-)")

matSpecular   = [1.0,1.0,1.0,0.0]
matShininess  = [50.0]
lightPosition = [8.0,8.0,8.0,0.0]
lightRadiance = [1.0,1.0,1.0,1.0]
glMaterialfv(GL_FRONT,GL_SPECULAR, matSpecular)
glMaterialfv(GL_FRONT,GL_SHININESS,matShininess)
glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE)
glLightfv(GL_LIGHT0,GL_POSITION,lightPosition)
glLightfv(GL_LIGHT0,GL_DIFFUSE, lightRadiance)
glLightfv(GL_LIGHT0,GL_SPECULAR,lightRadiance)
glEnable(GL_LIGHTING)
glEnable(GL_LIGHT0)
glDepthFunc(GL_LEQUAL)
glEnable(GL_DEPTH_TEST)
glEnable(GL_COLOR_MATERIAL)
glEnable(GL_NORMALIZE)	

glutDisplayFunc(draw)
glutKeyboardFunc(keyboard)
glutSpecialFunc(special)
glutReshapeFunc(reshape)
glutIdleFunc(idle)

glClearColor( 0.9, 0.9, 0.8, 0.0 );
glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
glMatrixMode(GL_PROJECTION);
glLoadIdentity();
gluPerspective(65.0,1.0,1.0,100.0);

if "-info" in sys.argv:
  print("GL_RENDERER   = ",glGetString(GL_RENDERER).decode())
  print("GL_VERSION    = ",glGetString(GL_VERSION).decode())
  print("GL_VENDOR     = ",glGetString(GL_VENDOR).decode())

print('===========================================')  
print(' b           : show bathymetry ')
print(' e           : show elevation (by default) ')
print(' UP/DOWN     \_ rotate the Earth ')
print(' RIGHT/LEFT  / ')
print(' r           : restart ')
print(' p           : pause ')
print(' s           : show/hide edges')
print(' ESC         : exit ')
print('===========================================')
 
[nNode,X,Y,H,nElem,elem] = tsunami.readMesh(theMeshFile)
print(nElem);                                                                   #TODO: Delete this.
try :
  E = tsunami.readResult(theResultFiles,0,nElem)
except FileNotFoundError:
  E = np.zeros([nElem,3])

glutMainLoop()

# -------------------------------------------------------------------------
