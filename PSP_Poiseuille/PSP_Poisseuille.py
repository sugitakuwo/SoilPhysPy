#PSP_Poisseulle.py
from vpython import scene, curve, arrow, vector
import numpy as np

scene.background=vector(1,1,1)
scene.title = ""

l = 5.     
R = 1.     
D_P = 40.   
eta = 1.   

angles=arange(1.1*pi,2.1*pi,pi/20.)
n = 500
for i in range(n):
 spring=curve(color=vector(0,1,1), radius=0.06)
 for phi in angles:
     spring.append(pos=vector(l*(float(i)/float(n)-0.55), R*cos(phi), R*sin(phi)))

for i in range(11):
  x = 0
  y = (float(i)/11.-0.5)*2.*R
  r = np.sqrt(x*x+y*y)
  if (r < R):
   arrow(pos=vector(0,y,x),axis = vector(D_P * (R*R-r*r)/(eta*l*4.),0,0),
                       shaftwidth = 0.035, color=vector(0,0,0))