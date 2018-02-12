from sensordata import *
from math import *
import numpy as np
from numpy import dot,diag

"""
# get nodes position and node velocity from the sensordata
N0 = getnodepos()
Nd0 = getnodevel()
"""
# control the height of the prism only.
L = np.array([[0,0,1]])
R = np.array([[0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1]])

# tensegrity prism
Cb = np.array([[-1,0,0,0,1,0],[0,-1,0,0,0,1],[0,0,-1,1,0,0]])
Cs = np.array([[0,0,0,0,-1,1],[0,0,0,1,0,-1],[0,0,0,-1,1,0],[0,-1,0,0,1,0],[0,0,-1,0,0,1],[-1,0,0,1,0,0],[0,0,-1,0,1,0],[-1,0,0,0,0,1],[0,-1,0,1,0,0]])
bar_len = 1.4980      # each bar is 60 cm in length
b0 = bar_len*np.eye(3)

m = 1
nb = 3
ns = 9
th = 30

# example node position and node velocity to test the control code.

N = np.array([[ -0.5000,0,0.5000,-0.5774,0.2887,0.2887],
   [-0.2887,0.5774,-0.2887,0,0.5000,-0.5000],
         [0,0,0,1.0000,1.0000,1.0000]])
Nd = np.zeros((3,6))
W = np.zeros((3,6))
B = dot(N,Cb.T)
l_hat = diag(diag(dot(B.T,B)))**0.5
gamma_hat0 = np.zeros((9,9))
S = dot(N,Cs.T)
s0 = np.sqrt(diag(dot(S.T,S)))
P = np.concatenate((np.eye(3),np.zeros((3,3))),0)
D = N[:,0:3]
Yt = 1.2*np.ones((1,3))
a = 2
b=2
