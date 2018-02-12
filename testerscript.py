from scipy.optimize import lsq_linear as lsqlin
import numpy as np
from numpy import eye
from scipy import sparse

C = np.array([[0.9501,0.7620,0.6153,0.4057],[0.2311,0.4564,0.7919,0.9354],[0.6068,0.0185,0.9218,0.9169],[0.4859,0.8214,0.7382,0.4102],[0.8912,0.4447,0.1762,0.8936]])

d = np.array([[0.0578],[0.3528],[0.8131],[0.0098],[0.1388]])
A = np.array([[0.2027,0.2721,0.7467,0.4659],[0.1987,0.1988,0.4450,0.4186],[0.6037,0.0152,0.9318,0.8462]])
b = np.array([0.5251,0.2026,0.6721])
Aeq = np.array([[3,5,7,9]])
beq = np.array([[4]])

lb = -0.1*np.ones(4)
ub = 2*np.ones(4)

res = lsqlin(A,b,(lb,ub),tol=1e-10,lsmr_tol='auto',verbose=0)

print res.x
