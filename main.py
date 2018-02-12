from numpy import *
from initialize_struct import *
from numpy.linalg import multi_dot as multidot
from numpy.linalg import *
from scipy.optimize import lsq_linear as lsqlin

def diagz(N):
    # stupid python doesnt let u diagonalize peacefully
    if N.shape[0]==1:
        N = np.reshape(N,(1,N.shape[1]))[0]
    else:
        pass
    return diag(N)
def BIGTHETA(i,j,P):
    thetamatf = np.zeros((3,P.shape[1]))
    thetamatf[i][j] = 1
    return thetamatf

def calclagmat(N, Cb, B, Bd, Cs, gamma_hat, Minv, P, W, l_hat, m_hat):
    # function to calculate the lagrange multiplier Omega
    C = dot(Cb, P)
    A = np.ones((N.shape[0]*P.shape[1],1))
    for i in range(0, P.shape[1]):
        for j in range(0, 3):
            lambdaochat = np.ones((1,1))
            for k in range(0, Cb.shape[0]):
                lambdaochat=np.concatenate((lambdaochat,np.array([[-0.5 * B[j][k]*C[k][i]*l_hat[k][k]**-2]])),axis=1)
            lambdaochat = np.delete(lambdaochat,0,1)
            lambdaochat = diagz(lambdaochat)
            phiij = multidot([N, Cb.T, lambdaochat, Cb, Minv, P]) + multidot([BIGTHETA(j, i, P), P.T, Minv, P])
            A = np.concatenate((A,phiij.reshape((N.shape[0]*P.shape[1],1))),axis=1)
    A = np.delete(A,0,1)
    lambdaoohat = diagz(diagz(0.5*multidot([l_hat**0.5,B.T,multidot([N,Cs.T,gamma_hat,Cs])-W,Cb.T])-(1.0/12.0)*multidot([l_hat**0.5,m_hat,Bd.T,Bd])))
    phi0 = multidot([multidot([N,Cs.T,gamma_hat,Cs])-multidot([N,Cb.T,lambdaoohat,Cb])-W,Minv,P])
    B = phi0.reshape((N.shape[0]*P.shape[1],1))
    Omg = dot(linalg.pinv(A),B)

    Omg = Omg.reshape((3,Omg.size/3))
    return Omg

def calclagmat2(N, Cb, B, Bd, Cs, gamma_hat, Minv, P, W, l_hat, m_hat):
    kC = multidot([P.T,Cb.T])
    kD = multidot([Cb,Minv,P])
    kE = multidot([P.T,Minv,P])
    E = eye(3)
    S = multidot([N,Cs.T])
    kA = multidot([-1*S,gamma_hat,Cs,Minv,P])+multidot([B,diagz(diagz(
        multidot([0.5*l_hat^-2,B.T,multidot([S,gamma_hat,Cs])-W,Cb.T])
        -(1.0/12.0)*multidot([l_hat^-2,m_hat,Bd,Bd.T]))),Cb,Minv,P])\
         +multidot([W,Minv,P])
    Omg = 0
    for i in range(0,kC.shape[1]):
        temp = kC[:][i]
        Omg = Omg+(1.0/(2.0*l_hat[i][i]^2))*kron(temp.T)

def control(N, Nd, W, s0, P, D, Yt, a, b,gamma_hat0,l_hat):
    """
    :param W: External force
    :param s0: initial string length
    :param P: Constraint matrix
    :param D: Constraint matrix values
    :param Yt: Target position matrix
    :param a: control coefficient parameter
    :param b: control constant paramater
    :return: 0
    """

    B = dot(N, Cb.T)
    S = dot(N, Cs.T)
    Bd = dot(Nd, Cb.T)
    Sd = dot(Nd, Cs.T)
    Cr = 0.5 * abs(Cb)
    m_hat = m * eye(nb)
    M = (1.0 / 12.0) * multidot([Cb.T, m_hat, Cb]) + multidot([Cr.T, m_hat, Cr])
    Minv = 3.0 * multidot([Cb.T, inv(m_hat), Cb]) + 4.0 * multidot([Cr.T, inv(m_hat), Cr])
    # control to find the gamma
    EYE = eye(Cb.shape[0])
    biglambda=np.ones((1,Cs.shape[0]))
    tau=np.ones((1,1))
    Omg = calclagmat(N, Cb, B, Bd, Cs, gamma_hat0, Minv, P, W, l_hat, m_hat)
    for i in range(0,Cb.shape[0]):
        BTch = B[:,i].reshape((1,B.shape[0]))
        EYch = EYE[:,i].reshape((EYE.shape[0],1))
        biglambda=np.concatenate((biglambda,(-1.0/(2.0*l_hat[i,i]**2))*multidot([BTch,S,diagz(multidot([Cs,Cb.T,EYch]).T)])),axis=0)
        tau=np.concatenate((tau,(1.0/(2.0*l_hat[i,i]**2))*multidot([BTch,(W+dot(Omg,P.T)),Cb.T,EYch])+(1.0/(12.0*pow(l_hat[i,i],2)))*m*pow(norm(Bd[:,i]),2)),axis=0)
    # delete the first (dummy) row created
    biglambda = np.delete(biglambda,0,0)
    tau = np.delete(tau,0,0)
    Acon = a*eye(R.shape[1]) # control parameter a
    Bcon = b*eye(R.shape[1]) # control parameter b
    BTu = multidot([L,(W+dot(Omg,P.T)),Minv,R])+multidot([L,Nd,R,Acon])+dot(multidot([L,N,R])-Yt,Bcon)
    BIGTAU = np.ones((1,Cs.shape[0]))
    meu = np.ones((1,1))
    EYE = eye(R.shape[1])
    for i in range(0,R.shape[1]):
        EYch = EYE[:,i].reshape((EYE.shape[0],1))
        BIGTAU=np.concatenate((BIGTAU,dot(L,dot(S,diagz(multidot([Cs,Minv,R,EYch]).T))+multidot([B,diagz(multidot([Cb,Minv,R,EYch]).T),biglambda]))),axis=0)
        meu=np.concatenate((meu,dot(BTu,EYE[:,i])-multidot([L,B,diagz(multidot([Cb,Minv,R,EYE[:,i]])),tau])),axis=0)
    BIGTAU = np.delete(BIGTAU,0,0)
    meu = np.delete(meu,0,0)
    meu = np.reshape(meu,(1,meu.shape[0]))[0]
    # now run the optimization - lsqlin preferably.
    res = lsqlin(BIGTAU,meu,bounds=(np.zeros(ns),np.ones(ns)),method='trf',tol=1e-10,verbose=0)
    print dot(BIGTAU,res.x)-meu
    return res


gamma = control(N, Nd, W, s0, P, D, Yt, a, b,gamma_hat0,l_hat)
print gamma.x