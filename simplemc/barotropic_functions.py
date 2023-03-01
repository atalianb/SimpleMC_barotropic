import numpy as np
#from scipy import interpolate
###
##Constants
####
Gm=2.5188027*10**(-3)
Gp=0.499246*10**(-6)
G=4.258*10**(4)
Gv=1.24*10**(5)
y0=np.zeros(2)
h=0.01
x0=0.00001
###
##Initial condition
y0[0]=0.#no cambia,
######
def Integrate(func,x0,y0,x,h):
    # Finds value of y for a given x using step size h 
    # and initial value y0 at x0.
    def RK4(func,x0,y0,h):
        K0 = h*func(x0,y0)
        K1 = h*func(x0 + 0.5*h, y0 + 0.5*K0)
        K2 = h*func(x0 + 0.5*h, y0 + 0.5*K1)
        K3 = h*func(x0 + h, y0 + K2)
        return (K0 + 2.*K1 + 2.*K2 + K3)/6.
    X = []
    Y = []
    X.append(x0)
    Y.append(y0)
    while y0[0] < 4.90000000000:
    #while y0[0]>0.99:
        y0 = y0 + RK4(func,x0,y0,h) #update next value of y
        x0 = x0 + h #update next value of x
        X.append(x0)
        Y.append(y0)
    return np.array(X),np.array(Y)
######
def Vc_(kappa,rho0):
    #ec de estado P=kap*rho
    K = kappa
    y0[1] = K*rho0
    def rhs(x,y):
        dy=np.zeros(2)
        dy[0]=4*np.pi*Gm*x**2*(y[1]/K)
        dy[1]=-Gp*y[0]/x**2*(y[1]/K)
        return dy
    X,Y=Integrate(rhs,x0,y0,10.,h)
    v=np.sqrt(Y[:,0]/X)
    RC=v*np.sqrt(G)
    return X,RC
