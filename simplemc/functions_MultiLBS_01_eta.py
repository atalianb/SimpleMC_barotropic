import numpy as np
from scipy import optimize
from scipy import integrate
from scipy import interpolate
########
G_kpc = 4.302e-6#kpc/SolarMass(km/s)^2
####
##Integration values
###
x0_0 = 2.#first integration
x0_0v1 = 0.0# Start of integration. Use 0.01 for continuity in l=3 solution, for 0,1,2 0.0 it's ok
xf_0v1 = 10.5# End of integration
step_0 = 0.2#Step to integrate from x0_0 until xf_0v1 is reached
h = 0.01
######
#Shooting values
####
u0 = -3.90026683
u1 = -2.81979045#-0.5#1.2# 1st trial value of unknown init. cond. gamma_100
u2 = -2.13270971#-0.4#2.7# 2nd trial value of unknown init. cond. gamma_210
#u3 = 4.65#-1.9#4.9# 3th trial value of unknown init. cond. gamma_320
u4 = 1.
u5 = 1.
u_012_V = np.array([u5,u0,u1,u4,u2])#,u3])
#####
#data_path = "simplemc/data/Blok_McGaugh_&_Rubin_(2001)/"
#data = np.loadtxt(data_path+'U11557.dat')
#vecRp_data = np.array([row[1] for row in data])# galactocentric distance [kpc]
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
    while x0 < x:
         # Count number of iterations using step size or
        h = min(h,x - x0)# step height h
        y0 = y0 + RK4(func,x0,y0,h)#update next value of y
        x0 = x0 + h#update next value of x
        X.append(x0)
        Y.append(y0)
    return np.array(X),np.array(Y)
def shooting01_V(func,u,x0,x,xf,step,eta,h=0.01):
    def IC_V(u):
        return np.array([u[0],0.,u[1],0.,u[2],0.,u[3],0.,u[4],0.])
    def res(u):# Boundary condition residual
        X,Y = Integrate(func,x0,IC_V(u),x,h)
        y = Y[len(Y) - 1]#last value of Y
        r = np.zeros(len(u))
        r[0] = y[0]#y0(inf)=0
        r[1] = y[2]/x + y[3]#y_3(inf)/r + y_4(inf)=0
        r[2] = y[6]#
        r[3] = 1.5 - y[5] - y[9]
        r[4] = y[9] - eta*y[5]
        #r[3] = y[10]#y6(inf)=0
        #r[1] = y[1] + np.sqrt(2.*y[4]**2.)*y[0]#y_2(inf)+sqrt(2y_5^2)y_1(inf)=0
        return r
    while x<=xf:
        root = optimize.root(res,u)
        #print(root.x)
        u = root.x
        root_temp = optimize.root(res,root.x)
        X,Y = Integrate(func,x0,IC_V(root_temp.x),x,h)
        x = x+step
    return X,Y,root_temp
########
def rhs(x,y):
    dy = np.zeros(10)
    dy[0] = y[1]
    dy[2] = y[3]
    dy[4] = 0.
    dy[5] = y[0]**2.*x**2.
    dy[6] = y[7]
    dy[8] = 0.
    dy[9] = y[6]**2.*x**2.
    if(x==0.):
        dy[1] = 2.*(y[2]-y[4])*y[0]/3.
        dy[3] = y[0]**2./3.+ 3.*x**2.*y[6]**2./5.
        dy[7] = 2.*(y[2] - y[8])*y[6]/5.
    else:
        dy[1] = 2.*(y[2]-y[4])*y[0] - 2.*y[1]/x
        dy[3] = y[0]**2. + 3.*x*2.*y[6]**2. - 2.*y[3]/x
        dy[7] = 2.*(y[2] -y[8])*y[6] -4.*y[7]/x
    return dy
#######
def Mass_func(r,phi,l):
    Int = np.zeros(len(r))
    dr = np.diff(r)[0]
    if l==0.:
        phi_array = np.array(phi[:,0])
    if l==1.:
        phi_array = np.array(phi[:,6])
    for i in range(0,len(r)-1):
        Int[i+1] = dr*(phi_array[i+1]**2.*r[i+1]**(2.*l+2.)) + Int[i]
    return Int
#####
#Units
def Vc2_cir(r,eps,M):
    units =8.95e10*eps**2.
    return (units*M)/r
####
##Units for r in kpc
###
def r_units(r,eps,m_a):
    return (6.39e-27*r)/(eps*m_a)
#####
def Vc_(r,m_a,eps,eta):
    X,Y,root_f=shooting01_V(rhs,u_012_V,x0_0v1,x0_0,xf_0v1,step_0,eta)
    M_r0 = Mass_func(X,Y,l=0.)+Mass_func(X,Y,l=1.)#Integrates rho(r) to obtain M(r)
    Vc2_r0 = Vc2_cir(X,eps,M_r0)#Vc^2[km/s]^2 theoretical
    X0_units = r_units(X,eps,m_a)#r[kpc] theoretical
    M_r0_units = M_r0*eps*1.34e-10/m_a
    if X0_units[-1]<r[-1]:
        #array from last element of the r[kpc] theoretical to the last element of the data array,
        # with 80 elements. It can be replaced by np.arange(X0_units[-1],vecRp_data[-1],0.1) 
        #but you have to be careful in the next function with interpolate
        r_array = np.linspace(X0_units[-1],r[-1],80)
        Vc2_rmayor = G_kpc*M_r0_units[-1]/r_array#Computes Vc^2 with with the last result from M(r)
        Vc2_total = np.append(Vc2_r0,Vc2_rmayor)#creates an array of Vc^2 with Vc2_r0 and Vc2_rmayor
        r_total = np.append(X0_units,r_array)
        return r_total,Vc2_total
    else:
        return X0_units,Vc2_r0
def Vc_inter(r,m_a,eps,eta):
    Vc = Vc_(r,m_a,eps,eta)
    #If you want to use np.arange in the previous function, It is recommended to use extrapolate
    f = interpolate.interp1d(Vc[0],Vc[1],fill_value='extrapolate')
    Vc_new = f(r)
    return Vc_new
