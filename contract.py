from scipy.sparse import lil_matrix,eye
from scipy.sparse.linalg import spsolve as solve
from numpy import *

##define default parameters
ngp=256

def setup_grid():
    """Setup finite difference grid in alpha space, i.e. equispaced within the
    contracting band.
    """
    global x,delta,dx,dx2,one
    x,delta=linspace(0.,1.,ngp,retstep=True)
    one=eye(ngp)
    dx=lil_matrix((ngp,ngp))
    dx2=lil_matrix((ngp,ngp))
    for i in xrange(ngp):
        dx[i,(i+1)%ngp]=1.
        dx[i,i-1]=-1.
        #####
        dx2[i,(i+1)%ngp]=1.
        dx2[i,i-1]=1.
        dx2[i,i]=-2.
    dx=dx/delta/2.
    dx2=dx2/delta**2
setup_grid()

def lhs(w):
    """determine left hand side of force balance equation"""
    #bulk
    out=4.*l**2*dx2/w**2-one
    #bc at center
    out[0,0]=1.
    out[0,1]=0.
    out[0,-1]=0.
    ##force balance at x=+-w/2
    out[-1,-1]=3./(delta*w)*l**2
    out[-1,-2]=-4./(delta*w)*l**2
    out[-1,-3]=1./(delta*w)*l**2
    out[-1,0]=0
    return out

#right hand side force balance
zeta=lambda rho : rho*(rho-rho0)

def rhs(w,rho):
    """determine right hand side of force balance equation"""
    #bulk
    zz=zeta(rho)
    out=-s0*l*dx*zz*2/w
    #bc
    out[0]=0.
    out[-1]=-s0*l*zz[-1]
    return out

#compute velocities
def get_v(w,rho):
    """compute velocities"""
    ll=lhs(w)
    rr=rhs(w,rho)
    return solve(ll,rr)

#compute time derivatives
def get_dt(w,rho):
    """compute time derivatives"""
    v=get_v(w,rho)
    dtw=2.*v[-1]

    dtrho=-2./w*dx*((v-dtw*x/2.)*rho)-dtw/w*rho
    #flux is odd at center! bc
    dtrho[0]=-2*(v[1]-dtw*x[1]/2.)*rho[1]/(w*delta)-dtw/w*rho[0]
    ##2nd order approx
    dtrho[-1]=-(3.*(v[-1]-dtw*x[-1]/2.)*rho[-1]-4.*(v[-2]
-dtw*x[-2]/2.)*rho[-2]+(v[-3]-dtw*x[-3]/2.)*rho[-3])/(delta*w)

    dtrho[-1]=dtrho[-1]-dtw/w*rho[-1]
    return dtw,dtrho

def dt_packed(y,t):
    """return time derivatives formatted for odeint"""
    w=y[0]
    rho=y[1 : ngp+1]
    dtw,dtrho=get_dt(w,rho)
    out=zeros(ngp+1)
    out[0]=dtw
    out[1 : ngp+1]=dtrho
    return out

def init_packed():
    """return initial state formatted for odeint"""
    init=zeros(ngp+1)
    init[0]=w_init
    init[1 : ngp+1]=rho_init
    return init

def compute(_l,_s0,_r0,_w_init,t=linspace(0,80.,2000)):
    """compute a time course using odeint"""
    from scipy.integrate import odeint
    global l,s0,rho0,rho_init,w_init
    #set parameters
    l=_l
    s0=_s0
    rho0=1.
    rho_init=_r0*ones(ngp)
    w_init=_w_init
    return odeint(dt_packed,init_packed(),t)

if __name__=='__main__':
    #set default values
    print "running"
    a,b,r0=array([  2.26413429e+00,   7.09505347e-05,   3.26194452e-01]) #
    a=a/3.
    b=b/3.
    b=b/4.
    l=sqrt(a/b)*pi/r0
    s0=-1./sqrt(a*b)*pi/r0
    w_init=2100
    print("starting computation")
    o=compute(l,s0,r0,w_init)
    print("done")
    print o[:,0]
