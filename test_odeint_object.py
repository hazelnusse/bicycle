from scipy import *
from pylab import *
from scipy import integrate

# define a class with a member function that returns the derivative
class test:
    def deriv(self,y,t):
        return array([y[1],-y[0]])

obj = test()

# integration parameters
start=0
end=10
numsteps=10000
time=linspace(start,end,numsteps)
y0=array([0.0005,0.2])

# integrate the system
y=integrate.odeint(obj.deriv,y0,time)
plot(time,y[:,0])
show()
