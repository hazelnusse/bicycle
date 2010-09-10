
from numpy import zeros, sin, pi

def pendulum_torque(t):
    return 0 #10*sin(2*pi*t + pi/6)

def pendulum_f(x, t, params):
    g, l, m, z = params

    theta = x[0]
    omega = x[1]

    z[0] = sin(theta)
    
    thetap = omega
    omegap = -g/l*z[0] + torque/(m*l*l)
    
    f = zeros(2)
    f[0] = thetap
    f[1] = omegap 

    params[3] = z[0]
    return f


def pendulum_outputs(x, g, l, m, z):
    
    ke = 1/2*m*(l*thetap)**2
    return ke

