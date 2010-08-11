
from numpy import zeros, sin, pi

def pendulum_torque(t):
    return 0 #10*sin(2*pi*t + pi/6)

def pendulum_f(x, t, params):
    g, l, m, z = params
    theta = x[0]
    thetap = x[1]

    f = zeros(2)
    z = sin(theta)
    f[0] = thetap
    f[1] = -g/l*z + pendulum_torque(t) / (m*l*l)
    #f[1] = -g/l*sin(theta) + pendulum_torque(t) / (m*l*l)
    params[3] = z
    return f


def pendulum_outputs(x, g, l, m, z):
    
    ke = 1/2*m*(l*thetap)**2
    return ke

