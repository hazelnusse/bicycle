from numpy import zeros, sin, pi

class DynamicSystem:
    """Dynamic System class.

    """

    def __init__(self):

        # parameter names and their values
        self.parameters = {'g':9.81,
                           'm':1.,
                           'l':1.}
        # state names and their initial conditions
        self.states = {'theta':0.,
                       'omega':0.,
                       'thetap':0.,
                       'omegap':0.}
        # numerical integration parameters
        self.numint = {'ti':0.,
                       'tf':10.}

    def f(self, x, t, params):
        g, l, m, z = params

        theta = x[0]
        omega = x[1]

        torque = inputs(t)[0]

        z[0] = sin(theta)

        thetap = omega
        omegap = -g/l*z[0] + torque/(m*l*l)

        f = zeros(2)
        f[0] = thetap
        f[1] = omegap

        params[3] = z[0]
        return f

    def inputs(self, t):
        print pi
        print t
        torque = 10*sin(2*pi*t + pi/6)
        print torque
        inputs = zeros(1)
        inputs[0] = torque
        return inputs
