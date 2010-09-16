from numpy import zeros, sin, pi
from numpy import linspace, array
from scipy.integrate import odeint
from matplotlib.pyplot import figure, plot, show, legend, xlabel

class DynamicSystem:
    """Dynamic System class.

    """

    def __init__(self):

        # model name
        self.name = 'pendulum'

        # parameter names and their values
        self.parameters = {'g':9.81,
                           'l':1.,
                           'm':1.}

        # state names and their initial conditions
        self.states = ['theta',
                       'omega']

        self.units = ['radians',
                      'radians per second']

        # sets the initial conditions of the states
        self.x = array([0.,
                        0.])

        # sets the time to the initial time
        self.t = 0.

        # initializes the zees
        self.z = zeros(1)

        # numerical integration parameters
        self.numint = {'ti':0.,
                       'tf':10.,
                       'steps':100}

    def f(self, x, t):
        '''Returns the derivative of the states'''

        # defines the parameters from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

        # sets the current state
        theta = x[0]
        omega = x[1]

        # sets the zees
        self.z[0] = sin(theta)

        # calculates inputs
        torque = self.inputs(t)

        # calculates the derivatives of the states
        thetap = omega
        omegap = -g/l*self.z[0] + torque/(m*l*l)

        # plug in the derivatives for returning
        f = zeros(2)
        f[0] = thetap
        f[1] = omegap

        return f

    def inputs(self, t):
        torque = 0 #10*sin(2*pi*t + pi/6)
        return torque

    def simulate(self):
        self.t = linspace(self.numint['ti'], self.numint['tf'],
                          self.numint['steps'])

        self.y = odeint(self.f, self.x, self.t)

    def plot(self):
        '''Makes a plot of the simulation'''
        plot(self.t, self.y)
        legend(self.states)
        xlabel('Time [sec]')
        show()
