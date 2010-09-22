from numpy import zeros, sin, pi
from numpy import linspace, array
from scipy.integrate import odeint
from matplotlib.pyplot import figure, plot, show, legend, xlabel
import pickle

class DynamicSystem:
    """Dynamic System class.

    """

    # model name
    name = 'pendulum'

    # parameter names and their values
    parameters = {'g':9.81,
                       'l':1.,
                       'm':1.}

    # state names and their initial conditions
    states = ['theta',
                   'omega']

    units = ['radians',
                  'radians per second']

    # sets the initial conditions of the states
    x_init = array([0.,
                         0.])

    # intialize state vector
    x = zeros(len(states))

    # sets the time to the initial time
    t = 0.

    # initializes the zees
    z = zeros(1)

    # numerical integration parameters
    numint = {'ti':0.,
                   'tf':10.,
                   'steps':100}

    def __init__(self):
        '''Does nothing'''

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
        torque = 10*sin(2*pi*t + pi/6)
        return torque

    def simulate(self):
        # time vector
        t = linspace(self.numint['ti'],
                     self.numint['tf'],
                     self.numint['steps'])

        x = zeros((len(t), len(self.states)))

        x[0] = self.x_init

        for i in range(len(t)-1):
            # set the small interval
            t_int = [t[i], t[i+1]]
            # return the next state
            x[i + 1] = odeint(self.f, x[i], t_int)[1, :]
            self.t = t[i+1]
            self.x = x[i + 1]

        # make a dictionary of the intergration and save it to file
        intDict = {'t':t,'x':x,'model':self.name, 'params':self.parameters}
        pickle.dump(intDict, open(self.name + '.p', 'w'))

    def plot(self):
        '''Makes a plot of the simulation'''
        intDict = pickle.load(open(self.name + '.p'))
        plot(intDict['t'], intDict['x'])
        legend(self.states)
        xlabel('Time [sec]')
        show()

class LinearDynamicSystem(DynamicSystem):
    def __init__(self, equi_points):
        '''This function should take the equilibrium points and calculate the
        linear system: A, B, C, D either numerically or with analytic
        expressions'''
        self.name = "Linear Pendulum"
        self.equib = equi_points

    def f(self, x, t):
        '''Returns the derivative of the states'''

        # defines the parameters from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

        # sets the current state
        theta = x[0]
        omega = x[1]

        # sets the zees
        self.z[0] = theta

        # calculates inputs
        torque = self.inputs(t)

        self.A = array([[0.,   1.],
                        [-g/l, 0.]])

        self.B = array([0., 1./(m*l*l)])

        u = torque

        #xp = self.A*x + self.B*u

        #print xp
        #thetap = xp[0]
        #omegap = xp[1]

        # calculates the derivatives of the states
        thetap = omega
        omegap = -g/l*self.z[0] + torque/(m*l*l)

        # plug in the derivatives for returning
        f = zeros(2)
        f[0] = thetap
        f[1] = omegap

        return f
