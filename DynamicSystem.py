from numpy import zeros, sin, cos, pi, dot
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
                  'm':1.,
                  'i':1.}

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
    z = zeros(17)

    # numerical integration parameters
    numint = {'ti':0.,
              'tf':10.,
              'steps':100}

    def __init__(self):
        '''Does nothing, but will parse the input text file'''

    def f(self, x, t):
        '''Returns the derivative of the states'''

        # defines the parameters from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

        # sets the current state
        theta = x[0]
        omega = x[1]

        # calculates inputs
        torque = self.inputs(t)

        # sets the zees
        self.z[1] = cos(theta)
        self.z[2] = sin(theta)
        self.z[3] = self.z[1]**2 + self.z[2]**2
        self.z[4] = l*self.z[3]
        self.z[8] = g*m
        self.z[9] = torque*self.z[3] - 0.5*self.z[8]*self.z[2]*self.z[4]
        self.z[10] = i*self.z[3]
        self.z[11] = self.z[3]*self.z[10] + 0.25*m*self.z[4]**2
        self.z[12] = self.z[9]/self.z[11]
        self.z[13] = self.z[8]*self.z[1]*self.z[4]/self.z[11]
        self.z[14] = (m*self.z[4]**2+4*i*self.z[3]**2)*omega
        self.z[15] = g*l*m
        self.z[16] = self.z[15]*self.z[2]

        # calculates the derivatives of the states
        thetap = omega
        omegap = self.z[12]

        # calculate the outputs
        k = 0.125*(m*self.z[4]**2+4*i*self.z[3]**2)*omega**2
        p = -0.5*g*l*m*self.z[1]
        th2 = 2*theta

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

        # initialize the state vector
        x = zeros((len(t), len(self.states)))

        # set the initial conditions
        x[0] = self.x_init

        for i in range(len(t)-1):
            # set the interval
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

    name = "Linear Pendulum"
    
    def __init__(self, equi_points):
        '''This function should take the equilibrium points and calculate the
        linear system: A, B, C, D either numerically or with analytic
        expressions'''
        self.equib = equi_points
        # sets the zees for the equilbrium points
        DynamicSystem.f(self,self.equib,0.)
        # defines the A, B, C, D matrices
        self.A = zeros((2,2))
        self.B = zeros(2)
        self.C = zeros((5,2))
        self.D = zeros(5)
        self.A[0,0] = 0
        self.A[0,1] = 1
        self.A[1,0] = -0.5*self.z[13]
        self.A[1,1] = 0
        self.B[0] = 0
        self.B[1] = self.z[3]/self.z[11]
        self.C[0,0] = 1
        self.C[0,1] = 0
        self.C[1,0] = 0
        self.C[1,1] = 1
        self.C[2,0] = 0
        self.C[2,1] = 0.25*self.z[14]
        self.C[3,0] = 0.5*self.z[16]
        self.C[3,1] = 0
        self.C[4,0] = 2
        self.C[4,1] = 0
        self.D[0] = 0
        self.D[1] = 0
        self.D[2] = 0
        self.D[3] = 0
        self.D[4] = 0

    def f(self, x, t):
        '''Returns the derivative of the states'''

        # defines the parameters from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

        # calculates inputs
        u = self.inputs(t)

        xp = dot(self.A,x) + dot(self.B,u)

        thetap = xp[0]
        omegap = xp[1]

        # plug in the derivatives for returning
        f = zeros(2)
        f[0] = thetap
        f[1] = omegap

        return f
