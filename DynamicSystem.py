from numpy import zeros, sin, cos, pi, dot
from numpy import linspace, array, rank
from numpy.linalg import eig
from scipy.integrate import odeint
from matplotlib.pyplot import figure, plot, show, legend, xlabel, title
from matplotlib.pyplot import scatter, colorbar, cm, grid, axis
import pickle
import os
import string

class DynamicSystem:
    """
    Dynamic System class.


    """

    # model name
    name = 'Dynamic System'

    # parameter names and their values
    parameters = {}

    # state names
    state_names = ['x1','x2']

    # initialize state vector
    x = zeros(len(state_names))

    # output names
    output_names = ['y1', 'y2']

    # initialize output vector
    y = zeros(len(output_names))

    # input names
    input_names = ['u1']

    # initialize input vector
    u = zeros(len(input_names))

    # sets the time to the initial time
    t = 0.

    # initializes the zees
    z = zeros(1)

    # sets the initial conditions of the states
    x_init = array([0., 0.])

    # numerical integration parameters
    numint = {'ti':0.,
              'tf':1.,
              'steps':10}

    def __init__(self, model=None, parameters=None):
        '''Will parse the input text file'''
        # set parameter for saving files related to this system
        self.filename = string.join(string.split(self.name), "")
        self.directory = os.path.join('models', self.filename)

    def f(self, x, t):
        '''
        Returns the derivative of the states.

        Parameters:
        -----------
        x : ndarray
            State vector
        t : ndarray
            Time

        Returns:
        --------
        f : ndarray
            dx/dt

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

        '''

        # defines the parameters from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))

        # sets the current state
        x1 = x[0]
        x2 = x[1]

        # calculates inputs
        u = self.inputs(t)

        # sets the zees
        self.z[0] = 0.

        # calculates the derivatives of the states
        x1p = x2
        x2p = 1.

        # plug in the derivatives for returning
        f = zeros(2)
        f[0] = x1p
        f[1] = x2p

        return f

    def inputs(self, t):
        '''
        Returns the inputs to the system.

        Parameters:
        -----------
        t : ndarray
            Time

        Returns:
        --------
        u : ndarray
            u(t)

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

        '''
        u = 1.
        return u

    def outputs(self, x):
        '''
        Returns the outputs of the system.

        Parameters:
        -----------
        x : ndarray
            Current state

        Returns:
        --------
        y : ndarray
            y(t)

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

        '''
        y = zeros(len(self.output_names))
        y[0] = x[0]
        y[1] = x[1]

        return y

    def simulate(self):
        '''
        Simulates the system.

        Parameters:
        -----------

        Returns:
        --------

        Raises:
        -------

        See also:
        ---------

        Examples:
        ---------

        '''
        # time vector
        t = linspace(self.numint['ti'],
                     self.numint['tf'],
                     self.numint['steps'])

        # initialize the state vector
        x = zeros((len(t), len(self.state_names)))
        y = zeros((len(t), len(self.output_names)))
        u = zeros((len(t), len(self.input_names)))

        # set the initial conditions
        x[0] = self.x_init
        y[0] = self.outputs(x[0])
        u[0] = self.inputs(t[0])

        for i in range(len(t)-1):
            # set the interval
            t_int = [t[i], t[i+1]]
            #print "self.t before int = ", self.t
            #print "self.u before int = ", self.u
            #print "self.x before int = ", self.x
            #print "self.z before int = ", self.z
            #print "self.y before int = ", self.y
            # return the next state
            x[i + 1] = odeint(self.f, x[i], t_int)[1, :]
            # calculate the outputs and store them
            y[i + 1] = self.outputs(x[i + 1])
            u[i + 1] = self.inputs(t[i + 1])
            # update all the attributes
            self.t = t[i + 1]
            self.x = x[i + 1]
            self.y = y[i + 1]
            self.u = u[i + 1]
            #print "self.t after int = ", self.t
            #print "self.u after int = ", self.u
            #print "self.x after int = ", self.x
            #print "self.z after int = ", self.z
            #print "self.y after int = ", self.y

        # make a dictionary of the integration and save it to file
        intDict = {'t':t,
                   'x':x,
                   'y':y,
                   'model':self.name,
                   'params':self.parameters}

        # save the simulation to file
        self.save_sim(intDict)

    def save_sim(self, intDict):
        '''
        Save simulation to file

        '''
        if os.path.isdir(self.directory):
            pass
        else:
            os.system('mkdir ' + self.directory)
        pickle.dump(intDict, open(self.directory + self.filename + '.p', 'w'))

    def plot(self):
        '''
        Makes a plot of the simulation

        '''
        intDict = pickle.load(open(self.directory + self.filename + '.p'))
        plot(intDict['t'], intDict['y'])
        legend(self.output_names)
        xlabel('Time [sec]')
        show()

class Pendulum(DynamicSystem):
    """
    A simple one degree of freedom pendulum with length, mass and moment of
    inertia and input torque.

    """

    # model name
    name = 'Pendulum'

    filename = string.join(string.split(name), "")
    directory = 'models/' + filename + '/'

    # parameter names and their values
    parameters = {'g':9.81,
                  'l':1.,
                  'm':1.,
                  'i':1.}

    # state names and their initial conditions
    state_names = ['theta',
                   'omega']

    # sets the initial conditions of the states
    x_init = array([0.,
                    0.])

    # intialize state vector
    x = zeros(len(state_names))

    # output names
    output_names = ['theta',
                    'omega',
                    'kinetic energy',
                    'potential energy',
                    '2*theta']

    # initialize output vector
    y = zeros(len(output_names))

    # input names
    input_names = ['torque']

    # initialize input vector
    u = zeros(len(input_names))

    # sets the time to the initial time
    t = 0.

    # initializes the zees
    z = zeros(17)

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

        # plug in the derivatives for returning
        f = zeros(2)
        f[0] = thetap
        f[1] = omegap

        return f

    def inputs(self, t):
        torque = 10*sin(2*pi*t + pi/6)
        u = torque
        return u

    def outputs(self, x):
        # defines the parameters locally from the attribute
        for parameter, value in self.parameters.items():
            exec(parameter + ' = ' + str(value))
        theta = x[0]
        omega = x[1]
        # calculate the outputs
        k = 0.125*(m*self.z[4]**2+4*i*self.z[3]**2)*omega**2
        p = -0.5*g*l*m*self.z[1]
        th2 = 2*theta
        y = array([theta, omega, k, p, th2])
        return y

class LinearPendulum(Pendulum):

    name = "Linear Pendulum"

    def __init__(self, equi_points):
        '''This function should take the equilibrium points and calculate the
        linear system: A, B, C, D either numerically or with analytic
        expressions'''
        self.equib = equi_points
        self.linear(self.equib)

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

    def linear(self, equi_points):
        '''
        Sets the A, B, C, D matrices based on the equi_points.

        '''
        # sets the zees for the equilbrium points
        Pendulum.f(self, equi_points, 0.)
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

    def eig(self, *args, **kwargs):
        '''
        Calculates the eigenvalues of the system.
        '''

        # is the first arg a parameter?
        if args[0] in self.parameters.keys():
            par_range = linspace(kwargs['range'][0],
                                 kwargs['range'][1],
                                 100)
            w = zeros((len(par_range), rank(self.A)), dtype=complex)
            for i, val in enumerate(par_range):
                # set the parameter
                exec("self.parameters['" + args[0] + "']=" + str(val))
                # calculate the A matrix
                self.linear(self.equib)
                # calculate the eigenvalues
                w[i] = eig(self.A)[0]
            return w, par_range
        else:
            return eig(self.A)

    def plot(self, typ=None, *args, **kwargs):
        '''Makes a plot of the simulation'''
        # plot a graph with all the outputs
        if typ == None:
            intDict = pickle.load(open(self.name + '.p'))
            plot(intDict['t'], intDict['x'])
            legend(self.state_names)
            xlabel('Time [sec]')
        elif typ == 'loci':
            par = kwargs['param']
            par_range = kwargs['range']
            exec("w, p = self.eig(par, range=" + str(par_range) + ')')
            for j in range(w.shape[1]):
                scatter(w[:, j].real, w[:, j].imag, s=2, c=p,
                                    cmap=cm.gist_rainbow,
                                    edgecolors='none')
                colorbar()
                grid()
                axis('equal')
                title('Roci loci wrt to {param}'.format(param=par))
        show()
