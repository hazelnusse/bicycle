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
    name = 'DynamicSystem'

    filename = ''.join(name.split())
    directory = 'models/' + filename + '/'

    # numerical integration parameters
    intOpts = {'ti':0.0,
               'tf':1.0,
               'ts':0.1}

    # parameter names and their values
    parameters = {}

    # state names
    stateNames = ['x1',
                  'x2']

    # sets the initial conditions of the states
    initialConditions = [0.0,
                         0.0]

    # input names
    inputNames = ['u1']

    # output names
    outputNames = ['y1',
                  'y2']

    # initialize state vector
    x = zeros(len(stateNames))

    # initialize output vector
    y = zeros(len(outputNames))

    # initialize input vector
    u = zeros(len(inputNames))

    # initializes the zees
    z = zeros(1)

    # sets the time to the initial time
    t = 0.0


    def __init__(self, model=None, parameters=None):
        '''Nothing'''

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
        for i, name in enumerate(self.stateNames):
            exec(name + ' = ' + 'x[' + str(i) + ']')

        # calculates inputs
        u = self.inputs(t)
        for i, name in enumerate(self.inputNames):
            exec(name + ' = ' + 'self.u[' + str(i) + ']')

        # sets the zees
        self.z[0] = 0.

        # calculates the derivatives of the states
        x1p = x2
        x2p = 1.

        # plug in the derivatives for returning
        f = zeros(len(self.stateNames))
        for i, name in enumerate(self.stateNames):
            exec('f[' + str(i) + '] = ' + name + 'p')

        return f

    def inputs(self, t):
        '''Returns the inputs to the system.

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
        y = zeros(len(self.outputNames))
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
        t = linspace(self.intOpts['ti'],
                     self.intOpts['tf'],
                     (self.intOpts['tf']-self.intOpts['ti'])/self.intOpts['ts'])

        print self.stateNames
        print self.outputNames
        print self.inputNames

        # initialize the vectors
        x = zeros((len(t), len(self.stateNames)))
        y = zeros((len(t), len(self.outputNames)))
        u = zeros((len(t), len(self.inputNames)))

        print u.shape

        # set the initial conditions
        x[0] = self.initialConditions
        y[0] = self.outputs(x[0])
        u[0] = self.inputs(t[0])

        for i in range(len(t)-1):
            print 't =', t[i]
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
        legend(self.outputNames)
        xlabel('Time [sec]')
        show()

class LinearDynamicSystem(DynamicSystem):

    name = "LinearDynamicSystem"

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
