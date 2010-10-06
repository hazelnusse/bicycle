import os
import pickle
from uncertainties import ufloat, unumpy
from numpy import pi, zeros, vstack, dot, mean, array, shape, ones, zeros_like

class Bicycle(object):
    # these are the various parameter sets
    ptypes = ['Benchmark', 'Sharp', 'Moore', 'Peterson']

    def __new__(cls, shortname):
        '''Returns a NoneType object if there is no directory'''
        # is there a data directory for this bicycle? if not, tell to put some
        # fucking data in the folder so we have something to work with!
        try:
            if os.path.isdir('bicycles/' + shortname) == True:
                print "We have foundeth a directory named: bicycles/" + shortname
                return super(Bicycle, cls).__new__(cls)
            else:
                raise ValueError
        except:
            a = "Are you nuts?! Make a directory with basic data for your "
            b = "bicycle in bicycles/shortname, where 'shortname' is the "
            c = "capitalized one word name of your bicycle. Then I can "
            d = "actually created a bicycle object."
            print a+b+c+d
            return None

    def __init__(self, shortname):
        '''
        Sets the parameters if there any that are already saved.

        shortname: string
            shortname of your bicicleta, one word, first letter is capped and
            should match a directory under bicycles/
        '''

        self.shortname = shortname
        self.directory = ('bicycles/' + shortname + '/')
        self.params = {}

        # are there any parameters already listed? grab any parameters and
        # store them
        match = False
        for typ in self.ptypes:
            pfile = False
            for fname in os.listdir(self.directory):
                # name of parameter file
                fnamep = self.shortname + typ + '.p'
                fnametxt = self.shortname + typ + '.txt'
                # check for a pickle file first
                if fname == fnamep:
                    print "Found a pickle file", fname
                    # grab the .p file first if there is one
                    f = open(self.directory + fname, 'r')
                    self.params[typ] = pickle.load(f)
                    # set this flag so that
                    pfile = True
                    match = True
                    f.close()
                # then look for the .txt files, but only if there wasn't a
                # pickled version
                elif fname == fnametxt and pfile == False:
                    print "found a txt file", fname
                    match = True
                    f = open(self.directory + fname, 'r')
                    self.params[typ] = {}
                    # parse the text file
                    for i, line in enumerate(f):
                        list1 = line[:-1].split(',')
                        # if there is an uncertainty value try to make a ufloat
                        try:
                            self.params[typ][list1[0]] = ufloat((eval(list1[1]),
                                                                 eval(list1[2])))
                        # else keep it as a float
                        except:
                            self.params[typ][list1[0]] = eval(list1[1])
                else:
                    pass

        if match == False:
            print "There are no parameters, try calculate_from_measured"

    def save(self, filetype='pickle'):
        '''
        Saves all the parameters to file.

        filetype : string
            'pickle' : python pickled dictionary
            'matlab' : matlab .mat file
            'text' : comma delimited text file

        '''

        if filetype == 'pickle':
            for k, v in self.params.items():
                thefile = self.directory + self.shortname + k + '.p'
                f = open(thefile, 'w')
                pickle.dump(v, f)
                f.close()
        elif filetype == 'matlab':
            # this should handle the uncertainties properly
            print "Doesn't work yet"

        elif filetype == 'text':
            print "Doesn't work yet"

    def calculate_from_measured(self):
        '''
        Calculates the parameters from measured data.

        '''
        # is there a ____Measured.p? if so, load it in
        if os.path.isfile(self.directory + self.shortname + 'Measured.p'):
            # load the measured data file
            f = open(self.directory + self.shortname + 'Measured.p', 'r')
            ddU = pickle.load(f)
            f.close()
        # else get the data from the original text file and save a pickled
        # version
        else:
            f = open(self.directory + self.shortname + 'Measured.txt')
            ddU = {}
            for line in f:
                list1 = line[:-1].split('=')
                list2 = list1[1].split(',')
                # 1. it is a string (name and shortname)
                # 4. it is an array with a sigma array
                # 5. if is an array with no sigma array
                # 6. it is a variable name with no value
                # 2. it is a float and a sigma
                # 3. it is a float with no sigma
                if list1[0] == 'name' or list1[0] == 'shortname':
                    ddU[list1[0]] = list1[1]
                elif list1[1] == '':
                    # then there is no value for this one
                    ddU[list1[0]] = None
                elif list1[1][0] == '[':
                    list2 = list1[1].split('],[')
                    # then this one is an array
                    noms = [float(a) for a in list2[0][1:].split(',')]
                    if list2[1] == '':
                        stds = zeros_like(noms)
                    else:
                        stds = [float(a) for a in list2[1][:-1].split(',')]
                    ddU[list1[0]] = unumpy.uarray((noms, stds))
                else:
                    nom = float(list2[0])
                    if list2[1] == '':
                        std = 0.
                    else:
                        std = float(list2[1])
                    ddU[list1[0]] = ufloat((nom, std))

        # calculate all the benchmark parameters
        par = {}

        # calculate the wheel radii
        par['rR'] = ddU['rearWheelDist']/2./pi/ddU['rearWheelRot']
        par['rF'] = ddU['frontWheelDist']/2./pi/ddU['frontWheelRot']

        # steer axis tilt in radians
        par['lambda'] = pi/180.*(90. - ddU['headTubeAngle'])

        # calculate the front wheel trail
        forkOffset = ddU['forkOffset']
        par['c'] = (par['rF']*unumpy.sin(par['lambda'])
                      - forkOffset)/unumpy.cos(par['lambda'])

        # wheelbase
        par['w'] = ddU['wheelbase']

        # calculate the frame rotation angle
        # alpha is the angle between the negative z pendulum (horizontal) and the
        # positive (up) steer axis, rotation about positive y
        alphaFrame = ddU['frameAngle']
        # beta is the angle between the x bike frame and the x pendulum frame, rotation
        # about positive y
        betaFrame = par['lambda'] - alphaFrame*pi/180

        # calculate the slope of the CoM line
        frameM = -unumpy.tan(betaFrame)

        # calculate the z-intercept of the CoM line
        # frameMassDist is positive according to the pendulum ref frame
        frameMassDist = ddU['frameMassDist']
        cb = unumpy.cos(betaFrame)
        frameB = -frameMassDist/cb - par['rR']

        # calculate the fork rotation angle
        betaFork = par['lambda'] - ddU['forkAngle']*pi/180.

        # calculate the slope of the fork CoM line
        forkM = -unumpy.tan(betaFork)

        # calculate the z-intercept of the CoM line
        forkMassDist = ddU['forkMassDist']
        cb = unumpy.cos(betaFork)
        tb = unumpy.tan(betaFork)
        forkB = - par['rF'] - forkMassDist/cb + par['w']*tb

        # intialize the matrices for the center of mass locations
        frameCoM = zeros((2), dtype='object')
        forkCoM = zeros((2), dtype='object')

        comb = array([[0, 1], [0, 2], [1, 2]])
        # calculate the frame center of mass position
        # initialize the matrix to store the line intersections
        lineX = zeros((3, 2), dtype='object')
        # for each line intersection...
        for j, row in enumerate(comb):
            a = unumpy.matrix(vstack([-frameM[row], ones((2))]).T)
            b = frameB[row]
            lineX[j] = dot(a.I, b)
        frameCoM[:] = mean(lineX, axis=0)
        # calculate the fork center of mass position
        # reinitialize the matrix to store the line intersections
        lineX = zeros((3, 2), dtype='object')
        # for each line intersection...
        for j, row in enumerate(comb):
            a = unumpy.matrix(vstack([-forkM[row], ones((2))]).T)
            b = forkB[row]
            lineX[j] = dot(a.I, b)
        forkCoM[:] = mean(lineX, axis=0)

        par['xB'] = frameCoM[0]
        par['zB'] = frameCoM[1]
        par['xH'] = forkCoM[0]
        par['zH'] = forkCoM[1]

        self.params['Benchmark'] = par
