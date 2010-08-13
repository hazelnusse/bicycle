# States
# The state order and the default initial conditions
omega = 0.0
theta = 0.0
# This gets mapped to an ordered dictionary of the states, description and
# intial conditions
# x = {}
# x['0'] = ['omega', 0.0]
# x['1'] = ['theta', 0.0]

# Parameters (default numerical values)
m = 1.0
g = 9.81
l = 1.0

# Equations of Motion
omegap = -g/l*sin(theta)
thetap = omega
# this gets mapped to the following in the integration function
# f[0] = -g/l*sin(theta)
# f[1] = omega

# Model outputs
ke = 1/2*m*(l*thetap)**2
