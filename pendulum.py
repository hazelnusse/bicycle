# State ordering
# x[0] = theta
# x[1] = theatap

# Parameters (default numerical values)
m = 1.0
g = 9.81
l = 1.0

# Equations of Motion
f[0] = thetap
f[1] = -g/l*sin(theta)

# Model outputs
ke = 1/2*m*(l*thetap)**2
