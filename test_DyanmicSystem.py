import DynamicSystem as ds

pen = ds.DynamicSystem()

print "Name of the system: ", pen.name

print "Parameters: ", pen.parameters

print "List of state names: ", pen.states

print "Initial conditions: ", pen.x_init

print "Integration parameters: ", pen.numint

print "Evalute f: ", pen.f([1., 2.], 5.0)

print "Evalute the input: ", pen.inputs(5.)

print "Simulate and plot"
pen.simulate()

pen.plot()

print pen.z

linpen = ds.LinearDynamicSystem([0., 0.])

linpen.simulate()

linpen.plot('loci', param='g',range=(0.,5.))
