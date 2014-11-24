import revreaddy as r2
import numpy as np
import matplotlib.pyplot as plt

numberOfParticles = 100
simLength = 1000

msds = np.zeros((numberOfParticles,simLength))
for i in range(0,numberOfParticles):
	msds[i] = r2.start(simLength)

means = np.zeros(simLength)
for i in range(0,simLength):
	mean = 0.
	for j in range(0,numberOfParticles):
		mean += msds[j,i]
	means[i] = mean / float(numberOfParticles)

plt.plot(means)
plt.show()
