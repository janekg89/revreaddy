import revreaddy as r2
import numpy as np
import matplotlib.pyplot as plt

msds = np.zeros((1000,1000))
for i in range(0,1000):
	msds[i] = r2.start(1000)

means = np.zeros(1000)
for i in range(0,1000):
	mean = 0.
	for j in range(0,1000):
		mean += msds[j,i]
	means[i] = mean / 100.

plt.plot(means)
plt.show()
