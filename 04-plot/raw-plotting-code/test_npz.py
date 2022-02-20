import numpy as np

data = np.load('Test.npz')
dos = data['dos']
print(np.sum(dos))