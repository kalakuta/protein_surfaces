import numpy as np
from sklearn import manifold

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter

test_data = np.empty((0,3))
colours = []


for i in range(500):
    theta = np.random.uniform(high = 2 * np.pi)
    if theta < np.pi:
    	colours.append([(np.cos(theta / 2)) ** 2, (np.sin(theta / 2)) ** 2, 0])
    else:
        colours.append([(np.cos(theta / 2)) ** 2, 0,  (np.sin(theta / 2)) ** 2])
    phi = np.random.uniform(high = np.pi / 2)
    x, y, z = np.cos(theta) * np.cos(phi) + 0 * np.random.normal(scale = 0.1), np.sin(theta) * np.cos(phi) + 0 * np.random.normal(scale = 0.1), np.sin(phi) + 0 * np.random.normal(scale = 0.1)
    test_data = np.vstack([test_data, [x, y, z]])

    
fig = plt.figure(figsize=(20, 20))
#plt.suptitle("Manifold Learning with %i points, %i neighbors"
#             % (1000, n_neighbors), fontsize=14)





ax = fig.add_subplot(231, projection='3d')
ax.scatter(test_data[:, 0], test_data[:, 1], test_data[:, 2], c=colours, cmap=plt.cm.Spectral)
ax.view_init(4, -72)

n_neighbors = 10
n_components = 2


mds = manifold.MDS(n_components, max_iter=200, n_init=1, n_jobs=-1)
flattened_data = mds.fit_transform(test_data)

ax = fig.add_subplot(233)
plt.scatter(flattened_data[:, 0], flattened_data[:, 1], c=colours, cmap=plt.cm.Spectral)
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

flattened_data = manifold.LocallyLinearEmbedding(n_neighbors, n_components, eigen_solver='auto', method='standard').fit_transform(test_data)

ax = fig.add_subplot(234)
plt.scatter(flattened_data[:, 0], flattened_data[:, 1], c=colours, cmap=plt.cm.Spectral)
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

flattened_data = manifold.LocallyLinearEmbedding(n_neighbors, n_components, eigen_solver='auto', method='modified').fit_transform(test_data)

ax = fig.add_subplot(235)
plt.scatter(flattened_data[:, 0], flattened_data[:, 1], c=colours, cmap=plt.cm.Spectral)
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

flattened_data = manifold.Isomap(n_neighbors, n_components).fit_transform(test_data)

ax = fig.add_subplot(236)
plt.scatter(flattened_data[:, 0], flattened_data[:, 1], c=colours, cmap=plt.cm.Spectral)
ax.xaxis.set_major_formatter(NullFormatter())
ax.yaxis.set_major_formatter(NullFormatter())
plt.axis('tight')

plt.show()
