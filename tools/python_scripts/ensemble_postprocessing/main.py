from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def f(x, y):
    z = np.log(x ** 2) + y ** 2
    z1 = 1.5*np.log(4*x ** 2) + 4*y ** 2
    z2 = 2*np.log(3*x ** 2) + y ** 2
    #z3 = 4*np.log(x ** 2 + 5*y ** 2)
    z4 = np.log(2*x ** 2) + 2*y ** 2+0.01
    return [z, z1, z2, z4]


x = np.outer(np.linspace(-0.4, 0.4, 120), np.ones(120))
y = x.copy().T # transpose
[z, z1, z2, z4] = f(x, y)


# Legend configuration  ----------------------------------------------------------------
fake2Dline = mpl.lines.Line2D([0],[0], linestyle="none", c='k', marker = 'o')
fake2Dline2 = mpl.lines.Line2D([1],[0], linestyle="none", c='r', marker = 'o')
fake2Dline3 = mpl.lines.Line2D([2],[0], linestyle="none", c='g', marker = 'o')
fake2Dline4 = mpl.lines.Line2D([3],[0], linestyle="none", c='c', marker = 'o')
line = [fake2Dline,fake2Dline2,fake2Dline3,fake2Dline4]
labels = ['log(x^2) + y^2','1.5log(4x^2) + 4y^2','2log(3x^2) + y^2','log(2x^2) + 2y^2 + 0.01']

# Surfaces of Ensemble functions  --------------------------------------------------------
fig1 = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(x, y, z, cmap='viridis', edgecolor='black')
ax.plot_surface(x, y, z1, cmap='viridis', edgecolor='red')
ax.plot_surface(x, y, z2, cmap='viridis', edgecolor='green')
ax.plot_surface(x, y, z4, cmap='viridis', edgecolor='cyan')

ax.set_title('Ensemble Functions Surface')
ax.legend(line, labels, numpoints = 1)
fig1.show()

# Level curves of Ensemble functions  -----------------------------------------------------
fig2, ax = plt.subplots()
ax.contour(x, y, z, 3, colors='black')
ax.contour(x, y, z1, 3, colors='red')
ax.contour(x, y, z2, 3, colors='green')
ax.contour(x, y, z4, 3, colors='c')

plt.title('Ensemble Functions')
plt.legend(line, labels)
fig2.show()

# 2D Ensemble functions   ------------------------------------------------------------------
y = 0
[z, z1, z2, z4] = f(x, y)

fig3 = plt.figure()

plt.plot(x, z2, lw=2, color='green', label=labels[0])
plt.plot(x, z, lw=2, color='black', label=labels[1])
plt.plot(x, z1, lw=2, color='red', label=labels[2])
plt.plot(x, z4, lw=2, color='c', label=labels[3])


plt.xlabel('x')
plt.ylabel('f(x,y), y=0')

plt.title('Ensemble Functions')
plt.legend(line, labels)
fig3.show()

# Print all  ------------------------------------------------------------------------------
plt.show()