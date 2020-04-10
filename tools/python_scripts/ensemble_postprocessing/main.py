import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import



# def f(x, y):
#     z = np.log(x ** 2) + y ** 2
#     z1 = 1.5*np.log(4*x ** 2) + 4*y ** 2
#     z2 = 2*np.log(3*x ** 2) + y ** 2
#     #z3 = 4*np.log(x ** 2 + 5*y ** 2)
#     z4 = np.log(2*x ** 2) + 2*y ** 2+0.01
#     return [z, z1, z2, z4]

# def f(x, y):
#     z = 0.01*(x - 1)**2 + (y - x**2)**2
#     z1 = 0.01*(x+10)**2 + ((y-10) - (x+10)**2)**2
#     z2 = 0.01*(x - 10)**2 + (10*y - (x-10)**2)**2
#     #z3 = 4*np.log(x ** 2 + 5*y ** 2)
#     z4 = 0.01*(x - 1)**2 + (y - (x+4)**2)**2
#     return [z, z1, z2, z4]

def f(x, y):
    z = 100*(y-x**2)**2 + (x-1)**2
    z1 = 95*((y+4)-(x-0.8)**2)**2 + (x-0.8)**2
    z2 = 97 * (y + 0.3 - (x-0.4) ** 2) ** 2 + (x - 1) ** 2
    z3 = 103 * ((y+0.3) - (x+0.4) ** 2) ** 2 + (x - 1.2) ** 2 + 1
    z4 = 94 * ((y-1.8) - (x-0.3) ** 2) ** 2 + (x + 0.8) ** 2
    z5 = 98 * (y - (x+0.7) ** 2) ** 2 + (x + 0.3) ** 2
    z6 = 95 * ((y-1.8) - (x-0.5) ** 2) ** 2 + (x - 1) ** 2
    z7 = 106*(y-(x-0.7)**2)**2 + (x-0.2)**2
    z8 = 96*((y+4)-x**2)**2 + (x-1.3)**2
    z9 = 105*((y-2)-x**2)**2 + (x+0.7)**2
    z10 = 90* ((y + 0.6) - (x-0.2) ** 2) ** 2 + (x - 1) ** 2
    expValue = (1/11)*(z+z1+z2+z3+z4+z5+z6+z7+z8+z9+z10)
    return [z, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10,expValue]

# def f(x, y):
#     z   = 10*2 + (x**2 - 10*np.cos(2* np.pi*x)) + (y**2 - 10* np.cos(2*np.pi*y))
#     z1  = 15*2 + (x**2 - 15*np.cos(2* np.pi*x)) + (y**2 - 15* np.cos(2*np.pi*y))
#     z2  = 10*2 + ((x-0.2)**2 - 10*np.cos(2* np.pi*(x-0.2))) + (y**2 - 10* np.cos(2*np.pi*y))
#     z3  = 10*2 + (x**2 - 10*np.cos(2* np.pi*x+0.4)) + ((y+0.2)**2 - 10* np.cos(2*np.pi*(y+0.2)+0.5))
#     z4  = 11*2 + (x**2 - 10*np.cos(2* np.pi*x+0.6)) + ((y+0.15)**2 - 10* np.cos(2*np.pi*(y+0.15)))
#     z5  = 15*2 + ((x-0.2)**2 - 10*np.cos(2* np.pi*(x-0.2))) + (y**2 - 0.5 - 10* np.cos(2*np.pi*y))
#     z6  = 14*2 + ((x+0.3)**2 - 10*np.cos(2* np.pi*(x+0.3))) + ((y-0.1)**2 - 10* np.cos(2*np.pi*y))
#     z7  = 5*2 + ((x-0.4)**2 + 10*np.cos(2* np.pi*(x-0.4))) + ((y-0.4)**2 - 10* np.cos(2*np.pi*(y-0.4)))
#     z8  = 14*2 + ((x+0.2)**2 - 15*np.cos(2* np.pi*x)) + ((y-0.2)**2 - 0.2 - 15* np.cos(2*np.pi*y))
#     z9  = 13*2 + (x**2 - 10*np.cos(2* np.pi*x+0.6)) + ((y+0.1)**2 - 10* np.cos(2*np.pi*(y-0.2)))
#     z10 = 12*2 + ((x-0.6)**2 - 10*np.cos(2* np.pi*x+1)) + ((y-0.5)**2 - 10* np.cos(2*np.pi*y))
#     expValue = (1 / 11) * (z + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10)
#     return [z, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, expValue]

## Rosenbrock Ensembles Settings --------------------------------------------------------
x = np.outer(np.linspace(-4, 4, 80), np.ones(80))
y = np.outer(np.linspace(-5, 15, 80), np.ones(80)).T # transpose
[z, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, ExpValue] = f(x, y)

## Rastrigin Ensembles Settings --------------------------------------------------------
# x = np.outer(np.linspace(-4, 4, 80), np.ones(80))
# y = np.outer(np.linspace(-4, 4, 80), np.ones(80)).T # transpose
# [z, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, ExpValue] = f(x, y)


# Legend configuration  ----------------------------------------------------------------
fake2Dline = mpl.lines.Line2D([0],[0], linestyle="none", c='purple', marker = 'o')
fake2Dline1 = mpl.lines.Line2D([1],[0], linestyle="none", c='darkblue', marker = 'o')
fake2Dline2 = mpl.lines.Line2D([2],[0], linestyle="none", c='r', marker = 'o')
fake2Dline3 = mpl.lines.Line2D([3],[0], linestyle="none", c='coral', marker = 'o')
fake2Dline4 = mpl.lines.Line2D([4],[0], linestyle="none", c='sienna', marker = 'o')
fake2Dline5 = mpl.lines.Line2D([5],[0], linestyle="none", c='orange', marker = 'o')
fake2Dline6 = mpl.lines.Line2D([6],[0], linestyle="none", c='gold', marker = 'o')
fake2Dline7 = mpl.lines.Line2D([7],[0], linestyle="none", c='green', marker = 'o')
fake2Dline8 = mpl.lines.Line2D([8],[0], linestyle="none", c='c', marker = 'o')
fake2Dline9 = mpl.lines.Line2D([9],[0], linestyle="none", c='blue', marker = 'o')
fake2Dline10 = mpl.lines.Line2D([10],[0], linestyle="none", c='violet', marker = 'o')

line = [fake2Dline, fake2Dline1,fake2Dline2,fake2Dline3,fake2Dline4,fake2Dline5,fake2Dline6,fake2Dline7,fake2Dline8,fake2Dline9,fake2Dline10]
labels = ['z', 'z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'z9', 'z10']

# Surfaces of Ensemble functions  --------------------------------------------------------
fig1 = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(x, y, z, cmap='viridis', edgecolor='purple')
ax.plot_surface(x, y, z1, cmap='viridis', edgecolor='darkblue')
ax.plot_surface(x, y, z2, cmap='viridis', edgecolor='red')
ax.plot_surface(x, y, z3, cmap='viridis', edgecolor='coral')
ax.plot_surface(x, y, z4, cmap='viridis', edgecolor='sienna')
ax.plot_surface(x, y, z5, cmap='viridis', edgecolor='orange')
ax.plot_surface(x, y, z6, cmap='viridis', edgecolor='gold')
ax.plot_surface(x, y, z7, cmap='viridis', edgecolor='green')
ax.plot_surface(x, y, z8, cmap='viridis', edgecolor='c')
ax.plot_surface(x, y, z9, cmap='viridis', edgecolor='blue')
ax.plot_surface(x, y, z10, cmap='viridis', edgecolor='violet')


ax.set_title('Ensemble Functions Surface')
ax.legend(line, labels, numpoints = 1)
fig1.show()

# Individual surface plot   ---------------------------------------------------------------

fig = plt.figure()

# set up the axes for the first plot
ax = fig.add_subplot(2, 6, 1, projection='3d')
surf = ax.plot_surface(x, y, z, cmap='viridis', edgecolor='purple')
ax.set_title('Main Function')

# set up the axes for the second plot
ax = fig.add_subplot(2, 6, 2, projection='3d')
ax.plot_surface(x, y, z1, cmap='viridis', edgecolor='darkblue')
ax.set_title('Variation z1')

# set up the axes for the tirdy plot
ax = fig.add_subplot(2, 6, 3, projection='3d')
ax.plot_surface(x, y, z2, cmap='viridis', edgecolor='red')
ax.set_title('Variation z2')

# set up the axes for the fourth plot
ax = fig.add_subplot(2, 6, 4, projection='3d')
ax.plot_surface(x, y, z3, cmap='viridis', edgecolor='coral')
ax.set_title('Variation z3')

# set up the axes for the fiveth plot
ax = fig.add_subplot(2, 6, 5, projection='3d')
ax.plot_surface(x, y, z4, cmap='viridis', edgecolor='sienna')
ax.set_title('Variation z4')

# set up the axes for the sixth plot
ax = fig.add_subplot(2, 6, 6, projection='3d')
ax.plot_surface(x, y, z5, cmap='viridis', edgecolor='orange')
ax.set_title('Variation z5')

# set up the axes for the seventh plot
ax = fig.add_subplot(2, 6, 7, projection='3d')
ax.plot_surface(x, y, z6, cmap='viridis', edgecolor='gold')
ax.set_title('Variation z6')

# set up the axes for the eighth plot
ax = fig.add_subplot(2, 6, 8, projection='3d')
ax.plot_surface(x, y, z7, cmap='viridis', edgecolor='green')
ax.set_title('Variation z7')

# set up the axes for the nineth plot
ax = fig.add_subplot(2, 6, 9, projection='3d')
ax.plot_surface(x, y, z8, cmap='viridis', edgecolor='c')
ax.set_title('Variation z8')

# set up the axes for the tenth plot
ax = fig.add_subplot(2, 6, 10, projection='3d')
ax.plot_surface(x, y, z9, cmap='viridis', edgecolor='blue')
ax.set_title('Variation z9')

# set up the axes for the eleventh plot
ax = fig.add_subplot(2, 6, 11, projection='3d')
ax.plot_surface(x, y, z10, cmap='viridis', edgecolor='violet')
ax.set_title('Variation z10')

# set up the axes for the twenth plot
ax = fig.add_subplot(2, 6, 12, projection='3d')

ax.plot_surface(x, y, z, cmap='viridis', edgecolor='purple')
ax.plot_surface(x, y, z1, cmap='viridis', edgecolor='darkblue')
ax.plot_surface(x, y, z2, cmap='viridis', edgecolor='red')
ax.plot_surface(x, y, z3, cmap='viridis', edgecolor='coral')
ax.plot_surface(x, y, z4, cmap='viridis', edgecolor='sienna')
ax.plot_surface(x, y, z5, cmap='viridis', edgecolor='orange')
ax.plot_surface(x, y, z6, cmap='viridis', edgecolor='gold')
ax.plot_surface(x, y, z7, cmap='viridis', edgecolor='green')
ax.plot_surface(x, y, z8, cmap='viridis', edgecolor='c')
ax.plot_surface(x, y, z9, cmap='viridis', edgecolor='blue')
ax.plot_surface(x, y, z10, cmap='viridis', edgecolor='violet')


ax.set_title('Ensemble Functions Surface')
ax.legend(line, labels, numpoints = 1)

#
# # Level curves of Ensemble functions  -----------------------------------------------------
fig2, ax = plt.subplots()
ax.contour(x, y, z, 3, colors='purple')
ax.contour(x, y, z1, 3, colors='darkblue')
ax.contour(x, y, z2, 3, colors='red')
ax.contour(x, y, z3, 3, colors='coral')
ax.contour(x, y, z4, 3, colors='sienna')
ax.contour(x, y, z5, 3, colors='orange')
ax.contour(x, y, z6, 3, colors='gold')
ax.contour(x, y, z7, 3, colors='green')
ax.contour(x, y, z8, 3, colors='c')
ax.contour(x, y, z9, 3, colors='blue')
ax.contour(x, y, z10, 3, colors='violet')


plt.title('Ensemble Functions')
plt.legend(line, labels)
fig2.show()

# Expected Value prints ------------------------------------------------------------------
fig4 = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(x, y, ExpValue, cmap='viridis', edgecolor='black')
ax.set_title('Expected Value Surface')
fig4.show()

# Grid Points --------------------------------------------------------------------------
# Rastrigin
# fig6, ax = plt.subplots()
# ax.contour(x, y, z, 3, colors='purple')
# X = np.array([3.8,1.2,-3.2,-2.3])
# Y = np.array([3.5,-0.7,-3,2])
# ax.scatter(X,Y)
# ax.grid()
# ax.plot([-4,4], [0, 0], 'k-', lw=2)
ax.plot([0,0], [-4, 4], 'k-', lw=2)
# Rosenbrock
fig6, ax = plt.subplots()
X = np.array([3.9, 1.2, -4, -2.2, 0])
Y = np.array([14.5, -0.8, -5, 5, 15])
ax.scatter(X,Y)
ax.plot([0, 0], [-5,15], 'k-', lw=2)
ax.plot([-4, 4], [0,0], 'k-', lw=2)
ax.contour(x, y, z, 3, colors='purple')
ax.grid()
# 2D Ensemble functions   ------------------------------------------------------------------
y = 0
[z, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, ExpValue] = f(x, y)

fig3 = plt.figure()

plt.plot(x, z, lw=2, color='purple', label=labels[0])
plt.plot(x, z1, lw=2, color='darkblue', label=labels[1])
plt.plot(x, z2, lw=2, color='red', label=labels[2])
plt.plot(x, z3, lw=2, color='coral', label=labels[3])
plt.plot(x, z4, lw=2, color='sienna', label=labels[4])
plt.plot(x, z5, lw=2, color='orange', label=labels[5])
plt.plot(x, z6, lw=2, color='gold', label=labels[6])
plt.plot(x, z7, lw=2, color='green', label=labels[7])
plt.plot(x, z8, lw=2, color='c', label=labels[8])
plt.plot(x, z9, lw=2, color='blue', label=labels[9])
plt.plot(x, z10, lw=2, color='violet', label=labels[10])



plt.xlabel('x')
plt.ylabel('f(x,y), y=0')

plt.title('Ensemble Functions')
plt.legend(line, labels)
fig3.show()

# Expected Value prints ------------------------------------------------------------------
fig5 = plt.figure()
plt.plot(x, ExpValue, lw=2, color='black')
plt.title('Expected Value Function')
fig5.show()





# Print all  ------------------------------------------------------------------------------
plt.show()