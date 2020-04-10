from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from ast import literal_eval
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Choose your function ---------------------------
function = "RastriginNoEnsemble"

# Setiings from each function ---------------------
point_numbers = {
    'Rastrigin' : ['p1', 'p2', 'p3', 'p4'],
    'Rosenbrock' : ['p1', 'p2', 'p3', 'p4', 'p5'],
    'RastriginNoEnsemble' : ['p1', 'p2', 'p3', 'p4'],
    'RosenbrockNoEnsemble' : ['p1', 'p2', 'p3', 'p4', 'p5']
}

variations = {
    'Rastrigin': ['v1', 'v2'],
    'Rosenbrock': ['v1', 'v2', 'v3'],
    'RastriginNoEnsemble': ['v1', 'v2'],
    'RosenbrockNoEnsemble': ['v1', 'v2', 'v3']
}

if function == 'Rastrigin' or function == 'RastriginNoEnsemble':
    x = np.outer(np.linspace(-4, 4, 40), np.ones(40))
    y = np.outer(np.linspace(-4, 4, 40), np.ones(40)).T  # transpose
elif function == 'Rosenbrock' or function == 'RosenbrockNoEnsemble':
    x = np.outer(np.linspace(-4, 4, 160), np.ones(160))
    y = np.outer(np.linspace(-5, 15, 160), np.ones(160)).T  # transpose

functions = {
    'Rastrigin': 10*2 + (x**2 - 10*np.cos(2* np.pi*x)) + (y**2 - 10* np.cos(2*np.pi*y)),
    'Rosenbrock': 100*(y-x**2)**2 + (x-1)**2,
    'RastriginNoEnsemble': 10*2 + (x**2 - 10*np.cos(2* np.pi*x)) + (y**2 - 10* np.cos(2*np.pi*y)),
    'RosenbrockNoEnsemble': 100*(y-x**2)**2 + (x-1)**2
}

for variation in variations[function]:

    counter_num = 1
    counter_num_aux = 0
    for point_num in point_numbers[function]:
        # Settings ===============
        incunbemt_solution = []
        radius = []


        path1 = function + '/' + str(point_num) + '/data_iterations_' + function + '_' + str(point_num) + str(variation) + '.txt'
        with open(path1) as iterations:
            for line in iterations:
                vec = literal_eval(line)
                incunbemt_solution.append(vec[5])
                radius.append(vec[3])

        fig = plt.figure()
        col = int(len(radius)/6)

        for i in range(1,len(radius)):

            # set up the axes for the first plot
            if i < col*6:
                ax = fig.add_subplot(6, col, i )#, projection='3d')

                # Level curves of Ensemble functions  -----------------------------------------------------
                ax.contour(x, y, functions[function], 3, colors='hotpink')
                fig.suptitle('Point ' + str(point_num) + " Variation: " + str(variation))

                plt.scatter(incunbemt_solution[i][0], incunbemt_solution[i][1], s=radius[i]*10000, edgecolors='k', c='k')

plt.show()
#fig.savefig('plotcircles.png')