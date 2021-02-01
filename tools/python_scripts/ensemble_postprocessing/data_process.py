from ast import literal_eval
import numpy as np
import matplotlib.pyplot as plt
import ast

# Choose your function ---------------------------
function = "RosenbrockNoEnsemble"

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

for variation in variations[function]:

    fig, axs = plt.subplots(2, len(variations[function]))
    counter_num = 1
    counter_num_aux = 0
    for point_num in point_numbers[function]:
        # Settings ===============
        iteration = []
        fval = []
        incunbemt_solution = []
        alternative_fval1 = []
        alternative_fval2 = []
        alternative_fval3 = []
        alternative_fval4 = []
        alternative_fval5 = []
        alternative_fval6 = []
        alternative_fval7 = []
        alternative_fval8 = []
        alternative_fval9 = []
        alternative_fval10 = []
        alternative_fval11 = []


        path1 = function + '/' + str(point_num) + '/data_iterations_' + function + '_' + str(point_num) + str(variation) + '.txt'
        path2 = function + '/'+ str(point_num) + '/data_' + function + '_' + str(point_num) + str(variation) + '.txt'
        with open(path1) as iterations:
            counter = 0
            for line in iterations:
                vec = literal_eval(line)
                counter2 = 0
                with open(path2) as datas:
                    for data in datas:
                        dat = literal_eval(data)
                        if vec[5] == dat[0] and vec[1] == dat[2] and counter2 >= counter:
                            iteration.append(vec[0])
                            fval.append(vec[1])
                            incunbemt_solution.append(vec[5])
                            if function == "Rosenbrock" or function == "Rastrigin":
                                alternative_fval1.append(dat[1][0])
                                alternative_fval2.append(dat[1][1])
                                alternative_fval3.append(dat[1][2])
                                alternative_fval4.append(dat[1][3])
                                alternative_fval5.append(dat[1][4])
                                alternative_fval6.append(dat[1][5])
                                alternative_fval7.append(dat[1][6])
                                alternative_fval8.append(dat[1][7])
                                alternative_fval9.append(dat[1][8])
                                alternative_fval10.append(dat[1][9])
                                alternative_fval11.append(dat[1][10])

                        counter2 = counter2 +1
                counter = counter +1

        if counter_num == len(variations[function])+1:
            counter_num = 1
            counter_num_aux = 1
        if function == "Rosenbrock" or function == "Rastrigin":
            axs[counter_num_aux, counter_num-1].plot(iteration, alternative_fval1, lw=2, color='hotpink')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval2, lw=2, color='violet')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval3, lw=2, color='red')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval4, lw=2, color='coral')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval5, lw=2, color='sienna')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval6, lw=2, color='orange')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval7, lw=2, color='gold')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval8, lw=2, color='green')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval9, lw=2, color='c')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval10, lw=2, color='blue')
            axs[counter_num_aux,counter_num-1].plot(iteration, alternative_fval11, lw=2, color='pink')

        axs[counter_num_aux,counter_num-1].plot(iteration, fval, lw=2, color='black', label=str(fval[-1]))
        fig.suptitle('Variation ' + str(variation))
        axs[counter_num_aux, counter_num-1].set_title('Point ' + str(point_num))
        axs[counter_num_aux, counter_num - 1].legend(framealpha=1, frameon=True)
        counter_num = counter_num +1

plt.show()