"""
IMPORTS
"""

import datetime
import platform
import os
import matplotlib
from collections import Counter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation
import numpy as np
import tkinter as tk
import ttk as ttkstyle
from matplotlib import style
from mpl_toolkits.axes_grid1 import make_axes_locatable

Configurations_path = os.getcwd() + '\\Configurations\\'
print Configurations_path
plt.switch_backend('TkAgg')
version = '1.0.2 Alpha'
fig1 = plt.figure()
ax1 = fig1.add_subplot(121)
div = make_axes_locatable(ax1)
cax = div.append_axes('right', '5%', '5%')
ax2 = fig1.add_subplot(122, projection='3d')
mng = plt.get_current_fig_manager()
mng.window.state('zoomed')
q = 2


class Charge:
    """
    Defines the data type for a charge
    """

    def __init__(self, x, y, charge):
        self.x = x  # x position of the charge
        self.y = y  # y position of the charge
        self.charge = charge  # This is the charge in coulombs for the charge.


def factor_b_closest_a(a,
                       b):  # Deletes all factors greater than a, even if a factor of b is closer to a than to the next factor i < b.
    """
    Returns the factor of b closest to a
    :param a: int, looking for factor closest too
    :param b: int, factors of which are close to a
    :return: int, factor of b closest to a.
    """
    factors = sorted([i for i in range(1, b) if b % i == 0])  # Defines a list of all the factors of a.
    distance = [a - i for i in
                factors]  # Defines a list of the distance from a to the factor. (a-i) where i is a factors
    q = [i for i in distance if i >= 0]  # removes all factors which are larger than i
    return factors[q.index(min(q))]  # returns the factor of minimum distance.


def create_levels(potentials, n, rng, res):  # Contains unused parameter
    """
    This defines the location of the equipotential lines in the plot.
    :param potentials: This is the data set of the voltage.
    :param n: int, number of levels to return
    :param rng: tuple, unused
    :param res:
    :return: list of levels.
    """
    return sorted([potentials[((res / n)) * i][((res / n)) * i] for i in range(0, n)])  #


def create_meshgrid(resolution, data):
    """
    This creates the meshgrid using the data and resolution.
    :param resolution: Number of data points in each row
    :param data: input data
    :return: the meshgrid
    """
    if len(data) == 1:  # If the data contains only one charge
        """
        If data includes only one charge
        --------------------------------
        The grid, grid_x, grid_y are defined as a meshgrid with domain i.x-10, i.x+10
        --------------------------------
        """
        grid_x, grid_y = np.meshgrid(np.linspace(data[0].x - 10, data[0].x + 10, resolution),
                                     np.linspace(data[0].x - 10, data[0].x + 10, resolution))
    elif abs(max(i.x for i in data) - min(i.x for i in data)) > abs(max(i.y for i in data) - min(i.y for i in data)):
        """
        If the domain of x > domain of y, then
        --------------------------------------
        grid_x,grid_y are meshgrid with domain -1.1 * min(i.x) -> 1.1 max(i.x)
        --------------------------------------
        """
        grid_x, grid_y = np.meshgrid(
            np.linspace(-1.1 * abs(min(i.x for i in data)), 1.1 * abs(max(i.x for i in data)), resolution),
            np.linspace(-1.1 * abs(min(i.x for i in data)), 1.1 * abs(max(i.x for i in data)), resolution))
    else:
        grid_x, grid_y = np.meshgrid(
            np.linspace(-1.1 * abs(min(i.y for i in data)), 1.1 * abs(max(i.y for i in data)), resolution),
            np.linspace(-1.1 * abs(min(i.y for i in data)), 1.1 * abs(max(i.y for i in data)), resolution))
    return grid_x, grid_y


def generate_potential(data, grid_x, grid_y):
    try:
        potential = sum([(9e9 * i.charge) / np.sqrt((i.x - grid_x) ** 2 + (i.y - grid_y) ** 2) for i in data])
    except RuntimeWarning:
        print 'Warning 2001: Encountered a division by 0, passing.'
    return potential


def generate_vectors(grid_x, grid_y, data):
    """
    Returns the vector details of the electric field
    :param grid_x: input x
    :param grid_y: input y
    :param data: List containing Charges
    :return: The output vectors
    """
    try:
        grid_u = sum(
            [9e9 * (i.charge / (np.sqrt((grid_x - i.x) ** 2 + (grid_y - i.y) ** 2) ** 2)) * np.cos(
                np.arctan2(grid_y - i.y, grid_x - i.x)) for
             i in data]) / np.sqrt(grid_x ** 2 + grid_y ** 2)
        grid_v = sum(
            [9e9 * (i.charge / (np.sqrt((grid_x - i.x) ** 2 + (grid_y - i.y) ** 2) ** 2)) * np.sin(
                np.arctan2(grid_y - i.y, grid_x - i.x)) for
             i in data]) / np.sqrt(grid_x ** 2 + grid_y ** 2)
    except RuntimeWarning():
        print 'Warning 2002: Encountered division by 0, passing.'
        pass
    return grid_u, grid_v


def print_charges(data):
    print str(datetime.datetime.now())
    if len(data) == 0:
        print "No Charges"
        return False
    else:
        pass
    row_comp = '+' + (9 * '-') + (len(str(len(data))) * '-') + \
               '+' + ((max(len(str(i.x)) for i in data) + 2) * '-') + \
               '+' + ((max(len(str(i.y)) for i in data) + 2) * '-') + \
               '+' + ((max(len(str(i.charge)) for i in data) + 2) * '-') + '+'
    blank_comp = '|' + (9 * ' ') + (len(str(len(data))) * ' ') + \
                 '|' + ((max(len(str(i.x)) for i in data) + 2) * ' ') + \
                 '|' + ((max(len(str(i.y)) for i in data) + 2) * ' ') + \
                 '|' + ((max(len(str(i.charge)) for i in data) + 2) * ' ') + '|'
    data_comp = []
    for i in range(len(data)):
        data_comp.append('|' + ' ' + 'Charge ' + str(i + 1) + ' ' +
                         '|' + ' ' + str(data[i].x) + (
                             (max(len(str(i.x)) for i in data)) - len(str(data[i].x))) * ' ' + ' ' +
                         '|' + ' ' + str(data[i].y) + (
                             (max(len(str(i.y)) for i in data)) - len(str(data[i].y))) * ' ' + ' ' +
                         '|' + ' ' + str(data[i].charge) + (
                             (max(len(str(i.charge)) for i in data)) - len(str(data[i].charge))) * ' ' + ' |')

    for i in range(len(data)):
        print row_comp
        print blank_comp
        print data_comp[i]
        print blank_comp
    print row_comp

    return True


def get_inputs(base_charges, base_resolution):
    """
    Gets the necessary inputs for the program
    :return: (data),resolution
    """
    global Configurations_path
    global resolution
    input_finished = False
    charges = base_charges
    resolution = base_resolution
    while input_finished is False:
        os.system('cls')
        print_charges(charges)
        correct_input = False
        while correct_input is False:
            correct_input = True
            selection = raw_input("Add, Advanced, Import, exit, Delete, Continue: ")
            if selection == 'add' or selection == 'Add' or selection == 'new':
                try:
                    charges.append(Charge(float(raw_input('What is the x position of the charge? ')),
                                          float(raw_input("What is the y position of the charge? ")),
                                          float(raw_input('What is the particle\'s Charge? '))))
                except ValueError:
                    print "Error 0001: Input must be a floating point number."
                    correct_input = False
            elif selection == "Advanced" or selection == 'advanced' or selection == 'options':
                advanced_loop = False
                while advanced_loop == False:
                    os.system('cls')
                    print '[0] Change Configuration Path\n[1] System Info\n[2] Change Resolution\n[Enter] Back...'
                    adv_selection = raw_input('Which option to you wish to configure? ')
                    if adv_selection == '0' or adv_selection == 'Change Configuration Path':
                        os.system('cls')
                        config_path_loop = False
                        while config_path_loop == False:
                            print "Current Configuration Path: (%s) Current Working Directory: (%s) " % (
                            Configurations_path, os.getcwd())
                            new_path = raw_input("New Configuration Path (Must be from CWD): ")
                            if os.path.isdir(os.getcwd() + new_path) == True:
                                Configurations_path = new_path
                                config_path_loop = True
                            elif new_path == 'back':
                                config_path_loop = True
                            else:
                                print "Input failed, the path was not existant."
                                config_path_loop = False
                    elif adv_selection == '1' or adv_selection == "System Info":
                        os.system('cls')
                        print "Publisher: Boolean Designs"
                        print "Build Lead: Nathan Diggins"
                        print "EFAPI Version: %s" % version
                        print "Computer Name: %s" % platform.node()
                        print "OS type: %s" % platform.system()
                        waiter = raw_input('Press any key to continue: ')
                    elif adv_selection == '2' or adv_selection == "Change Resolution":
                        try:
                            print "Current Resolution: %s" % (str(resolution))
                            resolution = int(raw_input('New Resolution: '))
                        except ValueError:
                            print "Must be a number."
                    elif adv_selection == '' or adv_selection == 'back':
                        advanced_loop = True
                    else:
                        print "Invalid Input"
                        advanced_loop = False
            elif selection == 'Import' or selection == 'import' or selection == 'i':
                """

                This allows users to upload large data sets to the program without manual input

                """
                data_sets = os.listdir(Configurations_path)  # Lists all of the potential data files
                print "Found %s data sets in current configurations folder (%s). Configuration file can be changed in " \
                      "advanced settings. " % (str(len(data_sets)), Configurations_path)
                for file in data_sets:
                    print '[%s] Name: %s | Ext: %s' % (str(data_sets.index(file)),
                                                       file.split('.')[0],
                                                       file.split('.')[1])
                """

                Starts the file selection while loop

                """
                file_selection_loop = False
                while file_selection_loop == False:
                    selected_file = raw_input("Which file number do you wish to import? Enter 'exit' to go back: ")
                    try:
                        if selected_file == 'exit':
                            file_selection_loop = True
                        elif int(selected_file) <= len(data_sets) - 1 and int(selected_file) > -1:
                            """
                            File number exists and can be accessed

                            """
                            print Configurations_path + data_sets[int(selected_file)]
                            data_set = open(Configurations_path + data_sets[int(selected_file)], 'r+')
                            """
                            Has found and opened file

                            """
                            print 'Selected file %s.' % (selected_file)
                            file_data = data_set.read().split('\n')
                            print "Found %s charges in data set." % (str(len(file_data)))
                            for charge in file_data[:-1]:
                                temp_charge = charge.split(',')
                                if len(temp_charge) != 3:
                                    raise SyntaxError("Dataset improperly formatted")
                                else:
                                    pass
                                charges.append(Charge(int(temp_charge[0]), int(temp_charge[1]), int(temp_charge[2])))
                                data_set.close()
                            print "File successfully imported"
                            file_selection_loop = True
                        else:
                            print "Selection was invalid, must correspond to a listed number."
                            file_selection_loop = False
                    except ValueError:
                        print "Input must be a number."
                        file_selection_loop = False
                    except IOError:
                        print "File has been moved or deleted."
                        file_selection_loop = False
                    except SyntaxError:
                        print "Please review the data structure in the readme.txt"
                        file_selection_loop = False
                    except WindowsError:
                        print "Configuration Path is incorrect, please reconfigure."
                        file_selection_loop = False
            elif selection == 'Delete' or selection == 'delete' or selection == 'del':
                try:
                    del (charges[int(raw_input("Which charge do you wish to delete? ")) - 1])
                except ValueError:
                    print 'Error 0002: Input must be an integer.'
                    correct_input = False
            elif selection == 'Exit' or selection == 'Done' or selection == 'finish' or selection == 'exit':
                print "Exiting the Electrical Field API"
                exit()
            elif selection == 'Continue' or selection == '' or selection == 'go':
                return charges, resolution
            elif selection == 'res' or selection == 'resolution':
                try:
                    resolution = int(raw_input('What is the new resolution? '))
                except ValueError:
                    print 'Error 0003: Input must be an integer.'
                    correct_input = False
            else:
                print 'Error 0004: Input was improper. '
                correct_input = False


def create_graph(data, detail):
    ax1.clear()
    m_grid = create_meshgrid(detail, data)
    e_vectors = generate_vectors(m_grid[0], m_grid[1], data)
    # Plotting the charges
    color = np.log(np.sqrt(e_vectors[0] ** 2 + e_vectors[1] ** 2))
    stream_plot = ax1.streamplot(m_grid[0], m_grid[1], e_vectors[0], e_vectors[1], density=q, color=color,
                                 cmap=plt.cm.inferno)
    cont_plot_data = generate_potential(data, m_grid[0], m_grid[1])
    levels = create_levels(cont_plot_data, factor_b_closest_a(10, detail), [np.amin(m_grid[0]), np.amax(m_grid[0])],
                           detail)
    p = Counter(levels)  # This allows us to make sure that our levels list doesn't fall on an equipotential line
    if p.values().count(1) < len(
            levels) / 2:  # This is relatively arbitrary, we want to make sure its not occuring more than twice, but if it is it might be a point field
        """
        This means that there is an equipotential on y=x
        """
        contour_plot = ax1.contour(m_grid[0], m_grid[1], generate_potential(data, m_grid[0], m_grid[1]), 10, colors='r')
    elif p.values().count(1) >= len(levels) / 2:
        contour_plot = ax1.contour(m_grid[0], m_grid[1], generate_potential(data, m_grid[0], m_grid[1]), levels=levels,
                                   colors='r')
    else:
        raise ArithmeticError("Error 3001: Failed to determine equipotentiality")
    cax.cla()
    plt.colorbar(stream_plot.lines, cax=cax)
    # Plotting the Particles
    for i in data:
        if i.charge > 0:
            ax1.plot(i.x, i.y, 'ro')
        else:
            ax1.plot(i.x, i.y, 'bo')

    # Customizing the graph

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$y$')
    ax1.axis([np.amin(m_grid[0]), np.amax(m_grid[0]), np.amin(m_grid[0]), np.amax(m_grid[0])])
    ax1.set_title(r'$\vec E = \sum_{i=0}^{i=n} k \frac{Q_i}{R^2}$')
    ax1.set_aspect('equal', adjustable='box')
    ax1.spines['bottom'].set_position('zero')
    ax1.spines['left'].set_position('zero')
    ax1.spines['top'].set_color('none')
    ax1.spines['right'].set_color('none')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    ax1.grid()

    cont_plot_data = generate_potential(data, m_grid[0], m_grid[1])
    surface = ax2.contour3D(m_grid[0], m_grid[1], cont_plot_data, 150, cmap=plt.cm.jet)
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$y$')
    ax2.set_zlabel(r'$ln(z)$')
    ax2.set_xlim(np.amin(m_grid[0]), np.amax(m_grid[1]))
    ax2.set_ylim(np.amin(m_grid[0]), np.amax(m_grid[1]))
    ax2.set_zscale('log', basez=10)
    ax2.set_title(r'$\vec E = \sum_{i=0}^{i=n} k \frac{Q_i}{R^2}$')
    ax2.set_aspect('equal', adjustable='box')
    ax2.mouse_init()
    ax2.grid()


def animate(i):
    """
    This function animates the graph
    :param i:
    :return:
    """
    inputs = get_inputs(charges, resolution)
    ax1.cla()
    ax2.cla()
    create_graph(inputs[0], inputs[1])


animations = []
charges = [Charge(0.0, 0.0, 1.0)]
resolution = 50
ani = animation.FuncAnimation(fig1, animate, interval=1000)
plt.show()
