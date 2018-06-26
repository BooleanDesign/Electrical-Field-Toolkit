# Imports
import os
import matplotlib
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy as np
import tkinter as tk
import ttk as style

matplotlib.use("TkAgg")

# Classes
class App(tk.Tk):
    """
    This is the TK base code
    """

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        container = tk.Frame(self)  # This is the application window. Frame is the window
        """
        Styling
        """
        # tk.Tk.iconbitmap(self,default='file path') Changes the image, must be logo file .ico
        tk.Tk.wm_title(self, 'Electric Field Tool Kit')
        """
        Expand tells things to expand as large as possible in the window
        Fill is how we tell it to fill up the window.
        """
        container.pack(side="top", fill="both", expand=True)
        """
        Creating the grid
        -----------------
        * 0 is the minimum possible size of the grid_row or grid_column.
        * weight is the importance of one over the other.
        """
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        frame = Page(container, self)
        self.frames[Page] = frame

        frame.grid(row=0, column=0, sticky='nsew')  # This creates a grid structure for TKinter.
        self.show_frame(Page)

    def show_frame(self, cont):
        frame = self.frames[cont]  # Cont is all of the frames in the container

class Page(tk.Frame):
    def __init__(self, parent, controller, data, detail):
        tk.Frame.__init__(self, parent)
        title = tk.Label(self, text="Electrical Field Toolkit")
        title.pack()
class Charge:
    """
    Defines the data type for a charge
    """

    def __init__(self, x, y, charge):
        self.x = x  # x position of the charge
        self.y = y  # y position of the charge
        self.charge = charge  # This is the charge in coulombs for the charge.


# Low Level Definitions
def factor_b_closest_a(a,b):  # Deletes all factors greater than a, even if a factor of b is closer to a than to the next factor i < b.
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


def create_graph(data, detail):
    m_grid = create_meshgrid(detail, data)
    e_vectors = generate_vectors(m_grid[0], m_grid[1], data)
    ax = plt.subplot()
    # Plotting the charges
    color = np.log(np.sqrt(e_vectors[0] ** 2 + e_vectors[1] ** 2))
    stream_plot = ax.streamplot(m_grid[0], m_grid[1], e_vectors[0], e_vectors[1], density=q, color=color,
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
        contour_plot = ax.contour(m_grid[0], m_grid[1], generate_potential(data, m_grid[0], m_grid[1]), 10, colors='r')
    elif p.values().count(1) >= len(levels) / 2:
        contour_plot = ax.contour(m_grid[0], m_grid[1], generate_potential(data, m_grid[0], m_grid[1]), levels=levels,
                                  colors='r')
    else:
        raise ArithmeticError("Error 3001: Failed to determine equipotentiality")
    plt.colorbar(stream_plot.lines)
    # Plotting the Particles
    for i in data:
        if i.charge > 0:
            plt.plot(i.x, i.y, 'ro')
        else:
            plt.plot(i.x, i.y, 'bo')

    # Customizing the graph

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.axis([np.amin(m_grid[0]), np.amax(m_grid[0]), np.amin(m_grid[0]), np.amax(m_grid[0])])
    ax.set_title(r'$\vec E = \sum_{i=0}^{i=n} k \frac{Q_i}{R^2}$')
    ax.set_aspect('equal', adjustable='box')
    ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_position('zero')
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.grid()

    """
    Creating potential grid
    """
    plt.figure()
    ax2 = plt.subplot(projection='3d')
    surface = ax2.contour3D(m_grid[0], m_grid[1], cont_plot_data, 150, cmap=plt.cm.jet)
    plt.colorbar(surface)
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$y$')
    ax2.set_zlabel(r'$ln(z)$')
    ax2.set_xlim(np.amin(m_grid[0]), np.amax(m_grid[1]))
    ax2.set_ylim(np.amin(m_grid[0]), np.amax(m_grid[1]))
    ax2.set_zscale('log', basez=10)
    ax2.set_title(r'$\vec E = \sum_{i=0}^{i=n} k \frac{Q_i}{R^2}$')
    ax2.set_aspect('equal', adjustable='box')
    ax2.grid()

    plt.show()


def print_charges(data):
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


# Higher Level Definitions
def get_inputs():
    """
    Gets the necessary inputs for the program
    :return: (data),resolution
    """
    input_finished = False
    charges = []
    resolution = 50
    while input_finished is False:
        os.system('cls')
        print_charges(charges)
        correct_input = False
        while correct_input is False:
            correct_input = True
            selection = raw_input("Add,Delete,Continue,Resolution: ")
            if selection == 'add' or selection == 'Add':
                try:
                    charges.append(Charge(float(raw_input('What is the x position of the charge? ')),
                                          float(raw_input("What is the y position of the charge? ")),
                                          float(raw_input('What is the particle\'s Charge? '))))
                except ValueError:
                    print "Error 0001: Input must be a floating point number."
                    correct_input = False
            elif selection == 'Delete' or selection == 'delete':
                try:
                    del (charges[int(raw_input("Which charge do you wish to delete? ")) - 1])
                except ValueError:
                    print 'Error 0002: Input must be an integer.'
                    correct_input = False
            elif selection == 'Continue' or selection == 'Cont' or selection == 'go':
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


# Main Function
def main():
    main_loop = False
    while main_loop is False:
        inputs = get_inputs()
        create_graph(inputs[0], inputs[1])
        status = raw_input('Would you like to quit or continue? ')
        if status == 'quit' or status == 'yes':
            exit()
        elif status == 'no' or status == 'continue':
            main_loop = False
        else:
            print "Input failed"
            exit()



app = App()
app.mainloop()
