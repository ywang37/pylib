"""
Created on March 10, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np

#
#-----------------------------------------------------------------------------
#
def read_iter(filename, verbose=True):
    """ Read gctm.iteration file.
    """

    if verbose:
        print(' - read_iter: read ' + filename)

    data = np.loadtxt(filename, skiprows=13, 
            dtype = { 
            'names': ('IT', 'A', 'F', 'FdF0', 'PG', 'NG', 'PS', 'NS'),
            'formats': ('i', 'S3', 'e', 'f', 'f', 'f', 'f', 'f')
            })

    return data
#
#-----------------------------------------------------------------------------
#
def plot_cost(filename, ax=None, n_iter=None, select_iter=None,
        legend=True, verbose=True,
        ylabel='Normalized cost function',
        ):
    """ Plot reduction of cost function.

    Parameters
    ----------
    filename : str
        gctm.iteration file.
    ax : axes
        If ax is none, the function will create a figure and axes
    n_iter : int or None
        Only plot the first *n_iter* iterations. If n_iter is none,
        all iterterations are plotted.
    select_iter : in or None
        If *select_iter* is int, the selected iteration is marked.

    """

    # read gctm.iteration file
    data = read_iter(filename, verbose=verbose)
    if n_iter is not None:
        data = data[0:n_iter]

    if verbose:
        print(' - plot_cost: plot reuction of cost function')

    # plot
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.2, top=0.8, left=0.25, right=0.75)

    # Acccepted and unaccepted iterations
    IT   = data['IT']
    FdF0 = data['FdF0']
    ax.plot(IT, FdF0, marker='o', color='C0', markerfacecolor='w')

    # Only accepted iterations
    A = data['A']
    flag = []
    for i in range(len(A)):
        if (A[i] == b'T'):
            flag.append(True)
        else:
            flag.append(False)
    flag   = np.array(flag)
    A_IT   = IT[flag]
    A_FdF0 = FdF0[flag]
    ax.scatter(A_IT, A_FdF0, marker='o', c='C0', zorder=10)

    # Mark the selected iteration
    if select_iter is not None:
        ax.scatter([select_iter], [FdF0[select_iter-1]], marker='x', 
                c='r', zorder=20)

    # Legend
    if legend:
        handles = []
        handles.append(ax.scatter([], [], marker='o', 
            label='Not Accepted', c='w', edgecolors='C0'))
        handles.append(ax.scatter([], [], marker='o', 
            label='Accepted', c='C0'))
        handles.append(ax.scatter([], [], marker='x',
            label='Selected', c='r'))
    ax.add_artist(plt.legend(handles=handles))

    ax.set_xlabel('Number of iterations')
    ax.set_ylabel(ylabel)
#
#-----------------------------------------------------------------------------
#
def calc_fd_grad(cost1, cost2, perb, log_flag=True,
        method='center', ref=1.0):
    """ Calculate finite difference gradient.
    """

    if method == 'center':
        x2 = ref + perb
        x1 = ref - perb
    elif method == 'forward':
        x2 = ref + perb
        x1 = ref
    elif method == 'backward':
        x2 = ref
        x1 = ref - perb
    else:
        print(' - calc_fd_grad: method error')
        exit()

    if log_flag:
        x_diff =  np.log(x2) - np.log(x1)
    else:
        x_diff = x2 - x1

    grad = (cost2 - cost1) / x_diff

    return grad
#
#-----------------------------------------------------------------------------
#

#
#-----------------------------------------------------------------------------
#

#
#-----------------------------------------------------------------------------
#

