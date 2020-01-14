"""
Created on January 12, 2020

@author: Yi Wang
"""

def linear_interp_weight(x, x1, x2):
    """ Calculate the weight of linear interpolation.

    """

    if (x1 > x2):
        print(' - linear_interp_weight: x1 > x2.')
        print('x1 = {}'.format(x1))
        print('x2 = {}'.format(x2))
        exit()

    if ((x < x1) or (x > x2)):
        print(' - linear_interp_weight: (x < x1) or (x > x2).')
        print('x = {}'.format(x))
        print('x1 = {}'.format(x1))
        print('x2 = {}'.format(x2))
        exit()

    dist = x2 - x1
    w1 = (x2 - x) / dist
    w2 = (x - x1) / dist

    return (w1, w2)
