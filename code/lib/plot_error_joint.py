"""This function plots a fitted function along with std of parameters, and a
scatter of the original data
"""
##----------------------------PREAMBLE----------------------------------##
import matplotlib.pyplot as plt
import numpy as np
#external packages

##-------------------------FUNCTION DEFINITION-------------------------------##

def plot_error_joint(x, y, dy,y_0,y_1, y_2, ax = None, color = 'g', label = None, fill = False, xlabel = 'pick xlabel', ylabel = 'pick ylabel', title = "pick title"):
    """
    Parameters
    ----------
    x : array-like
        dependent variable for data

    y : array-like
        original data

    dy : array-like
        error on original data

    y_0 : array_like
        fitted function

    y_1 : array-like
        lower bound on fitted function values

    y_2 : array-like
        upper bound on fitted function values

    fill : bool
        if True, fill between upper and lower bounds
    """


    #prepare figure
    if ax is None:
        plt.rcParams.update({'font.size': 16})
        fig, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))
        plt.grid()

    
    if fill == True:
        #add filling
        ax.fill_between(x, y_1, y_2, color = color, alpha = 0.2)
    else:
        #plot std functions
        ax.plot(x, y_1, ':', c = 'r')
        ax.plot(x, y_2, ':', c='r')


    #scatter original data
    ax.scatter(x, y, c = color)
    ax.errorbar(x, y,dy, ls = 'none', capsize = 3, c = color)



    #plot fitted function
    if label is None:
        ax.plot(x, y_0, c = color, linewidth = 3)
    else:
        ax.plot(x, y_0, c = color, linewidth = 3, label = label)

    #figure properties
    plt.legend()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return
