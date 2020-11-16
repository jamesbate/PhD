import matplotlib.pyplot as plt  
import numpy as np 

def plot_error_joint(x, y, dy,y_0,y_1, y_2, ax = None, color = 'g', label = None):
    if ax is None:
        plt.rcParams.update({'font.size': 16})
        fig, ax = plt.subplots(1	,1, constrained_layout=True, figsize=(18, 9))
        plt.grid()
    ax.scatter(x, y, c = color)
    ax.errorbar(x, y,dy, ls = 'none', capsize = 3, c = color)
    ax.plot(x, y_1, ':', c = 'r')
    ax.plot(x, y_2, ':', c='r')
    if label is None:
        ax.plot(x, y_0, c = color, linewidth = 3)
    else:
        ax.plot(x, y_0, c = color, linewidth = 3, label = label)
    plt.legend()
    ax.set_title('Rabi Flops')
    ax.set_xlabel('Pulse Length (us)')
    ax.set_ylabel('Probability')

    return 