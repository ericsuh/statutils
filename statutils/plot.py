import numpy as np
import statsmodels
import scipy.stats
import matplotlib.pyplot as plt
import itertools

def ecdf(data, **kwargs):
    '''
    Plot empirical cumulative distribution function (eCDF) of data. Additional
    keyword arguments are passed on to matplotlib's "step" command. Returns the
    matplotlib figure object containing the eCDF.
    '''
    data = np.sort(np.asanyarray(data))
    return plt.step(data, statsmodels.distributions.ECDF(data)(data),
            **kwargs)

def density(data, **kwargs):
    '''
    Plot kernel density estimation of data. Additional keyword arguments are
    passed on to matplotlib's "plot" command. Returns matplotlib figure.
    '''
    data = np.sort(np.asanyarray(data))
    xs = np.linspace(data.min(),data.max(),100)
    plt.plot(xs, scipy.stats.gaussian_kde(data)(xs))

def pairs(data, names, **kwargs):
    """
    Plots a scatterplot matrix of subplots.  Each row of "data" is plotted
    against other rows, resulting in a nrows by nrows grid of subplots with the
    diagonal subplots labeled with "names".  Additional keyword arguments are
    passed on to matplotlib's "plot" command. Returns the matplotlib figure
    object containing the subplot grid.

    From http://stackoverflow.com/questions/7941207/is-there-a-function-to-make-scatterplot-matrices-in-matplotlib
    """

    numdata, numvars = data.shape
    fig, axes = plt.subplots(nrows=numvars, ncols=numvars, figsize=(8,8))
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    for ax in axes.flat:
        # Hide all ticks and labels
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        # Set up ticks only on one side for the "edge" subplots...
        if ax.is_first_col():
            ax.yaxis.set_ticks_position('left')
        if ax.is_last_col():
            ax.yaxis.set_ticks_position('right')
        if ax.is_first_row():
            ax.xaxis.set_ticks_position('top')
        if ax.is_last_row():
            ax.xaxis.set_ticks_position('bottom')

    # Plot the data.
    for i, j in zip(*np.triu_indices_from(axes, k=1)):
        for x, y in [(i,j), (j,i)]:
            axes[x,y].plot(data.ix[:,x], data.ix[:,y], **kwargs)

    # Label the diagonal subplots...
    for i, label in enumerate(names):
        axes[i,i].annotate(label, (0.5, 0.5), xycoords='axes fraction',
                ha='center', va='center')

    # Turn on the proper x or y axes ticks.
    for i, j in zip(range(numvars), itertools.cycle((-1, 0))):
        axes[j,i].xaxis.set_visible(True)
        axes[i,j].yaxis.set_visible(True)

    return fig
