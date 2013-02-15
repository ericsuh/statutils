import numpy as np
import scipy.interpolate

def qvalue(pvalues,
        tuning=None,
        pi0method='smoother',
        robust=False,
        smoothdf=3,
        smoothlog=False):
    '''
    Compute q-values given a list of p-values, which can then be used to
    classify p-values based on particular false discovery rate.

    This function is a Python reimplementation of the ``qvalue`` function in
    John Storey's and Alan Dabney's ``qvalue`` R package, which is an
    implementation of his q-value calculation algorithm. [1]_ This Python
    function aims to reproduce the R implementation reasonably closely, save
    some Pythonic/NumPythonic idioms, and renaming/dropping a few parameters.

    Parameters
    ----------

    pvalues : array_like, 1D
        p-values from multiple hypothesis tests for which to calculate
        corresponding q-values. Values must be within [0,1].
    tuning : float or array_like, optional
        Tuning parameters for the ``pi0`` estimation algorithm. Corresponds to
        ``lambda`` in the R package function. Must be either a float, an
        array_like of size 1, or a 1D array_like of size > 4. Values must fall
        within [0,1).
    pi0method : ``'smoother'`` or ``'bootstrap'``, optional
        Method for automatically choosing the tuning parameter in the
        ``pi0`` estimation, the proportion of true null hypotheses.
    robust : bool, optional
        An indicator of whether it is desired to make the estimate more robust
        for small p-values and a direct finite sample estimate of pFDR.
    smoothdf : int, optional
        Degrees of freedom to use in the spline smoother.
    smoothlog : bool, optional
        Whether spline smoothing should be done on a log scale.

    Returns
    -------
    qvalues : ndarray
        Array of q-values corresponding to each p-value.
    pi0 : float
        Estimated proportion of true null hypotheses in total data set.

    References
    ----------
    .. [1] Storey JD. (2002) A direct approach to false discovery rates.
        Journal of the Royal Statistical Society, Series B, 64: 479-498.
    '''
    # Parameter checking
    pvalues = np.asanyarray(pvalues)
    if pvalues.min() < 0 or pvalues.max() > 1 or np.isnan(pvalues).any():
        raise ValueError('p-values must fall in interval [0,1] and not be NaN')
    if tuning is None:
        tuning = np.linspace(0, 0.9, 19) # separated by 0.05
    else:
        tuning = np.asanyarray(tuning)
        tuning = tuning.reshape(tuning.size) # making 1-d array
    if tuning.size < 1 or (tuning.size > 1 and tuning.size < 4):
        raise ValueError('Number of tuning parameters must be either 1 or '
            'more than 4')
    if ((tuning < 0) | (tuning >= 1)).any():
        if tuning.size == 1:
            raise ValueError('Tuning parameter must be within [0,1).')
        else:
            raise ValueError('Tuning parameters must be within [0,1).')
    if pi0method not in ['smoother', 'bootstrap']:
        raise ValueError('pi0method must be one of "smoother" or "bootstrap"')

    # Compute pi0
    pi0 = (pvalues[:, np.newaxis] >= tuning).mean(axis=0) / (1 - tuning)
    if tuning.size > 1:
        if pi0method == 'smoother':
            if smoothlog:
                pi0 = np.log(pi0)
            spi0 = scipy.interpolate.splrep(tuning, pi0, k=smoothdf)
            pi0 = scipy.interpolate.splev(tuning.max(), spi0)
            if smoothlog:
                pi0 = np.exp(pi0)
        elif pi0method == 'bootstrap':
            pboot = np.random.choice(pvalues, (pvalues.size, 100))
            pi0boot = (pboot[:, :, np.newaxis] >= tuning).mean(axis=0) \
                    / (1-tuning)
            mse = ((pi0boot - pi0.min())**2).sum(axis=0)
            pi0 = pi0[mse.argmin()]
    pi0 = min(pi0, 1)
    if pi0 <= 0:
        raise Exception('Estimated pi0 <= 0. Check that you have valid '
                'p-values or use another pi0method')

    # Calculate q-values
    ranks = (pvalues[:, np.newaxis] <= pvalues).sum(axis=0)
    if robust:
        qvalues = pi0 * pvalues.size * pvalues / \
                (ranks * (1-(1-pvalues)**pvalues.size))
    else:
        qvalues = pi0 * pvalues.size * pvalues / ranks
    ordering = pvalues.argsort()
    maxidx = pvalues.size - 1
    qvalues[ordering[maxidx]] = min(
        qvalues[ordering[maxidx]],
        1
    )
    for i in xrange(0, maxidx).__reversed__():
        qvalues[ordering[i]] = np.min([
            qvalues[ordering[i]],
            qvalues[ordering[i+1]],
            1])

    return (qvalues, np.asarray(pi0).item(0))
