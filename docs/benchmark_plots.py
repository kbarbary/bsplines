import pickle

import matplotlib.pyplot as plt


def benchmark_figure(results, title, unit):
    """
    Parameters
    ----------
    results : dict
        Dictionary where key is spline name and value is a dictionary of
        timing results. The dictionary contains keys ``sizes``, ``times``,
        which are both arrays.
    title : str
        Axes title.
    unit : str
        Unit of ``sizes`` (e.g., knots or points).
    """
    plt.style.use('bmh')
    fig, ax = plt.subplots()

    for key, result in results.items():
        throughputs = result['sizes'] / (1e6 * result['times'])
        ax.semilogx(result['sizes'], throughputs, ls='-', label=key)

    ax.set_xlabel(unit)
    ax.set_ylabel("throughput ({} / us)".format(unit))
    ax.set_title(title)
    ax.legend(loc='upper left')

    return fig


def plot(pikname, title, unit):
    with open(pikname, 'rb') as f:
        results = pickle.load(f)
    return benchmark_figure(results, title, unit)


