"""
plot.py
=======

Convenience plotting functions for UnitTransformDist objects.
Requires matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_pdf(dist, ax=None, show=True, x_points=1000, **kwargs):
    """
    Plot the probability density function (PDF) of a UnitTransformDist.

    Parameters
    ----------
    dist : UnitTransformDist
        The distribution object to plot.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new figure is created.
    show : bool, default True
        If True, call plt.show() at the end.
    x_points : int, default 1000
        Number of points between 0 and 1 to evaluate the PDF.
    **kwargs : additional keyword arguments passed to plt.plot().

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes used for plotting.
    """
    if ax is None:
        fig, ax = plt.subplots()
    x = np.linspace(0.001, 0.999, x_points)  # avoid singularities at boundaries
    y = dist.pdf(x)
    ax.plot(x, y, **kwargs)
    ax.set_xlabel('y')
    ax.set_ylabel('PDF')
    ax.set_title(f'PDF of {dist.name}')
    ax.grid(True, alpha=0.3)
    if show:
        plt.show()
    return ax


def plot_cdf(dist, ax=None, show=True, x_points=1000, **kwargs):
    """
    Plot the cumulative distribution function (CDF) of a UnitTransformDist.

    Parameters
    ----------
    dist : UnitTransformDist
        The distribution object to plot.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new figure is created.
    show : bool, default True
        If True, call plt.show() at the end.
    x_points : int, default 1000
        Number of points between 0 and 1 to evaluate the CDF.
    **kwargs : additional keyword arguments passed to plt.plot().

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes used for plotting.
    """
    if ax is None:
        fig, ax = plt.subplots()
    x = np.linspace(0.001, 0.999, x_points)
    y = dist.cdf(x)
    ax.plot(x, y, **kwargs)
    ax.set_xlabel('y')
    ax.set_ylabel('CDF')
    ax.set_title(f'CDF of {dist.name}')
    ax.grid(True, alpha=0.3)
    if show:
        plt.show()
    return ax


def plot_sf(dist, ax=None, show=True, x_points=1000, **kwargs):
    """
    Plot the survival function (SF) of a UnitTransformDist.

    Parameters
    ----------
    dist : UnitTransformDist
        The distribution object to plot.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new figure is created.
    show : bool, default True
        If True, call plt.show() at the end.
    x_points : int, default 1000
        Number of points between 0 and 1 to evaluate the SF.
    **kwargs : additional keyword arguments passed to plt.plot().

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes used for plotting.
    """
    if ax is None:
        fig, ax = plt.subplots()
    x = np.linspace(0.001, 0.999, x_points)
    y = dist.sf(x)
    ax.plot(x, y, **kwargs)
    ax.set_xlabel('y')
    ax.set_ylabel('Survival')
    ax.set_title(f'Survival function of {dist.name}')
    ax.grid(True, alpha=0.3)
    if show:
        plt.show()
    return ax


def plot_all(dist, x_points=1000, figsize=(12, 8)):
    """
    Create a 2x2 figure showing PDF, CDF, SF, and a sample histogram (with overlay PDF).

    Parameters
    ----------
    dist : UnitTransformDist
        The distribution object to plot.
    x_points : int, default 1000
        Number of points for the smooth curves.
    figsize : tuple, default (12, 8)
        Figure size (width, height) in inches.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    axes : array of matplotlib.axes.Axes
        The four axes.
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    ax1, ax2, ax3, ax4 = axes.ravel()

    # PDF
    x = np.linspace(0.001, 0.999, x_points)
    y = dist.pdf(x)
    ax1.plot(x, y)
    ax1.set_title('PDF')
    ax1.grid(True, alpha=0.3)

    # CDF
    y = dist.cdf(x)
    ax2.plot(x, y)
    ax2.set_title('CDF')
    ax2.grid(True, alpha=0.3)

    # SF
    y = dist.sf(x)
    ax3.plot(x, y)
    ax3.set_title('Survival')
    ax3.grid(True, alpha=0.3)

    # Histogram with PDF overlay
    sample = dist.rvs(size=1000)
    ax4.hist(sample, bins=50, density=True, alpha=0.6, label='Sample histogram')
    ax4.plot(x, dist.pdf(x), 'r-', label='PDF')
    ax4.set_title('Histogram + PDF')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    fig.suptitle(f'Summary plots for {dist.name}', fontsize=14)
    plt.tight_layout()
    plt.show()
    return fig, axes