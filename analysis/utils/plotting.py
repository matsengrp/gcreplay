r"""Plotting utilities."""

import numpy as np
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess


def regplot_lowess_ci(
    data, x, y, ci_level, n_boot, ax=None, scatter_kwargs={}, line_kwargs={}
):
    """
    Create a LOWESS regression plot with confidence intervals.

    https://github.com/mwaskom/seaborn/issues/552

    Args:
        data (pd.DataFrame): DataFrame containing the data.
        x (str): Column name for the x-axis variable.
        y (str): Column name for the y-axis variable.
        ci_level (float): Confidence interval level (e.g., 95 for 95% CI).
        n_boot (int): Number of bootstrap samples for confidence interval estimation.
        ax (matplotlib.axes._axes.Axes, optional): Axes to plot on. If None, a new figure and axes will be created.
        scatter_kwargs (dict, optional): Additional keyword arguments for the scatter plot.
        line_kwargs (dict, optional): Additional keyword arguments for the line plot.

    Returns:
        matplotlib.axes._axes.Axes: The axes with the plot.
    """
    x_ = data[x].to_numpy()
    y_ = data[y].to_numpy()
    x_grid = np.linspace(start=x_.min(), stop=x_.max(), num=1000)

    def reg_func(_x, _y):
        return lowess(exog=_x, endog=_y, xvals=x_grid)

    beta_boots = sns.algorithms.bootstrap(
        x_,
        y_,
        func=reg_func,
        n_boot=n_boot,
    )
    err_bands = sns.utils.ci(beta_boots, ci_level, axis=0)
    y_plt = reg_func(x_, y_)

    # If no axis is provided, create one
    if ax is None:
        _, ax = plt.subplots()

    sns.scatterplot(x=x_, y=y_, ax=ax, **scatter_kwargs)
    sns.lineplot(x=x_grid, y=y_plt, ax=ax, **line_kwargs)
    ax.fill_between(x_grid, *err_bands, alpha=0.15, **line_kwargs)

    return ax
