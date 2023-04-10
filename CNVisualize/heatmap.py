#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Plot reads number in small bins for multiple samples as a heatmap

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def heatmap(df, windows_count, figsize=(25, 4)):
    """
    Plots a heatmap with black vertical lines separating the chromosomes.

    Args:
        df (pd.DataFrame): The DataFrame containing the data to be plotted.
        windows_count (dictionary): Number of windows in each chromosome.
        figsize (tuple of int, optional): The size of the figure. Defaults to (20, 10).

    Returns:
        matplotlib.figure.Figure: The resulting heatmap figure.
    """

    fig = plt.figure(figsize=figsize)
    sns.heatmap(df, cmap='YlOrRd', yticklabels=False)

    region_sizes = list(windows_count.values())
    line_position = 0
    for region_size in region_sizes[:-1]:
        line_position += region_size
        plt.axvline(x=line_position, color='black', lw=1.5)

    xticks = np.cumsum(region_sizes) - np.array(region_sizes) / 2
    xtick_labels = ['chr' + key for key in windows_count.keys()]
    plt.xticks(xticks, xtick_labels)

    return fig
