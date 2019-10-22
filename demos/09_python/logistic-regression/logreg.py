#!/usr/bin/env python3

# Import our libraries.
# Things would fail here if we hadn't set up docker properly.

import sys

import numpy as np

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa
from mpl_toolkits.mplot3d import Axes3D  # noqa


def np_groupby(arr, groups):
    for i in np.unique(groups):
        yield arr[groups == i]


def iter_groups(names, arr, groups, colours=plt.cm.Set1.colors):
    grouped_arr = np_groupby(arr, groups)
    return zip(names, grouped_arr, colours)


def add_scatter_3d(ax, label, arr, colour):
    ax.scatter(
        arr[:, 0],
        arr[:, 1],
        arr[:, 2],
        c=[colour],
        edgecolor='k',
        s=40,
        depthshade=False,
        label=label,
    )
    return


def plot_scatter_3d(data, labels, actual, predicted):
    # set up a figure twice as wide as it is tall
    fig = plt.figure(figsize=plt.figaspect(0.5))

    # set up the axes
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')

    # Plot the data
    # Left
    ax1.set_title("Actual labels")
    for label, arr, colour in iter_groups(labels, data, actual):
        add_scatter_3d(ax1, label, arr, colour)

    # Right
    ax2.set_title("Predicted labels")
    for label, arr, colour in iter_groups(labels, data, predicted):
        add_scatter_3d(ax2, label, arr, colour)

    # Remove the axis ticks
    for ax in [ax1, ax2]:
        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])

    # Add a legend
    legend = ax2.legend(bbox_to_anchor=[1.1, 1])
    legend.get_frame().set_edgecolor("1.0")

    fig.tight_layout()
    return fig, (ax1, ax2)


def main():

    if len(sys.argv) > 1:
        outfile = sys.argv[1]
    else:
        outfile = "python_docker_test.png"

    # Load an example dataset.
    iris = datasets.load_iris()
    # Plot the darn thing.

    # Perform some simpler machine learning tasks.
    # Predict the species of the iris based on sepal and petal features.
    model = LogisticRegression(
        solver="lbfgs",
        multi_class="auto",
        max_iter=200,
        C=0.1
    )
    model.fit(iris.data, iris.target)

    # Run PCA on the data so we can plot in 3 dimensions instead of 4.
    predicted = model.predict(iris.data)
    X_reduced = PCA(n_components=3).fit_transform(iris.data)
    fig, _ = plot_scatter_3d(
        X_reduced,
        iris.target_names,
        iris.target,
        predicted
    )
    fig.suptitle("First three principal components", fontsize=14)
    fig.savefig(outfile)
    return


if __name__ == "__main__":
    main()
