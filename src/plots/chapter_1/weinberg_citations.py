#!/usr/bin/env python3

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

SMALL_SIZE = 12
MEDIUM_SIZE = 16
sns.set_color_codes("muted")

plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize

plt.tight_layout()
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    }
)

years = [
    1961,
    1962,
    1963,
    1964,
    1965,
    1966,
    1967,
    1968,
    1969,
    1970,
    1971,
    1972,
    1973,
    1974,
    1975,
    1976,
    1977,
    1978,
    1979,
]

author_data = {
    "Weinberg": [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        1,
        1,
        9,
        81,
        150,
        260,
        235,
        268,
        339,
        354,
        406,
    ],
    "Glashow": [1, 2, 0, 1, 1, 1, 2, 0, 1, 2, 3, 10, 8, 18, 13, 14, 29, 27, 69],
    "Salam and Ward": [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 2, 18, 29, 59, 37, 43, 41, 30, 31],
    "Glashow, Iliopoulos, Maiani": [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        2,
        10,
        35,
        49,
        153,
        288,
        337,
        284,
        227,
        212,
    ],
    "'t Hooft": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 19, 35, 63, 34, 46, 32, 33, 49],
    "Gross and Wilczek": [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        34,
        176,
        103,
        124,
        121,
        144,
        141,
    ],
}


def cumsum(lst):
    return [sum(lst[0:i:1]) for i in range(len(lst) + 1)][1:]


author_data = {k: cumsum(v) for k, v in author_data.items()}

data = {"Year": [], "Citations": [], "Author": []}
for i, y in enumerate(years):
    for k in author_data.keys():
        data["Citations"].append(author_data[k][i])
        data["Author"].append(k)
        data["Year"].append(y)


cites = pd.DataFrame(data=data)
fig, ax = plt.subplots(figsize=[8, 5])
g = sns.lineplot(x="Year", y="Citations", hue="Author", data=cites, ax=ax)
ax.legend(ncol=2, loc="upper left", frameon=True)
xlim, ylim = (1965, 1979), (0, 2000)
ax.set_xlim(*xlim)
ax.set_ylim(*ylim)
ax.set_yticks(range(ylim[0], ylim[1] + 1, 200))
ax.set_xticks(range(xlim[0], xlim[1] + 1))

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[1:], labels=labels[1:])

axins = inset_axes(
    ax,
    3.2,
    1.2,
    loc=3,
    bbox_to_anchor=(0.17, 0.23),
    bbox_transform=ax.figure.transFigure,
)
gz = sns.lineplot(x="Year", y="Citations", hue="Author", data=cites, ax=axins)
axins.xaxis.set_ticks_position("top")
xlim_zoom, ylim_zoom = (1969, 1973), (0, 20)
axins.set_xlim(*xlim_zoom)  # apply the x-limits
axins.set_ylim(*ylim_zoom)  # apply the y-limits
axins.set_ylabel("")
axins.set_xlabel("")

mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")
plt.legend([], [], frameon=False)

img_path = "../../img/chapter_1/"
plt.savefig(img_path + "weinberg_citations.pdf", bbox_inches="tight")
