import numpy as np
import seaborn as sns
import itertools
import matplotlib.pyplot as plt

# -------------------------------------------------- #
# --- general fig,ax env

def get_axes(L, max_col=3, fig_frame=(5, 4), res=200):
    cols = L if L <= max_col else max_col
    rows = int(L / max_col) + int(L % max_col != 0)
    fig, axes = plt.subplots(rows,
                             cols,
                             figsize=(cols * fig_frame[0], rows * fig_frame[1]),
                             dpi=res)
    if L > 1:
        axes = axes.flatten()
        for s in range(L, max_col*rows):
            for side in ['bottom', 'right', 'top', 'left']:
                axes[s].spines[side].set_visible(False)
            axes[s].set_yticks([])
            axes[s].set_xticks([])
            axes[s].xaxis.set_ticks_position('none')
            axes[s].yaxis.set_ticks_position('none')
    return fig, axes


def remove_frame(axes):
    for side in ['bottom', 'right', 'top', 'left']:
        axes.spines[side].set_visible(False)
    axes.set_yticks([])
    axes.set_xticks([])
    axes.xaxis.set_ticks_position('none')
    axes.yaxis.set_ticks_position('none')
    

# -------------------------------------------------- #
# --- FES

def setFES1D(X, bins, kb='kcal', temp=300,
             interval=None, fill_empty=True):
    kb_mode = {
        'kJ': 0.00831446261815324,
        'kcal': 0.00198720425864083,
        'unit': 1.0 / temp
    }
    KBT = kb_mode[kb] * temp
    hist, edges = np.histogram(X,
                               bins=bins, range=interval,
                               density=True)
    FES = -1 * KBT * np.log(hist)
    FES = FES - np.min(FES)
    if fill_empty:
        max_ = np.max(FES[FES != np.inf])
        FES[FES == np.inf] = max_
    # get min value
    edges_ = edges[:-1]
    min_value = edges_[FES == np.min(FES)]
    return FES, edges, min_value


def convtoFES2D(X, Y, bins, kb='kcal', temp=300,
                fill_empty=True, interval=None):
    kb_mode = {
        'kJ': 0.00831446261815324,
        'kcal': 0.00198720425864083,
        'unit': 1.0 / temp
    }
    KBT = kb_mode[kb] * temp
    H, xedges, yedges = np.histogram2d(X, Y, bins=bins,
                                       range=interval, density=True)
    FES = -1 * KBT * np.log(H)
    FES = FES - np.min(FES)
    if fill_empty:
        max_ = np.max(FES[FES != np.inf])
        FES[FES == np.inf] = max_
    return FES, xedges, yedges


def setFES2D(x, y, bins, kb='kcal', temp=300,
             fill_empty=True, range_fes=None):
    X = np.asarray(x)
    Y = np.asarray(y)

    FES, xed, yed = convtoFES2D(X.flatten(),
                                Y.flatten(),
                                bins,
                                kb=kb, temp=temp,
                                fill_empty=fill_empty,
                                interval=range_fes)

    Xmin = xed.min()
    Xmax = xed.max()
    Ymin = yed.min()
    Ymax = yed.max()
    xx = np.arange(Xmin, Xmax, ((Xmax - Xmin) / bins))
    yy = np.arange(Ymin, Ymax, ((Ymax - Ymin) / bins))
    XX, YY = np.meshgrid(xx, yy)
    return XX, YY, FES.T


def plotFES2D(X, Y, Z, levels, figure, axes,
              colorbar=True, contlabels=True, ghost=False):
    cont = axes.contour(X, Y, Z, levels,
                        colors='k', linewidths=0.5, zorder=2)
    if not ghost:
        surf = axes.contourf(X, Y, Z, levels,
                             cmap='coolwarm_r', zorder=1)
        if colorbar:
            figure.colorbar(surf)
    if contlabels:
        axes.clabel(cont, inline=True, colors='k',
                    fontsize=8, fmt='%1.1f', zorder=3)