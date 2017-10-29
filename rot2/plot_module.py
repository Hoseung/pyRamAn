import numpy as np


def get_lower_upper(x,y,ax,nbins=10,
                     std=False,
                     percentile=0,
                     xpos_type=None,
                     bintype = "uniform_size",
                     linedata = "mean",
                     alpha=0.5,
                     color='blue'):
    if xpos_type is None:
        if bintype == "uniform_size":
            xpos_type = "median"
        elif bintype == "uniform_width":
            xpos_type = "mean"

    if percentile == 0:
        std = True
    if std:
        n, bins = np.histogram(x, bins=nbins)
        sy, bins = np.histogram(x, bins=nbins, weights=y)
        sy2, bins = np.histogram(x, bins=nbins, weights=y*y)
        mean = sy / n
        std = np.sqrt(sy2/n - mean*mean)
        y_upper = mean + std
        y_lower = mean - std
        xpos = (bins[1:] + bins[:-1])/2
    elif percentile > 0:
        if bintype == "uniform_width":
            x = x[np.argsort(x)]
            y = y[np.argsort(x)]
            n, bins = np.histogram(x, bins=nbins)
            bins[-1] += 0.0001 # last element in the last bin.
            binned_x = x[np.digitize(x, bins)]
            dgt = np.digitize(x, bins)
            binned_x = [x[dgt == i] for i in range(1,nbins+1)]
            binned_y = [y[dgt == i] for i in range(1,nbins+1)]
            xpos = (bins[1:] + bins[:-1])/2
        elif bintype == "uniform_size":
            binned_x = np.array_split(x[np.argsort(x)], nbins)
            binned_y = np.array_split(y[np.argsort(x)], nbins)

        mean, y_upper, y_lower=[], [], []
        if bintype != "uniform_width":
            xpos = []

        if percentile > 1:
            percentile /=100

        #print(binned_y)
        for xb, yb in zip(binned_x, binned_y):
            if len(xb) < 1 :
                print("Too small sized sample. aborting.. ")
                return
            if bintype != "uniform_width":
                xpos.append((min(xb) + max(xb))/2)
            if linedata == "median":
                print("Taking median")
                mean.append(np.median(yb))# it's actually mean.
            elif linedata == "mean":
                print("Taking mean")
                mean.append(np.mean(yb))
            ybsrt = np.argsort(yb)
            #print(ybsrt, np.floor(len(yb)*(1-percentile)))
            y_upper.append(yb[ybsrt[int(np.floor(len(yb)*(0.5-percentile/2)))]])
            y_lower.append(yb[ybsrt[int(np.floor(len(yb)*(0.5+percentile/2)))]])

    return xpos, np.array(y_lower), np.array(y_upper)

def mean_std(x,y,ax, nbins=10, **kwargs):
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)
    ax.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, **kwargs)

def mean_alone(x,y,ax, nbins=10, **kwargs):
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    #sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    #std = np.sqrt(sy2/n - mean*mean)
    ax.plot((_[1:] + _[:-1])/2, mean, **kwargs)


def plot_scatter_mean_binned(ax, xx, yy,
                             nbins = 10,
                             do_fill_between=True,
                             do_mean_std = False,
                             do_mean_alone = False,
                             do_scatter = True,
							 percentile = 0,
                             label=None,
                             scatter_args=dict(color="gray"),
                             std_args=dict(color="orange", linewidth=3),
                             percentile_args=dict(color="orange", linewidth=3)):
    if do_scatter:
        #scatter_args = {"edgecolor":"none", "alpha":1.0}
        ax.scatter(xx, yy, **scatter_args, label="isolated (H)")

    if do_fill_between:
        fill_between_args ={"bintype":"uniform_size",
                 "linedata":"median",
                 "nbins":nbins, "percentile":25}
        xpos, y_lower, y_upper = get_lower_upper(xx, yy, ax,  **fill_between_args)
        ax.plot(xpos, y_lower, lw=2, color="navy")
        ax.plot(xpos, y_upper, lw=2, color="navy")
        #ax.fill_between(xpos, y_lower, y_upper, color = "r", alpha=0.2)
        #ax.plot(xpos, mean, color="blue", linewidth=2)

        fill_between_args ={"bintype":"uniform_size",
                 "linedata":"median",
                 "nbins":nbins, "percentile":75}
        xpos, y_lower, y_upper = get_lower_upper(xx, yy, ax,  **fill_between_args)
        ax.plot(xpos, y_lower, lw=2, color="blue")
        ax.plot(xpos, y_upper, lw=2, color="blue")

    if do_mean_std:
        mean_std(xx, yy, ax, nbins=nbins, **err_bar_args)

    if do_mean_alone:
        mean_alone(xx, yy, ax, nbins=nbins, color='b', linewidth=3)

    if percentile != 0:
        fill_between_args ={"bintype":"uniform_size",
                 "linedata":"median",
                 "nbins":nbins, "percentile":percentile}
        xpos, y_upper, y_lower = get_lower_upper(xx, yy, ax,  **fill_between_args)
        ax.errorbar(xpos, 0.5*(y_lower+y_upper), [y_lower, y_upper], **percentile_args)#, fmt="none")


    #line = ax.scatter(range(1),range(1),edgecolor='none',marker='o', facecolor="b",alpha=1.0, label=label)
    if label is not None: ax.legend()
    #ax.plot()
    #ax.legend(line,"name",numpoints=3,loc=1)
    #ax.set_ylabel(r"$\Delta \lambda_{R_{e}}$ per merger", fontsize=16)
    #ax.set_xlabel("Merger mass ratio", fontsize=16)
