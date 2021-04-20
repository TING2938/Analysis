import matplotlib

def configPlot(fontsize=16, lw=2, axeslw=1.2, figsize=(4.7, 3.07)):
    matplotlib.rcParams['font.size'] = fontsize
    matplotlib.rcParams['axes.titlesize'] = fontsize
    matplotlib.rcParams['axes.linewidth'] = axeslw
    matplotlib.rcParams['xtick.major.width'] = axeslw
    matplotlib.rcParams['ytick.major.width'] = axeslw
    matplotlib.rcParams['xtick.minor.width'] = axeslw * 0.85
    matplotlib.rcParams['ytick.minor.width'] = axeslw * 0.85

    matplotlib.rcParams['lines.linewidth'] = lw
    matplotlib.rcParams['legend.fontsize'] = fontsize
    matplotlib.rcParams['figure.figsize'] = figsize

