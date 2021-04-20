#!/usr/bin/env python.exe
import argparse
import matplotlib.pyplot as plt
import numpy as np

list_type = lambda s: [int(i) for i in s.split()]

parse = argparse.ArgumentParser(description="plot")
parse.add_argument("fnm", type=str, help="file name")
parse.add_argument("-c", "--column", type=list_type,
                   default=[2], help="column to plot")
parse.add_argument("-x", type=int, default=1, help="x data")
parse.add_argument("-t", "--type", type=int, default=1,
                   help="1, separately; 2, together; 3, sum")
parse.add_argument("-xlabel", type=str, default=None, help="x label")
parse.add_argument("-ylabel", type=str, default=None, help="y label")
parse.add_argument("-xlim", type=list_type, default=None, help="x limit")
parse.add_argument("-ylim", type=list_type, default=None, help="y limit")
args = parse.parse_args()
#print(args)

def xvgPlot(fnm, column=[2], plot_type=1, x=1, xlabel=None, ylabel=None, xlim=None, ylim=None):
    """ 
    fnm:
        file name to draw
    column (default: 2):
        column to draw or columns to draw (Put it in "")
    plot_type (default: 1):
        1) Draw each column separately
        2) Draw each column together
        3) Sum over each column, and then draw it
    x (default: 1):
        x data list, x = 0 means plot without x.
    """
    d = np.loadtxt(fnm, comments=['#', '@'])
    if x <= 0:
        if plot_type == 1:
            for i in column:
                fig, ax = plt.subplots()
                ax.plot(d[:, i-1], label=f"column {i}")
                ax.legend()
                if xlabel:
                    ax.set_xlabel(xlabel)
                if ylabel:
                    ax.set_ylabel(ylabel)
                if xlim:
                    ax.set_xlim(xlim)
                if ylim:
                    ax.set_ylim(ylim)
        elif plot_type == 2:
            fig, ax = plt.subplots()
            for i in column:
                ax.plot(d[:, i-1], label=f"column {i}")
            if xlabel:
                ax.set_xlabel(xlabel)
            if ylabel:
                ax.set_ylabel(ylabel)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)
            ax.legend()
        else:
            fig, ax = plt.subplots()
            y = np.zeros(len(d))
            for i in column:
                y += d[:, i-1]
            ax.plot(x, y, label=f"column {column}")
            if xlabel:
                ax.set_xlabel(xlabel)
            if ylabel:
                ax.set_ylabel(ylabel)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)
            ax.legend()
    else:
        x = d[:, x-1]
        if plot_type == 1:
            for i in column:
                fig, ax = plt.subplots()
                ax.plot(x, d[:, i-1], label=f"column {i}")
                if xlabel:
                    ax.set_xlabel(xlabel)
                if ylabel:
                    ax.set_ylabel(ylabel)
                if xlim:
                    ax.set_xlim(xlim)
                if ylim:
                    ax.set_ylim(ylim)
                ax.legend()
        elif plot_type == 2:
            fig, ax = plt.subplots()
            for i in column:
                ax.plot(x, d[:, i-1], label=f"column {i}")
            if xlabel:
                ax.set_xlabel(xlabel)
            if ylabel:
                ax.set_ylabel(ylabel)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)
            ax.legend()
        else:
            fig, ax = plt.subplots()
            y = np.zeros(len(x))
            for i in column:
                y += d[:, i-1]
            ax.plot(x, y, label=f"column {column}")
            if xlabel:
                ax.set_xlabel(xlabel)
            if ylabel:
                ax.set_ylabel(ylabel)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)
            ax.legend()
    plt.show()


if __name__ == "__main__":
    xvgPlot(fnm=args.fnm, column=args.column, plot_type=args.type, x=args.x,
            xlabel=args.xlabel, ylabel=args.ylabel, xlim=args.xlim, ylim=args.ylim)
