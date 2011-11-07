import sys
import os
from pylab import *
from matplotlib import rc
from matplotlib.pyplot import figure, axes, semilogx, plot, xlabel, ylabel, title, grid, savefig, show

rc('text', usetex=True)

plotter = semilogy

def read_file(path):
    handle = open(path, 'r')
    lines = [line.strip() for line in handle.readlines()]
    handle.close()

    n = int(lines[1].split()[1])
    time = float(lines[4].split()[1])

    return (n, time)

def read_folder(path):
    files = os.listdir(path)
    data = [read_file(path + '/' + f) for f in files]
    data.sort(key = lambda x : x[0])
    x, y = zip(*data)

    return (list(x), list(y))

def plot_folder(path, labels):
    x, y = read_folder(path)

    print 'Label for %s?' % path
    lbl = texify(raw_input())
    labels.append(lbl)

    plotter(x, y, label = lbl)

def make_plot(labels, outpath):
    xlabel(texify('Number of Graph Nodes'), fontsize=14)
    ylabel(texify('Time elapsed'), fontsize=14)
    title(texify('Floyd-Warshall OpenMP Speedup Plots'), fontsize=18)
    legend(tuple(labels))
    grid(True)
    savefig(outpath, transparent=True)

def texify(string):
    return r"$\mathrm{%s}$" % string.replace(' ', '\ ')

if __name__ == '__main__':
    labels = []

    if len(sys.argv) < 3:
        print 'Arguments, please: folder1 folder2 ... outpath'
    else:
        for exp in sys.argv[1:-1]:
            plot_folder(exp, labels)

        make_plot(labels, sys.argv[-1])

