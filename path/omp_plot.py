import sys
import os
from pylab import *
from matplotlib import rc
from matplotlib.pyplot import figure, axes, semilogx, plot, xlabel, ylabel, title, grid, savefig, show

rc('text', usetex=True)

plotter = loglog

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
    title(texify('Floyd-Warshall OpenMP Weak Scaling'), fontsize=18)
    legend(tuple(labels), loc=4)
    grid(True)
    savefig(outpath, transparent=True)

def strong_plot(outpath, x, y):
    plot(x, y)

    xlabel(texify('Number of Threads'), fontsize=14)
    ylabel(texify('Time elapsed'), fontsize=14)
    title(texify('Floyd-Warshall OpenMP Strong Scaling'), fontsize=18)
    grid(True)
    savefig(outpath, transparent=True)

def texify(string):
    return r"$\mathrm{%s}$" % string.replace(' ', '\ ')

if __name__ == '__main__':
    labels = []

    if len(sys.argv) < 4:
        print 'Arguments, please: [weak|strong] folder1 folder2 ... outpath'
    elif sys.argv[1] == 'weak':
        for exp in sys.argv[2:-1]:
            plot_folder(exp, labels)

        make_plot(labels, sys.argv[-1])
    elif sys.argv[1] == 'strong':
        data = []
        for exp in sys.argv[2:-1]:
            _, y = read_folder(exp)
            threads = int(exp[-1])
            time = y[-1]

            data.append((threads, time))

        data.sort(key = lambda x : x[0])
        x, y = zip(*data)

        strong_plot(sys.argv[-1], x, y)
