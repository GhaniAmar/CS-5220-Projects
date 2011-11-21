import sys
import os
from pylab import *
from matplotlib import rc
from matplotlib.pyplot import figure, axes, semilogx, plot, xlabel, ylabel, title, grid, savefig, show

rc('text', usetex=True)

plotter = plot

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
    ylabel(texify('Time elapsed (s)'), fontsize=14)
    title(texify('Floyd-Warshall MPI Weak Scaling'), fontsize=18)
    legend(tuple(labels), loc=2)
    grid(True)
    savefig(outpath, transparent=True)

def strong_plot(outpath, x, y):
    single = y[0]

    plot(x, [single/i for i in y])

    xlabel(texify('Number of Threads'), fontsize=14)
    ylabel(texify('Time taken (s)'), fontsize=14)
    title(texify('Floyd-Warshall MPI Strong Scaling'), fontsize=18)
    grid(True)
    savefig(outpath, transparent=True)

def timing_plot(folder, data, outpath):
    (n, computation, communication, total) = data

    threads = int(folder[-1])

    plotter(n, computation, label = texify('Computation'))
    plotter(n, communication, label = texify('Communication'))
    plotter(n, total, label = texify('Total time'))

    xlabel(texify('Number of Graph Nodes'), fontsize=14)
    ylabel(texify('Time elapsed (s)'), fontsize=14)
    title(texify('Time Breakdown for MPI Ring Implementation with %d Processors' % threads),
          fontsize=18)
    legend((texify('Computation'), texify('Communication'), texify('Total time')), loc=2)
    grid(True)
    savefig(outpath, transparent=True)

def single_plot(folder, data, outpath):
    (n, computation, communication, total) = data

    threads = int(folder[-1])

    plotter(n, communication, label = texify('Communication'))

    xlabel(texify('Number of Graph Nodes'), fontsize=14)
    ylabel(texify('Time elapsed (s)'), fontsize=14)
    title(texify('Communication Time for MPI Ring Implementation with %d Processors' % threads),
          fontsize=18)
    grid(True)
    savefig(outpath, transparent=True)


def texify(string):
    return r"$\mathrm{%s}$" % string.replace(' ', '\ ')

def read_timing(path):
    file = open(path, 'r')
    lines = [x.strip() for x in file.readlines()]
    file.close()

    computation = float(lines[0].split()[1])
    communication = float(lines[1].split()[1])
    total = float(lines[2].split()[1])

    return (computation, communication, total)

def read_timings(folder):
    files = os.listdir(folder)
    files.sort(key = lambda x : int(x[:-4]))

    n = []
    computation = []
    communication = []
    total = []

    acc = []

    for f in files:
        n.append(int(f[:-4]))
        comp, comm, totes = read_timing(folder + '/' + f)
        computation.append(comp)
        communication.append(comm)
        total.append(totes)

    return (n, computation, communication, total)

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

    elif sys.argv[1] == 'timing':
        folder = sys.argv[2]
        outpath = sys.argv[3]
        data = read_timings(folder)
        timing_plot(folder, data, outpath)

    elif sys.argv[1] == 'single':
        folder = sys.argv[2]
        outpath = sys.argv[3]
        data = read_timings(folder)
        single_plot(folder, data, outpath)
v
