from make_qsub import make_qsub
import subprocess
import os
import sys

def node_qsub(app, nodelist, threads, results, scriptpath):
    buffer = []
    buffer.append('#!/bin/bash\n')
    buffer.append('#\n')
    buffer.append('#$ -cwd\n')
    buffer.append('#$ -j y\n')
    buffer.append('#$ -S /bin/bash\n')
    buffer.append('export OMP_NUM_THREADS=%d\n\n' % threads)

    for nodes in nodelist:
        buffer.append('./%s -n %d > %s/%d.out\n' % (app, nodes, results, nodes))

    file = open(scriptpath, 'w')
    file.writelines(buffer)
    file.close()

def node_experiment(nodes, threads, results):
    qsub = 'run_omp_%d_%d.qsub' % (threads, nodes)

    make_qsub('path-omp.x -n %d' % nodes, threads,
              '%s/%d.out' % (results, nodes), qsub)
    os.system('qsub %s' % qsub)

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'Arguments please: start stop step threads'
    else:
        _, start, stop, step, threads = sys.argv
        results = 'exp_omp_' + threads
        qsub = 'run_omp_%s_nodes.qsub' % threads
        os.system('rm -r -f %s' % results)
        os.system('mkdir %s' % results)
        os.system('make')

        node_qsub('path-omp.x', range(int(start), int(stop), int(step)),
                  int(threads), results, qsub)

        os.system('qsub %s' % qsub)
