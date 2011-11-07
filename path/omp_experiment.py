from make_qsub import make_qsub
import subprocess
import os
import sys

def node_experiment(nodes, threads, results):
    qsub = 'run_omp_%s_%d.qsub' % (threads, nodes)

    make_qsub('path-omp.x -n %d' % nodes, threads,
              '%s/%d.out' % (results, nodes), qsub)
    os.system('qsub %s' % qsub)

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'Arguments please: start stop step threads'
    else:
        _, start, stop, step, threads = sys.argv
        results = 'exp_omp_' + threads
        os.system('rm -r -f %s' % results)
        os.system('mkdir %s' % results)
        os.system('make')

        for nodes in range(int(start), int(stop), int(step)):
            node_experiment(nodes, threads, results)
