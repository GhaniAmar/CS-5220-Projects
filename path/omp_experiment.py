from make_qsub import make_qsub
import subprocess
import sys

def node_experiment(nodes, results):
    threads = 8
    qsub = 'run_omp_%d.qsub' % nodes

    make_qsub('path-omp.x -n %d' % nodes, str(threads),
              '%s/%d.out' % (results, nodes), qsub)
    subprocess.call('qsub %s' % qsub)

if __name__ == '__main__':
    if len(sys.argv < 4):
        print 'Arguments please: start stop step'
    else:
        results = 'omp_exp'
        start, stop, step, _ = sys.argv
        subprocess.call('rm -r -f %s' results)
        subprocess.call('mkdir %s' results)
        subprocess.call('make')

        for nodes in range(int(start), int(stop), int(step)):
            node_experiment(nodes, results)
