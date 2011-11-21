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
    buffer.append('#$ -pe orte %d\n\n' % threads)

    for nodes in nodelist:
        buffer.append('/opt/openmpi/bin/mpirun -np $NSLOTS ./%s -n %d > %s/%d.out\n' %
                      (app, nodes, results, nodes))

    file = open(scriptpath, 'w')
    file.writelines(buffer)
    file.close()

def thread_qsub(app, nodes, threads, results, scriptpath):
    buffer = []
    buffer.append('#!/bin/bash\n')
    buffer.append('#\n')
    buffer.append('#$ -cwd\n')
    buffer.append('#$ -j y\n')
    buffer.append('#$ -S /bin/bash\n')
    buffer.append('#$ -pe orte %d\n\n' % threads)

    buffer.append('/opt/openmpi/bin/mpirun -np $NSLOTS ./%s -n %d > %s/%d.out\n' %
                  (app, nodes, results, threads))

    file = open(scriptpath, 'w')
    file.writelines(buffer)
    file.close()

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'Arguments please: [weak|strong] start stop step threads'
    else:
        _, experiment, start, stop, step, threads = sys.argv

        if experiment == 'weak':
            results = 'exp_mpi_' + threads
            qsub = 'run_mpi_%s_nodes.qsub' % threads
            os.system('rm -r -f %s' % results)
            os.system('mkdir %s' % results)
            os.system('make')

            node_qsub('path-mpi.x', range(int(start), int(stop), int(step)),
                      int(threads), results, qsub)

            os.system('qsub %s' % qsub)

        elif experiment == 'strong':
            results = 'exp_mpi_strong'
            qsub = 'run_mpi_strong.qsub'
            os.system('rm -r -f %s' % results)
            os.system('mkdir %s' % results)
            os.system('make')

            for i in range(1, int(threads) + 1):
                thread_qsub('path-mpi.x', int(start), i, results, qsub)
                os.system('qsub %s' % qsub)

        elif experiment == 'timing':
            results = 'exp_mpi_timing_' + threads
            qsub = 'run_mpi_timing_%s_nodes.qsub' % threads
            os.system('rm -r -f %s' % results)
            os.system('mkdir %s' % results)
            os.system('make')

            node_qsub('mpi-timing.x', range(int(start), int(stop), int(step)),
                      int(threads), results, qsub)

            os.system('qsub %s' % qsub)

        elif experiment == 'naive':
            results = 'exp_naive_timing_' + threads
            qsub = 'run_naive_timing_%s_nodes.qsub' % threads
            os.system('rm -r -f %s' % results)
            os.system('mkdir %s' % results)
            os.system('make')

            node_qsub('naive-paths.x', range(int(start), int(stop), int(step)),
                      int(threads), results, qsub)

            os.system('qsub %s' % qsub)
