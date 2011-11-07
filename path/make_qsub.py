import sys

def make_qsub(path, threads, outpath, qsubpath):
    buffer = []
    buffer.append('#!/bin/bash')
    buffer.append('#')
    buffer.append('#$ -cwd')
    buffer.append('#$ -j y')
    buffer.append('#$ -S /bin/bash')
    buffer.append('')
    buffer.append('export OMP_NUM_THREADS=%s' % threads)
    buffer.append('./%s > %s' % (path, outpath))
    buffer.append('exit 0')

    file = open(qsubpath, 'w')
    for line in buffer:
        file.write(line + '\n')
    file.close

# threads outpath qsubpath path
# path passes along remainder of command line arguments to function
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Arguments: threads outpath qsubpath path'
    else:
        print sys.argv
        threads = sys.argv[1]
        outpath = sys.argv[2]
        qsubpath = sys.argv[3]

        path = ' '.join(sys.argv[4:])

        make_qsub(path, threads, outpath, qsubpath)


