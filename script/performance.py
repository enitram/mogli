#! /bin/env python

import argparse
import functools

import itertools
import os

import subprocess

from multiprocessing import Pool


def run(ms, molid1, in_dir, gen_int, min_core, runs, timeout):
    molid2 = ms[0]
    shell = ms[1]

    values = []
    gen_type = {
            0: 'NOPT',
            1: 'UC',
            2: 'DEG1',
            3: 'UC_DEG1'
        }[gen_int]

    for run in xrange(runs):
	try:
            timings = subprocess.check_output(['/usr/bin/timeout', str(timeout)] +
                                              ['../build/report',
                                              '-x', str(molid1),
                                              '-y', str(molid2),
                                              '-i', in_dir,
                                              '-b', str(gen_int),
                                              '-m', str(min_core),
                                              '-s', str(shell)])
            if not timings:
                return
        except subprocess.CalledProcessError as e:
            if e.returncode != 124:
                traceback.print_exc()
            return ['{0}x{1}\t{2}\t{3}\t{4}\n'.format(molid1, molid2, run, shell, gen_type)]

        line = '{0}x{1}\t{2}\t{3}\t{4}\t{5}'.format(molid1, molid2, run, shell, gen_type, timings)
        values.append(line)
    return values


def main():
    parser = argparse.ArgumentParser(description='Measure performance.')
    parser.add_argument('-b', choices=['0','1','2','3'], required=True,
                        help='Product graph generation type\n'
                             '     0 - No opt\n'
                             '     1 - Uncon rule\n'
                             '     2 - Deg-1 rule\n'
                             '     3 - Uncon & Deg-1 rule')
    parser.add_argument('-i', type=str, required=True, help='Input directory')
    parser.add_argument('-o', type=str, required=True, help='Output directory')
    parser.add_argument('-f', type=int, required=False, default=0,
                        help='Start from this index')
    parser.add_argument('-m', type=int, required=False, default=3,
                        help='Minimal core size')
    parser.add_argument('-s', type=int, required=False, default=3,
                        help='Maximal shell size')
    parser.add_argument('-r', type=int, required=False, default=5,
                        help='Number of runs')
    parser.add_argument('-p', type=int, required=False, default=-1,
                        help='Number of processes')
    parser.add_argument('-t', type=int, required=False, default=600,
                        help='Timeout for one match')

    args = parser.parse_args()
    gen_int = int(args.b)
    in_dir = args.i
    out_dir = args.o
    proc = args.p

    min_shell = {
        0: 0,
        1: 0,
        2: 1,
        3: 1,
    }[gen_int]

    ext = {
        0: '_nopt.csv',
        1: '_uc.csv',
        2: '_deg1.csv',
        3: '_uc_deg1.csv',
    }[gen_int]

    molids = sorted([f[:-4] for f in os.listdir(in_dir) if f.endswith('.lgf')])

    with open(os.path.join(out_dir, 'results'+ext), 'w') as fi:
        fi.write('MATCH\tRUN\tSHELL\tGEN-TYPE\tNODES\tFRAGMENTS\tPG-TIME\tBK-TIME\tGEN-TIME\tRUNTIME\n')

    pool = Pool() if proc < 1 else Pool(proc)
    for i in xrange(min(args.f, len(molids)-2), len(molids)-1):
        print('%d/%d' % (i+1,len(molids)))
        lines = pool.map(functools.partial(run,
                                           molid1=molids[i],
                                           in_dir=in_dir,
                                           gen_int=gen_int,
                                           min_core=args.m,
                                           runs=args.r,
                                           timeout=args.t),
                         itertools.product(molids[i+1:], xrange(min_shell, args.s+1)))
        lines = filter(lambda x: x, lines)
        if len(lines) > 0:
            with open(os.path.join(out_dir, 'results'+ext), 'a') as fi:
                fi.writelines(itertools.chain.from_iterable(lines))

    pool.close()
    pool.join()

if __name__ == '__main__':
    main()
