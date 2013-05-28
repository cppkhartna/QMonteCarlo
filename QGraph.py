#!/usr/bin/env python2
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt
import pylab as p
import os
import sys
from ctypes import *
import numpy as np
import scipy as sp
import argparse

class replica(Structure):
    pass
replica._fields_ = [
        ("x", POINTER(c_double)),
        ("next", POINTER(replica))
        ]

def parabola(xs, ys, err):
    plt.figure()
    plt.errorbar(xs, ys, yerr=err, fmt='o')
    plt.show()

def draw2d(energies):
    ax = plt.subplot(111)
    xs = []
    ys = []
    i = 0
    for e in energies:
        xs.append(i)
        ys.append(e)
        i += 1
    ax.plot(xs, ys)
    plt.show()

def draw3d(replicas, mol):
    fig=p.figure()
    ax = p3.Axes3D(fig)
    ax.autoscale(True)

    xs = []
    ys = []
    zs = []

    rep = replicas.contents

    try:
        while True:
            i = rep.x
            xs.append(i[0])
            ys.append(i[1])
            zs.append(i[2])
            rep = rep.next.contents
    except ValueError:
        pass

    if (mol == 3):
        rep = replicas.contents
        try:
            while True:
                i = rep.x
                xs.append(i[3])
                ys.append(i[4])
                zs.append(i[5])
                rep = rep.next.contents
        except ValueError:
            pass

    ax.scatter3D(xs, ys, zs)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    p.show()

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest = 'aim', type=int, nargs=1,
            help='1 - atom; 2 - ion; 3 - molecule')
    parser.add_argument('-c', dest = 'rep', type=int, nargs=1,
            help='initial number of replicas')
    parser.add_argument('-i', dest = 'itr', type=int, nargs=1,
            help='iterations')
    parser.add_argument('-r', dest = 'rad', type=float, nargs=1,
            help='ion/molecule radius')
    parser.add_argument('-g', dest = 'gra', type=int, nargs=1,
            help='graph type 2d/3d (0/1)')
    parser.add_argument('-p', dest = 'par', type=int, nargs=1,
            help='graph of parabole')
    args = parser.parse_args()

    ROOT = os.path.dirname(os.path.abspath(__file__))
    qlib = CDLL('%s/qlib.so' % ROOT)
    run  = qlib.run
    run.restype = POINTER(replica)

    atom = 1
    ion = 2
    molecule = 3

    if (args.aim != None and args.aim[0] >= 1 and args.aim[0] <= 3
            and args.rep != None 
            and args.itr != None 
            and args.rad != None
            and args.gra != None):
        # create array for doubles {{{
        es_type = c_double * args.itr[0]
        es = es_type()
        for i in range(args.itr[0]):
            es[i] = i
        # }}}


        if (args.par == None):
            replicas = run(c_int(args.rep[0]), c_int(args.itr[0]),
                    c_int(args.aim[0]), c_double(args.rad[0]), es)
            if (args.gra[0] == 0):
                draw3d(replicas, args.aim[0])
            else:
                draw2d(es)
        else:
            energies = []
            varz = []
            means = []
            magic = 3000
            if (args.itr[0] < magic):
                print "There are not enough iterations (%s)" % magic
                sys.exit()
            xs = np.arange(1.3, 3.0, 0.1)
            for rad in xs:
                replicas = run(c_int(args.rep[0]), c_int(args.itr[0]),
                        c_int(args.aim[0]), c_double(rad), es)
                varz.append(np.var(es[magic:]))
                means.append(np.mean(es[magic:]))
            parabola(xs, means, varz)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
