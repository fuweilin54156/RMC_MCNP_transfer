# -*- coding:utf-8 -*-

import argparse


def parse(desc=""):
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('inp', metavar='inp', type=str, help="the path to the input file",
                        default="workspace/inp")
    parser.add_argument('--platform', metavar='PLATFORM', dest="platform", default="linux", type=str,
                        help="platform on which to run the job")
    parser.add_argument('--mpi', metavar="MPI_NUM", dest='mpi', default=4, type=int,
                        help="the number of MPI processes to use")
    parser.add_argument('--omp', metavar="OMP_NUM", dest='omp', default=1, type=int,
                        help="the number of OpenMP threads to use in each process")
    parser.add_argument('--assem', metavar="ASSEM_NUM", dest='assem', default=0, type=int,
                        help="the number of assemblies in the model, default not to calculate CTF")
    parser.add_argument('--continue-inp', metavar="CONTINUE_INP", dest='conti_inp', default="", type=str,
                        help="the name of the input file for continuous calculation")
    parser.add_argument('--continue', metavar="CONTINUE", dest="conti", action="store_const", const=True, default=False,
                        help="whether this calculation is a continuous one or not")

    return parser.parse_args()
