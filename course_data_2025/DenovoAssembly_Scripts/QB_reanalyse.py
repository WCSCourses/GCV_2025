#!/usr/bin/env python

"A Python re-factoring of the perl script qb_post_process.pl"

import argparse
import sys

import QuasiBAM

import itertools as it

################################################################################
"Parses the command line, overwriting args from config/QuasiBAM.yaml."


def CommandLine_args():

    ap = argparse.ArgumentParser()

    ap.add_argument(dest="tabular", help="/path/to/input/tabular/file")

    ap.add_argument(
        "-c",
        dest="cons",
        default=20,
        help="Consensus inclusion threshold frequency. [20]",
    )
    ap.add_argument(
        "-d",
        dest="depth",
        default=100,
        help="Consensus inclusion threshold depth. [100]",
    )
    # ap.add_argument('-n', dest='N_to_gap', action='store_true',
    #                 help='Report consensus Ns as gaps (-)')

    ap.add_argument(
        "-mg",
        dest="mask_gaps",
        action="store_false",
        help="Mask gaps where base/gap both >cons [True]",
    )

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit()

    args = ap.parse_args()

    name = args.tabular.rsplit('.', 1)[0]
    args.out_file = f"{name}.c{args.cons}.d{args.depth}.fas"

    args.cons = [*map(float, str(args.cons).split(","))]
    args.depth = [*map(int, str(args.depth).split(","))]

    if len(args.cons) == 1 or len(args.depth) == 1:
        args.cd = [*it.product(args.cons, args.depth)]
    elif len(args.cons) == len(args.depth):
        args.cd = [*zip(args.cons, args.depth)]

    return args


################################################################################


def QB_reanalyse():

    # Command Line args
    args = CommandLine_args()

    # Main body
    QuasiBAM.make_fasta(
        tabular=args.tabular,
        cd=args.cd,
        consolidate=True,
        mask_gaps=args.mask_gaps,
        out_file=args.out_file,
    )

    return 0


################################################################################

if __name__ == "__main__":

    sys.exit(QB_reanalyse())

################################################################################
