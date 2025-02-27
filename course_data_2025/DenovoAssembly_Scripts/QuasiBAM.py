#!/usr/bin/env python

"""
A Python re-factoring of Richard Myers' C++ quasi_bam application. Most legacy
functionality has been retained, including directional analysis, error-checking
for read-terminal variant frequencies as well as read-derived codon counts. Mean
quality scores for each base included in the analysis are omitted. The ability
to produce frame-specific outputs based upon sub-target sequences is currently
unavailable. In its place is the facility to perform analysis on BAM files with
multiple targets. Tabular and sequence outputs are generated for all and only
those targets in the passed FASTA file - other RNAME targets in the BAM file are
ignored.

In this version, an expanded format for amino acid, codon and insert sequence
data is employed - rather than comma-separated identity:frequency pairs, the new
output format comprises space-separated identity:count:frequency triplets. No
data on the consensus is given in the tabular file as the command line options
preclude a single parameter for consensus frequency and depth.

In a further change, the multiprocessing library is used to facilitate
threading. A local <gU.threaded> function controls multiprocessing pooling using
map or itertools.starmap to pass lists of argument(s) to a function.
"""

import argparse
import functools
import logging
import os
import sys

import itertools as it
import numpy as np
import pandas as pd
import subprocess as sp

from collections import Counter, defaultdict
from glob import glob
from types import SimpleNamespace

import generalUtilities as gU

__version__ = "2.0.1gcv"    #   WCSC-specific version
__date__ = "240212"
__author__ = "David Bibby"


"functools.partials for 'pure' base functions and 'gap=0' indexes."
local_arr2seq = functools.partial(gU.arr2seq, _bases="-ACGT")
local_seq2arr = functools.partial(gU.seq2arr, _bases="-ACGT")
local_translate = functools.partial(gU.translate_pure, offset=1)


def sum_final_axis(arr):
    "Function to summate array data."

    return [np.sum(arr, axis=-1)]


def df_merge(left, right):
    "Local pd.merge for joining DataFrames"

    if right.empty:
        return left
    return pd.merge(left, right, on="locus")


def get_str(df, cols, name, f_thr, spacer=" "):
    "Triplet string pd.Series generator for codons, amino acids and insertions"

    d = df[df.freq >= f_thr].sort_values("freq", ascending=False)
    if d.empty:
        ser = pd.Series(data={k: "" for k in np.unique(df.index)}, name=name, dtype=str)
        ser.index.name = "locus"
        return ser
    d.freq = round(d.freq, 3)
    # Converts each entry into a colon-separated triplet
    d["strings"] = d.apply(lambda x: ":".join(map(str, x[cols])), axis=1)

    # Groups triplets and space-separates them in a grouped series
    return d.groupby("locus").strings.agg(spacer.join).rename(name)


# -------------------------------------------------------------------------------
def stack_offset(arr, n=3):
    """
    Takes 1D <arr> and stacks it <n> times, offsetting by 1 each time.
    Used to generate 3-fold stacked arrays for translating base-by-base.
    """

    return np.stack([arr[i : arr.size - n + 1 + i] for i in range(n)])


################################################################################
class local_read_array:
    "Object class for adding functionality to gU.read_array class."

    def __init__(self, sam_line, qual):

        # Get base gU.read_array data
        r_arr = gU.read_array(sam_line, bases="NACGT")
        self.p_arr, self.s_arr, q_arr, c_arr = r_arr()
        self.orientation = (int(sam_line.FLAG) & 16) >> 4  # 0 = Fwd, 1 = Rev

        # Erase soft-clips
        self.p_arr[c_arr == 4] = 0

        # Parse inserts
        self.i_list = list()
        i_arr = (c_arr == 8).nonzero()[0]
        for pos in np.unique(self.p_arr[i_arr]):
            indexes = self.p_arr[i_arr] == pos
            if (q_arr[i_arr][indexes] < qual).any():
                continue
            # Location within read (for error-checking)
            read_pos = min(*(self.p_arr == pos).nonzero())
            # Insertion sequence as string
            ins = local_arr2seq(self.s_arr[i_arr][indexes])
            self.i_list.append((read_pos, pos, ins))
            # Remove the insertion from the array once stored.
            self.p_arr[i_arr[indexes]] = 0

        # Erase low-qual bases
        self.p_arr[(q_arr < qual) & (c_arr != 1)] = 0

        # Find where loci are consecutive (i.e. non-zero). Then find where the
        # next locus is also non-zero to get positions of reportable codons.
        d_arr = np.diff(self.p_arr) == 1
        self.d_arr = d_arr[1:] * (~np.diff(d_arr))

    def __call__(self):
        "Returns POS, SEQ, codon indexes, insertions, and orientation"

        return (self.p_arr, self.s_arr, self.d_arr, self.i_list, self.orientation)


class counting_array:
    """
    Base class for counting nucleotides / codons. Orientation and terminal
    position are the last two dimensions, if present.
    """

    def __init__(self, dimensions):

        "Initialise a zeros array with <dimensions> shape"

        self.arr = np.zeros([*dimensions], dtype=np.int32)

    def inc(self, *coordinates):

        "Increment array"

        self.arr[coordinates] += 1

    def get_array(self):

        "Return array"

        return self.arr


class count_arrays:
    """
    Each instance nolds a nucleotide counter and a codon counter.
    - nt_arr    A 4D array of loci:base/gap:orientation:read-terminal counts
    - cd_arr    A 5D array of loci:*triplet bases:orientation counts

    __add__ parses the passed read_array and increments the counters appropriately.

    input
    -----
    seq_len     int     the length of the reference sequence
    read_len    int     the max length of a read
    terminal    float   the fraction of read ends considered terminal
    """

    def __init__(self, seq_len, read_len, terminal=0.0):

        self.nt_arr = counting_array(dimensions=(seq_len + 1, 5, 2, 2))
        self.cd_arr = counting_array(dimensions=(seq_len + 1, 5, 5, 5, 2))

        # Dimensionality of the insert array must account for orientation and
        # error-checking
        self.i_list = []

        self.er_arr = np.zeros(read_len, dtype=np.int8)
        self.tail = int(terminal * read_len)
        self.er_arr[: self.tail] = 1
        self.er_arr[read_len - self.tail :] = 1

    @gU.memo
    def error_array(self, size):

        "A cache of index arrays into the last dimension of the nt_arr"

        excess = size - self.er_arr.size
        if excess <= 0:
            return self.er_arr[:size]
        return np.insert(self.er_arr, self.er_arr.size // 2, [0] * excess)

    def get_arrays(self):

        "Simple data return (arrays)"

        return (self.nt_arr.get_array(), self.cd_arr.get_array())

    def get_inserts(self):

        "Simple data return (dict)"

        return self.i_list

    def add_read(self, read):

        "Add a read_array to the nt_ and cd_ counting arrays"

        # Get read data
        p_arr, s_arr, d_arr, i_list, ori = read()

        # Get error_array
        er_arr = self.error_array(p_arr.size)

        # Increment nucleotide array
        self.nt_arr.inc(p_arr, s_arr, er_arr, [ori])

        # Parse amino acid counts and increment codon array
        c1, c2, c3 = stack_offset(s_arr)[:, d_arr]
        self.cd_arr.inc(p_arr[:-2][d_arr], c1, c2, c3, [ori])

        # Increment the appropriate insert dictionary
        for (read_pos, locus, ins) in i_list:
            self.i_list.append((locus, ins, er_arr[read_pos], ori))

    def add_reads(self, read_list):

        for read in read_list:
            self.add_read(read)

        return self


################################################################################
def CommandLine_args():
    """
    Parses the command line, overwriting args from config/QuasiBAM.yaml.
    """

    ap = argparse.ArgumentParser()

    # Positional non-optional arguments
    ap.add_argument("bam", type=str, help="BAM file")
    ap.add_argument("fas", type=str, help="FASTA file containing references")

    # Flags
    # Analysis flags

    ap.add_argument(
        "-st",
        dest="strand",
        action="store_true",
        help="Turn ON independent strand analysis",
    )

    # Output FASTA formatting flags

    ap.add_argument(
        "-mg",
        dest="mask_gaps",
        action="store_false",
        help="Mask gaps where base/gap both >cons",
    )

    # Variant position error checking
    ap.add_argument(
        "-e", dest="error_check", action="store_false", help="Turn OFF error-checking"
    )

    # Parameters
    # Variant position error checking
    ap.add_argument(
        "-e1",
        dest="read_end",
        type=int,
        default=6,
        help="Number of terminal bases to check for terminal errors [6]",
    )
    ap.add_argument(
        "-e2",
        dest="read_freq_thr",
        type=float,
        default=0.9,
        help="Error threshold frequency [0.9]"
    )

    # Analysis inclusion thresholds
    ap.add_argument(
        "-m",
        dest="mapq",
        type=int,
        default=30,
        help="Analysis inclusion threshold MAPQ (read) [30]"
    )
    ap.add_argument(
        "-mc",
        dest="mapl",
        type=int,
        default=50,
        help="Analysis inclusion threshold LEN (read) [50]"
    )
    ap.add_argument(
        "-q",
        dest="qual",
        type=int,
        default=34,
        help="Analysis inclusion threshold QUAL (base) [34]"
    )

    # Reporting and consensus inclusion thresholds, and combination flag
    ap.add_argument(
        "-f",
        dest="freq",
        type=float,
        default=1.0,
        help="Reporting threshold frequency [1.0]"
    )
    
    ap.add_argument(
        "-c",
        dest="cons",
        type=str,
        default="20",
        help="Consensus inclusion threshold frequency. To apply "
        "multiple values, separate by comma [20.0]",
    )
    
    ap.add_argument(
        "-d",
        dest="depth",
        type=str,
        default="100",
        help="Consensus inclusion threshold depth. To apply "
        "multiple values, separate by comma [100]",
    )
    ap.add_argument(
        "-cd",
        dest="product",
        action="store_true",
        help="If set, all combinations of CONS and DEPTH are "
        "applied. Otherwise, paired entries only are taken.",
    )

    # New arguments
    # Output filename basename
    ap.add_argument(
        "-s",
        dest="stem",
        help="Output filestem"
    )

    ap.add_argument(
        "-t",
        dest="threads",
        type=int,
        default=4,
        help="Threads to use when parallelising [4].",
    )

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit()

    args = ap.parse_args()

    return args


# -------------------------------------------------------------------------------
def QuasiBAM_args():
    """
    Adds module-level variables to the <args> SimpleNamespace. Performs basic
    validation of passed arguments and input files.
    """

    ######################################
    # Validate the CONS and DEPTH values #
    ######################################

    args.cons = [*map(float, str(args.cons).split(","))]
    args.depth = [*map(int, str(args.depth).split(","))]

    if len(args.cons) == 1 or len(args.depth) == 1 or args.product:
        args.cd = [*it.product(args.cons, args.depth)]
    elif len(args.cons) == len(args.depth):
        args.cd = [*zip(args.cons, args.depth)]
    else:
        args.log.error(f"{args.fas}: CONS, DEPTH and -cd flag are incompatible.")
        args.log.error(f"\tCONS: {args.cons}.")
        args.log.error(f"\tDEPTH: {args.depth}.")
        args.log.error(f"\tCD FLAG: {args.product}.")
        return 1

    args.log.debug(f'{args.fas}: CONS and DEPTH: {" ".join(map(str, args.cd))}')
    ###################################################
    # Test for co-ordinate-sorted BAM file and index. #
    ###################################################

    # Confirm BAM file is a properly-headered BAM file.
    # try:
    pc = gU.samtools("view", H=args.bam, func=sp.run)
    # except Exception:
    #     args.log.error(f"{args.bam} cannot be read by samtools. Exiting")
    #     return 1

    # Test for co-ordinate ordering of the BAM file. Sort if necessary.
    if not [*filter(lambda x: "SO:coordinate" in x, pc.stdout.split("\n"))]:
        args.log.debug(f"{args.bam}: Not sorted by co-ordinate. Sorting.")
        tmp = f"{args.bam}.tmp"
        gU.samtools("sort", args.bam, func=sp.run, stdout=open(tmp, "w"))
        sp.run(["mv", tmp, args.bam])
        args.log.debug(f"{args.bam} now sorted by co-ordinate.")

    # Index the BAM.
    gU.samtools("index", args.bam, func=sp.run)

    # Get the read length for error-checking.
    args.read_len = gU.get_read_length(args.bam)
    if not args.read_len:
        args.log.error(f"No reads in {args.bam}. Exiting.")
        return 1
    if not args.error_check:
        args.read_end = 0
    elif args.read_end > 1:
        args.read_end = args.read_end / args.read_len

    # Find the mapped targets from the BAM header.
    pc = gU.samtools("view", "-H", args.bam, func=sp.run)
    RNAMEs = set(
        map(
            lambda x: x.split(":")[-1],
            filter(lambda x: x.startswith("SN:"), pc.stdout.split()),
        )
    )

    # Filter the FASTA by mapped targets only, and identify subsequence ORFs
    args.fastas = []
    args.orfs = defaultdict(list)
    for n, s in gU.fasta_parser(args.fas):
        if n in RNAMEs:
            args.fastas.append((n, s))
        else:
            for name, seq in args.fastas:
                if s in seq:
                    args.orfs[name].append((n, seq.index(s) - 1, len(s)))

    if len(args.fastas) == 0:
        args.log.error(f"None of {args.fas} mapped in {args.bam}. Exiting")
        return 1

    if len(args.fastas) > 1 and args.mapq > 0:
        # When mapping to multiple similar references, MAPQ is often set to 0
        # for reads with strong secondary mapping.
        args.log.info(
            f"{args.fas}: Mutliple references detected. Specifying "
            "MAPQ > 0 may result in substantial loss of data."
        )

    args.log.debug(f"BAM file: {args.bam}.")
    args.log.debug(f"FASTA reference(s) file: {args.fas}.")
    args.log.debug(f"{len(args.fastas)} reference(s) to be parsed.")
    args.log.debug(
        f"{sum(len(j) for i, j in args.orfs.items())} sub-sequence(s) to be parsed."
    )
    args.log.debug(f"Max read length: {args.read_len}.")
    args.log.debug(
        f"Terminal bases to error-check: {int(args.read_end * args.read_len)}"
    )


################################################################################
def process_ref(name, seq):
    """
    Manages the stages of processing a reference. Iterates over parsing BAM,
    counting nt/codons, error-checking and output generation.
    """

    args.log.info(f"{name}: Processing")

    # Establish the base DataFrame from the reference sequence.
    RefN = pd.Series(list(seq), index=np.arange(len(seq)) + 1, name="Ref_N")
    RefAA = pd.Series(
        map(
            lambda tr: "X" if len(tr) > 4 else "".join(sorted(tr)),
            map(gU.translate_mixed, stack_offset(gU.seq2arr(seq)).T),
        ),
        index=np.arange(len(seq) - 2) + 1,
        name="Ref_AA",
    )

    base_df = pd.concat((RefN, RefAA), axis=1).fillna("")

    ###################################################################
    #     PARSE SAM-LINES INTO READ_ARRAY INSTANCES (parallelised)    #
    ###################################################################

    # Get the lines from <args.bam> with RNAME==<name> as a generator
    pc = gU.samtools("view", args.bam, name, func=sp.run)
    SAM_lines = pc.stdout.split("\n")
    # if len(SAM_lines) < args.min_reads:
    #     args.log.error(f"{name}: Insufficient reads mapped")
    #     return 1

    r_arrs = gU.threaded(
        func=local_read_array,
        data=[
            *zip(filter(None, map(gU.get_sam_line, SAM_lines)), it.repeat(args.qual))
        ],
        procs=4,
    )

    args.log.debug(f"{name}: {len(r_arrs)} read arrays created.")

    ######################################################################
    #    PARSE READ_ARRAYS INTO COUNT_ARRAYS INSTANCES (parallelised)    #
    ######################################################################

    # Initialise the count_arrays instances (one for each thread), and
    # parallel process the read_arrays.
    sublists = [
        (
            count_arrays(len(seq), args.read_len, args.read_end),
            r_arrs[i :: args.threads],
        )
        for i in range(args.threads)
    ]
    c_arrs = gU.threaded(func=count_arrays.add_reads, data=sublists, procs=4)

    # Summate both the parallel nt and codon arrays from the count_arrays
    arrs = [
        *map(lambda x: np.sum(x, axis=0), zip(*map(count_arrays.get_arrays, c_arrs)))
    ]

    args.log.debug(f"{name}: Nucleotide and codon count arrays merged.")

    #######################################################################
    #    PROCESS THE COUNTING ARRAYS + ERROR CHECK AND STRAND ANALYSIS    #
    #######################################################################

    # Sum the orientation axes (-1 in both cases) in any case.
    # If stranded analysis is required, also split each array and append.
    def process_arrs(arr):
        arrs = sum_final_axis(arr)
        if args.strand:
            arrs.extend([*map(np.squeeze, np.split(arr, 2, axis=-1))])
        return arrs

    arrs = [*map(process_arrs, arrs)]

    # Parse each set of arrays and perform error-checking if required. Output
    # contains two arrays: depth/freq & cd_arr.
    arrs = zip(*gU.threaded(func=parse_counts, data=[*zip(*arrs)], procs=4))

    args.log.debug(f"{name}: Error-checking complete.")

    #######################################################################
    #    PROCESS THE INSERTION LISTS - ERROR CHECK AND STRAND ANALYSIS    #
    #######################################################################

    # Convert the gU.threaded insertion dictionaries to a Counter/
    fqd_arrs = next(arrs)
    ins_C = Counter(gU.chain(map(count_arrays.get_inserts, c_arrs)))

    # Process the Counter, using the Depth column(s) of fqd_arrs.
    ins_dfs = parse_ins(ins_C, np.dstack(fqd_arrs)[:, -1])
    args.log.debug(f"{name}: Insertions processed.")

    #################################################################
    #    PARSE CODON COUNT_ARRAYS INTO DATAFRAMES (parallelised)    #
    #################################################################

    def parse_cd_arr(arr):
        "Splits a codon array into chunks and converts to strings."
        loci_arrs = np.array_split(np.unique(arr.nonzero()[0]), args.threads)
        chunks = [(loci_arr, arr[loci_arr]) for loci_arr in loci_arrs]
        parsed_arr = pd.concat(gU.threaded(func=parse_codons, data=chunks, procs=4))
        return parsed_arr

    # Create DataFrames of codon/amino acid triplet strings from codon array(s).
    codon_dfs = [*map(parse_cd_arr, next(arrs))]
    args.log.debug(f"{name}: Codon data processed.")

    #########################################################
    #    SIMPLE ANNOTATION OF THE FREQUENCY/DEPTH ARRAYS    #
    #########################################################

    freq_dfs = [*map(freq_parse, fqd_arrs)]
    args.log.debug(f"{name}: Frequency/depth data processed.")

    #############################################
    #    CREATE OUTPUT TABULAR AND FAS FILES    #
    #############################################

    # Create the tabular file(s)
    args.log.debug(f"{name}: Producing tabular file(s)")

    # If args.strand, each of the *_dfs will be length=2, otherwise 1.
    stem = args.stem if args.stem else name
    data = [
        (e, pd.concat((base_df, *dfs), axis=1).fillna(""), stem)
        for e, dfs in enumerate(zip(freq_dfs, ins_dfs, codon_dfs))
    ]
    gU.threaded(func=make_tabular, data=data, procs=4)

    # Create the FASTA file(s)
    args.log.debug(f"{name}: Producing FASTA file(s)")
    gU.threaded(func=export_fasta, data=glob(f"{stem}*.tabular"), procs=4)

    args.log.debug(f"{name}: QuasiBAM analysis complete")


# -------------------------------------------------------------------------------
def parse_counts(nt_arr, cd_arr):
    "Generates frequency, depth & codon arrays. Performs error-checking."
    # Combine error-checked and non-checked counts
    ntd_arr = sum_final_axis(nt_arr)[0]

    # Sum the combined counts along the nt axis - depth array
    d_arr = sum_final_axis(ntd_arr)[0][:, None]

    # Create a frequency array from the nt array
    fq_arr = np.zeros_like(ntd_arr, dtype=np.float64)
    np.true_divide(ntd_arr * 100, d_arr, out=fq_arr, where=(d_arr > 0))
    fq_arr[fq_arr < args.freq] = 0

    # Replace the nt array with its sum along the error axis (-1)
    if args.error_check:

        er_arr = np.zeros_like(fq_arr)
        np.true_divide(nt_arr[..., 1], ntd_arr, out=er_arr, where=(fq_arr > 0))

        # Find non-zero coverage of the mapping (may not be full reference!)
        (cv_loci,) = np.where(np.squeeze(d_arr))
        tail = int(args.read_end * args.read_len)

        # Set the error array to zero at the terminal X bases of the coverage
        er_arr[: cv_loci[1] + tail] = 0
        er_arr[cv_loci[-1] - tail :] = 0

        # Exclude both ends, where *all* nucleotides will be read-terminal
        loci, bases = (er_arr >= args.read_freq_thr).nonzero()

        # Set the nt frequencies to zero at the erroring loci
        fq_arr[(loci, bases)] = 0

        # Set the three affected codons to zero
        cd_arr[(loci, bases)] = 0
        cd_arr[(loci - 1), :, bases] = 0
        cd_arr[(loci - 2), ..., bases] = 0

    return (np.hstack((fq_arr, d_arr)), cd_arr)


# -------------------------------------------------------------------------------
def parse_ins(C, arr):
    "Generates insert DataFrame(s). Performs error-checking."

    def no_insertions():
        df = pd.DataFrame(index=np.arange(arr.shape[0])[1:])
        df["Ins"] = ""
        df["I_Desc"] = ""
        return df

    # If there are no insertions at all!
    if not C:
        return [no_insertions()] * arr.shape[1]

    index = ["locus", "ins", "tail", "ori"]

    # Prepare a DataFrame from the insert Counter with a four-level index
    df = pd.DataFrame.from_dict(C, orient="index", columns=["n"])
    df.set_index(pd.MultiIndex.from_tuples(df.index, names=index), inplace=True)
    df.reset_index(inplace=True, level="ori")

    def parse_inserts(e):
        try:
            d = df[~(df.ori.isin([2 - e]))].groupby(index[:-1]).sum().reset_index()
            if args.error_check:
                # Pivot on the error-checking flag ('tail')
                d = d.pivot(index=index[:-2], columns="tail", values="n").fillna(0)
                d["n"] = d.sum(axis=1).astype(np.int32)
                # Eliminate the erroring counts
                d = d[d[1] / d.n < args.read_freq_thr].n.reset_index()
            d["freq"] = 100 * d.n / arr[d.locus, e]
            I_Desc = get_str(d, ["ins", "n", "freq"], "I_Desc", args.freq)
            d = d.groupby("locus").sum()
            d = df_merge(d[d.freq >= args.freq], I_Desc)[["freq", "I_Desc"]]
            d.rename(columns={"freq": "Ins"}, inplace=True)
            return d
        except Exception:
            return no_insertions()

    return [*map(parse_inserts, range(arr.shape[1]))]


# -------------------------------------------------------------------------------
def parse_codons(indices, cd_arr):
    """
    Parses a codon array into a DataFrame of triplet & amino acid counts &
    frequencies.
    """

    nz_arr = cd_arr.nonzero()

    # Create "codon" DataFrame of loci, triplet base IDs, and counts
    cd_arr = np.stack((indices[nz_arr[0]], *nz_arr[1:], cd_arr[nz_arr])).T
    codon = pd.DataFrame(cd_arr, columns=["locus", "c1", "c2", "c3", "n"])

    # Create "depth" Series and join back to "codon"
    depth = codon.groupby("locus").n.sum().rename("AA_depth")
    codon = codon.join(depth, on="locus")

    # Add relevant columns to "codon" (seqs, amino acids, freq)
    codon["triplet"] = codon.apply(lambda x: [x.c1, x.c2, x.c3], axis=1)
    codon["Cod"] = codon.triplet.apply(local_arr2seq)
    codon["AA"] = codon.triplet.apply(local_translate)
    codon["freq"] = 100 * codon.n / codon.AA_depth

    # Create "amino acid" DataFrame by reducing redundant triplets in "codon"
    amino_acid = codon.groupby(by=["locus", "AA", "AA_depth"]).sum().reset_index()

    # Convert "codon" to Series of spaced formatted frequency data strings
    codon = get_str(codon, ["Cod", "n", "freq"], "Cod", args.freq)

    # Create "consensus" Series from "amino_acid" comprising strings of code(s)
    # cons = get_str(amino_acid, ['AA'], 'Cons_AA', args.cons, spacer='')

    # Convert "amino_acid" to Series of spaced formatted frequency data strings
    amino_acid = get_str(amino_acid, ["AA", "n", "freq"], "AA", args.freq)

    # All series are of the same length and iteratively merged into a DataFrame
    # df = reduce(df_merge, (depth, cons, codon, amino_acid))
    df = functools.reduce(df_merge, (depth, codon, amino_acid))

    return df


# -------------------------------------------------------------------------------
def freq_parse(arr):
    "Converts a frequency/depth array into a DataFrame."

    df = pd.DataFrame(arr, columns=["Gap", "A", "C", "G", "T", "Depth"])
    df.index.name = "locus"

    return df.astype({"Depth": np.int32}).loc[1:]


# -------------------------------------------------------------------------------
def make_tabular(e, df, name):
    "Generates a tabular file from a combined DataFrame. Strand-aware."

    columns = [
        "Ref_N",
        "Depth",
        "A",
        "C",
        "G",
        "T",
        "Gap",
        "Ins",
        "I_Desc",
        "Ref_AA",
        "AA_depth",
        "Cod",
        "AA",
    ]

    def export_tabular(df, filename):
        args.log.debug(f"\t{filename}")
        df.to_csv(filename, sep="\t", columns=columns)

    # Integer-ise the depth fields
    df.loc[df.AA_depth == "", "AA_depth"] = 0
    df.loc[df.Depth == "", "Depth"] = 0
    df = df.astype({"AA_depth": np.int32, "Depth": np.int32})

    # Export basic tabular files
    suffix = ("", "fwd.", "rev.")[e]
    tabular = f"{name}.{suffix}tabular"
    export_tabular(df.loc[1:], tabular)

    # Export ORF tabular file(s) if any ORFs are specified.
    columns[9:9] = ["AA_pos"]
    for (orf, start, length) in args.orfs[name]:
        sub_df = df.iloc[start + 1 : start + length + 1].copy()
        arr = np.arange(length)
        sub_df.index = arr + 1
        sub_df["AA_pos"] = (arr // 3) + 1
        sub_df.loc[
            (sub_df.index % 3) != 1, ["AA_pos", "Ref_AA", "AA_depth", "Cod", "AA"]
        ] = ""
        export_tabular(sub_df, tabular.replace("tabular", f"{orf}.tabular"))

    return tabular


# -------------------------------------------------------------------------------
def export_fasta(tabular):
    "Internal wrapper for <make_fastas>. (See QB_reanalyse.py)"

    make_fasta(
        tabular=tabular,
        cd=args.cd,
        log=args.log,
        mask_gaps=args.mask_gaps,
        consolidate=not args.product,
    )


# -------------------------------------------------------------------------------
def make_fasta(tabular, cd, consolidate=True, log=None, mask_gaps=False, out_file=None):
    """
    Generates FASTA(s) from a tabular file and cons/depth parameters.
    Consolidates data into a single FASTA if requested.
    """

    # <cons_depth> is a list of (cons, depth) tuples
    if not log:
        log = gU.dummy_log()

    def apply_inserts(arr, inserts):
        for locus, seq in sorted(inserts, reverse=True):
            arr = np.insert(arr, locus, seq)
        return arr

    bases = gU.bases + gU.bases.lower()
    df = pd.read_csv(tabular, sep="\t", index_col=0)
    name = gU.filestem(out_file if out_file else tabular)

    for col in df.columns:
        fill_val = 0.0
        if col in ("Ref_N", "I_Desc", "Ref_AA", "Cod", "AA"):
            fill_val = ""
        elif col in ("Depth", "AA_depth"):
            fill_val = 0
        df[col].fillna(fill_val, inplace=True)

    arrs = []
    inserts = set()

    for cons, depth in cd:
        # Establish basic sequence from base/gap frequencies
        arr = np.full(len(df), 15, dtype=np.int8)
        b_arr = np.array((df.loc[:, "A":"Gap"] >= cons) * [2, 4, 8, 16, 1]).sum(axis=1)
        arr[~(b_arr == 0)] = b_arr[~(b_arr == 0)] >> 1

        # Return regions of low depth to Ns
        arr[df.Depth < depth] = 15

        # If mixed gap/base is flagged, report in lower case
        arr[(df.Gap >= cons) & (arr > 0)] += 16 * mask_gaps

        # Check Ins/I_Desc and include those satisfying cons/depth thresholds
        for locus, insertions in df.loc[df.Ins > 0, "I_Desc"].iteritems():
            ins, count, freq = insertions.split()[0].split(":")
            if float(freq) >= cons and int(count) >= depth:
                inserts.add((locus, tuple(gU.seq2arr(ins))))

        if not consolidate:
            # Write individual cons/depth fastas as they are derived
            if out_file:
                log.error(
                    f"Multiple cons/depth permutations passed - {out_file} not written."
                )
            fas = f"{name}.c{cons}.d{depth}.fas"
            arr = apply_inserts(arr, inserts)
            gU.write_fastas(
                filename=fas, fastas=[(fas[:-4], gU.arr2seq(arr, _bases=bases))]
            )
            log.debug(f"\t- {fas}")
            inserts = set()
        else:
            arrs.append(arr)

    if consolidate:
        # Combine FASTAs to generate a composite FASTA.
        out_file = out_file if out_file else f"{name}.fas"
        arr = np.stack(arrs).astype(np.int16)
        arr[arr == 15] = 256
        arr = np.bitwise_or.reduce(arr)
        arr[arr == 256] = 15
        arr = apply_inserts(arr, inserts).astype(np.int8)
        gU.write_fastas(
            filename=out_file, fastas=[(name, gU.arr2seq(arr, _bases=bases))]
        )
        log.info(f"\t{out_file}")


################################################################################
def QuasiBAM():
    "Main function."

    # CommandLine args - these override the workflow defaults
    _args = CommandLine_args()
    for k, v in vars(_args).items():
        args.__dict__[k] = v
    args.__dict__["log"] = logging.getLogger("QuasiBAM")
    
    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(logging.INFO)
    streamhandler.setFormatter(formatter)
    args.log.addHandler(streamhandler)
    args.log.setLevel(5)
    
    
    # Script-specific args
    if QuasiBAM_args():
        args.log.error(f"{args.fas}: Exiting. Check error log.")
        return 1

    # Main body
    for name, seq in sorted(args.fastas):
        process_ref(name, seq)

    return 0


################################################################################
if __name__ == "__main__":

    args = SimpleNamespace()
    sys.exit(QuasiBAM())

################################################################################
