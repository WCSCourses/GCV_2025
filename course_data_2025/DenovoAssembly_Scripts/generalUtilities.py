"A set of functions used throughout Genomancer"

import logging
import os
import re
import sys

from collections import namedtuple
from glob import glob
from random import shuffle

import itertools as it
import multiprocessing as mp
import numpy as np
import pandas as pd
import subprocess as sp

__author__ = "David F Bibby"
__version__ = "2.2.2gcv"    # WCSC-specific
__date__ = "240212"

#################################
#     SHORT FUNCS & ALIASES     #
#################################

chain = it.chain.from_iterable

def get_scriptPath(fname):
    return os.path.abspath(os.path.dirname(fname))

def iter_zip(n, s):
    return zip(*[iter(s)] * n)

def dummy_log(name="dummy"):
    return logging.getLogger(name)

#######################
#      CONSTANTS      #
#######################

# Nucleotides, amino acids, FASTQs, colors

bases = "-ACMGRSVTWYHKDBN"
amino_acids = "ACDEFGHIKLMNPQRSTVWY*"
codon_arr = np.array(
    [
        [[8, 11, 8, 11], [16, 16, 16, 16], [14, 15, 14, 15], [7, 7, 10, 7]],
        [[13, 6, 13, 6], [12, 12, 12, 12], [14, 14, 14, 14], [9, 9, 9, 9]],
        [[3, 2, 3, 2], [0, 0, 0, 0], [5, 5, 5, 5], [17, 17, 17, 17]],
        [[20, 19, 20, 19], [15, 15, 15, 15], [20, 1, 18, 1], [9, 4, 9, 4]],
    ],
    dtype=np.int8,
)

cigar_codes = {"D": 1, "M": 2, "S": 4, "I": 8}

class read_array:
    """Parses the sam_line into arrays, and stores the QUAL value"""

    def __init__(self, sam_line, bases=bases, qual=0):

        self.sam_line = sam_line

        # c_arr: CIGAR
        self.c_arr = cigar2arr(sam_line.CIGAR)

        # p_arr: POS
        self.p_arr = np.full_like(self.c_arr, int(sam_line.POS) - 1, dtype=np.int32)
        self.p_arr += np.cumsum((self.c_arr != 8))
        self.p_arr -= (self.c_arr != 4).nonzero()[0][0]

        # s_arr: SEQ
        self.s_arr = np.zeros_like(self.c_arr)
        self.s_arr[~(self.c_arr == 1)] = seq2arr(sam_line.SEQ, _bases=bases)

        # q_arr: QUAL
        self.q_arr = np.zeros_like(self.c_arr)
        self.q_arr[~(self.c_arr == 1)] = qual2arr(sam_line.QUAL)

        # other
        self.orientation = (int(sam_line.FLAG) & 16) >> 4

    def __call__(self):

        return (self.p_arr, self.s_arr, self.q_arr, self.c_arr)

def memo(func):
    """A memoisation decorator"""

    func_dict = {}

    def outFunc(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in func_dict:
            func_dict[key] = func(*args, **kwargs)
        return func_dict[key]

    return outFunc


def threaded(func, data, procs=1):
    """
    Parallelises tasks using a multiprocessing pool

    func:  The function to be parallelised
    data:  A list of data items to be passed to each process.
    procs: The number of parallel processes to run concurrently [1]

    Outputs are the return values from func for each item in groups. Because the
    pool runs asynchronously, the order of outputs will not correspond to that
    of data
    """

    with mp.Pool(processes=procs) as pool:

        # Check for iterable group members - single values need <map>
        pool_func = pool.map
        if data and type(data[0]) in (list, tuple):
            pool_func = pool.starmap
        # Shuffles the data order to avoid uneven queues
        shuffle(data)
        outputs = pool_func(func, data)

    return outputs

################################################################################
#                         FASTA MANIPULATION TOOLS                             #
################################################################################


def fasta_parser(filename):
    """
    An alternative to Bio.SeqIO.FastaIO.SimpleFastaParser.

    filename: Can be either a filepath of type string or a file handle.
    """

    def getStream(handle):
        if type(handle) == str:
            return open(handle)
        return handle

    def is_not_header(x):
        return not x.startswith(">")

    with getStream(filename) as iterLines:
        for line in it.dropwhile(is_not_header, iterLines):
            header = line[1:]
            break
        else:
            return

        sequence = ""
        for line in iterLines:
            if is_not_header(line):
                sequence += line.strip().replace("\r", "")
            else:
                yield (header.rstrip(), sequence)
                header = line[1:]
                sequence = ""

    yield (header.rstrip(), sequence)


def write_fastas(filename, fastas, mode="w"):
    """
    Writes FASTAs to a named file.

    filename: The destination file
    fastas:   The FASTAs. These can be either a list of (header, sequence)
              tuples/lists, or a dictionary of header:sequence key:value pairs.
    mode:     Set to 'a' to append FASTAs to a file. (Does not check for newline
              at the end of the destination file)
    """

    if type(fastas) == dict:
        fastas = list(fastas.items())
    with open(filename, mode) as out:
        out.write("".join(f">{n}\n{s}\n" for n, s in fastas))

    return filename


def seq2arr(seq, _bases=bases):
    """
    Transforms a sequence into a numpy array with bases encoded into 8-bit
    integers.

    seq:   The sequence (str). Terminal whitespaces are stripped.
    bases: A string of bases. The index in this string of each character in seq
           (upper case) is the encoded integer.

    Outputs a numpy array of length len(seq) and type np.int8, unless a
    character in seq is not present in bases, upon which None is re
    """

    if _bases.upper() == _bases:
        seq = seq.upper()
    seq = seq.strip()
    if re.search(f"[^{_bases}]", seq):
        return None
    arr = np.fromiter(map(_bases.index, seq), dtype=np.int8)

    return arr


def arr2seq(arr, _bases=bases):
    """
    The inverse function to seq2arr - converts a numpy 1D array of type np.int8
    into a sequence string. This function cannot be memoised as numpy arrays are
    represented incompletely as str (the args/kwargs key for the func_dict).

    arr:    The numpy array of np.int8 values.
    bases:  The string into which the array values index.

    Outputs a string of length arr.size.
    """

    seq = "".join(_bases[b] for b in arr)

    return seq


@memo
def translate_mixed(codon):
    """
    Takes a codon in seq2arr-style encoding, and translates the amino acid(s).

    codon: A numpy array of length 3, representing the encoded bases of a codon
           triplet *with the default bases string as the index*.

    Outputs a set of single-letter amino acid codes comprising all possible
    translations of the codon. If any member of the codon triplet is a gap, i.e.
    encoded by a zero (or less), the set condenses to {"X"}
    """

    if (codon < 1).any():
        return set("X")

    codons = it.product(*(
        np.nonzero(x[:3:-1])[0]
        for x in map(
            np.unpackbits,
            codon.astype(np.uint8)
        )
    ))

    translation = set(map(translate_pure, codons))

    return translation


@memo
def translate_pure(arr, offset=0):
    """
    Translates a known pure-base codon in the form of an array with values 0-3
    representing ACGT.

    arr:    A pure-base codon in 'ACGT'.index form.
    offset: If the array was set up indexing '-ACGT' to account for the
            possibility of gaps, use offset=1 to compensate.

    Outputs the single-letter code amino acid corresponding to the codon. If any
    member of the triplet is not in ACGT (e.g. gap, 'N'), returns "X".
    """

    arr = np.array(arr) - offset
    if (arr < 0).any():
        return "X"
    aa = amino_acids[codon_arr[tuple(arr)]]

    return aa


################################################################################
#                          FASTQ MANIPULATION TOOLS                            #
################################################################################


def interleave_fastq(fastqs, *args, **kwargs):
    """Interleaves (possibly gzipped) paired end FASTQs"""

    cmd = ["interleave-reads.py"]
    
    pc = modify_and_run(cmd, *args, final=fastqs, **kwargs)


################################################################################
#                   BWA AND SAMTOOLS TYPE & FUNCTIONS                          #
################################################################################


def get_read_length(bam):
    """Returns the max length of the reads in bam. None if no reads."""

    pc = samtools("view", bam)
    pc = shell("head", n=1000, stdin=pc.stdout, func=sp.Popen)
    pc = awk("{if(length($10)>m)m=length($10)}END{print m}", stdin=pc.stdout)

    read_length = pc.stdout.strip()
    if read_length:
        read_length = int(read_length)

    return read_length


sam_line = namedtuple(
    "sam_line",
    [
        "QNAME",
        "FLAG",
        "RNAME",
        "POS",
        "MAPQ",
        "CIGAR",
        "RNEXT",
        "PNEXT",
        "TLEN",
        "SEQ",
        "QUAL",
    ],
)


def get_sam_line(in_line, filtered=True):
    """
    Creates a named tuple from the first 11 fields of a SAM file line. Returns
    None if the line is empty, has fewer than 11 fields or filtered is True and
    the RNAME or CIGAR values are "*"."""

    if not in_line:
        return
    line = in_line.split("\t")
    if len(line) < 11:
        return
    out_line = sam_line._make(line[:11])
    if filtered:
        if out_line.RNAME == "*" or out_line.CIGAR == "*":
            return None

    return out_line


@memo
def cigar_regex(omit):
    """Returns (n, letter) tuples from a CIGAR string, omitting "omit" codes."""
    codes = "".join(set("SMIDH") - set(omit))
    re_cigar = re.compile(rf"((\d+)([{codes}]))")

    return re_cigar


@memo
def parse_cigar(cigar, omit="H"):
    return [(int(x[1]), x[2]) for x in cigar_regex(omit).findall(cigar)]


@memo
def cigar2arr(cigar):
    """
    Uses parse_cigar (and cigar_regex) to returns a numpy array of indexes into
    'SMIDH' from the CIGAR string. For example - "2S4M1I2D2M1S" would generate:
        [0 0 1 1 1 1 2 3 3 1 1 0], corresponding to
        [S S M M M M I D D M M S]
    """

    c_arr = np.fromiter(
        chain(it.repeat(cigar_codes[j], int(i)) for i, j in parse_cigar(cigar)),
        dtype=np.int8,
    )

    return c_arr


# -------------------------------------------------------------------------------
@memo
def qual2arr(qual):
    """Returns an array of integers corresponding to quality string elements"""
    q_arr = np.fromiter(map(ord, qual), dtype=np.int8) - 32

    return q_arr


################################################################################
#                               FILE UTILITIES                                 #
################################################################################


def filestem(fileobject):
    return os.path.basename(fileobject).rsplit(".", 1)[0]


def exists(filename):
    """Returns true only if a file exists AND it is of non-zero size."""

    return os.path.exists(filename) and os.path.getsize(filename) > 0


################################################################################
#                              TOOL WRAPPERS                                   #
################################################################################


def modify_and_run(cmd, *args, **kwargs):
    """Appends args/kwargs to subprocess commands and runs them"""

    out_file = kwargs.pop("out_file", None)

    if out_file and exists(out_file):
        return modify_and_run(cmd=["true"])

    final = kwargs.pop("final", [])
    func = kwargs.pop("func", sp.run)
    stdin = kwargs.pop("stdin", None)
    stdout = kwargs.pop("stdout", sp.PIPE)
    stderr = kwargs.pop("stderr", sp.PIPE)
    text = kwargs.pop("text", True)

    cmd.extend([*map(str, args)])
    for k, v in kwargs.items():
        if v is None:
            continue
        cmd.extend([f"-{k}", str(v)])

    if type(final) == str:
        cmd.append(final)
    else:
        cmd.extend(final)

    # if verbose:
    #     log.cmdline(f'\t{" ".join(cmd)}')
    # 
    pc = func(cmd, stdout=stdout, stderr=stderr, stdin=stdin, text=text)

    if pc.returncode and pc.stderr:
        log.error(f'\t{" ".join(cmd)}')
        log.error(pc.stderr)

    return pc


def log_subprocess(log, stream):
    stream = stream.strip().split("\n") if stream else []
    for line in stream:
        log.error(line)


################################################################################
"SHELL"


def shell(cmd, *args, **kwargs):
    """All-purpose shell command parser - unlike subprocess, it allows logging"""
    if not type(cmd) == list:
        cmd = [cmd]

    pc = modify_and_run(cmd, *args, **kwargs)
    return pc


def remove(*args, **kwargs):
    """filename(s) must be the last of the *args"""

    cmd = get_app(real_name="rm") + ["-rf"]

    pc = modify_and_run(cmd, *args, **kwargs)

    return pc



def awk(*args, **kwargs):
    """Pass the awk command within *args (allows pre-and post-command arguments)"""

    cmd = get_app()

    pc = modify_and_run(cmd, *args, **kwargs)
    return pc


################################################################################


def get_app(name=None, real_name=None):
    """Finds application command in $PATH unless specified in an environment variable"""

    caller = real_name if real_name else sys._getframe(1).f_code.co_name
    name = name if name else caller

    return [os.environ.get(caller, name)]


################################################################################
"MAPPING, DE NOVO ETC."

def samtools(*args, **kwargs):

    cmd = get_app()

    kwargs.setdefault("func", sp.Popen)

    pc = modify_and_run(cmd, *args, **kwargs)
    return pc

################################################################################


def jf_count(*args, fastq=None, **kwargs):
    """Jellyfish count"""

    if not fastq:
        raise GenomancerError()

    cmd = get_app(real_name="jellyfish")

    cmd.append("count")

    kwargs.setdefault("o", f"{filestem(fastq)}.jf")
    kwargs.setdefault("m", 31)
    kwargs.setdefault("t", 4)

    pc = modify_and_run(cmd, *args, out_file=kwargs["o"], final=fastq, **kwargs)
    return pc


def jf_dump(*args, jf=None, **kwargs):
    """Jellyfish dump"""
    if not jf:
        raise GenomancerError()

    cmd = get_app(real_name="jellyfish")

    cmd.append("dump")

    if "out_file" in kwargs:
        kwargs["o"] = kwargs.pop("out_file")

    pc = modify_and_run(cmd, *args, final=jf, **kwargs)
    return pc
