#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from warnings import warn
from pathlib import Path
import multiprocessing as multi
import subprocess
import sys

from src import utils
from src import classes


def main():
    ### DEFINE AND CHECK ARGS
    default_prefix=Path(__file__).name.split(sep='.')[0]
    default_fuzzy_aln_percent=0.9
    default_buffer=0

    parser = ArgumentParser(description="Takes a set of ONT (or illumina) ddRAD reads and a set of barcodes and demultiplexes them.")
    parser.add_argument('-r1', '--read1', help='Path to a fastq file. If -r2 unset, interpreted as a single ONT readfile. Interpreted as READ1 of paired illumina read files otherwise. Required.', type=str, required=True)
    parser.add_argument('-r2', '--read2', help='Path to READ2 of paired illumina read files. Should only be set if processing illumina reads..', type=str)
    parser.add_argument('-d', '--demux', help='Path to the demux file. Required.', type=str, required=True)
    parser.add_argument('-b', '--buffer', help='The integer number of base pairs outside of the long element the short element can align to. Defaults to 0 (internal matching only).', type=int, default=default_buffer)
    parser.add_argument('-t', '--threads', help='The number of threads (actually CPU cores) this script can use to parallelize processes. Defaults to half the capacity of the host machine.', type=int, default=multi.cpu_count() // 2)    
    parser.add_argument('-p', '--prefix', help='Output prefix. Defaults to \'speeddemONT_out\' if not set.', default=f'{default_prefix}_out', type=str)
    parser.add_argument('-fa', '--fuzzy_aln_percent', help='The minimum percent identity needed to fuzzy-match a full index to a sequence. Not compatible with -rt=illumina.', default=default_fuzzy_aln_percent, type=float)
    parser.add_argument('-ns', '--num_short_mismatch', help='The number of mismatches allowed in an illumina barcode region. Defaults to 0, allowed options are [0,1,2].', default=0, type=int, choices=[0,1,2])
    parser.add_argument('-mc', '--max_comparisons', help='The maximum number of \'valid\' alignments the program can make in \'fuzzy\' alignments before picking the best match. Defaults to an arbitrarily high number otherwise. Lowering the value increases speed, but may produce misclassifications.', type=int)
    args = parser.parse_args()

    ### COLORS!
    BRIGHT_CYAN = '\033[96m'
    RED = '\033[91m'
    RESET = '\033[0m' # called to return to standard terminal text color


    ### Check inputs
    # Ensure that our alignment percent arguments are a value between 0 and 1
    if not (min(0, 1) < args.fuzzy_aln_percent <= max(0, 1)):
        error_message=f'''
        {RED}ERROR:{RESET} Alignment arguments (--fuzzy_aln_percent) must be percentages: eg, -fa .9
        \t--fuzzy_aln_percent={args.fuzzy_aln_percent}
        '''
        raise ValueError(error_message)

    ## if the user set args.max_comparisons, use their value
    # otherwise, set to 1 if 'exact', and an arbitrarily high number otherwise
    if (args.max_comparisons is None):
        if (args.fuzzy_aln_percent == float(1)):
            args.max_comparisons = 1
        else:
            args.max_comparisons = 80085 # this is far more than the number of barcodes in any realistic demux scenario

    # Warn the user if they've set inappropriate/inapplicable params for an Illumina run
    # buffer and fuzzy alignment percent don't do anything
    if (args.read2):
        warning_template = """
    {RED}WARNING:{RESET} --{name} set by user to non-default value
    \tDefault: --{name}={default}
    \tCurrent: --{name}={current}
    --{name} does not operate on illumina data.
    {RED}This does not affect your results{RESET}, but the argument is useless.
    It is recommended that you remove it.
    """
        warnings = []
        if (args.fuzzy_aln_percent != default_fuzzy_aln_percent):
            warnings.append(
                warning_template.format(
                    name='fuzzy_aln_percent', default=default_fuzzy_aln_percent, current=args.fuzzy_aln_percent, RED=RED,RESET=RESET
                )
            )
        if (args.buffer != default_buffer):
            warnings.append(
                warning_template.format(
                    name='buffer', default=default_buffer, current=args.buffer, RED=RED,RESET=RESET
                )
            )
        if warnings:
            warn("".join(warnings))


    ### ONCE ALL ARGS ARE VALID, PRINT USER INFO AND CALL A WORKER SCRIPT
    utils.print_args(args)
    
    # FIXME - provide a list of args to go into ONT, and one to go into illumina.
    # Use that list to set warnings and provide a set of args for str() to work on to DRY code
    demux_dir = Path(__file__).parent / "demux"
    if (not args.read2):
        script_path = demux_dir / "demux_ONT.py"
        result = subprocess.run(
            [sys.executable, str(script_path), 
                "--read1", str(args.read1),
                "--demux", str(args.demux),
                "--buffer", str(args.buffer),
                "--threads", str(args.threads),
                "--prefix", str(args.prefix),
                "--fuzzy_aln_percent", str(args.fuzzy_aln_percent),
                "--num_short_mismatch", str(args.num_short_mismatch),
                "--max_comparisons", str(args.max_comparisons),
                ],
            text=True
        )
    else:
        script_path = demux_dir / "demux_illumina.py"
        result = subprocess.run(
            [sys.executable, str(script_path), 
                "--read1", str(args.read1),
                "--read2", str(args.read2),
                "--demux", str(args.demux),
                "--threads", str(args.threads),
                "--prefix", str(args.prefix),
                "--num_short_mismatch", str(args.num_short_mismatch),
                ],
            text=True
        )

    if result.stdout:
        print(result.stdout) 

if __name__ == "__main__":
    sys.exit(main())