#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../speeddemONT_illumina') # FIXME - brittle hardcoded link to the directory above to resolve `src`. Fix it later
import argparse
from src import utils
from src import classes


def main() -> int:
    # We're not repeating help definitions or defaults here.
    read_type="illumina"
    parser = argparse.ArgumentParser(description=f"Demux {read_type} data.")
    parser.add_argument('-r1', '--read1', type=str, required=True)
    parser.add_argument('-r2', '--read2', type=str, required=True)
    parser.add_argument('-d', '--demux', type=str, required=True)
    parser.add_argument('-t', '--threads', type=int, required=True)    
    parser.add_argument('-p', '--prefix', type=str, required=True)
    parser.add_argument('-ns', '--num_short_mismatch', type=int, required=True)
    args = parser.parse_args()

    #### READ IN FILES
    utils.print_user_info(f'Reading in {read_type} sequences:')
    SimpleSeqRecord_lst1 = utils.parse_seqfile(args.read1) # uses FastqGeneralIterator to read big FAs cheaply
    SimpleSeqRecord_lst2 = utils.parse_seqfile(args.read2)
    print(len(SimpleSeqRecord_lst1))

    return 0

if __name__ == "__main__":
    sys.exit(main())
