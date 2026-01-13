#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../speeddemONT_illumina') # FIXME - brittle hardcoded link to the directory above to resolve `src`. Fix it later
import argparse
from src import utils
from src import DemuxClasses as dc


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
    parser.add_argument('-na', '--no_adapters', type=bool, required=True)

    args = parser.parse_args()

    #### READ IN FILES
    utils.print_user_info(f'Reading in {read_type} sequences:')
    SimpleSeqRecord_lst1 = utils.parse_seqfile(args.read1) # uses FastqGeneralIterator to read big FAs cheaply
    SimpleSeqRecord_lst2 = utils.parse_seqfile(args.read2)
    
    # should this be a different format???
    print(f'{len(SimpleSeqRecord_lst1)} reads in READ1')
    print(f'{len(SimpleSeqRecord_lst2)} reads in READ2')


    error_code_not_implemented=100
    print("Congratulations! You have tried to demux Illumina data. The programmer hasn't finished implementing this yet.")
    print(f"We\'re going to return with error code {error_code_not_implemented} now.")

    # check to see if the user's demux groups all have the same index values.
    # if they do have the same, print a warning, but don't override

    return error_code_not_implemented

if __name__ == "__main__":
    sys.exit(main())
