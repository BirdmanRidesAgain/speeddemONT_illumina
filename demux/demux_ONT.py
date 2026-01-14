#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../speeddemONT_illumina') # FIXME - brittle hardcoded link to the directory above to resolve `src`. Fix it later
import argparse
from src import utils
from src import DemuxClasses as dc
from itertools import product
import multiprocessing as multi
from tqdm import tqdm

def main() -> int:
    # We're not repeating help definitions or defaults here.
    read_type="ONT"
    parser = argparse.ArgumentParser(description=f"Demux {read_type} data.")
    parser.add_argument('-r1', '--read1', type=str, required=True)
    parser.add_argument('-d', '--demux', type=str, required=True)
    parser.add_argument('-b', '--buffer', type=int, required=True)
    parser.add_argument('-t', '--threads', type=int, required=True)    
    parser.add_argument('-p', '--prefix', type=str, required=True)
    parser.add_argument('-fa', '--fuzzy_aln_percent', type=float, required=True)
    parser.add_argument('-ns', '--num_short_mismatch', type=int, required=True)
    parser.add_argument('-mc', '--max_comparisons', type=int, required=True)
    args = parser.parse_args()

    #### READ IN FILES
    utils.print_user_info(f'Reading in {read_type} sequences:')
    SimpleSeqRecord_lst = utils.parse_seqfile(args.read1) # uses FastqGeneralIterator to read big FAs cheaply

    SampleID_dict = utils.make_SampleID_dict(args.demux, args.fuzzy_aln_percent, args.num_short_mismatch, args.buffer)

    utils.print_user_info('Making alignments:')
    DCA_lst=[]
    chunk_size = 10000 # this is an arbitrary number
    if len(SimpleSeqRecord_lst) > chunk_size:
        utils.print_user_info(f'\tLarge input. Running a burnin of {int(chunk_size/10)} replicates to optimize alignment order', bold=False)
        SampleID_dict=utils.optimize_SampleID_dict_order(SimpleSeqRecord_lst, int(chunk_size/10), SampleID_dict)

    # if the user set args.max_comparisons, use their value
    # otherwise, set to 1 if 'exact', and 'number of DCs' if fuzzy
    if (args.max_comparisons is None):
        if (args.fuzzy_aln_percent == float(1)):
            args.max_comparisons = 1
        else:
            args.max_comparisons = sum(len(SampleID_dict[key][0]) for key in SampleID_dict)

    input_lst = list(product(SimpleSeqRecord_lst, [SampleID_dict], [args.max_comparisons]))
    DCA_lst_valid = []

    utils.print_user_info('\tBeginning main loop:', bold=False)
    input_lst_of_lsts = utils.chunk_input_lst(input_lst, chunk_size)
    pool = multi.Pool(processes = args.threads)
    for tranche in tqdm(input_lst_of_lsts):
        DCA_sublst = pool.map(utils.make_DCA, tranche)
        DCA_lst.extend(DCA_sublst)
    pool.close()
    pool.join()

    utils.print_user_info('Checking alignment validity:')
    sample_id_dict, invalid_dict = utils.split_DCA_lst(SampleID_dict, DCA_lst)

    # Begin writing plots and FQs to outdir
    # Create one DemuxxedSample for each unique sample_id

    utils.print_user_info(f'Writing fastq files to {utils.BRIGHT_CYAN}{args.prefix}/reads{utils.RESET}:')
    fq_lst = []
    for SampleID, SampleID_info in SampleID_dict.items():
        DS = dc.DemuxxedSample(SampleID, SampleID_info[1])
        fq_lst.append(DS.init_FastqFile_from_Demuxxed_Sample(outdir=f'{args.prefix}/reads'))

    num_written_in_parallel=args.threads//10 # set this to use the computer's memory, not threads
    print(num_written_in_parallel)
    if num_written_in_parallel < 1:
        num_written_in_parallel=1
    pool = multi.Pool(processes = num_written_in_parallel) # FIXME - please improve the number of threads that can be run by writing without BioPython
    pool.map(utils.write_fastq, fq_lst)
    pool.close()
    pool.join()    

    # Plot results and print percent
    utils.print_user_info(f'\nWriting summary statistics to {utils.BRIGHT_CYAN}{args.prefix}{utils.RESET}:')
    df_success = utils.make_df_from_SampleID_dict(SampleID_dict)
    plot_success = utils.plot_number_of_SimpleSeqRecords_per_SampleID(df_success)
    plot_success.savefig(f'{args.prefix}/{args.prefix}_demult_success.png', dpi=300)
    df_success.to_csv(f'{args.prefix}/{args.prefix}_seqs_per_SampleID_stats.tsv', sep = "\t", index=False)
    utils.print_demultiplexing_summary(df_success)

    #plot2, df2 = plot_reasons_for_SimpleSeqRecord_invalidity(invalid_dict)
    #plot2.savefig(f'{args.prefix}/{args.prefix}_demult_failure.png', dpi=300)
    #df2.to_csv(f'{args.prefix}/{args.prefix}_failed_seqs_stats.tsv', sep = "\t", index=False)


    return 0

if __name__ == "__main__":
    sys.exit(main())