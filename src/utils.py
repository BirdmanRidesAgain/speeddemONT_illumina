#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('../speeddemONT_illumina') # FIXME - brittle hardcoded link to the directory above to resolve `src`. Fix it later

from termcolor import cprint
from sys import stderr
from tqdm import tqdm
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from mimetypes import guess_type
from functools import partial
import pandas as pd
from collections import defaultdict
import random
import matplotlib.pyplot as plt
import seaborn as sns


from src import DemuxClasses as dc

BRIGHT_CYAN = '\033[96m'
RED = '\033[91m'
RESET = '\033[0m' # called to return to standard terminal text color


warning_template = """
    {RED}WARNING:{RESET} --{name} set by user to non-default value
    \tDefault: --{name}={default}
    \tCurrent: --{name}={current}
    --{name} does not operate on {read_type} data.
    {RED}This does not affect your results{RESET}, but the argument is useless.
    It is recommended that you remove it.
    """

def print_user_info(message: str, bold: bool=True):
    '''
    Formats user-important strings as boldface magenta, and writes them to stderr.
    Writing to stderr keeps them from being picked up by redirection.
    '''
    if bold:
        cprint(message, "light_magenta", attrs=["bold"], file=stderr)
    else:
        cprint(message, "light_magenta", file=stderr)

def print_logo():
    '''
    Prints an obnoxious ASCII-logo in rainbow lettering to stdout.
    '''
    cprint('                                 _________________               ________________________              ', "light_magenta", attrs=["bold"], file=stderr)
    cprint('        _______________________________  /_____  /___________ _____// __ \__// | / //__  __/         ', "light_magenta", attrs=["bold"], file=stderr)
    cprint('       ___  ___/__  __ \  _ \  _ \  __  /_  __  /_  _ \_  __ \__ \// / / /_//  |/ /__// /          ', "light_magenta", attrs=["bold"], file=stderr)
    cprint('     ____(__  )__  /_/ /  __/  __/ /_/ / / /_/ / /  __/  / / / / // /_/ /_// /|  / _// /         ', "light_magenta", attrs=["bold"], file=stderr)
    cprint('    ____/____/ _  ____/\___/\___/\__,_/  \__,_/  \___//_/ /_/ /_/\_____/ //_/ |_/  //_/        ', "light_magenta", attrs=["bold"], file=stderr)
    cprint('               /_/                                                                           ', "light_magenta", attrs=["bold"], file=stderr)

def print_args(args):
    '''
    Prints user-defined arguments.
    '''
    BOLD='\033[1m'
    RESET='\033[0m'
    print_user_info('User-defined arguments:')
    for key, value in vars(args).items():
        print(f"{BOLD}\t{key}{RESET}: {value}")

def parse_seqfile(filepath: str):
    '''
    Takes in a path to a sequence file (compressed or uncompressed).
    Returns a list of all records in the file as SimpleSeqRecord objects.
    '''
    encoding = guess_type(filepath)[1]
    _open = partial(gzip.open, mode = 'rt') if encoding == 'gzip' else open
    seqfile_lst = []
    with _open(filepath) as f:
        for id, seq, qual in FastqGeneralIterator(f):
            seqfile_lst.append(dc.SimpleSeqRecord(id, seq, qual))
    return seqfile_lst

def parse_demux_file_without_adapters(filepath: str):
    '''
    Parser and checker for a file that designed to demux WITHOUT adapters (ie, index_full, index, barcode_full, barcode).
    This is only intended for use in demuxxing Illumina data.

    FIXME:
    We need to enforce the following standards:
        - 2-column TSV file expected to contain a specific header: [sample_id,barcode] 
        - All cols (except sample_id) must be lowercase DNA nucleotides
        - `barcode` must be 6-9 bp long
        - Barcodes cannot contain dups.
    
    A minimal example looks like this:
    sample_id	barcode
    R10N00251	ACACCT
    R20N00088	ACAGCA

    Raises a ValueError if the barcode is not 6 or 9 nucleotides long or if barcodes are repeated.
    Returns a data frame.
    '''
    df = pd.read_csv(filepath, sep='\t', header='infer')
    header=['sample_id','barcode']

    if list(df.columns) != header:
        raise ValueError(
            "Invalid demux header."
            f"Expected columns: {header}. "
            f"Found: {list(df.columns)}"
        )
    
    num_cols=2
    if df.shape[1] != num_cols:
        raise ValueError(f"Expected {num_cols} columns, found {df.shape[1]}")
    
    # Check columns 3 and 5 (index 2 and 4) for string length 6 or 9
    for col in [df.columns[2], df.columns[4]]:
        invalid_rows = df[~df[col].apply(lambda x: isinstance(x, str) and len(x) in (6, 9))]
        if not invalid_rows.empty:
            raise ValueError(
            f"Barcodes and indices '{col}' must be 6 or 9 nucleotides long. "
            f"Invalid rows:\n{invalid_rows[[col]].to_string(index=True)}"
            )

    #return df

def parse_demux_file_with_adapters(filepath: str):
    '''
    Parser and checker for a file that designed to demux with adapters (ie, index_full, index, barcode_full, barcode).
    This is the default, and 'safest' parser- it's the ONLY demux parser for ONT data, and is used for Illumina by default.
    The only way to NOT use it is to set `no_adapters` on an illumina run.

    FIXME - please actually enforce these
    We enforce the following standards:
        - 5-column TSV file expected to contain a specific header: [sample_id,index_full,index,barcode_full,barcode]
        - All cols (except sample_id) must be lowercase DNA nucleotides
        - `index` and `barcode` must be 6-9 bp long
    
    A minimal example looks like this:
    sample_id	index_full	index	barcode_full	barcode
    R10N00251	CAA...AATT	CGTGAT	AAT...GCA	ACACCT
    R20N00088	CAA...ATTT	CGTGAT	AAT...GCA	ACAGCA

    Raises a ValueError if either the barcode or index are not 6 or 9 nucleotides long.
    Returns a data frame.
    '''
    df = pd.read_csv(filepath, sep='\t', header='infer')
    header=['sample_id','index_full','index','barcode_full','barcode']

    if list(df.columns) != header:
        raise ValueError(
            "Invalid demux header."
            f"Expected columns: {header}. "
            f"Found: {list(df.columns)}"
        )
    num_cols=5
    if df.shape[1] != num_cols:
        raise ValueError(f"Expected {num_cols} columns, found {df.shape[1]}")
    
    # Check columns 3 and 5 (index 2 and 4) for string length 6 or 9
    for col in [df.columns[2], df.columns[4]]:
        invalid_rows = df[~df[col].apply(lambda x: isinstance(x, str) and len(x) in (6, 9))]
        if not invalid_rows.empty:
            raise ValueError(
            f"Barcodes and indices '{col}' must be 6 or 9 nucleotides long. "
            f"Invalid rows:\n{invalid_rows[[col]].to_string(index=True)}"
            )

    # FIXME - please implement checking for lowercase DNA strings

    return df

def convert_demux_df_to_SampleID_dict(df: pd.DataFrame, fuzzy_aln_percent: float, num_short_mismatch: float, buffer: int):
    '''
    Ingests a data frame and converts it to a dictionary object organized by sample_id.
    The values are a list of DCs, allowing multiple DCs per sample_id.
    '''
    # Build mapping from sample_id to associated barcode/index values
    sample_id_dict = defaultdict(list)
    # Group the DataFrame by sample_id, and for each group, iterate over all rows (values)
    for sample_id, group_df in df.groupby('sample_id'):
        DCs_in_SampleID_lst = [] # these are both initialized as empty
        demuxxed_seq_ids = []
        for _, row in group_df.iterrows():
            index_full_CE = dc.ConstructElement(row['index_full'], 'long', fuzzy_aln_percent, buffer)
            index_CE = dc.ConstructElement(row['index'], 'short', num_short_mismatch)
            barcode_full_CE = dc.ConstructElement(row['barcode_full'], 'long', fuzzy_aln_percent, buffer)
            barcode_CE = dc.ConstructElement(row['barcode'], 'short', num_short_mismatch)
            DC = dc.DemuxConstruct(sample_id, index_full_CE, index_CE, barcode_full_CE, barcode_CE)
            DCs_in_SampleID_lst.append(DC) 
        # initialize the seqs demuxxed as an empty list
        sample_id_dict[row['sample_id']] = [DCs_in_SampleID_lst, demuxxed_seq_ids]
    return sample_id_dict

def make_SampleID_dict(filepath: str, adapter_aln_percent: float, num_short_mismatch: int, buffer: int, no_adapters: bool=False):
    '''
    Parses our demux file. and `convert_demux_df_to_DemuxConstruct_lst1
    '''
    if no_adapters:
        demux_df = parse_demux_file_without_adapters(filepath)

    else:
        demux_df = parse_demux_file_with_adapters(filepath) # this is the safer, slower, and ONT-only option
    SampleID_dict = convert_demux_df_to_SampleID_dict(demux_df, adapter_aln_percent, num_short_mismatch, buffer)
    return SampleID_dict

def chunk_input_lst(lst: list, n_rds: int = 10000):
    '''
    Splits the incoming input_lst into groups of n iterations.
    We define 'n' arbitrarily as 10000 - I am unsure what the optimal value is.

    We are doing this to prevent the pooled processes from stalling.
    '''
    input_list_of_lists = [lst[i:i + n_rds] for i in range(0, len(lst), n_rds)]
    return input_list_of_lists

def optimize_SampleID_dict_order(SimpleSeqRecord_lst: list, n_samples: int, id_dict: dict):
    '''
    Gets the probability distribution of sample_ids in the `DC_dict` by classifying `n_samples`.
    Returns an ordered dict.
    '''
    random.seed(a=n_samples)
    for SimpleSeqRecord in tqdm(random.sample(SimpleSeqRecord_lst, n_samples)):
        for sample_id, sample_id_info in id_dict.items():
            DC_lst=sample_id_info[0]
            for DC in DC_lst:
                DCA=dc.DemuxConstructAlignment(SimpleSeqRecord, DC)
                DCA.check_DemuxConstructAlignment_validity()
                if DCA.valid:
                    id_dict[sample_id][1].append(DCA.SimpleSeqRecord)
    # code to sort the DC_dict
    # Sort based on reverse of Values
    sample_id_dict_sorted = {k: v for k, v in sorted(id_dict.items(), key=lambda item: len(item[1][1]), reverse=True)}    
    
    # Overwrite value[1] with blank list for every key
    # This prevents us from running the same samples twice, and means the dict is only reordered
    for key in sample_id_dict_sorted:
        sample_id_dict_sorted[key][1] = []
    return sample_id_dict_sorted

def write_fastq(fq):
    '''
    Helper function to parallelize writing fastq files to the drive by providing a mappable function.
    '''
    fq.write_FastqFile_to_outdir()


def make_DCA(input_lst: list):
    '''
    A wrapper around the default constructor for `DemuxConstructAlignment`.
    Takes a list of tuples as input, making it more amenable to multiprocessing.
    '''
    # renaming items so humans can interpret this
    SimpleSeqRecord = input_lst[0]
    SampleID_dict = input_lst[1]
    max_valid_comparisons = input_lst[2]

    valid_DCA_lst = [] # when alignments are fuzzy, you can get many DCAs
    for SampleID_info in SampleID_dict.values():
        DC_lst = SampleID_info[0]
        for DC in DC_lst:
            if len(valid_DCA_lst) >= max_valid_comparisons:
                break
            DCA=dc.DemuxConstructAlignment(SimpleSeqRecord, DC)
            DCA.check_DemuxConstructAlignment_validity()
            if DCA.valid:
                valid_DCA_lst.append(DCA)

    if (valid_DCA_lst):
        lowest_DCA = min(valid_DCA_lst, key=lambda dca: dca.editDistance) # find the DCA with the lowest edit distance
        lowest_DCA.trim_ConstructElements_out_of_SimpleSeqRecord()
        return(lowest_DCA)
    else:
        return(DCA) # you have to return the invalid ones to see why they failed

def split_DCA_lst(SampleID_dict: dict, DCA_lst: list):
    '''
    Takes in the sample_id_dict and divides it into successes (which can be assigned a sample_id) and failures (which cannot).
    
    Returns a 'valid' and a 'invalid' dict.
    Which is weird, but...
    '''
    SampleID_dict['not_assigned'] = [[], []] # We add this key in to account for failed sequences.
    invalid_dict = {} # empty dict; will contain all seqid of failed reads, plus their filter values. we want to know why they failed
    
    for DCA in tqdm(DCA_lst):
        if DCA.valid:
            SampleID_dict[DCA.DemuxConstruct.SampleID][1].append(DCA.SimpleSeqRecord)
        else:
            # this records all the info we have for the failed seqs.
            SampleID_dict['not_assigned'][1].append(DCA.SimpleSeqRecord)
            invalid_dict[DCA.SimpleSeqRecord.id] = DCA.valid_dict # this is a separate dict b/c it corresponds to plot 2, which has a different schema

    return [SampleID_dict, invalid_dict]

def make_df_from_SampleID_dict(SampleID_dict: dict):
    '''
    Converts the SampleID_dict into a data frame for plotting+printing.
    '''
    data = []
    for SampleID, (dc_list, seq_list) in SampleID_dict.items():
        data.append({
            'sample_id': SampleID,
            'count': len(seq_list),
            'num DemuxConstructs': len(dc_list)
        })
    df = pd.DataFrame(data)
    print(df)
    return(df)

def print_demultiplexing_summary(df: pd.DataFrame):
    '''
    Calculates and prints summary statistics for demultiplexing results.
    
    Takes a DataFrame with columns 'sample_id' and 'count' (as returned by 
    plot_number_of_SimpleSeqRecords_per_sample_id) and prints:
    - Number of successfully categorized samples (sample_id != 'not_assigned')
    - Number of successfully categorized reads
    - Total number of reads processed
    
    All values are reported as both absolute numbers and percentages.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing demultiplexing results with columns:
        - 'sample_id': str, sample identifiers (may include 'not_assigned' for uncategorized reads)
        - 'count': int, number of reads per sample_id
    
    Returns
    -------
    None
        Prints summary statistics to stdout.
    
    '''
    total_rds = df['count'].sum()
    binned_rds = df[df['sample_id'] != 'not_assigned']['count'].sum()
    num_samples = len(df[df['sample_id'] != 'not_assigned'])
    
    # Print summary statistics
    print(f"{binned_rds:,}/{total_rds:,} ({100 * binned_rds / total_rds:.2f}%) in {num_samples} sample IDs successfully demuxxed")
    print('\n')
    print(df.to_string(index=False))

def plot_number_of_SimpleSeqRecords_per_SampleID(df: pd.DataFrame):
    '''
    Plots the dataframe from sample_id_dict using seaborn.
    '''
    fig = plt.figure(figsize=(14,8), layout='constrained')
    ax = sns.barplot(data=df, y='sample_id', x='count', hue='num DemuxConstructs')
    for container in ax.containers:
        ax.bar_label(container, fmt='%d', padding=3, fontsize=12)
    ax.grid(True, axis='x', linestyle='--', alpha=0.5)

    # add labels
    label_font = {'fontsize':14,'fontweight':'bold'}
    ax.set_xlabel('Number of reads',fontdict=label_font)
    ax.set_ylabel('Sample IDs (bins)',fontdict=label_font)
    ax.tick_params(axis='y', labelrotation=0)
    ax.set_title('Successfully demultiplexed reads per sample',loc='left',fontsize=16,fontweight='bold',)
    fig.add_axes(ax)
    return fig

def plot_reasons_for_SimpleSeqRecord_invalidity(invalid_dict: dict):
    '''
    Takes in the invalid_dict to plot where SimpleSeqRecords are primarily lost and record sumstats.
    Counts the number of sequences that failed at each filter step.
    '''
    # Initialize counters for each filter step
    num_reads=len(invalid_dict)
    filter_counts = {
        'all_CEAs_valid': 0,
        'no_long_CEA_concatamers_valid': 0,
        'CEAs_in_CEAP_same_orientation': 0,
        'CEAs_short_inside_CEAs_long': 0,
        'CEAPs_opposite_orientation': 0,
    }
    
    # Count failures at each filter step
    for seq_id, valid_dict in invalid_dict.items():
        for filter_step, is_valid in valid_dict.items():
            if filter_step in filter_counts and not is_valid:
                filter_counts[filter_step] += 1
    
    # Convert to DataFrame for plotting
    data = []
    data.append({'filter_step': 'total num failed reads', 'count': num_reads})
    for filter_step, count in filter_counts.items():
        data.append({
            'filter_step': filter_step,
            'count': (num_reads - count)
        })
    df = pd.DataFrame(data)
    
    fig = plt.figure(figsize=(14,8), layout='constrained')
    ax = sns.barplot(
        data=df,
        x='filter_step', y='count'
    )
    for container in ax.containers:
        ax.bar_label(container, fmt='%d', padding=3, fontsize=12)
    ax.grid(True, axis='y', linestyle='--', alpha=0.5)

    # add labels
    label_font = {'fontsize':14,'fontweight':'bold'}
    ax.set_xlabel('Filter step',fontdict=label_font)
    ax.set_ylabel('Number of reads lost',fontdict=label_font)
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_title(
        'Number of SimpleSeqRecords after each filter',
        loc='left',
        fontsize=16,
        fontweight='bold',
    )

    fig.add_axes(ax)
    return [fig, df]
