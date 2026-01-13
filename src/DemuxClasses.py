#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import subprocess
from io import TextIOWrapper
from shutil import which

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import edlib
from itertools import chain
import os

class Boundary:
    '''
    Simple class representing the zero-indexed start/end boundaries of two sequences.
    Alignments are made with BioPython's Align module.
    '''
    def __init__(self, seq_len: int, start_idx=np.nan, end_idx=np.nan, buffer: int = 0, editDistance:int = np.nan, valid: bool = False):
        self.seq_len = seq_len
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.buffer = buffer
        self.editDistance = editDistance 
        self.span = [np.nan, np.nan]
        self.valid = valid#

    def __str__(self):
        str=f'''
        Full_seq_len\tAln_start_idx\tAln_end_index\tBuffer\tScore\tValidity
        {self.seq_len}\t{self.start_idx}\t{self.end_idx}\t{self.buffer}\t{self.editDistance}\t{self.valid}
        '''
        return(str)#

    def get_Boundary_span(self):
        '''
        Calculates the span of the boundary, taking the buffer into account when present.
        If the buffer puts the span past the start of the seq (ie, if it's less than 0), defaults to 0.
        If the buffer puts it over the length of the seq, defaults to the seq length.#

        If the start or end idx are np.nan, then just leave this unset.
        '''
        if not self.valid:
            raise ValueError("Cannot get span of an invalid alignment.")
        span_start_idx = self.start_idx - self.buffer
        if (span_start_idx < 0):
            span_start_idx=np.int64(0)
        span_end_idx = self.end_idx + self.buffer
        if (span_end_idx > self.seq_len):
            span_end_idx = self.seq_len
        self.span=[span_start_idx, span_end_idx]
        
    def check_Boundary_validity(self):
        '''Boundaries are valid if both the start and end indices are valid.'''
        if np.isnan(self.start_idx) or np.isnan(self.end_idx):
            self.valid = False
        else:
            self.valid = True
    
    def set_Boundary(self, start_idx: int, end_idx: int, buffer: int, editDistance: int):
        '''
        Fills in a boundary and checks its validity.
        '''
        self.start_idx=start_idx
        self.end_idx=end_idx
        self.buffer=buffer
        self.editDistance=editDistance
        self.check_Boundary_validity()
        if self.valid:
            self.get_Boundary_span()#

class SimpleSeqRecord():
    '''
    A lighter-weight implementation of Biopython's SimpleSeqRecord object, designed to work with `FastqGeneralIterator`.
    It contains the title, seq, and qual strings of the original SimpleSeqRecord object, and otherwise acts as a drop-in replacement.
    '''
    def __init__(self, id: str, seq: str, qual: str):
        self.id = id
        self.seq = Seq(seq)
        self.qual = qual
    
    def __str__(self):
        str = f'''
        id\t{self.id}
        seq\t{self.seq}
        qual\t{self.qual}
        '''
        return(str)

    def make_SeqRecord(self):
        '''
        Converts a SimpleSeqRecord to a SeqRecord.
        '''
        seq_record = SeqRecord(self.seq, id=self.id)
        seq_record.letter_annotations["phred_quality"] = [ord(q) - 33 for q in self.qual]#
        return(seq_record)#

class ConstructElement:
    '''
    Individual elements of a larger DNA construct, implemented as a Seq object with metadata.
    In the context of ONTddRADparse, they correspond to the index, barcode, index_full and barcode_full elements.
    We include the 'length' of construct, (long, short, etc.), the alignment specificity (`aln_percent`), and the raw sequence.
    '''
    def __init__(self, seq: Seq, CE_type: str, aln_percent: float, buffer: int=0):
        self.seq = seq
        self.CE_type = CE_type
        self.aln_percent = aln_percent
        self.buffer = buffer#

    def __str__(self):
        str = f'''
        seq\tCE_ype\taln_percent\tbuffer
        {self.seq}\t{self.CE_type}\t{self.aln_percent}\t{self.buffer}
        '''
        return(str)#

class ConstructElementAlignment:
    '''
    Represents a ConstructElement aligned to a SimpleSeqRecord.
    Generates/houses boundary objects and stores the validity of their alignment.
    A ConstructElementAlignment is valid if XORing RBoundary and RBoundary is true.
    '''

    def __init__(self, SimpleSeqRecord: 'SimpleSeqRecord', ConstructElement: 'ConstructElement'):
        self.SimpleSeqRecord = SimpleSeqRecord
        self.ConstructElement = ConstructElement
        self.orientation = []
        seq_len=len(self.SimpleSeqRecord.seq)
        self.FBoundary = Boundary(seq_len) # you could probably set the boundaries right now
        self.RBoundary = Boundary(seq_len)
        self.editDistance = np.nan
        self.valid = False

    def __str__(self):
        str=f'''
        SimpleSeqRecord\tConstructElement_type\Orientation\taln_valid\tFBoundary_start_idx\tFBoundary_end_idx\tRBoundary_start_idx\tRBoundary_end_idx
        {self.SimpleSeqRecord.id}\t{self.ConstructElement.CE_type}\t{self.ConstructElement.orientation}\t{self.valid}\t{self.FBoundary.start_idx}\t{self.FBoundary.end_idx}\t{self.RBoundary.start_idx}\t{self.RBoundary.end_idx}
        '''
        return(str)

    def align_ConstructElement(self):
        '''
        Aligns the subset to the seq in the 5->3 ('F') and 3->5 ('R') direction and fills in boundaries.
        Also updates the "orientation" and "editDistance" tags of the ConstructElementAlignment.#

        '''
        self.set_ConstructElement_Boundary('F')
        self.set_ConstructElement_Boundary('R')
        
        # We allow both to be appended to account for short sequences possibly appearing more than once by chance
        if self.FBoundary.valid: # aligned in F
            self.orientation.append('F')
        if self.RBoundary.valid: # aligned in R
            self.orientation.append('R')

        # we don't check for editDistances with short CEs.
        # this is because short CEs should match exactly within long CEs
        # also because they can match multiple times. there's no reason to score them
        if (self.ConstructElement.CE_type == 'long'):
            if self.orientation == ['F']:
                self.editDistance = self.FBoundary.editDistance
            if self.orientation == ['R']:
                self.editDistance = self.RBoundary.editDistance

    def set_ConstructElement_Boundary(self, orientation: str):
        '''
        Align in the specified direction ('F' for forward, 'R' for reverse).
        Sets the appropriate boundary (FBoundary or RBoundary) based on orientation.
        '''
        if orientation == 'F':
            target_seq = self.SimpleSeqRecord.seq
            boundary = self.FBoundary
        elif orientation == 'R':
            target_seq = self.SimpleSeqRecord.seq.reverse_complement()
            boundary = self.RBoundary
        else:
            raise ValueError("Orientation must be either 'F' or 'R'.")
        
        # We will only use the alignment method if we have to. Exact matches can get string methods.
        if self.ConstructElement.CE_type == 'long':
            idxes, editDistance = align_target(target_seq, self.ConstructElement.seq, self.ConstructElement.aln_percent)
        elif self.ConstructElement.CE_type == 'short':
            subseq=self.ConstructElement.seq.lower()
            start_idx = target_seq.lower().find(subseq) # 'lower' handles any mismatches due to softmasking
            if start_idx != -1:
                idxes = [start_idx, start_idx+(len(subseq)-1)]
                editDistance = 0 # string-substrings are an exact search, so it must have an edit dist of 0
            else:
                idxes = [np.nan, np.nan]
                editDistance = len(self.ConstructElement.seq) # it still has to be 0 so we can appropriately add them together.
        else:
            raise ValueError("CE_type must be either 'short' or 'long'.")

        # set the boundary
        boundary.set_Boundary(idxes[0], idxes[1], self.ConstructElement.buffer, editDistance)#

    def check_ConstructElementAlignment_validity(self):
        '''
        Checks the orientation of each boundary, as well as the length of the construct to determine validity.
        Valid either when XORed or when a short element is aligned in both directions.
        '''
        # Base rule: valid when exactly one of aligned_in_orientation_F / aligned_in_orientation_R is True (logical XOR)
        orientation_F = bool(('F' in self.orientation) and ('R' not in self.orientation))
        orientation_R = bool(('F' not in self.orientation) and ('R' in self.orientation))
        xor_valid = bool(orientation_F ^ orientation_R)#

        # Additional rule: for short constructs, allow both boundaries aligned
        short_orientation_FR_valid = bool((self.ConstructElement.CE_type == 'short') and ('F' in self.orientation) and ('R' in self.orientation))#

        self.valid = xor_valid or short_orientation_FR_valid

    def check_ConstructElementAlignment_concatamer_validity(self):
        '''
        Checks for concatamers in binding in cis to each other.
        Concatamers resulting from the same seq binding in trans are already removed by checkConstructElementAlignment_validity.
        '''
        subseq = self.ConstructElement.seq
        max_edit_dist = len(subseq) - len(subseq)*self.ConstructElement.aln_percent

        if ('F' in self.orientation):
            seq = self.SimpleSeqRecord.seq
        elif ('R' in self.orientation):
            seq = self.SimpleSeqRecord.seq.reverse_complement()
        aln=edlib.align(query=subseq, target=seq, k=max_edit_dist, mode='HW', task='locations')
        
        # We discard all concatamers; there is no trimming.
        if (len(aln.get('locations')) > 1):
            self.valid = False
            return False
        self.valid=True
        return True

class ConstructElementAlignmentPair:
    '''
    Two complementary `ConstructElements` - one long, and one short. In the context of ONTddRADparser, they consist of:
        - `index` and `index_full`
        - `barcode` and `barcode_full`
    The short one is around 6-9 bp, and is allowed to exist more than once in the SimpleSeqRecord.
    The long one should only be found once - otherwise, it is a concatamer.
    '''
    def __init__(self, CEA_long: 'ConstructElementAlignment', CEA_short: 'ConstructElementAlignment'):
        if CEA_long.SimpleSeqRecord.id != CEA_short.SimpleSeqRecord.id:
            raise ValueError("ConstructElementAlignments must be based on the same sequence.")
        self.valid = False
        self.CEA_long = CEA_long
        self.CEA_short = CEA_short
        self.orientation = [] # can be ['F','R','invalid']
        self.editDistance = np.nansum([CEA_long.editDistance, CEA_short.editDistance])
        self.CEA_short_in_CEA_long = False # can be T or F
    
    def __str__(self):
        str = f'''
        overall_valid\torientation\tshort_CE_inside_long_CE
        {self.valid}\t{self.orientation}\t{self.CEA_short_in_CEA_long}
        '''
        return(str)

    def get_ConstructElementAlignmentPair_orientation(self):
        '''
        Determines whether the pair is `F`, `R`, or [], and updates `self.orientation`.
        '''
        if len(self.orientation) == 2:
            raise ValueError("You are trying to get the orientation of something that already has been oriented.")
        for i in ['F','R']:
            if (i in self.CEA_long.orientation):
                if (i not in self.CEA_short.orientation):
                    self.orientation = [] # this may not be necessary. But it'll return false in any checks
                    return False
                self.orientation.append(i)
                return True

        # It is valid to have more than one orientation in a CEA. However, a CEAP should never have more than one alignment
        if len(self.orientation) > 1:
            raise ValueError("Something has gone wrong here.")#

    def check_short_ConstructElementAlignment_in_long_ConstructElementAlignment(self):
        '''
        Performs position checking to ensure that CEA_short is inside CEA_long.
        '''
        if not self.orientation: # if you didn't have an orientation before, get one first.
            self.get_ConstructElementAlignmentPair_orientation()#

        # now we need to check that the short element is found inside of the long element
        if ('F' in self.orientation):
            short_span = self.CEA_short.FBoundary.span
            long_span = self.CEA_long.FBoundary.span
        elif ('R' in self.orientation): # now we need to check that the short element is found inside of the long element
            short_span=self.CEA_short.RBoundary.span
            long_span=self.CEA_long.RBoundary.span
        else:
            # failing here means you have an invalid alignment.
            self.CEA_short_in_CEA_long = False
            return False
            
        if ((short_span[0] < long_span[0]) or (short_span[1] > long_span[1])):
            self.CEA_short_in_CEA_long = False
            return False
        self.CEA_short_in_CEA_long = True
        return True

    def check_ConstructElementAlignmentPair_validity(self):
        '''
        Checks if the CEAP is valid, based on self.orientation and self.CEA_short_in_CEA_long, and updates the self.valid tag.
        '''
        if not self.orientation:
            #self.get_ConstructElementAlignmentPair_orientation()
            if self.CEA_short_in_CEA_long:
                self.valid = True
                return True
        self.valid = False
        return False

class DemuxConstruct:
    '''
    Represents the construct sequence data we add to our sequences to demux them.
    Consists of a sample ID and four ConstructElement objects: two 6-9bp ones which must be exact matches, and two longer ones where error is allowed.
    Canonically, should be read in from your demux file.
    '''
    def __init__(self, SampleID, index_full: 'ConstructElement', index:'ConstructElement', barcode_full:'ConstructElement', barcode: 'ConstructElement'):
        self.SampleID = SampleID
        self.index_full = index_full
        self.index = index
        self.barcode_full = barcode_full
        self.barcode = barcode#

    def __str__(self):
        str = f'''
        SimpleSeqRecord\t{self.SampleID}
        index_full\t{self.index_full}
        index\t\t{self.index}
        barcode_full\t{self.barcode_full}
        barcode\t{self.barcode}'''
        return(str)#

class DemuxConstructAlignment:
    '''
    Represents the alignments between a SecRecord object and a DemuxConstruct.
    Contains the SimpleSeqRecord and DemuxConstruct.
    Also contains two ConstructElementPairAlignments, (idx and barcode).
    Finally, contains a 'valid' and 'reason' tag.
    '''#

    def __init__(self, SimpleSeqRecord: 'SimpleSeqRecord', DemuxConstruct: 'DemuxConstruct'):
        self.valid_dict = {
            'all_CEAs_valid': False,
            'no_long_CEA_concatamers_valid': False,
            'CEAs_in_CEAP_same_orientation': False,
            'CEAs_short_inside_CEAs_long': False,
            'CEAPs_opposite_orientation': False,
        }

        self.valid = False     # we initialize validity as false until proven otherwise
        self.SimpleSeqRecord = SimpleSeqRecord
        self.DemuxConstruct = DemuxConstruct
        self.orientation = []
    
        # construct and validate CEAs
        CEA_index_full=ConstructElementAlignment(SimpleSeqRecord, DemuxConstruct.index_full)
        CEA_index=ConstructElementAlignment(SimpleSeqRecord, DemuxConstruct.index)
        CEA_barcode_full=ConstructElementAlignment(SimpleSeqRecord, DemuxConstruct.barcode_full)
        CEA_barcode=ConstructElementAlignment(SimpleSeqRecord, DemuxConstruct.barcode)#

        self.align_all_ConstructElements([CEA_index_full, CEA_index, CEA_barcode_full, CEA_barcode])#

        self.CEAP_index = ConstructElementAlignmentPair(CEA_long=CEA_index_full, CEA_short=CEA_index)
        self.CEAP_barcode = ConstructElementAlignmentPair(CEA_long=CEA_barcode_full, CEA_short=CEA_barcode)
        self.editDistance = self.CEAP_index.editDistance + self.CEAP_barcode.editDistance
        self.CEAP_index.get_ConstructElementAlignmentPair_orientation()
        self.CEAP_barcode.get_ConstructElementAlignmentPair_orientation()
        self.CEAP_index.check_short_ConstructElementAlignment_in_long_ConstructElementAlignment()
        self.CEAP_barcode.check_short_ConstructElementAlignment_in_long_ConstructElementAlignment()#

    def __str__(self):
        str = f'''
        SimpleSeqRecord\t{self.SimpleSeqRecord.id}
        CEAP_index\t{self.CEAP_index}
        CEAP_barcode\t{self.CEAP_barcode}
        '''
        return(str)#

    def align_all_ConstructElements(self, CEA_lst):
        for CEA in CEA_lst:
            CEA.align_ConstructElement()
            CEA.check_ConstructElementAlignment_validity()
            if not CEA.valid:
                self.valid_dict['all_CEAs_valid'] = False
                return False
        self.valid_dict['all_CEAs_valid'] = True
        return True

    def check_all_ConstructElementAlignments_validity(self):
        '''
        Wrapper around ConstructElementAlignment.check_ConstructElementAlignment_validity().
        Ports the logic into the main script so we can weed out SeqQbjects that fail the first line of validation.
        '''
        CEA_lst = [self.CEAP_index.CEA_long, self.CEAP_index.CEA_short, self.CEAP_barcode.CEA_long, self.CEAP_barcode.CEA_short]
        for CEA in CEA_lst:
            CEA.check_ConstructElementAlignment_validity()
            if not CEA.valid:
                self.valid_dict['all_CEAs_valid'] = False
                return False
        self.valid_dict['all_CEAs_valid'] = True
        return True

    def check_all_ConstructElementAlignments_concatamer_validity(self):
        '''
        Checks validated CEAs for concatamers.
        '''
        CEA_longs = [self.CEAP_index.CEA_long, self.CEAP_barcode.CEA_long]
        for CEA_long in CEA_longs:
            CEA_long.check_ConstructElementAlignment_concatamer_validity()
            if not CEA_long.valid:
                self.valid_dict['no_long_CEA_concatamers_valid'] = False
                return False
        self.valid_dict['no_long_CEA_concatamers_valid'] = True
        return True

    def check_all_ConstructElementAlignmentPairs_orientation_validity(self):
        '''
        Checks that both CEAPs have an orientation. (If there is no orientation, one couldn't be assigned and its false.
        '''
        CEAPs = [self.CEAP_index, self.CEAP_barcode]
        for CEAP in CEAPs:
            if not CEAP.orientation:
                self.valid_dict['CEAs_in_CEAP_same_orientation'] = False
                return False
        self.valid_dict['CEAs_in_CEAP_same_orientation'] = True
        return True

    def check_all_ConstructElementAlignmentPairs_short_in_long_validity(self):
        '''
        Checks that the indices of the short CEA are inside the boundaries of the long one.
        '''
        CEAPs = [self.CEAP_index, self.CEAP_barcode]
        for CEAP in CEAPs:
            if not CEAP.CEA_short_in_CEA_long:
                self.valid_dict['CEAs_short_inside_CEAs_long'] = False
                return False
        self.valid_dict['CEAs_short_inside_CEAs_long'] = True        
        return True

    def check_ConstructElementAlignments_opposite_orientation_validity(self):
        '''
        Checks if there are any violations in the interaction of the two CEAPs and updates 'self.valid' as necessary.
        Barcode and index must have alternating orientations.
        '''
        for i in [self.CEAP_index, self.CEAP_barcode]:
            if self.CEAP_index.orientation == self.CEAP_barcode.orientation:
                self.valid_dict['CEAPs_opposite_orientation'] = False
                return False
        self.orientation.append(self.CEAP_index.orientation)
        self.orientation.append(self.CEAP_barcode.orientation)
        self.orientation = list(chain(*self.orientation))
        self.valid_dict['CEAPs_opposite_orientation'] = True
        return True

    def check_DemuxConstructAlignment_validity(self):
        '''
        Checks to see if all items in 'valid_dict' are True, and sets the valid tag appropriately.
        '''
        filter_lst = [
        'check_all_ConstructElementAlignments_validity',
        'check_all_ConstructElementAlignments_concatamer_validity',
        'check_all_ConstructElementAlignmentPairs_orientation_validity',
        'check_all_ConstructElementAlignmentPairs_short_in_long_validity',
        'check_ConstructElementAlignments_opposite_orientation_validity',
        ]

        for filter in filter_lst:
            if not getattr(self, filter)():
                break

        if all(filter_value for filter_value in self.valid_dict.values()):
            self.valid = True
            return True

    def trim_ConstructElements_out_of_SimpleSeqRecord(self):
        '''
        Use the coordinates from the long CEAs to remove them from the SimpleSeqRecord.
        '''
        if self.orientation == ['F','R']:
            FSpan = self.CEAP_index.CEA_long.FBoundary.span
            RSpan = self.CEAP_barcode.CEA_long.RBoundary.span
        elif self.orientation == ['R','F']:
            FSpan = self.CEAP_barcode.CEA_long.FBoundary.span
            RSpan = self.CEAP_index.CEA_long.RBoundary.span
        else:
            raise ValueError("Orientation is invalid for DCA")

        # Flip RSpan to compensate for reverse complement: convert to forward coordinates and swap
        seq_len = np.int64(len(self.SimpleSeqRecord.seq))
        RSpan = [seq_len - RSpan[1], seq_len - RSpan[0]]

        # Slice around both spans and concatenate
        self.SimpleSeqRecord.seq = self.SimpleSeqRecord.seq[0:FSpan[0]] + \
                                    self.SimpleSeqRecord.seq[FSpan[1]:RSpan[0]] + \
                                    self.SimpleSeqRecord.seq[RSpan[1]:seq_len]
        self.SimpleSeqRecord.qual = self.SimpleSeqRecord.qual[0:FSpan[0]] + \
                                    self.SimpleSeqRecord.qual[FSpan[1]:RSpan[0]] + \
                                    self.SimpleSeqRecord.qual[RSpan[1]:seq_len]
        return True

class DemuxxedSample:
    '''
    Represents each individual sample with that has sequence data associated with it.
    Aggregates SimpleSeqRecords from DemuxConstructAlignments that share the same sample_id.
    Uses that information to create a FastqFile by calling that class's method.
    '''
    def __init__(self, SampleID, SimpleSeqRecord_lst: list = []):
        self.SampleID = SampleID
        self.SimpleSeqRecord_lst = SimpleSeqRecord_lst#

    def init_FastqFile_from_Demuxxed_Sample(self, outdir='.'):
        f = FastqFile(filename = self.SampleID, outdir=outdir, SimpleSeqRecord_lst = self.SimpleSeqRecord_lst)
        return(f)
    
class FastqFile:
    '''
    Represents a name and a set of associated sequences.
    Written to the working directory by default, optionally takes an output directory.
    
    This class might be redundant, and its methods should be transferred to DemuxxedSample.
    '''
    format = 'fastq'

    def __init__(self, filename, outdir='.', SimpleSeqRecord_lst = []):
        if (not filename.endswith('.fq.gz')):
            filename = filename + '.fq.gz'
        self.filename = filename
        if (not outdir.endswith('/')):
            outdir = outdir + '/'
        self.outdir = outdir
        self.filepath = f'{self.outdir}{self.filename}'
        self.SimpleSeqRecord_lst = SimpleSeqRecord_lst

    def write_FastqFile_to_outdir(self):
        '''
        Thin wrapper around SeqIO.write.
        Uses `subprocess` to compress files - looks for `pigz` by default; otherwise uses `gzip`.
        Also creates an outdir if it does not exist.
        '''
        os.makedirs(self.outdir, exist_ok=True)

        pigz_path = which("pigz")
        if pigz_path:
            with open(self.filepath, "wb") as outfile, \
                    subprocess.Popen([pigz_path, "-c"], stdin=subprocess.PIPE, stdout=outfile) as proc:
                with TextIOWrapper(proc.stdin, encoding="utf-8") as handle:
                    for i in self.SimpleSeqRecord_lst:
                        new_SeqRecord=i.make_SeqRecord()
                        SeqIO.write(new_SeqRecord, handle, self.format)
                    handle.flush()
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    raise RuntimeError("pigz failed while writing FASTQ output.")
        else:
            with gzip.open(self.filepath, "wt") as handle:
                SeqIO.write(self.SimpleSeqRecord_lst, handle, self.format)#

def align_target(seq: Seq, subseq: Seq, aln_percent: float = 1.0):
    '''
    Takes two Bio.Seq objects (seq/target and query/subseq), and aligns them.
    Returns a tuple of the start and end indices of 'subseq' to 'seq'.
    Returns a list of two indices which can be formatted to a Boundary object.
    '''
    max_edit_dist = len(subseq) - len(subseq)*aln_percent#

    aln = edlib.align(query=subseq, target=seq, mode='HW', k=max_edit_dist, task='locations')
    if (aln.get('editDistance')==-1):
        idxes=[np.nan, np.nan]
        editDistance=len(subseq)
    else:
        aln_boundaries=list(aln.get('locations')[0])
        idxes=[min(aln_boundaries), max(aln_boundaries)]
        editDistance=aln.get('editDistance')
    return [idxes, editDistance] # renaming it to keep it consistent with the CE_short return boundary methods#

#def main():#

#    seq_name='test_seq'
#    input_seq='TTTTCTGTCCTGTACTTCGTTCAGTTAGGTATTGTTCAAGCAGAAGACAGAATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATTCTTTGGTTTATTTAAGATTAAACGCATACATGTTTGCAGAATTGCTGCCACAATAAGGTACTTTATCTTTTATACACAGGGAGCATTTACATTTTATACATATGAAAATATACCTCTAAGGACTTTTTTTTTTGCAACAAAACACAGCAGTTACCGACGCCTGATTCCCAGCTGGGGTAAGTCAGCTTGCAGACACTGTAGGAGCTGTGATGGTTGTAGCAGCTGAGATCTAGATGTACAGTCATGTCGAGTTTCCTTAAATATATGACACAAATGTACTACTCTCACACTCAGAGTGGCAACTTAGCACAACTATCTGCCAGCGCATGAGCATCTCTCAGTCCCAGAGAAAACCTTTATCCCTTACCTACACAGGTAATCTTTGAAACACTGACAGCACAACAACTAAACAAAATCTTGTATCATAAAGGCATAGAATTTAGGTTCTTTTTGATGAGAAAATCATCAAATACAGCTCTGAGCCAAAGCCCAGTGATGTTCCAGCCCTTTCTGCTGCCACCACCAGGCATGTCCCATGAAGGCCTGCAGAAGCTAGATAGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAGCAATACGT'
#    index_full='CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT'
#    barcode_full='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA'
#    index='CGTGAT'
#    barcode='ACACCT'#

#    #seq, DC, DCA=create_sample_data(seq_name, input_seq, index_full, index, barcode_full, barcode)

#if __name__=="__main__":
#    main()