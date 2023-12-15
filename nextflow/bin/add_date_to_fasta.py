#!/usr/bin/env python

from Bio import SeqIO
import warnings
import os
import sys

def add_date_to_fasta(original_file, updated_file=None, delim='@', keyword='naive', keyword_time=0, other_seq_time=20):
    """ 
    Add sequence time to fasta headers.
    It allows to add a specific time to one sequence, and set another time to the rest of sequences.
    
    CAUTION: If no sequence match keyword exactly, then 'other_seq_time' is added to all sequences. It will print a warning. 
    Parameters
    ----------
    original_file : str
        The fasta file path to be modified
    updated_file : str
        The modified fasta file will be written into. 
        Default is None, which will create a new file with '_with_time' appended to the original file name.
    delim : str
        A delimiter that seperates the orginal record and the added time
    keyword : str
        A particular sequence name to be set at a give time
    keyword_time : float
        Time of the keyword sequence
    other_seq_time : float
        Time of the rest of the sequences
        
    """
    keyword_found = False

    if updated_file is None:
        path = os.path.dirname(original_file)
        input_fn = os.path.basename(original_file)
        stem = os.path.splitext(input_fn)[0] # Note this is a hack. fn with multiple dots will not work.
        updated_file = os.path.join(path, stem + '_with_time.fasta')
    
    with open(original_file) as original, open(updated_file, 'w') as updated:
        records = SeqIO.parse(original_file, 'fasta')
        for record in records:
            if record.id == keyword:
                keyword_found = True
                record.id = record.id + delim + str(keyword_time)
                record.description = '' # Need this line; strips the orginal header
            else:
                record.id = record.id + delim + str(other_seq_time)
                record.description = ''
            SeqIO.write(record, updated, 'fasta-2line')
    if keyword_found == False:
        warnings.warn("No keyword found. All sequences are assigned with the same time.")
    return

if __name__ == '__main__':
    add_date_to_fasta(sys.argv[1], keyword='naive', keyword_time=0, other_seq_time=20)
