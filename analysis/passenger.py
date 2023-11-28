"""
Common functions for analyzing passenger data.
"""

import glob
import os
import re
import gzip

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Align import PairwiseAligner

def umi_histogram(fastq_gz_path, max_count):
    conscount_values = []
    
    # Regular expression to find CONSCOUNT field and extract the number
    pattern = re.compile(r"CONSCOUNT=(\d+)")

    # Parse the fastq file
    with gzip.open(fastq_gz_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            match = pattern.search(record.description)
            if match:
                conscount = int(match.group(1))
                # All values greater than max_count are set to max_count
                conscount = min(conscount, max_count)
                conscount_values.append(conscount)
    
    # Plot the histogram
    # bins should be one more than max_count to include the last bin
    plt.hist(conscount_values, bins=max_count+1, range=(0, max_count+1), edgecolor='black')
    
    fastq_basename = fastq_gz_path.split('/')[-1]
    plt.title(f'Distribution of CONSCOUNT for {fastq_basename}')
    plt.xlabel('CONSCOUNT')
    plt.ylabel('Frequency')
    
    ticks = list(range(0, max_count + 1))
    labels = [str(i) for i in range(0, max_count)] + [f"=>{max_count}"]
    plt.xticks(ticks, labels, rotation=90)
    
    plt.show()


def oneline_print_alignment(alignment):
    aligned_seq1_str = str(alignment[0])
    aligned_seq2_str = str(alignment[1])
    match_str = ''.join('|' if a == b and a != '-' else '.' if a != '-' and b != '-' else ' ' for a, b in zip(aligned_seq1_str, aligned_seq2_str))

    print(aligned_seq1_str)
    print(match_str)
    print(aligned_seq2_str)


def blast_df_of_blast_files(blast_paths):

    dfs = []

    for blast_file in blast_paths:
        blast_results = pd.read_csv(blast_file, sep="\t", header=None, names=["query", "subject", "identity", "length", "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"])
        # Add a "dataset" column to the dataframe
        blast_results["dataset"] = os.path.basename(blast_file).split('.')[0]

        fasta_file = blast_file.replace(".blast.tsv", ".fasta")
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        # Note that the sequences in the FASTA file are on the reverse strand, so we're doing reverse complement here.
        matching_seqs = {name: str(seq_dict[name].seq) for name in blast_results['subject'] if name in seq_dict}
        seq_df = pd.DataFrame(list(matching_seqs.items()), columns=['subject', 'sequence'])
        merged_df = pd.merge(blast_results, seq_df, on='subject', how='inner')
        dfs.append(merged_df)

        # Concatenate the dataframes together
        blast_df = pd.concat(dfs, ignore_index=True)

        blast_df = blast_df[blast_df["length"] == 20]

        assert len(blast_df) == len(set(blast_df["subject"]))

        # assert that s_start is greater than s_end for every row, meaning that the subject is on the positive strand
        assert np.all(blast_df["s_start"] < blast_df["s_end"])

        blast_df = blast_df.drop(columns=['query'])
        return blast_df

        
def correct_alignment(aligned_seq):
    # For some reason, the aligner adds a bunch of gaps at the end of the sequence
    # so that the last base of the sequence is aligned to the last base of the template.
    # Here we shift those gaps over. See below for demonstrative tests of what
    # we want to fix.
    aligned_seq = str(aligned_seq)
    pattern = r'(-+)([ACGTN])$'
    def replace_pattern(match):
        gaps, base = match.groups()
        return base + '-' * len(gaps)
    corrected_seq = re.sub(pattern, replace_pattern, aligned_seq)
    return corrected_seq


def test_correct_alignment():
    aligned_seq1 = "...AGAAATAAAA----------C"
    correct_seq1 = "...AGAAATAAAAC----------"
    aligned_seq2 = "...AGAAATAAAA-----A----C"
    correct_seq2 = "...AGAAATAAAA-----AC----"
    aligned_seq3 = "...AGAAATAAAA---------AC"
    correct_seq3 = "...AGAAATAAAA---------AC"
    aligned_seq4 = "...AGAAATA--AAA---------"
    correct_seq4 = "...AGAAATA--AAA---------"
    assert correct_alignment(aligned_seq1) == correct_seq1
    assert correct_alignment(aligned_seq2) == correct_seq2
    assert correct_alignment(aligned_seq3) == correct_seq3
    assert correct_alignment(aligned_seq4) == correct_seq4


def oneline_print_alignment(alignment):
    aligned_seq1_str = str(alignment[0])
    aligned_seq2_str = str(alignment[1])

    # Print sequence 1
    print(aligned_seq1_str)
    
    match_str = ''.join('|' if a == b and a != '-' else '.' if a != '-' and b != '-' else ' ' for a, b in zip(aligned_seq1_str, aligned_seq2_str))

    print(match_str)
    
    # Print sequence 2
    print(aligned_seq2_str)


def summarize_alignment(alignment):
    """
    This function takes in an alignment object where the first sequence is the
    template, and returns a list with
    the sequence of gap lengths in the template and the positions of mutations,
    indexed by non-gap sites in the template sequence.
    We do not consider read gaps as mutations because most of the time they are
    just from sequence length variation.
    
    :param alignment: Biopython alignment object
    :return: Tuple (list of gap lengths in the template, list of mutation positions)
    """
    gap_lengths = []
    mutation_positions = []
    n_positions = []
    mutation_bases = []

    # Extracting aligned sequences from alignment object
    template_seq = alignment[0]
    read_seq = alignment[1]

    gap_count = 0
    non_gap_index = 0  # This is used to index positions by non-gap sites in the template sequence

    for template_base, read_base in zip(template_seq, read_seq):
        if template_base == '-':
            gap_count += 1
        else:
            if gap_count > 0:
                gap_lengths.append(gap_count)
                gap_count = 0  # Reset the gap count after a gap sequence has ended

            if template_base != read_base and read_base != '-' and read_base != 'N':
                mutation_positions.append(non_gap_index)
                mutation_bases.append(read_base)
            
            if read_base == 'N':
                n_positions.append(non_gap_index)
            
            non_gap_index += 1  # Increment the non-gap index whenever a non-gap character is encountered in the template
    
    # In case the sequence ends with gaps
    if gap_count > 0:
        gap_lengths.append(gap_count)
    
    # Because we're aligning to the trimmed template we should start and end with a gap.
    bookended_by_gaps = template_seq[0] == '-' and template_seq[-1] == '-'
    
    return gap_lengths, bookended_by_gaps, mutation_positions, mutation_bases, n_positions


def test_summarize_alignment():
    """
    original alignment:
    --AAA---AA
    A-C-AAAACN

    after deleting columns with gaps in template:
    AAAAA
    C-ACN

    So we have mutations at positions 0, 3 (site 1 is a gap)
    """
    alignment = ['--AAA---AA', 'A-C-AAAACN']
    
    # Call the function with this manual alignment
    gap_lengths, bookended_by_gaps, mutation_positions, mutation_bases, n_positions = summarize_alignment(alignment)

    # Expected outputs
    expected_gap_lengths = [2, 3]
    expected_mutation_positions = [0, 3]
    expected_mutation_bases = ["C", "C"]
    expected_n_positions = [4]
    
    assert gap_lengths == expected_gap_lengths, f"Expected {expected_gap_lengths}, but got {gap_lengths}"
    assert bookended_by_gaps == False
    assert mutation_positions == expected_mutation_positions, f"Expected {expected_mutation_positions}, but got {mutation_positions}"
    assert mutation_bases == expected_mutation_bases, f"Expected {expected_mutation_bases}, but got {mutation_bases}"
    assert n_positions == expected_n_positions, f"Expected {expected_n_positions}, but got {n_positions}"
    

def make_mutation_histogram(chigy_believable):
    # Find the maximum mutation count to set the bin edges
    max_mutation_count = chigy_believable['mutation_count'].max()
    bin_edges = range(0, max_mutation_count + 2)  # +1 for inclusive last bin, +1 for the rightmost edge

    # Set up the FacetGrid
    g = sns.FacetGrid(chigy_believable, row="dataset", height=4, aspect=2, sharex=True, sharey=False)

    # Map a histogram to each subplot
    g.map(plt.hist, "mutation_count", bins=bin_edges, align='left', rwidth=0.8)

    # Adjust the subplot parameters for better layout
    g.fig.subplots_adjust(hspace=0.4)  # adjust the space between plots

    # Set y-axis label and make the layout tight
    g.set_axis_labels(y_var="Frequency")
    plt.tight_layout()

    # Show the plots
    plt.show()


def compute_thing_counts(df, thing_name):
    """
    Compute mutation counts by position from a filtered dataframe.
    """
    max_position = df[thing_name].apply(lambda x: max(x) if x else 0).max()
    mutation_counts_by_position = np.zeros(max_position + 1, dtype=int)
    
    for positions_list in df[thing_name]:
        if positions_list:
            for position in positions_list:
                mutation_counts_by_position[position] += 1
    
    return mutation_counts_by_position


def mutation_frequency_by_position_of(filtered_df):
    mutation_counts_by_position = compute_thing_counts(filtered_df, "mutation_positions")
    total_read_count = len(filtered_df)
    return mutation_counts_by_position / total_read_count


def compute_mutation_counts_by_base(df, wt_sequence):
    """
    Compute mutation counts by position and by mutated base from a filtered dataframe.
    Entries corresponding to the wild-type base in wt_sequence will have a count of -1.
    Returns a DataFrame with columns 'to_A', 'to_C', 'to_G', 'to_T' and index as position.
    """
    max_position = df['mutation_positions'].apply(lambda x: max(x) if x else 0).max()
    bases = ['A', 'C', 'G', 'T']
    data = {f'to_{base}': [0] * (max_position + 1) for base in bases}

    for index, row in df.iterrows():
        for position, base in zip(row['mutation_positions'], row['mutation_bases']):
            if base in bases:
                data[f'to_{base}'][position] += 1
                
    for position, base in enumerate(wt_sequence): 
        if position <= max_position:
            data[f'to_{base}'][position] = -1

    mutation_counts_df = pd.DataFrame(data)
    return mutation_counts_df


def test_compute_mutation_counts_by_base():
    wt_sequence = "AGTTC"
    test_df = pd.DataFrame({
        'mutation_positions': [[ 1,   2],   [2,   3],   [3,   4], [  1]],
        'mutation_bases':     [['A', 'C'], ['G', 'A'], ['A', 'G'], ['C']]
    })
    expected_df = pd.DataFrame({
        'to_A': [-1, 1, 0, 2, 0],
        'to_C': [0, 1, 1, 0, -1],
        'to_G': [0, -1, 1, 0, 1],
        'to_T': [0, 0, -1, -1, 0],
    })

    result_df = compute_mutation_counts_by_base(test_df, wt_sequence)
    assert result_df.equals(expected_df), f"Expected:\n{expected_df}\nGot:\n{result_df}"


def create_mutation_heatmap(data_df):
    # Replace the column names
    data_df = data_df.rename(columns=lambda x: x.split('_')[-1])
    
    # Replace -1 with NaN
    data_df = data_df.replace(-1, np.nan)

    # Setting the dimensions of the plot
    fig, ax = plt.subplots(figsize=(len(data_df.columns) * 0.5, len(data_df) * 0.5))
    
    # Define the custom color map
    cmap = sns.light_palette("purple", as_cmap=True)
    cmap.set_bad('#FFFFE0')  # This will set NaN color to a less intense yellow

    # Plotting the heatmap
    sns.heatmap(data_df, cmap=cmap, linewidths=.5, ax=ax, cbar_kws={"label": "Mutation Count"},
                vmin=0, vmax=data_df.max().max(), square=True, yticklabels=True)

    ax.set_ylabel('Sites')
    ax.set_xlabel('Mutations to')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')

    plt.tight_layout()
    plt.show()


class Passenger:
    def __init__(self, chigy_stop_trimmed, region_dict):
        assert chigy_stop_trimmed == "".join(region_dict.values())
        self.chigy_stop_trimmed = chigy_stop_trimmed
        self.region_dict = region_dict
        self.aligner = PairwiseAligner()
        self.aligner.match_score = 1
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -5
        self.aligner.extend_gap_score = -2

    def run_alignment(self, template, sequence):
        # Run and correct the alignment, returning a tuple of strings.
        alignments = self.aligner.align(template, sequence)
        assert alignments is not None, "No alignment found."
        alignment = alignments[0]
        alignment = (correct_alignment(alignment[0]), str(alignment[1]))
        return alignment
    
    def pretty_print_alignments(self, alignment_df):
        for _, row in alignment_df.iterrows():
            sequence = row["sequence"]
            alignment = self.run_alignment(self.chigy_stop_trimmed, sequence)
            print(f"> {row['dataset']} {row['Sequence number']}")
            oneline_print_alignment(alignment)
            print()

    def perform_alignment_and_summary(self, row):
        """Note that the mutation_positions are with repsect to the trimmed template
        sequence, whereas the gap_lengths are with respect to the original template
        sequence."""
        sequence = row["sequence"]
        alignment = self.run_alignment(self.chigy_stop_trimmed, sequence)
        gap_lengths, bookended_by_gaps, mutation_positions, mutation_bases, n_positions = summarize_alignment(alignment)
        
        return pd.Series({
            'gap_segment_count': len(gap_lengths),
            'bookended_by_gaps': bookended_by_gaps,
            'mutation_positions': mutation_positions,
            'mutation_bases': mutation_bases,
            'n_positions': n_positions,
        })


    def processed_stop_df_of_blast_df(self, blast_df):
        alignment_summary_df = blast_df.apply(self.perform_alignment_and_summary, axis=1)

        # Concatenating the new DataFrame with the original one horizontally
        processed_stop_df = pd.concat([blast_df, alignment_summary_df], axis=1)

        processed_stop_df["mutation_count"] = processed_stop_df["mutation_positions"].apply(len)
        processed_stop_df["n_count"] = processed_stop_df["n_positions"].apply(len)
        processed_stop_df["Sequence number"] = processed_stop_df.index

        return processed_stop_df


    def make_mutation_rate_plot(self, chigy_believable):
        # Calculate the positions where each sequence ends in the full sequence
        labels = self.region_dict.keys()
        positions = np.cumsum([len(self.region_dict[key]) for key in labels])
        colors = ['blue', 'red', 'green', 'yellow', 'purple', 'orange', 'cyan']

        # Get the unique datasets
        unique_datasets = chigy_believable['dataset'].unique()

        # Set up the figure and the array of subplots
        fig, axs = plt.subplots(len(unique_datasets), 1, figsize=(12, 4 * len(unique_datasets)), dpi=150)

        # Iterate over unique datasets and create a subplot for each one
        for index, dataset in enumerate(unique_datasets):
            ax = axs[index] if len(unique_datasets) > 1 else axs
            filtered_df = chigy_believable[chigy_believable['dataset'] == dataset]
            
            mutation_frequency_by_position = mutation_frequency_by_position_of(chigy_believable)
            
            ax.scatter(range(len(mutation_frequency_by_position)), mutation_frequency_by_position, alpha=0.5, color="black")
            
            start_pos = 0
            for pos, color, label in zip(positions, colors, labels):
                ax.axvspan(start_pos, pos, facecolor=color, alpha=0.2)
                ax.text((start_pos + pos) / 2, ax.get_ylim()[1] * 0.95, label, horizontalalignment='center')
                start_pos = pos + 1e-9  # Add a tiny offset to avoid overlapping
            
            ax.set_title(f'Mutation Frequency by Position - Dataset: {dataset}')
            ax.set_xlabel('Position')
            ax.set_ylabel('Mutation Frequency')

        plt.tight_layout()
        plt.show()


CHIGY_LC_STOP_TRIMMED = "GACATTGTGATGACtCAGTCTCAAAAATTCATGTCCACATCAGTAGGAGACAGGGTCAGCGTCACCTGCAAGGCCAGTCAGAATGTGGGTACTAATGTAGCCTGGTATCAACAGAAACCAGGGCAATCTCCTAAAGCACTGATTTACTCGGCATCCTACAGGTACAGTGGAGTCCCTGATCGCTTCACAGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAATGTGCAGTCTGAAGACTTGGCAGAGTATTTCTGTCAGCAATATAACAGCTATCCTCTCACGTTCGGCTCGGGGACtAAGCTaGAAATAAAAC".upper()

LC_REGION_DICT = {
    "FW1": "GACATTGTGATGACTCAGTCTCAAAAATTCATGTCCACATCAGTAGGAGACAGGGTCAGCGTCACCTGC",
    "CDR1": "AAGGCCAGTCAGAATGTGGGTACTAATGTAGCC",
    "FW2": "TGGTATCAACAGAAACCAGGGCAATCTCCTAAAGCACTGATTTAC",
    "CDR2": "TCGGCATCCTACAGGTACAGT",
    "FW3": "GGAGTCCCTGATCGCTTCACAGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAATGTGCAGTCTGAAGACTTGGCAGAGTATTTCTGT",
    "CDR3": "CAGCAATATAACAGCTATCCTCTCACG",
    "FW4": "TTCGGCTCGGGGACTAAGCTAGAAATAAAAC"
}

CHIGY_LC = Passenger(CHIGY_LC_STOP_TRIMMED, LC_REGION_DICT)
