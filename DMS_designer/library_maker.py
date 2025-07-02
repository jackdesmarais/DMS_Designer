"""
Library Maker for Protein Mutagenesis Studies

This module provides tools for creating oligonucleotide libraries for deep mutational
scanning and protein mutagenesis studies. It includes functionality for generating
single and multiple amino acid mutants, creating synonymous nucleotide variants,
and analyzing library composition.

The module reads codon usage data from 'codon_usage.xlsx' to generate realistic
nucleotide sequences that encode the desired protein variants.

Classes
-------
LibraryMaker
    Main class for creating mutagenesis libraries

Functions
---------
aa_to_codon : Convert amino acid to codon using codon usage table
protein_to_orf : Convert protein sequence to open reading frame
orf_to_protein : Convert open reading frame to protein sequence
summarize_file : Generate analysis plots and statistics for library files

Examples
--------
>>> # Create a library with single and double mutants
>>> lm = LibraryMaker(wt_nt="ATGGCCGAA", single_mutants=True, doubles_to_make='equal')
>>> 
>>> # Analyze a library file
>>> summarize_file('my_library.csv', out_path='./analysis/')
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import logomaker
import mavenn
from os import path 
import itertools
from math import comb
import seaborn as sns

codon_df = pd.read_excel('codon_usage.xlsx')
codon_df.set_index('codon', inplace=True)

amino_acids = np.array([aa for aa in codon_df['aa'].unique()])
amino_acids.sort()

aa_to_codon_dict = {}
for aa in codon_df['aa'].unique():
    ix = (codon_df['aa']==aa) & (codon_df['use'])
    assert sum(ix)>0, f'Amino acid "{aa}" has no codon. Abort!'
    codons = codon_df.index[ix].values
    aa_to_codon_dict[aa] = codons
    
def aa_to_codon(aa, rng=None):
    """
    Convert an amino acid to a codon according to intrinsic probabilities.
    
    Parameters
    ----------
    aa : str
        Single amino acid character to convert to codon
    rng : numpy.random.Generator, optional
        Random number generator instance. If None, creates a new default generator.
        
    Returns
    -------
    str
        Three-nucleotide codon encoding the amino acid
        
    Raises
    ------
    AssertionError
        If the amino acid is not found in the codon usage table
        
    Notes
    -----
    Uses the codon usage probabilities defined in the loaded 'codon_usage.xlsx' file.
    Only codons marked as 'use=True' in the table are considered.
    """
    if rng is None:
        rng = np.random.default_rng()
    assert aa in aa_to_codon_dict, f'Amino acid "{aa}" has no codon. Abort!'
    codon = rng.choice(a=aa_to_codon_dict[aa], size=1)[0]
    return(codon)

def protein_to_orf(protein, rng=None):
    """
    Convert a protein sequence to an open reading frame using codon usage table.
    
    Parameters
    ----------
    protein : str
        Protein sequence as a string of amino acid characters
    rng : numpy.random.Generator, optional
        Random number generator instance. If None, creates a new default generator.
        
    Returns
    -------
    str
        Nucleotide sequence encoding the protein according to codon usage probabilities
        
    Notes
    -----
    Each amino acid is converted to a codon using the probabilities defined in
    the codon usage table. Multiple calls with the same protein may produce
    different nucleotide sequences due to synonymous codon variation.
    """
    if rng is None:
        rng = np.random.default_rng()
    return(''.join([aa_to_codon(aa,rng=rng) for aa in protein]))

def orf_to_protein(orf):
    """
    Convert an open reading frame to a protein sequence.
    
    Parameters
    ----------
    orf : str
        Nucleotide sequence to translate. Must be a multiple of 3.
        
    Returns
    -------
    str
        Protein sequence translated from the nucleotide sequence
        
    Raises
    ------
    AssertionError
        If the nucleotide sequence length is not a multiple of 3
        
    Notes
    -----
    Uses the codon-to-amino acid mapping defined in the loaded 'codon_usage.xlsx' file.
    """
    # Convert an ORF to protein sequence
    assert len(orf)%3==0
    num_codons = int(len(orf)/3)
    codons = [orf[(i*3):(i*3+3)] for i in range(num_codons)]
    return(''.join([codon_df['aa'][codon] for codon in codons]))

def summarize_file(file_name, out_path = './', fig_width=12, context='talk'):
    """
    Generate comprehensive analysis and visualization of a mutagenesis library file.
    
    This function reads a CSV file containing library variants and generates multiple
    plots and statistics including variant counts, mutant distributions, heatmaps,
    and position-by-position analysis.
    
    Parameters
    ----------
    file_name : str
        Path to the CSV file containing the library data. Expected columns:
        'name', 'protein', 'orf', 'oligo'
    out_path : str, default './'
        Directory path where output figures will be saved
    fig_width : int, default 12
        Base width for generated figures. Height is typically equal to width
    context : {'talk', 'paper', 'notebook', 'poster'}, default 'talk'
        Seaborn plotting context for figure styling
        
    Returns
    -------
    None
        
    Notes
    -----
    The function generates the following output files:
    - variant_count_[filename].pdf: Count of variants by mutation number
    - mutant_count_[filename].pdf: Count of unique mutants by mutation number  
    - mutant_fraction_[filename].pdf: Fraction of possible mutants generated
    - mutant_observations_[filename].pdf: Distribution of observations per mutant
    - [N]pt_heatmap_[filename].pdf: Heatmaps showing amino acid frequencies at each position
    - [N]pt_mutant_heatmap_[filename].pdf: Heatmaps for unique mutants only
    - 2pt_2d_heatmap_[filename].pdf: Position-by-position matrix for double mutants
        
    The function also prints summary statistics to the console including:
    - Wild-type sequence information
    - Counts of variants by mutation number
    - Total library size
    - Fraction of possible mutants generated
    """
    sns.set_context(context)
    sns.set_palette('colorblind')

    print(f'Processing {file_name}')
    df = pd.read_csv(file_name)
    df = df.join(df.name.str.extract(r'(?P<mut_num>\d+)pt_(?P<mutation_desc>.+)_v(?P<version>\d+)'))
    df.mut_num = df.mut_num.astype(int)
    df = df.set_index('name', drop=True)
    file_name = path.basename(file_name)

    wt_aa = df.loc['0pt_WT_v1','protein']
    print(f'WT aa seq (L={len(wt_aa)}):\n\t{wt_aa}')

    wt_orf = df.loc['0pt_WT_v1','orf']
    print(f'WT dna seq (L={len(wt_orf)}):\n\t{wt_orf}')

    wt_oligo = df.loc['0pt_WT_v1','oligo']
    print(f'WT oligo seq (L={len(wt_oligo)}):\n\t{wt_oligo}')

    fig, ax = plt.subplots(1,1,figsize=[fig_width,fig_width])
    mut_nums, counts = np.unique(df.mut_num, return_counts=True)
    plt.scatter(mut_nums, counts)
    plt.xlabel("Mutation number")
    plt.ylabel("Variant count")
    plt.xticks(mut_nums, mut_nums)
    plt.savefig(out_path+f'variant_count_{file_name[:-4]}.pdf')
    for mut_num, count in zip(mut_nums, counts):
        print(f'Number of {mut_num} variants: {count}')
    print(f'Total number of variants: {len(df)}')


    fig, ax = plt.subplots(1,1,figsize=[fig_width,fig_width])
    mut_nums, counts = np.unique(df.drop_duplicates(subset='mutation_desc').mut_num, return_counts=True)
    plt.scatter(mut_nums, counts)
    plt.xlabel("Mutation number")
    plt.ylabel("Mutant count")
    plt.xticks(mut_nums, mut_nums)
    plt.savefig(out_path+f'mutant_count_{file_name[:-4]}.pdf')
    for mut_num, count in zip(mut_nums, counts):
        print(f'Number of {mut_num} mutant variants: {count}')

    fig, ax = plt.subplots(1,1,figsize=[fig_width,fig_width])
    mut_nums, counts = np.unique(df.drop_duplicates(subset='mutation_desc').mut_num, return_counts=True)
    total_possible = [comb(len(wt_aa), mut_num)*(19**mut_num) if mut_num!=1 else comb(len(wt_aa), mut_num)*(20**mut_num) for mut_num in mut_nums]
    plt.scatter(mut_nums, counts/total_possible)
    plt.xlabel("Mutation number")
    plt.ylabel("Fraction of possible mutants")
    plt.xticks(mut_nums, mut_nums)
    plt.savefig(out_path+f'mutant_fraction_{file_name[:-4]}.pdf')
    for mut_num, count in zip(mut_nums, counts/total_possible):
        print(f'Fraction of possible {mut_num} mutants: {count}')
    
    fig, ax = plt.subplots(1,1,figsize=[fig_width,fig_width])
    counts = df[['mut_num','mutation_desc','version']].groupby(['mut_num','mutation_desc']).count().reset_index()
    p = sns.boxplot(data=counts,x='mut_num', y='version', showfliers=False, ax=ax)
    sns.stripplot(data=counts,x='mut_num', y='version', ax=p)
    plt.ylabel('Number of observations per mutant')
    plt.xlabel('Mutant number')
    plt.savefig(out_path+f'mutant_observations_{file_name[:-4]}.pdf')

    # Compute the number of variants with each amino acid at each position
    wt_mat = logomaker.alignment_to_matrix(df[df.mut_num==0]['protein'])
    wt_mat.index = wt_mat.index+1
    wt_mask=wt_mat>0
    for name, group in df.groupby('mut_num'):
        if name>0:
            ct_mat = logomaker.alignment_to_matrix(pd.concat([group['protein'],pd.Series([wt_aa])]))
            ct_mat.index = ct_mat.index+1
            ct_mat[wt_mask]=np.nan
            alphabet = ct_mat.columns.values
        
            n_pos = len(ct_mat.index)
            n_aas = len(alphabet)

            print(f'Total counts in {name} mutant heatmap: {ct_mat.sum().sum()}')


            fig, ax = plt.subplots(1,1,figsize=[fig_width*(3/2),fig_width/1.5])
            cmap = mpl.cm.magma
            cmap.set_bad('#516572')
            mavenn.heatmap(ct_mat, 
                        ax=ax,
                        alphabet=alphabet, 
                        clim=[ct_mat.min().min(),ct_mat.max().max()], 
                        cmap=cmap, 
                        seq=wt_aa, 
                        seq_kwargs={'color':'w', 's':5})
            ax.set_title(f'{name}pt variants');
            ax.set_xlabel('variant position')
            ax.set_ylabel('variant aa')

            # There must be an easier way to do this; should build into MAVE-NN
            xticks = np.array(range(4,len(wt_aa),5))
            ax.set_xticks(xticks)
            ax.set_xticklabels([f'{x+1}' for x in xticks])

            heatmap_name = f'{name}pt_heatmap_{file_name[:-4]}.pdf'
            fig.savefig(out_path+heatmap_name)
            print(f'Heatmap saved as {heatmap_name}')
            
            fig, ax = plt.subplots(1,1,figsize=[fig_width,fig_width])
            value, count = np.unique(ct_mat.stack(), return_counts=True)
            plt.scatter(value, count)
            plt.title(f'Number of {name}pt variants per single mutant', fontsize=18)
            plt.xlabel(f'Number of {name}pt variants', fontsize=18)
            plt.ylabel('Number of single mutants', fontsize=18)
            heatmap_name = f'2pt_doubles_p_singles_{file_name[:-4]}.pdf'
            fig.savefig(out_path+heatmap_name)


            #recalculate for mutants not variants
            ct_mat = logomaker.alignment_to_matrix(pd.concat([group.drop_duplicates(subset='mutation_desc')['protein'],
                                                              pd.Series([wt_aa])]))
            ct_mat.index = ct_mat.index+1
            ct_mat[wt_mask]=np.nan
            alphabet = ct_mat.columns.values
        
            n_pos = len(ct_mat.index)
            n_aas = len(alphabet)

            print(f'Mutant observances in {name} mutant heatmap: {ct_mat.sum().sum()}')


            fig, ax = plt.subplots(1,1,figsize=[fig_width*(3/2),fig_width/1.5])
            cmap = mpl.cm.magma
            cmap.set_bad('#516572')
            mavenn.heatmap(ct_mat, 
                        ax=ax,
                        alphabet=alphabet, 
                        clim=[ct_mat.min().min(),ct_mat.max().max()], 
                        cmap=cmap, 
                        seq=wt_aa, 
                        seq_kwargs={'color':'w', 's':5})
            ax.set_title(f'{name}pt mutants');
            ax.set_xlabel('mutant position')
            ax.set_ylabel('mutant aa')

            # There must be an easier way to do this; should build into MAVE-NN
            xticks = np.array(range(4,len(wt_aa),5))
            ax.set_xticks(xticks)
            ax.set_xticklabels([f'{x+1}' for x in xticks])

            heatmap_name = f'{name}pt_mutant_heatmap_{file_name[:-4]}.pdf'
            fig.savefig(out_path+heatmap_name)
            print(f'Heatmap saved as {heatmap_name}')
            
            fig, ax = plt.subplots(1,1,figsize=[fig_width,fig_width])
            value, count = np.unique(ct_mat.stack(), return_counts=True)
            plt.scatter(value, count)
            plt.title(f'Number of {name}pt mutants per single mutant', fontsize=18)
            plt.xlabel(f'Number of {name}pt mutants', fontsize=18)
            plt.ylabel('Number of single mutants', fontsize=18)
            heatmap_name = f'2pt_mutant_doubles_p_singles_{file_name[:-4]}.pdf'
            fig.savefig(out_path+heatmap_name)
    
    ct_2pt = sum(df.mut_num == 2)
    print(f'Number of 2pt variants: {ct_2pt}')
    if ct_2pt > 0:
        n_pos = len(wt_aa)
        n_aas = len(amino_acids)
        # Compute the number of variants with each amino acid at each position
        df_2pt = df[df.mut_num == 2]

        mut_aa = np.stack(df_2pt['protein'].apply(lambda x: np.array(list(x))))
        wt_aa_arr = np.array(list(wt_aa))
        seq, position = np.where(wt_aa_arr[None,:] != mut_aa)
        aas = mut_aa[seq, position]
        seqs = len(seq)//2
        seq = seq.reshape((seqs,2))
        assert (seq[:,0] == seq[:,0]).all()
        position = position.reshape((seqs,2))
        aas = aas.reshape((seqs,2))

        s_id, n_id, aa_id = np.where(amino_acids[None,None,:] == aas[:,:,None])
        aas = np.zeros(aas.shape, dtype=int)
        aas[s_id, n_id] = aa_id
        assert (aas!=0).all()

        pt2_mat = np.zeros((n_pos,n_pos,n_aas,n_aas))
        for (i, j), (aa_i, aa_j) in zip(position, aas):
            pt2_mat[i,j,aa_i, aa_j] += 1

        plt.figure(figsize=(fig_width,fig_width))
        plt.imshow(pt2_mat.sum(axis=-1).sum(axis=-1), cmap='magma')
        plt.title('Position by position 2pt counts matrix', fontsize=18)
        plt.xlabel('Position', fontsize=18)
        plt.ylabel('Position', fontsize=18)
        cbar = plt.colorbar(label='Total counts for 2pt pairs')
        cbar.ax.set_ylabel('Total counts for 2pt pairs',fontsize=18, rotation=270, labelpad=15)
        heatmap_name = f'2pt_2d_heatmap_{file_name[:-4]}.pdf'
        fig.savefig(out_path+heatmap_name)
        
    
    

class LibraryMaker():
    """
    Create oligonucleotide libraries for protein mutagenesis studies.
    
    This class generates comprehensive mutagenesis libraries including single and
    multiple amino acid mutants, synonymous nucleotide variants, and focused
    mutagenesis at specific positions. The library is designed for deep mutational
    scanning experiments.
    
    Parameters
    ----------
    wt_nt : str
        The nucleotide sequence of the wild-type variant. This sequence should
        encode an amino acid sequence, start on the first nucleotide of the first
        codon to be mutagenized, and end on the last nucleotide of the last
        codon to be mutagenized.
    single_mutants : bool, default True
        If True, generate every possible single amino acid mutant of the wild-type
        variant.
    doubles_to_make : str, float, or int, default 'equal'
        Controls the number of double mutants to generate:
        
        - 'equal': Generate random double mutants equal to the number of single mutants
        - 'max:N': Generate as many random double mutants as possible without exceeding
          N total oligos for synthesis
        - float: Generate that fraction of the number of single mutants as random
          double mutants
        - int: Generate exactly that number of random double mutants
    name : str, default 'Mutant_lib'
        Name of the library to be created. Used as prefix for output files.
    var_num : int, default 3
        Number of different synonymous nucleotide variants to generate for each mutant.
    wt_num : int, default 100
        Number of different synonymous nucleotide variants to generate for the
        wild-type sequence.
    max_tries : int, default 100
        Maximum number of attempts to generate a unique random mutant or variant
        before giving up.
    focus_positions : list of int, default []
        Amino acid positions (0-indexed) to be mutated in greater depth.
    focus_to_make : int, float, or str, default 0
        Controls focused mutagenesis when focus_positions is specified:
        
        - int/float: Generate that many random combinations of mutations at focus positions
        - 'singles': Generate all single mutants of the focus positions
        - 'doubles': Generate all double mutants of the focus positions  
        - 'singles_doubles': Generate all single and double mutants of the focus positions
        - 'all': Generate every possible combination of mutations at the focus positions
    focus_mutants_per : int, default 3
        Average number of mutations desired per focused mutant (used when focus_to_make
        is numeric).
    focus_mutants_per_spread : float, default 0
        Probability spread for mutation numbers other than focus_mutants_per.
        Probability decays arithmetically as you move away from the target.
    focus_mutants_per_min : int, default 0
        Minimum number of mutations per focused mutant.
    focus_mutants_per_max : int, default 0
        Maximum number of mutations per focused mutant.
    upstream_dna : str, default 'ATGGAAAGTGTCCCAGGAGACTACAGC'
        Upstream adapter sequence to append to the 5' end of oligonucleotides.
    downstream_dna : str, default 'CAGACAAGTAACTCAAGACTTAATGGC'
        Downstream adapter sequence to append to the 3' end of oligonucleotides.
    out_path : str, default './'
        Directory path for output files. Directory must already exist.
    seed : int, optional
        Random seed for reproducibility. If None, uses system default.
    fig_width : int, default 8
        Base width for generated analysis figures.
    context : {'talk', 'paper', 'notebook', 'poster'}, default 'talk'
        Seaborn plotting context for figure styling.
        
    Attributes
    ----------
    wt_nt : str
        Wild-type nucleotide sequence
    wt_protein : str
        Wild-type protein sequence
    L : int
        Length of the protein sequence
    variant_dict : dict
        Dictionary mapping variant names to protein sequences
    orf_dict : dict
        Dictionary mapping variant names to nucleotide sequences
    orf_df : pandas.DataFrame
        Final library data with columns: name, protein, orf, oligo
    file_name : str
        Name of the saved CSV file
    date_time : str
        Timestamp of library creation
        
    Notes
    -----
    The library generation process:
    1. Creates single mutants (if enabled)
    2. Creates random double mutants (if enabled)
    3. Creates focused mutants (if focus positions specified)
    4. Generates synonymous nucleotide variants for each protein variant
    5. Saves results to CSV file
    6. Generates analysis plots and statistics
    
    Examples
    --------
    >>> # Create a simple library with single mutants only
    >>> lm = LibraryMaker("ATGGCCGAA", single_mutants=True, doubles_to_make=0)
    >>> 
    >>> # Create a comprehensive library with focused mutagenesis
    >>> lm = LibraryMaker("ATGGCCGAA", 
    ...                   single_mutants=True,
    ...                   doubles_to_make='equal',
    ...                   focus_positions=[5, 10, 15],
    ...                   focus_to_make='all')
    """
    def __init__(self, wt_nt, single_mutants = True, doubles_to_make = 'equal', name='Mutant_lib',
                 var_num = 3, wt_num = 100, max_tries = 100,
                 focus_positions = [], focus_to_make=0, focus_mutants_per=3, focus_mutants_per_spread=0,focus_mutants_per_min=0,focus_mutants_per_max=0,
                upstream_dna = 'ATGGAAAGTGTCCCAGGAGACTACAGC',
                 downstream_dna = 'CAGACAAGTAACTCAAGACTTAATGGC',
                out_path = './', seed=None,
                fig_width=8, context='talk'):
        """
        This class creates oligo libraries for mutagenesis studies of protein function

        Parameters
        ----------
            wt_nt : string
                The nucleotide sequence of the WT variant. This sequence should encode an AA sequence, 
                start on the first NT of the first codon to be mutagenized, and end on the last nt of the last codon to be mutagenized
            single_mutants : bool Default is True
                If true, make every possible single AA mutant of the WT variant
            doubles_to_make : string, float, or int Default is 'equal
                If this is the string 'equal' make random double mutants equal to the number of single mutants made
                If this is a string starting with 'max:' follwed by an int, make as many random double mutants as possible 
                without allowing the number of oligos required to synthesize single and double mutants to exceed the number listed
                If this is a float, Make that fraction of the number of single mutants as random duble mutants
                If this is an int, Make that number of random double mutants
            name : string Default is 'Mutant_lib'
                name of the library to be created, will be appended as a prefix to the output files
            var_num : int default 3
                The number of different synonymous nt variants to be made for each mutant
            wt_num : int default 100
                The number of different synsonymous nt vatriants to be made for the WT sequence
            max_tries : int default 100
                The number of times to retry making a random mutant or variant if an already existing sequence was created
            focus_positions : list default []
                Amino acid positions of AAs to be mutated in greater depth
            focus_to_make : int, float, or {'singles', 'doubles', 'singles_doubles', 'all'} default 0
                This parameter is only used of focus_positions contains at least 1 element
                If this is an int or a float, make that many random combinations of mutations at the focus positions 
                If this is the string 'singles', make all single mutants of the focus positions
                If this is the string 'doubles', make all double mutants of the focus positions
                If this is the string 'singles_doubles', make all single and double muatants of the focus positions
                If this is the string 'all', make every combination of mutations at the focus positions
            focus_mutants_per : int default 0
                This parameter is only used of focus_positions contains at least 1 element and focus_to_make is an int or float >0
                this indicates the average number of mutations desired per mutant
            focus_mutants_per_spread : int default 0
                This parameter is only used of focus_positions contains at least 1 element and focus_to_make is an int or float >0
                This is the probability of producing a mutation number other than focus_mutants_per, mutation probability will be evenly distributed
                between higher and lower mutation numbers with arithmatic decay as you move away from focus_mutants_per
            focus_mutants_per_min
                This parameter is only used of focus_positions contains at least 1 element and focus_to_make is an int or float >0
                This is the minimum mutation number you want to see per mutant
            focus_mutants_per_max
                This parameter is only used of focus_positions contains at least 1 element and focus_to_make is an int or float >0
                This is the Maximum mutation number you want to see per mutant
            upstream_dna : string default 'ATGGAAAGTGTCCCAGGAGACTACAGC'
                This is the upstream adaptor sequence to append to the 5' end of oligos
            downstream_dna: string default 'CAGACAAGTAACTCAAGACTTAATGGC'
                This is the downstream adaptor sequence to append to the 3' end of oligos
            out_path : string default './'
                path to place output files. This directory must already exist
            seed : int Default None
                Random seed. Set to ensure reproducibility.
            fig_width : int default 8
                Size parameter for output figures
            context : {'talk', 'paper', 'notebook', 'poster'} default 'talk'
                Used to set plotting configuration. passed to seaborn.set_context
        """
        
        start = datetime.now()
        self.date_time = start.strftime("%y.%m.%d.%Hh.%Mm.%Ss")

        self.name = name

        self.rng = np.random.default_rng(seed=seed)
        
        self.wt_nt = wt_nt
        self.wt_protein = orf_to_protein(wt_nt)
        self.L = len(self.wt_protein)
        self.variant_dict = {'0pt_WT':self.wt_protein}
        
        self.out_path = out_path
        
        self.var_num = var_num
        self.wt_num = wt_num
        self.max_tries = max_tries
        
        self.upstream_dna = upstream_dna
        self.downstream_dna = downstream_dna

        self.focus_positions = focus_positions
        if len(self.focus_positions)>0:
            self.focus_to_make = focus_to_make
        
        if single_mutants:
            print('Making singles')
            self.num_singles = len(self.create_singles())
            print('%d single mutants made\n%f sec elapsed\n'%(self.num_singles,(datetime.now()-start).total_seconds()))
            
        if doubles_to_make == 'equal':
            self.num_doubles = self.num_singles
        elif isinstance(doubles_to_make, float):
            self.num_doubles = doubles_to_make*self.num_singles
        elif isinstance(doubles_to_make, int):
            self.num_doubles = doubles_to_make
        elif doubles_to_make.startswith('max:'):
            doubles_to_make = (int(doubles_to_make[len('max:'):])-self.var_num*self.num_singles-self.wt_num)//self.var_num
            self.num_doubles = max(0, doubles_to_make)
        else:
            assert False, "doubles_to_make has an invalid value"
        
        if self.num_doubles > 0:
            print('Making %d doubles'%self.num_doubles)
            self.doubles_made = len(self.create_doubles())
            print('%d doubles mutants made\n%f sec elapsed\n'%(self.doubles_made,(datetime.now()-start).total_seconds()))


        if len(self.focus_positions) > 0:
            print('Making fgocused mutations')
            if self.focus_to_make == 'singles':
                self.focus_made = len(self.create_focused_singles())
            elif self.focus_to_make == 'doubles':
                self.focus_made = len(self.create_focused_doubles())
            elif self.focus_to_make == 'singles_doubles':
                self.focus_made = len(self.create_focused_singles())
                self.focus_made = self.focus_made + len(self.create_focused_doubles())
            elif self.focus_to_make == 'all':
                self.focus_made = len(self.create_focused_all_combinations())
            elif isinstance(self.focus_to_make, float)|isinstance(self.focus_to_make, int):
                assert (focus_mutants_per>0)&(focus_mutants_per<=len(self.focus_positions)), 'focus_mutants_per must be a valid number'
                self.focus_mutants_per = focus_mutants_per

                assert (focus_mutants_per_spread>=0)&(focus_mutants_per_spread<=1), 'focus_mutants_per_spread must be a valid number'
                self.focus_mutants_per_spread = focus_mutants_per_spread

                assert (focus_mutants_per_min>0)&(focus_mutants_per_min<=focus_mutants_per), 'focus_mutants_per_min must be a valid number'
                self.focus_mutants_per_min = focus_mutants_per_min

                assert (focus_mutants_per_max>=focus_mutants_per)&(focus_mutants_per_max<=len(self.focus_positions)), 'focus_mutants_per_max must be a valid number'
                self.focus_mutants_per_max = focus_mutants_per_max
                self.num_focus= self.focus_to_make
                self.focus_made = len(self.create_focused_random_combinations())
            else:
                raise ValueError('Invalid focus_to_make option')
            print('%d focused mutants made\n%f sec elapsed\n'%(self.focus_made,(datetime.now()-start).total_seconds()))


        
        print('Making nt seqs')
        self.total_seqs = len(self.generate_nt_seqs())
        print('%d nt seqs made\n%f sec elapsed\n'%(self.total_seqs,(datetime.now()-start).total_seconds()))
        
        self.save_file()
        print('File saved\n%f sec elapsed'%((datetime.now()-start).total_seconds()))
        print('------------------------------------------------')
        print()
        print()
        summarize_file(self.file_name, out_path=self.out_path, fig_width=fig_width, context=context)
        print('Output checked\n%f sec elapsed'%((datetime.now()-start).total_seconds()))
        
    def create_focused_singles(self):
        """
        Create all single amino acid mutants at focused positions.
        
        Generates every possible single amino acid substitution at each position
        specified in focus_positions, excluding the wild-type amino acid.
        
        Returns
        -------
        dict
            Dictionary mapping variant names to protein sequences for focused
            single mutants. Names follow pattern '1pt_PPAA' where PP is the
            position (01-indexed) and AA is the amino acid.
            
        Notes
        -----
        Updates self.variant_dict with the new variants.
        """
        focused_dict = {}
        for pos in self.focus_positions:
            wt = self.wt_protein[pos]
            for aa in amino_acids[(amino_acids!=wt)]:
                var_protein = self.wt_protein[:pos] + aa + self.wt_protein[pos+1:]
                name = f'1pt_{pos+1:02d}{aa}'
                focused_dict[name] = var_protein

        # Add to variant dictionary
        self.variant_dict.update(focused_dict)
        return(focused_dict)
    
    def create_focused_doubles(self):
        """
        Create all double amino acid mutants at focused positions.
        
        Generates every possible combination of two amino acid substitutions
        at pairs of positions specified in focus_positions, excluding wild-type
        amino acids and stop codons.
        
        Returns
        -------
        dict
            Dictionary mapping variant names to protein sequences for focused
            double mutants. Names follow pattern '2pt_PP1AA1_PP2AA2' where PP
            is the position (01-indexed) and AA is the amino acid.
            
        Notes
        -----
        Updates self.variant_dict with the new variants.
        Excludes stop codons (*) from possible substitutions.
        """
        focused_dict = {}
        for pos1, pos2 in itertools.combinations(self.focus_positions, 2):
            wt1 = self.wt_protein[pos1]
            wt2 = self.wt_protein[pos2]
            for aa1, aa2 in itertools.product(amino_acids[(amino_acids!=wt1)&(amino_acids!='*')], amino_acids[(amino_acids!=wt2)&(amino_acids!='*')]):
                var_protein = self.wt_protein[:pos1] + aa1 + self.wt_protein[pos1+1:pos2] + aa2 + self.wt_protein[pos2+1:]
                name = f'2pt_{pos1+1:02d}{aa1}_{pos2+1:02d}{aa2}'
                focused_dict[name] = var_protein

        # Add to variant dictionary
        self.variant_dict.update(focused_dict)
        return(focused_dict)
    
    def create_focused_all_combinations(self):
        """
        Create all possible combinations of mutations at focused positions.
        
        Generates every possible combination of amino acid substitutions at the
        focused positions, from single mutations up to mutations at all focused
        positions simultaneously. Excludes wild-type amino acids and stop codons.
        
        Returns
        -------
        dict
            Dictionary mapping variant names to protein sequences for all focused
            mutant combinations. Names follow pattern 'Npt_PP1AA1_PP2AA2_...'
            where N is the number of mutations, PP is position (01-indexed),
            and AA is the amino acid.
            
        Notes
        -----
        Updates self.variant_dict with the new variants.
        Excludes stop codons (*) from possible substitutions.
        Generates combinations of all sizes from 1 to len(focus_positions).
        """
        focused_dict = {}
        for positions in itertools.chain.from_iterable(itertools.combinations(self.focus_positions, r+1) for r in range(len(self.focus_positions))):
            wts = (self.wt_protein[pos] for pos in positions)
            aa_options = (amino_acids[(amino_acids!=wt)&(amino_acids!='*')] for wt in wts)
            for aas in itertools.product(*aa_options):
                var_protein = self.wt_protein[:positions[0]]
                name = f'{len(aas)}pt_'
                cur_poses, next_poses = itertools.tee(positions)
                next(next_poses, None)
                for cur_pos, next_pos, aa in zip(cur_poses, next_poses, aas):
                    var_protein = var_protein+aa+self.wt_protein[cur_pos+1:next_pos]
                    name = name+f'{cur_pos+1:02d}{aa}_'
                var_protein = var_protein+aas[-1]+self.wt_protein[positions[-1]+1:]
                name = name+f'{positions[-1]+1:02d}{aas[-1]}'

                focused_dict[name] = var_protein

        # Add to variant dictionary
        self.variant_dict.update(focused_dict)
        return(focused_dict)
    
    def create_focused_random_combinations(self):
        """
        Create random combinations of mutations at focused positions.
        
        Generates random combinations of amino acid substitutions at focused
        positions according to specified probability distributions and constraints.
        
        Returns
        -------
        dict
            Dictionary mapping variant names to protein sequences for random
            focused mutant combinations. Names follow pattern 'Npt_PP1AA1_PP2AA2_...'
            where N is the number of mutations, PP is position (01-indexed),
            and AA is the amino acid.
            
        Notes
        -----
        Updates self.variant_dict with the new variants.
        
        The number of mutations per variant is controlled by:
        - focus_mutants_per: Target average number of mutations
        - focus_mutants_per_spread: Probability spread around the target
        - focus_mutants_per_min: Minimum mutations per variant
        - focus_mutants_per_max: Maximum mutations per variant
        
        If spread > 0, probability decays arithmetically as you move away
        from the target number of mutations.
        """
        low_num = self.focus_mutants_per-self.focus_mutants_per_min
        high_num = self.focus_mutants_per_max-self.focus_mutants_per

        if (low_num==0)&(high_num==0):
            spread = 0
        elif (low_num==0)|(high_num==0)&(self.focus_mutants_per_spread>0):
            raise ValueError('Mean cannot be at one edge of a range if spread is greater than 0')
        else:
            spread = self.focus_mutants_per_spread
        
        if spread>0:
            if low_num>0:
                low_max = spread/low_num
                low_d = low_max/(low_num+1)
                low = [(i)*low_d for i in range(low_num+1)]
            else:
                low=[0]

            if high_num > 0:
                high_max = spread/high_num
                high_d = high_max/(high_num+1)
                high = [(i)*high_d  for i in range(high_num+1)][::-1]
            else:
                high = [0]

            ps =  low+ [1-spread]+ high
            ps = ps[1:-1]
            choices = np.arange(self.focus_mutants_per_min,self.focus_mutants_per_max+1)
        else:
            ps=[1]
            choices = [self.focus_mutants_per]
        focused_dict = {}
        for var_num in range(self.num_focus):
            still_trying=True
            tries = 0
            while still_trying:
                tries = tries+1
                mut_number = self.rng.choice(choices, p=ps)
                mut_positions = self.rng.choice(a=self.focus_positions, size=mut_number, replace=False)
                mut_positions = np.sort(mut_positions)
                wts = (self.wt_protein[pos] for pos in mut_positions)
                aas = [self.rng.choice(a=amino_acids[(amino_acids!=wt)&(amino_acids!='*')]) for wt in wts]

                var_protein = self.wt_protein[:mut_positions[0]]
                name = f'{len(aas)}pt_'
                cur_poses, next_poses = itertools.tee(mut_positions)
                next(next_poses, None)
                for cur_pos, next_pos, aa in zip(cur_poses, next_poses, aas):
                    var_protein = var_protein+aa+self.wt_protein[cur_pos+1:next_pos]
                    name = name+f'{cur_pos+1:02d}{aa}_'
                var_protein = var_protein+aas[-1]+self.wt_protein[mut_positions[-1]+1:]
                name = name+f'{mut_positions[-1]+1:02d}{aas[-1]}'

                if name not in focused_dict.keys():
                    focused_dict[name] = var_protein
                    still_trying=False
                elif tries>=self.max_tries:
                    still_trying=False
                
        # Add to variant dictionary
        self.variant_dict.update(focused_dict)
        return(focused_dict)


    def create_singles(self):
        """
        Create all possible single amino acid mutants.
        
        Generates every possible single amino acid substitution at each position
        in the protein sequence, excluding the wild-type amino acid.
        
        Returns
        -------
        dict
            Dictionary mapping variant names to protein sequences for single
            mutants. Names follow pattern '1pt_PPAA' where PP is the position
            (01-indexed) and AA is the amino acid.
            
        Notes
        -----
        Updates self.variant_dict with the new variants.
        Creates (L * 20) single mutants where L is the protein length.
        The 20 includes the stop codons.
        """
        single_variant_dict = {}
        for i in range(self.L):
            wt_aa = self.wt_protein[i]
            for aa in amino_acids:
                if aa!=wt_aa:
                    var_protein = self.wt_protein[:i] + aa + self.wt_protein[i+1:]
                    name = f"1pt_{i+1:02d}{aa}"
                    single_variant_dict[name] = var_protein

        # Add to variant dictionary
        self.variant_dict.update(single_variant_dict)
        return(single_variant_dict)
    
    def create_doubles(self):
        """
        Create random double amino acid mutants.
        
        Generates random combinations of two amino acid substitutions at different
        positions in the protein sequence. Uses a tracking matrix to avoid
        generating duplicate combinations.
        
        Returns
        -------
        dict
            Dictionary mapping variant names to protein sequences for random
            double mutants. Names follow pattern '2pt_PP1AA1_PP2AA2' where PP
            is the position (01-indexed) and AA is the amino acid.
            
        Notes
        -----
        Updates self.variant_dict with the new variants.
        Uses self.already_seen matrix to track which combinations have been
        generated to avoid duplicates.
        Excludes stop codons (*) from possible substitutions.
        """
        already_seen = np.zeros((self.L, self.L, len(amino_acids), len(amino_acids)))
        double_variant_dict = {}
        for n in range(self.num_doubles):
            pair = self.rng.choice(a=self.L, size=2, replace=False)
            pair.sort()
            i,j = tuple(pair)

            # Chose aa for position i
            wt_aa_i = self.wt_protein[i]
            mask = (amino_acids!=wt_aa_i)&(amino_acids!='*')
            aa_i_idx = self.rng.choice(np.arange(len(amino_acids))[mask],size=1)[0]
            aa_i = amino_acids[aa_i_idx]
            

            # Chose aa for position j
            wt_aa_j = self.wt_protein[j]
            mask = (amino_acids!=wt_aa_j)&(amino_acids!='*')&(already_seen[i,j, aa_i_idx]==0)
            aa_j_idx = self.rng.choice(np.arange(len(amino_acids))[mask],size=1)[0]
            aa_j = amino_acids[aa_j_idx]
            
            already_seen[i,j, aa_i_idx, aa_j_idx] = already_seen[i,j, aa_i_idx, aa_j_idx]+1

            # Create protein sequence
            var_protein = list(self.wt_protein)
            var_protein[i] = aa_i
            var_protein[j] = aa_j
            var_protein = ''.join(var_protein)

            # Create name 
            name = f'2pt_{i+1:02d}{aa_i}_{j+1:02d}{aa_j}'

            double_variant_dict[name] = var_protein
            
        # Add to variant dictionary
        self.variant_dict.update(double_variant_dict)
        self.already_seen = already_seen
        return(double_variant_dict)
    
    
    def generate_nt_seqs(self):
        """
        Generate synonymous nucleotide sequences for all protein variants.
        
        Creates multiple synonymous nucleotide sequences for each protein variant
        using the codon usage table. The number of variants depends on whether
        it's the wild-type (wt_num variants) or a mutant (var_num variants).
        
        Returns
        -------
        dict
            Dictionary mapping variant names to nucleotide sequences.
            Names follow pattern 'PROTEIN_NAME_vN' where N is the variant number.
            
        Notes
        -----
        Updates self.orf_dict with the nucleotide sequences.
        
        For wild-type variants, generates wt_num synonymous sequences.
        For mutant variants, generates var_num synonymous sequences.
        
        The first wild-type variant uses the exact input nucleotide sequence.
        Subsequent variants are generated by randomly selecting one of the two most common codons for each aa.
        
        Attempts up to max_tries times to generate unique sequences before
        giving up on a particular variant.
        """
        protein_names = list(self.variant_dict.keys())
        protein_names.sort()

        self.orf_dict = {}
        
        for i, protein_name in enumerate(protein_names):
            if i%100 == 0:
                print(f'{i}: {protein_name}')

            # Get protein sequence
            protein_seq = self.variant_dict[protein_name]

            # Iterate over synonymous variants
            k=0
            t = 0
            if protein_name == '0pt_WT':
                K = self.wt_num
            else:
                K = self.var_num
            while k < K:

                # Draw coding sequence
                if (protein_name == '0pt_WT')&(k==0):
                    orf_seq = self.wt_nt
                else:
                    orf_seq = protein_to_orf(protein_seq, rng=self.rng)

                # Add if not already in dictionary
                if orf_seq not in self.orf_dict.values():
                    key = protein_name+f'_v{k+1}'
                    self.orf_dict[key] = orf_seq
                    k += 1
                else:
                    if t >= self.max_tries:
                        break
                    else:
                        t = t+1
        return(self.orf_dict)
    
    def save_file(self):
        """
        Save the library to a CSV file and generate oligonucleotide sequences.
        
        Creates a pandas DataFrame with all library variants and saves it to
        a CSV file. The DataFrame includes protein sequences, nucleotide
        sequences, and full oligonucleotide sequences with adapters.
        
        Notes
        -----
        Updates self.orf_df with the final library DataFrame.
        Updates self.file_name with the saved file path.
        
        The CSV file contains columns:
        - name: Variant identifier
        - protein: Amino acid sequence
        - orf: Nucleotide sequence (without adapters)
        - oligo: Full oligonucleotide sequence (with adapters)
        
        File naming convention:
        {name}_orf_library_N={total_variants}_time={timestamp}.csv
        """
        # Save variants to dataframe
        orf_df = pd.DataFrame()
        keys, values = zip(*self.orf_dict.items())
        orf_df['name'] = keys
        orf_df['protein'] = [orf_to_protein(orf) for orf in values]
        orf_df['orf'] = values
        orf_df['oligo'] = [self.upstream_dna + dna + self.downstream_dna for dna in values]
        orf_df.index.name = 'index'
        orf_df.sort_values(by='name', inplace=True)
        self.file_name = f'{self.name}_orf_library_N={len(orf_df)}_time={self.date_time}.csv'
        self.file_name = self.out_path+self.file_name
        orf_df.to_csv(self.file_name)
        print(self.file_name)
        self.orf_df = orf_df
    

if __name__ == '__main__':
    """
    Command-line interface for LibraryMaker.
    
    This script provides a command-line interface for creating mutagenesis
    libraries. It can be used to generate libraries from either nucleotide
    or amino acid sequences.
    
    Examples
    --------
    # Create library from nucleotide sequence
    python library_maker.py -w ATGGCCGAA -s True -d equal -l even_single_and_double_nt_library
    
    # Create library from amino acid sequence  
    python library_maker.py -w MKEL -a -s True -d 100 -l even_single_and_double_aa_library
    
    # Create library with custom parameters
    python library_maker.py -w ATGGCCGAA -s True -d max:10000 -n 5 -t 50 -l single_and_more_doubles_nt_library

    # Create library with focused mutations
    python library_maker.py -w ATGGCCGAA -s True -d max:10000 -f 1 3 5 -c singles -p 3 -e 0.5 -l single_and_more_doubles_with_all_focused_mutations_nt_library
    """
    from argparse import ArgumentParser
    desc = """
    Generate oligonucleotide libraries for protein mutagenesis studies.
    
    This script creates comprehensive mutagenesis libraries including single and
    multiple amino acid mutants, synonymous nucleotide variants, and focused
    mutagenesis. The output is a CSV file containing all variants with their
    protein sequences, nucleotide sequences, and full oligonucleotide sequences
    ready for synthesis.
    
    The library can be generated from either a nucleotide sequence (default)
    or an amino acid sequence (use -a flag). The script automatically handles
    codon usage optimization and generates multiple synonymous variants for
    each protein variant.
    """
    
    parser = ArgumentParser(
                        prog='library_maker',
                        description=desc,
                        epilog='For more information, see the LibraryMaker class documentation.')

    parser.add_argument('-w','--wt', type=str, required=True,
                       help='The nucleotide or amino acid sequence to mutagenize')
    parser.add_argument('-a','--is-aa', action='store_true',
                       help='Indicate that the input sequence is an amino acid sequence (default: nucleotide)')
    parser.add_argument('-s','--single-mutants',type=bool, default=True,
                       help='Include single mutants in the library (default: True)')
    parser.add_argument('-d','--doubles-to-make',default='equal',
                       help='Number of double mutants to generate. Options: "equal" (same as singles), "max:N" (up to N total oligos), float (fraction of singles), or int (exact number)')
    parser.add_argument('-n','--var-num',type=int, default=3,
                       help='Number of synonymous nucleotide variants per mutant (default: 3)')
    parser.add_argument('-t','--wt-num', type=int, default=100,
                       help='Number of synonymous nucleotide variants for wild-type (default: 100)')
    parser.add_argument('-j','--max_tries', type=int, default=100,
                       help='Maximum attempts to generate unique variants before giving up (default: 100)')
    parser.add_argument('-u', '--upstream_dna',type=str, default='ATGGAAAGTGTCCCAGGAGACTACAGC',
                       help='Upstream adapter sequence (default: ATGGAAAGTGTCCCAGGAGACTACAGC)')
    parser.add_argument('-r', '--downstream_dna',type=str, default='CAGACAAGTAACTCAAGACTTAATGGC',
                       help='Downstream adapter sequence (default: CAGACAAGTAACTCAAGACTTAATGGC)')
    parser.add_argument('-o','--out_path',type=str, default='./',
                       help='Output directory for library files (default: current directory)')
    parser.add_argument('-f','--focus_positions',type=int, nargs='*', default=[],
                       help='Amino acid positions to focus on for mutagenesis (e.g., "1 3 5")')
    parser.add_argument('-c','--focus_to_make',type=str, default='singles',
                       help='Type of focused mutations to generate. Options: "singles", "doubles", "singles_doubles", "all"')
    parser.add_argument('-p','--focus_mutants_per',type=int, default=3,
                       help='Number of focused mutations per variant (default: 3)')
    parser.add_argument('-e','--focus_mutants_per_spread',type=float, default=0,
                       help='Spread of focused mutations per variant (default: 0)')
    parser.add_argument('-i','--focus_mutants_per_min',type=int, default=0,
                       help='Minimum number of focused mutations per variant (default: 0)')
    parser.add_argument('-x','--focus_mutants_per_max',type=int, default=0,
                       help='Maximum number of focused mutations per variant (default: 0)')
    parser.add_argument('-v','--seed',type=int, default=None,
                       help='Random seed for reproducibility (default: None)')
    parser.add_argument('-g','--fig_width',type=int, default=8,
                       help='Width of output figures (default: 8)')
    parser.add_argument('-k','--context',type=str, default='talk',
                       help='Seaborn context for output figures (default: talk)')
    parser.add_argument('-l','--name',type=str, default='Mutant_lib',
                       help='Name of the library (default: Mutant_lib)')
    
    args = parser.parse_args()
    if args.is_aa:
        wt_nt = protein_to_orf(args.wt)
    else:
        wt_nt = args.wt
        
    LibraryMaker(wt_nt, 
                 single_mutants = args.single_mutants, 
                 doubles_to_make = args.doubles_to_make, 
                 name=args.name,
                 var_num = args.var_num,
                 wt_num = args.wt_num, 
                 max_tries = args.max_tries,
                 focus_positions = args.focus_positions, 
                 focus_to_make=args.focus_to_make, 
                 focus_mutants_per=args.focus_mutants_per, 
                 focus_mutants_per_spread=args.focus_mutants_per_spread,
                 focus_mutants_per_min=args.focus_mutants_per_min,
                 focus_mutants_per_max=args.focus_mutants_per_max,
                 upstream_dna = args.upstream_dna,
                 downstream_dna = args.downstream_dna,
                 out_path = args.out_path, 
                 seed=args.seed,
                 fig_width=args.fig_width, 
                 context=args.context)
    
    
        