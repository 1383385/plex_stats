from Bio.Seq import IUPAC
import os
import click
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from functools import partial
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd
import numpy as np
import itertools as it



# column names
PRIMER_SEQ = "Seq"
PRIMER_NAME = "Primer_Name"
PRIMER_LEN = "Primer_Length"
PRIMER_TM = "Tm"
PRIMER_DELTA_TM = "Delta_Tm"
PRIMER_GC = "GC_Percent"
PRIMER_LEN_HOMO = "Length_Longest_Homopolymer"
PRIMER_PERCENT_HOMO = "Percent_Homopolymer"
PRIMER_AMPLICON = "Amplicon"


def get_primer_stats(primer_fasta, outpath, prefix, tm_params):
    primer_table = []
    with open(primer_fasta, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            ambiguous_primers = expand_ambiguous_dna(record.seq)
            if len(ambiguous_primers) == 1:
                primer_table.append([record.id, record.seq])
            else:
                for i, primer in enumerate(ambiguous_primers):
                    primer_table.append(["{0}_{1}".format(
                        record.id, i+1), primer])
    primer_table = pd.DataFrame(primer_table, columns=[PRIMER_NAME, PRIMER_SEQ])
    primer_table = run_primer_stats(primer_table, tm_params)
    primer_table.to_csv(os.path.join(
        outpath, "{}.csv".format(prefix)), index=False)
    style_table = style_df(primer_table)

    write_to_html_file(
        style_table, title=prefix.upper(),
        filename=os.path.join(outpath, "{}.html".format(prefix)))
    
    style_table.to_excel(
        os.path.join(outpath, "{}.xlsx".format(prefix)),
        engine='xlsxwriter', index=False)


def expand_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = IUPAC.IUPACData.ambiguous_dna_values
    return tuple(map("".join, it.product(*map(d.get, seq))))

def highlight_cutoff(data, max_, bad_color="#D93A46", good_color="#C8D9E6"):
    """
    Color values based on whether they are greater than cutoff
    """
    is_bad = data > max_
    return [
        'background-color: {}'.format(bad_color) if v else
        'background-color: {}'.format(good_color) for v in is_bad]

def highlight_range(data, max_, min_, bad_color="#D93A46", good_color="#C8D9E6"):
    """
    Color values based on whether they are in or out of range
    """
    is_bad = np.logical_or(np.greater(data, max_), np.less(data, min_))
    return [
        'background-color: {}'.format(bad_color) if v else
        'background-color: {}'.format(good_color) for v in is_bad]

def highlight_pass(data, max_, good_color="#C8D9E6"):
    is_good = data < max_
    return ['background-color: {}; color: black'.format(good_color)
    if v else "" for v in is_good]


# def highlight_out_of_range(data, min_, max_, 
#                            bad_color="#D93A46", good_color="#C8D9E6"):
#     '''
#     Highlight values that are outside of a range
#     '''
#     attr_bad = 'background-color: {}'.format(bad_color)
#     attr_good = 'background-color: {}'.format(good_color)
#     if data.ndim == 1:  # Series from .apply(axis=0) or axis=1
#         out_of_range = np.logical_or(data > max_, data < min_)
#         return [attr_bad if v else attr_good for v in out_of_range]


# def highlight_cutoff(data, max_, bad_color="#D93A46", good_color="#C8D9E6"):
#     '''
#     Highlight values that are greater than cutoff
#     '''
#     attr_bad = 'background-color: {}'.format(bad_color)
#     attr_good = 'background-color: {}'.format(good_color)
#     if data.ndim == 1:  # Series from .apply(axis=0) or axis=1
#         out_of_range = data > max_
#         return [attr_bad if v else attr_good for v in out_of_range]


def parse_hybrid_columns(cols):
    stats_cols = []
    matrix_cols = []
    for c in cols:
        if any(c.startswith(stat) for stat in ["Mean", "Median", "Max"]):
            stats_cols.append(c)
        else:
            matrix_cols.append(c)
    return stats_cols, matrix_cols


def style_df(df):
    stat_cols, matrix_cols = parse_hybrid_columns(list(df.columns.values)[8:])
    hybrid_max_run_columns = [
        c for c in matrix_cols if c.endswith("Hybrid_Run")]
    hybrid_score_columns = [
        c for c in matrix_cols if c.endswith("Hybrid_Score")]
    hybrid_score_3prime_columns = [
        c for c in matrix_cols if c.endswith("3'_Score")]
    hybrid_percent_columns = [
        c for c in matrix_cols if c.endswith("Percent_Hybrid")]
    hybrid_mean = [c for c in stat_cols if c.startswith("Mean")]
    hybrid_median = [c for c in stat_cols if c.startswith("Median")]
    hybrid_max = [c for c in stat_cols if c.startswith("Max")]

    df = df[list(df.columns.values)[:8] + 
        hybrid_mean + hybrid_median + hybrid_max + 
        hybrid_score_columns + hybrid_score_3prime_columns +
        hybrid_max_run_columns + hybrid_percent_columns]
    gc_cutoff = partial(highlight_range, min_=0.4, max_=0.6)
    # df.style.apply(gc_cutoff, subset=[PRIMER_GC])
    dtm_cutoff = partial(highlight_cutoff, max_ = 5)
    # df.style.apply(dtm_cutoff, subset=[PRIMER_DELTA_TM])
    homopolymer_cutoff = partial(highlight_cutoff, max_=5)
    # df.style.apply(homopolymer_cutoff, subset=[PRIMER_LEN_HOMO])
    homopolymer_percent_cutoff = partial(highlight_cutoff, max_=.35)
    hybrid_run_pass = partial(highlight_pass, max_=6)
    hybrid_score_pass = partial(highlight_pass, max_ = 25)
    hybrid_percent_pass = partial(highlight_pass, max_ = 0.5)
    # df.style.apply(homopolymer_percent_cutoff, subset=[PRIMER_PERCENT_HOMO])
    # cm = sns.diverging_palette(240, 10, as_cmap=True)
    cm = sns.light_palette("red", as_cmap=True)

    # df.style.format(
    #     {c: '{:.1f}' for c in df.columns[2:] if df.dtypes[c] == float}).format(
    #     {c: '{:.1%}' for c in df.columns[2:] if "percent" in c.lower()}).apply(
    #         gc_cutoff, subset=[PRIMER_GC]).apply(
    #             dtm_cutoff, subset=[PRIMER_DELTA_TM]).apply(
    #                 homopolymer_cutoff, subset=[PRIMER_LEN_HOMO]).apply(
    #                     homopolymer_percent_cutoff,
    #                     subset=[PRIMER_PERCENT_HOMO]).background_gradient(
    #     subset=stat_cols + matrix_cols, cmap=cm).hide_index()
    s = df.style.format(
        {c: '{:.1f}' for c in df.columns[2:] if df.dtypes[c] == float}).format(
        {c: '{:.1%}' for c in df.columns[2:]
        if "percent" in c.lower()}).hide_index().apply(
        homopolymer_cutoff, subset=[PRIMER_LEN_HOMO]).apply(
        gc_cutoff, subset=[PRIMER_GC]).apply(
        dtm_cutoff, subset=[PRIMER_DELTA_TM]).apply(
                            homopolymer_percent_cutoff,
        subset=[PRIMER_PERCENT_HOMO]).background_gradient(
            subset=stat_cols, cmap=cm).background_gradient(
            subset=hybrid_max_run_columns,
            cmap=cm, axis=None).background_gradient(
            subset=hybrid_score_columns,
            cmap=cm, axis=None).background_gradient(
            subset=hybrid_score_3prime_columns,
            cmap=cm, axis=None).background_gradient(
                subset=hybrid_percent_columns,
                cmap=cm, axis=None).apply(
            hybrid_run_pass, subset=hybrid_max_run_columns).apply(
            hybrid_score_pass, subset=hybrid_score_columns).apply(
            hybrid_score_pass, subset=hybrid_score_3prime_columns).apply(
            hybrid_percent_pass, subset=hybrid_percent_columns)
    return s


def write_to_html_file(df, title='', filename='out.html'):
    '''
    Write an entire dataframe to an HTML file with nice formatting.
    '''

    result = '''
        <html>
        <head>
        <!-- Compiled and minified CSS -->
        <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/css/materialize.min.css">
        <style>

            h2 {
                text-align: center;
                font-family: Helvetica, Arial, sans-serif;
            }
            table { 
                margin-left: auto;
                margin-right: auto;
            }
            table, th, td {
                border: 1px solid black;
                border-collapse: collapse;
            }
            th, td {
                padding: 5px;
                text-align: center;
                font-family: Helvetica, Arial, sans-serif;
                font-size: 90%;
                background: white
            }

            table tbody tr:hover td{
                filter: brightness(85%);
            }
            .wide {
                width: 90%; 
            }

        </style>
        </head>
        <body>

        '''
    result += '<h2> %s </h2>\n' % title
    if type(df) == pd.io.formats.style.Styler:
        result += df.render()
    else:
        result += df.to_html(classes='wide', escape=False)
    result += '''

    </body>
    </html>
    '''
    with open(filename, 'w') as f:
        f.write(result)





def get_hybrid_matrix_stats(hybrid_matrix):
    stats = pd.DataFrame([])
    stats['Mean'] = hybrid_matrix.mean(axis=1)
    stats['Max'] = hybrid_matrix.max(axis=1)
    stats['Median'] = hybrid_matrix.median(axis=1)
    return hybrid_matrix.join(stats)


def run_primer_stats(df, tm_params):
    df = format_seq(
        get_homopolymer(
            get_gc(
                get_delta_tm(
                    get_tm(
                        get_length(
                            df), tm_params)))))

    score_matrix, score_matrix_3_prime, run_matrix, percent_matrix = \
        get_cross_hybrid(df, tm_params)
    score_matrix = get_hybrid_matrix_stats(score_matrix)
    score_matrix_3_prime = get_hybrid_matrix_stats(score_matrix_3_prime)
    run_matrix = get_hybrid_matrix_stats(run_matrix)
    percent_matrix = get_hybrid_matrix_stats(percent_matrix)
    run_score_matrix = score_matrix.join(
        run_matrix, rsuffix="_Hybrid_Run", lsuffix="_Hybrid_Score")
    percent_3prime_matrix = score_matrix_3_prime.join(
        percent_matrix, rsuffix="_Percent_Hybrid", lsuffix="_3'_Score")
    cross_hybrid = run_score_matrix.join(percent_3prime_matrix)     
    
    df = df.join(cross_hybrid, on=PRIMER_NAME)
    return df


def get_empty_cross_hybrid_matrix(index):
    return pd.DataFrame(
        np.zeros((len(index), len(index))),
        columns=index,
        index=index)

def get_cross_hybrid(df, tm_params):
    score_matrix = get_empty_cross_hybrid_matrix(df[PRIMER_NAME].values)
    score_matrix_3_prime = get_empty_cross_hybrid_matrix(df[PRIMER_NAME].values)
    run_matrix = get_empty_cross_hybrid_matrix(df[PRIMER_NAME].values)
    percent_matrix = get_empty_cross_hybrid_matrix(df[PRIMER_NAME].values)
    name_to_seq_map = {
        name: seq for name, seq in zip(df[PRIMER_NAME], df[PRIMER_SEQ])}

    for p1, p2 in it.combinations_with_replacement(df[PRIMER_NAME].values, 2):
        max_score, max_score_3_prime, max_complement, max_run = \
             calculate_cross_hybrid_score(
            name_to_seq_map[p1], name_to_seq_map[p2], tm_params)
        score_matrix.loc[p1, p2] = max_score
        score_matrix.loc[p2, p1] = max_score
        score_matrix_3_prime.loc[p1, p2] = max_score_3_prime
        score_matrix_3_prime.loc[p2, p1] = max_score_3_prime
        run_matrix.loc[p1, p2] = max_run
        run_matrix.loc[p2, p1] = max_run
        percent_matrix.loc[p1, p2] = max_complement
        percent_matrix.loc[p2, p1] = max_complement
    return score_matrix, score_matrix_3_prime, run_matrix, percent_matrix
    




def cross_hybrid_score_worker(primer_1, primer_2, tm_params):
    score_dict = {('G', 'C'): 4, ('A', 'T'): 2, ('C', 'G'): 4, ('C', 'A'): -0.6,
                  ('T',): -0.6, ('A',): -0.6, ('G',): -0.6, ('C',): -0.6,
                  ('A', 'C'): -0.6, ('C', 'T'): -0.6, ('A', 'G'): -0.4, ('G', 'A'): -0.4,
                  ('G', 'T'): -0.4, ('T', 'G'): -0.4, ('T', 'A'): 2, ('T', 'C'): -0.6}
    comp = [('C', 'G'), ('A', 'T'), ('T', 'A'), ('G', 'C')]
    weighted_comp = {('C', 'G'): 4, ('A', 'T'): 2, ('T', 'A'): 2, ('G', 'C'): 4}

    matches = list(map(lambda x: 1 if tuple(set(x)) in comp else 0,
                       zip(primer_1, primer_2)))
    # Calculate the number of bases that align
    max_comp = sum(matches)
    # weighted_matches = list(
    #     map(lambda x: weighted_comp.get(tuple(set(x)), 0),
    #     zip(primer_1, primer_2) ))
    # max_comp = (2 * sum(weighted_matches)) / \
    #     (4* float(len(primer_1) + len(primer_2))) * 100
    
    longest_run = 0
    complementary = [list(comp) for run, comp in it.groupby(matches)]
    for comp in complementary:
        if 1 in comp and len(comp) > longest_run:
            longest_run = len(comp)
    try:
        start = matches.index(1)
        stop = matches[::-1].index(1)
    except ValueError:
        # No hybridization
        return 0, 0, 0
    p1 = primer_1[start: -stop]
    p2 = primer_2[start: -stop]
    try: 
        score = mt.Tm_NN(primer_1, c_seq=primer_2, **tm_params)
    except ValueError:
        # super rough approximation when thermodynamic data is not available
        score = sum(map(lambda x: score_dict[tuple(set(x))], zip(p1, p2)))
    return score, max_comp, longest_run


def calculate_cross_hybrid_score(primer_1, primer_2, tm_params):
    p1 = primer_1 + ("*" * (len(primer_2) - 1))
    p2 = "*" * (len(primer_1) - 1) + primer_2[::-1]
    max_score = 0
    max_complement = 0
    max_run = 0
    max_score_3_prime = 0
    while len(p1):
        score, max_comp, longest_run = cross_hybrid_score_worker(
            p1, p2, tm_params)
        max_comp = (2 * max_comp) / (len(primer_1) + len(primer_2))
        if len(p1) > min(len(primer_1), len(primer_2)):
            if score > max_score_3_prime:
                max_score_3_prime = score
        if score > max_score:
            max_score = score

        if max_comp > max_complement:
            max_complement = max_comp
        if longest_run > max_run:
            max_run = longest_run
        p1 = p1[:-1]
        p2 = p2[1:]
    return max_score, max_score_3_prime, max_complement, max_run


def format_seq(df):
    df[PRIMER_SEQ] = df[PRIMER_SEQ].astype(str)
    return df

def get_amplicon(df):
    df[PRIMER_AMPLICON] = df[PRIMER_NAME].apply(
        lambda x: "_".join(x.split("_")[:-1]))
    return df

def get_length(df):
    df[PRIMER_LEN] = df[PRIMER_SEQ].apply(lambda x: len(x))
    return df

def get_tm(df, tm_params):
    df[PRIMER_TM] = df[PRIMER_SEQ].apply(lambda x: calculate_tm(x, **tm_params))
    return df

def get_delta_tm(df):
    avg_tm = df[PRIMER_TM].mean()
    df[PRIMER_DELTA_TM] = np.abs(df[PRIMER_TM] - avg_tm)
    return df

def get_gc(df):
    df[PRIMER_GC] = df[PRIMER_SEQ].apply(calculate_gc)
    return df

def get_homopolymer(df):
    df[PRIMER_LEN_HOMO], df[PRIMER_PERCENT_HOMO] = list(zip(*df[PRIMER_SEQ].apply(
        calculate_homopolymer)))
    return df

def calculate_tm(seq, Na=50, K=0, Tris=0, Mg=0, dNTPs=0):
    seq = str(seq)
    return mt.Tm_NN(seq, Na=Na, K=K,
                        Tris=Tris, Mg=Mg, dNTPs=dNTPs)
                    
def calculate_gc(seq):
    return ((seq.count("C") + seq.count('G')) / float(len(seq)))


def calculate_homopolymer(seq):
    homopolymers = [
        list(hp) for run, hp
        in it.groupby(seq)]
    n_homopolymer = 0
    for hp in homopolymers:
        if len(hp) >= 3:
            n_homopolymer += len(hp)
    longest_homopolymer = len(max(homopolymers, key=len))
    percent_homopolymer = n_homopolymer/float(len(seq))
    return longest_homopolymer, percent_homopolymer
