import os
import sys


import pandas as pd
import numpy as np


def gen_scatter(coverage_fp, reads_fp, output_dp):
    """ Generate scatter plot of preprocessed read positions

    Arguments
    reads_fp -- String file path to reads csv file created in coverage.py process
    """
    coverage_df = pd.read_csv(coverage_fp)
    for key in ['length', 'number_bp', 'number_reads', 'ratio_covered']:
        if key not in coverage_df.keys():
            sys.exit('Missing column(s):\n[ERROR]: {} not in {}'.format(key, coverage_fp))

    reads_df = pd.read_csv(reads_fp)
    for key in ['read_length', 'mapq', 'start', 'end', 'reference']:
        if key not in reads_df.keys():
            sys.exit('Missing column(s):\n[ERROR]: {} not in {}'.format(key, reads_fp))

    for ref in reads_df['reference'].unique():
        ref_reads = reads_df[reads_df['reference'] == ref]
        x, y, r = [], [], []
        for i, read in ref_reads.iterrows():
            for pos in xrange(read['start'], read['end']):
                x.append(pos)
                y.append(read['mapq'])
                r.append(ref)

        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns; sns.set(style='white', color_codes=True)

        df = pd.DataFrame({'Mapped Reads': x, 'Mapping Quality': y})
        g = sns.jointplot(x='Mapped Reads', y='Mapping Quality', data=df, stat_func=None, s=5,
                          joint_kws={'marker': 's'})
        g.fig.suptitle('{} Mapped Reads by Quality'.format(ref))
        g.savefig(os.path.join(output_dp, '{}.png'.format(ref).replace(' ', '_')))
