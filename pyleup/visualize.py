import os
import sys


import pandas as pd
import numpy as np

from bokeh.layouts import row, column
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer
from bokeh.plotting import figure, show, curdoc
from bokeh.io import export_png



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
    
        # create the scatter plot
        p = figure(plot_width=600, plot_height=600, min_border=10, min_border_left=50,
                    x_axis_location=None, y_axis_location=None, toolbar_location=None,
                    title='Experiment')
        r = p.scatter(x, y)

        # create the horizontal plot
        hhist, hedges = np.histogram(x, bins=20)
        hzeros = np.zeros(len(hedges)-1)
        hmax = max(hhist)*1.1

        ph = figure(toolbar_location=None, plot_width=p.plot_width, plot_height=200, x_range=p.x_range,
                    y_range=(-hmax, hmax), min_border=10, min_border_left=50, y_axis_location='right')



        # create the vertical plot
        vhist, vedges = np.histogram(y, bins=20)
        vzeros = np.zeros(len(vedges)-1)
        vmax = max(vhist)*1.1
        pv = figure(toolbar_location=None, plot_width=200, plot_height=p.plot_height, x_range=(-vmax, vmax),
                                                 y_range=p.y_range, min_border=10, y_axis_location="right")

        layout = column(row(p, Spacer(width=200, height=200)), row(ph, Spacer(width=200, height=200)))

        curdoc().add_root(layout)
        curdoc().title = 'Coverage Plot'

        #def update(attr, old, new):
        #    inds = np.array(new['1d']['indices'])
        #    if len(inds) == 0 or len(inds) == len(x):
        #        hhist1, hhist2 = hzeros, hzeros
        #        vhist1, vhist2 = vzeros, vzeros
        #    else:
        #        neg_inds = np.ones_like(x, dtype=np.bool)
        #        neg_inds[inds] = False
        #        hhist1, _ = np.histogram(x[inds], bins=hedges)
        #        vhist1, _ = np.histogram(y[inds], bins=vedges)
        #        hhist2, _ = np.histogram(x[neg_inds], bins=hedges)
        #        vhist2, _ = np.histogram(y[neg_inds], bins=vedges)

        #    hh1.data_source.data["top"]   =  hhist1
        #    hh2.data_source.data["top"]   = -hhist2
        #    vh1.data_source.data["right"] =  vhist1
        #    vh2.data_source.data["right"] = -vhist2

        #r.data_source.on_change('selected', update)

        export_png(curdoc(), filename=os.path.join(output_dp, '{}.png'.format(ref).replace(' ', '_')))
        #ref_df = pd.DataFrame({'reads':x, 'mapq':y, 'reference':r})
        #sns_plot = sns.jointplot(x='reads', y='mapq', data=ref_df, join_kws={'marker':'x'})
        #sns_plot.savefig(os.path.join(output_dp, '{}.png'.format(ref).replace(' ','_')))

    
    
