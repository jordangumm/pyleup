import logging
import fnmatch
import os, sys
import numpy as np
import glob
import pysam
import pandas
import click


def get_ref_lens(samfile):
    temp = []
    for i, dic in enumerate(samfile.header['SQ']):
        temp.append({'Seq':dic['SN'], 'Length':dic['LN']})
    return temp


def coverage_vectors(contigs_size):
    coverage = {}
    for i in contigs_size:
        temp = {}
        temp["positions"] = np.zeros(i["Length"])
        temp["nb_reads"] = 0
	temp["nb_bp"] = 0
	coverage[i["Seq"]] = temp
    return coverage


def process_bam(bam, output_dp):
    """ Create coverage and read data from bam file

    Returns file paths to generated coverage and read CSVs
    """
    bam_fn = os.path.basename(bam)
    coverage_fp = os.path.join(output_dp, bam_fn.replace('.bam', '_coverage.csv'))
    reads_fp = os.path.join(output_dp, bam_fn.replace('.bam', '_reads.csv'))

    samfile = pysam.AlignmentFile(bam, "rb")
    contigs_size = get_ref_lens(samfile)
    coverage = coverage_vectors(contigs_size)

    read_output = open(reads_fp, 'w+')
    read_output.write('read_length,mapq,start,end,reference')
    for l in samfile.fetch():
        if l.mapq < 10: continue
        if l.rlen < 50: continue
        read_output.write('\n{},{},{},{},{}'.format(l.rlen, l.mapq,
                    l.reference_start, l.reference_end, samfile.getrname(l.reference_id).split(',')[0]))
        coverage[samfile.getrname(l.tid)]["nb_reads"] += 1
        coverage[samfile.getrname(l.reference_id)]["positions"][l.reference_start:l.reference_end] = 1
        coverage[samfile.getrname(l.tid)]["nb_bp"] += l.rlen
    read_output.close()

    coverage_prop = {}
    for contig,vector in coverage.items():
        if vector['nb_bp'] == 0: # no reads, so output blank file
            output = pandas.DataFrame()
            output.to_csv(coverage_fp, index=False)
            continue
        temp = {}
        for i in contigs_size:
            if contig == i["Seq"]:
                temp["length"] = i["Length"]
                temp["ratio_covered"] = np.sum(vector["positions"])/float(len(vector["positions"]))
                temp["number_reads"] = vector["nb_reads"]
                temp["number_bp"] = vector["nb_bp"]
                if vector["nb_reads"] > 0 :
                    coverage_prop[contig] = temp

        output = pandas.DataFrame(coverage_prop).transpose()
        output = output.sort_values(['number_bp','ratio_covered'],ascending=[0,0])
        output.to_csv(coverage_fp, index=False)
    samfile.close()
    return coverage_fp, reads_fp

