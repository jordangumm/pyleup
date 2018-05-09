import logging
import fnmatch
import os
import numpy as np
import glob
import pysam
import pandas
import click

import concurrent.futures

def get_seq_len_from_bam(samfile):
    temp = []
    for i, dic in enumerate(samfile.header['SQ']):
        print 'dic: {}'.format(dic)
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

def parse(samfile,coverage):
    """ Completes read coverage depths over reference for sample """
    for l in samfile.fetch():
        coverage[samfile.getrname(l.tid)]["nb_reads"]+=1

        begin = l.reference_start+1
        end = begin + samfile.header['SQ'][0]['LN']
        print 'begin: {}'.format(begin)
        print 'end: {}'.format(end)
        #end = l.reference_start+l.reference_length
	coverage[samfile.getrname(l.reference_id)]["positions"][begin:end] = 1
	coverage[samfile.getrname(l.tid)]["nb_bp"] += abs(end-begin)

	return coverage


def process(inputfile,outputfile):
    samfile = pysam.AlignmentFile(inputfile, "rb")
    print 'samfile: {}'.format(samfile)
    contigs_size=get_seq_len_from_bam(samfile)
    coverage = coverage_vectors(contigs_size)
    coverage = parse(samfile,coverage)
    coverage_prop = {}
    for contig,vector in coverage.items():
        print 'contig: {}'.format(contig)
        print 'vector: {}'.format(vector)
        temp = {}
        for i in contigs_size:
            print 'contig size: {}'.format(i)
            print '{} == {}'.format(contig, i['Seq'])
            if contig == i["Seq"]:
                temp["length"] = i["Length"]
                temp["length_covered"] = np.sum(vector["positions"])/float(len(vector["positions"]))*100
                temp["number_reads"] = vector["nb_reads"]
                temp["number_bp"] = vector["nb_bp"]
                temp["coverage"] = vector["nb_bp"]/temp["length"]
                if vector["nb_reads"] > 0 :
                    coverage_prop[contig] = temp
        print 'coverage_prop: {}'.format(coverage_prop)

        output=pandas.DataFrame(coverage_prop).transpose()
        #output['name'] = output.index
        #output = output[['name', 'length', 'length_covered', 'coverage', 'number_bp', 'number_reads']]
        output=output.sort_values(['number_bp','length_covered'],ascending=[0,0])
        output.to_csv(outputfile, index=False)
    samfile.close()


@click.command()
@click.option('-b', '--bam-dp')
@click.option('-o', '--output-dp')
@click.option('-v', '--verbose')
def generate_coverage_file(bam_dp, output_dp, verbose):
    bam_list = glob.glob(bam_dp+'*_sorted.bam')
    print bam_list
    for bam_fp in bam_list:
        output_fp = bam_fp.replace('.bam', '_coverage.csv')

        print "we will process file {0} -> will go to {1}".format(bam_fp, output_fp)
        arguments = (bam_fp, output_fp)
        with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
            result = executor.submit(process, *arguments).result()

@click.command()
@click.option('-b', '--bam')
@click.option('-o', '--output')
@click.option('-v', '--verbose')
def generate_single_coverage_file(bam, output, verbose):
    output_fp = os.path.join(output, bam.split('/')[-1].replace('.bam', '_coverage.csv'))
    process(bam, output_fp)

if __name__ == "__main__":
    generate_single_coverage_file()
