import logging
import fnmatch
import os, sys
import numpy as np
import glob
import pysam
import pandas
import click


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
        if l.mapq < 20: continue # ignore reads (default = 0.9 correct mapping)
        if l.rlen < 50: continue # ignore extremely short reads
        print 'rlen: {}'.format(l.rlen)
        print 'mapping_quality: {}'.format(l.mapping_quality)
        print 'mapq: {}'.format(l.mapq)
        print 'reference_start: {}'.format(l.reference_start)
        print 'reference_end: {}'.format(l.reference_end)
        print

        #for elem in dir(l):
        #    if '__' in elem: continue
        #    print '{} : {}'.format(elem, getattr(l, elem))

        #print samfile.header
        #sys.exit()

        coverage[samfile.getrname(l.tid)]["nb_reads"]+=1

        #begin = l.reference_start+1
        begin = l.reference_start
        end = l.reference_end # begin + samfile.header['SQ'][0]['LN']
        print 'begin: {}'.format(begin)
        print 'end: {}'.format(end)
        #end = l.reference_start+l.reference_length
	coverage[samfile.getrname(l.reference_id)]["positions"][begin:end] += 1
	coverage[samfile.getrname(l.tid)]["nb_bp"] += l.rlen

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
        if vector['nb_bp'] == 0: continue
        temp = {}
        for i in contigs_size:
            print 'contig size: {}'.format(i)
            print '{} == {}'.format(contig, i['Seq'])
            if contig == i["Seq"]:
                temp["length"] = i["Length"]
                temp["ratio_covered"] = np.sum(vector["positions"])/float(len(vector["positions"]))
                temp["number_reads"] = vector["nb_reads"]
                temp["number_bp"] = vector["nb_bp"]
                temp["coverage"] = vector["nb_bp"]/temp["length"]
                if vector["nb_reads"] > 0 :
                    coverage_prop[contig] = temp
        print 'coverage_prop: {}'.format(coverage_prop)

        output=pandas.DataFrame(coverage_prop).transpose()
        #output['name'] = output.index
        #output = output[['name', 'length', 'length_covered', 'coverage', 'number_bp', 'number_reads']]
        output=output.sort_values(['number_bp','ratio_covered'],ascending=[0,0])
        output.to_csv(outputfile, index=False)
    samfile.close()


@click.command()
@click.option('-b', '--bam')
@click.option('-o', '--output')
@click.option('-v', '--verbose')
def gen_coverage_file(bam, output, verbose):
    output_fp = os.path.join(output, bam.split('/')[-1].replace('.bam', '_coverage.csv'))
    process(bam, output_fp)
    return output_fp

if __name__ == "__main__":
    gen_coverage_file()
