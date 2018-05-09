#!/usr/bin/perl

use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h='';
my $dir='';
my $genome_asked='';
my $out_file='';
my $cmd='';
my $out='';
GetOptions ('help' => \$h, 'h' => \$h, 'b=s'=>\$dir, 'g=s'=>\$genome_asked, 'o=s'=>\$out_file);
if ($h==1 || $dir eq "" || $out_file eq ""){ # If asked for help or did not set up any argument
    print "# Script to cluster new contigs

# Arguments:

-b : bam file directory
-o : out file
-g : (facultative) Genome Id to generate the RP for (leave blank for all genomes available)\n";
    die "\n";
}

my $default_th_TOV=0; # Th on the number of reads mapped for a sample to get out of the "Other TOV" category
my %th_reads_TOV=("PSAHP1_revcom_reperm_ORF00031"=>50); #  We used 50 by default, 500 for HP1

# my $palette_file="Color_palette.tsv";
my $code_other="Other_TOV90";

my $k=0;
my %counts;
my %check_seq;
my %store_col;
my @list=<$dir*.bam>;
open my $s1,">",$out_file;
print $s1 "Genome,Start,Stop,Id,Sample\n";
foreach my $bam_file (@list){
    $bam_file=~/(\w+)-bam\/(\w+)-(\w+)_sorted.bam/; # blastn 2.6.0 outfmt 6
    print "$bam_file\n";
    my $sample_name = (split("/", $bam_file))[-1];
    print "bamfile: $bam_file\tsample_name: $sample_name\n";
    my %store_hit=();
    my %count=();
    open my $bam,"samtools view $bam_file |";
    while(<$bam>){
        chomp($_);
        my @tab=split("\t",$_);
        if (($tab[2] eq $genome_asked) || ($genome_asked eq "")){
            print "tab: @tab\n";
            my $genome=$tab[2];
            print "genome: $genome\n";
            my $start=$tab[6];
            print "start: $start\n";
            my $length=$tab[3];
            print "length: $length\n";
            my $stop=$start+$length;
            print "stop: $stop\n";
            my $mismatch = $tab[4];
            print "mismatch: $mismatch\n";
            my $match = $length-$mismatch;
            print "\n";
            exit 0;

            my $ali_length = $match+$mismatch;
            my $id = $tab[0];  #int($match/($match+$mismatch)*10000)/100;
            $count{$genome}++;
            $store_hit{$genome}{"Hit_".$count{$genome}}="\"$genome\",$start,$stop,$id";
        }
    }
    close $bam;
    my $tag=0;
    foreach my $genome (keys %count){
        my $th=$default_th_TOV;
        my $sample_name_g=$sample_name;
        if (defined($th_reads_TOV{$genome})){$th=$th_reads_TOV{$genome}}
        if ($count{$genome}<$th){
            $sample_name_g=$code_other;
        }
        else{
            $tag=1;
        }
        if ($count{$genome}>0){
            foreach(keys %{$store_hit{$genome}}){
                print $s1 "$store_hit{$genome}{$_},$sample_name_g\n";
            }
        }
    }
        
}
close $s1;
