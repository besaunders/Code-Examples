#!/usr/bin/perl -w

use strict;

my $taxfile = shift @ARGV;
my %bltarget;
open TF, $taxfile or die "Cannot open taxonomy specification file\n";
while (<TF>) {
	chomp;
	my ($tax,$div,$subdiv) = split /\t/;
	if ($subdiv eq 'Arena') {
		$bltarget{$div}++;
	}
}
close TF;

foreach my $bldiv (sort {$a cmp $b} keys %bltarget) {
	my $aafile = 'Arena_' . $bldiv . '.aa';
	my $bldb = $bldiv . '_nr';
	my $outfile = $bldb . '.blnr';
	my $t0 = time;
	my $blcmd = "blastp -query $aafile -db $bldb -max_target_seqs 20 " .
		"-outfmt '6 qseqid pident nident qlen qstart qend slen sstart " .
		"send gapopen gaps evalue score sallseqid staxids sscinames " .
		"salltitles' -out $outfile";
	`$blcmd`;
	my $eltime = time - $t0;
	print STDERR "$bldiv: $eltime seconds\n";
}

