#!/usr/bin/perl -w

use strict;

my $taxfile = shift @ARGV;
my (%blinput,%bltarget);
open TF, $taxfile or die "Cannot open taxonomy specification file\n";
while (<TF>) {
	chomp;
	my ($tax,$div,$subdiv) = split /\t/;
	if ($subdiv eq 'Arena') {
		$blinput{$div}++;
	}
	else {
		$bltarget{$div}++;
	}
}
close TF;

my $divspec = shift @ARGV;

if (! defined $divspec) {
	$divspec = 'undefined';
}
elsif (! exists $blinput{$divspec}) {
	print STDERR "Warning: specified division $divspec undefined\n\n";
	$divspec = 'undefined';
}
else {
	print STDERR "Using $divspec Arena proteins only\n\n";
}

foreach my $blip (sort {$a cmp $b} keys %blinput) {
	next if ($blip ne $divspec);
	foreach my $bldiv (sort {$a cmp $b} keys %bltarget) {
		my $ipfile = 'Arena_' . $blip . '.aa';
		my $bldb = $bldiv . '_nr';
		my $outfile = $blip . '_' . $bldiv . '.blnr';
		my $t0 = time;
		my $blcmd = "blastp -query $ipfile -db $bldb -max_target_seqs 20 " .
			"-outfmt '6 qseqid pident nident qlen qstart qend slen sstart " .
			"send gapopen gaps evalue score sallseqid staxids sscinames " .
			"salltitles' -out $outfile";
		`$blcmd`;
		my $eltime = time - $t0;
		print STDERR "$outfile: $eltime seconds\n";
	}
}

