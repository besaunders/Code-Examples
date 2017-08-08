#!/usr/bin/perl -w

use strict;
use IO::File;

my $stime = time;
my $taxlin = 'C:\Users\bsaunders\db\nodes.dmp';
my %deftax;
my $taxfile = shift @ARGV;
my $blastdb = shift @ARGV;
unless (defined $blastdb) {
	print STDERR "Blast DB set to protein nr\n";
	$blastdb = 'nr';
}
die "Invalid taxonomy file\n" unless (-f $taxfile);

open TF, $taxfile;
while (<TF>) {
	chomp;
	my ($tax,$div,$subdiv) = split /\t/;
	$deftax{$tax} = $div;
}

my %lineage;
open LIN, $taxlin;
while (<LIN>) {
	my ($taxid,$sep,$parent) = split /\t/;
	push @{$lineage{$parent}}, $taxid;
}
close LIN;

sub getchildren {
my ($lineage,@parents) = @_;
	my @children;
	foreach my $parent (@parents) {
		if (exists $$lineage{$parent}) {
			push @children, @{$$lineage{$parent}};
		}
	}
	return @children;
}

my (%taxids,%divids);
my $maxdepth = 0;
foreach my $rootid (keys %deftax) {
	my $div = $deftax{$rootid};
	$taxids{$rootid} = $div;
	push @{$divids{$div}{'taxonomy'}}, $rootid;
	my @parentids = ($rootid);
	my $depth = 0;
	while ((my @childids = &getchildren(\%lineage,@parentids)) && ($depth < 100)) {
		foreach my $taxid (@childids) {
			$taxids{$taxid} = $div;
			push @{$divids{$div}{'taxonomy'}}, $taxid;
		}
	@parentids = @childids;
	$depth++;
	}
	($depth > $maxdepth) && ($maxdepth = $depth)
}

my %seqprct;
foreach my $div (keys %divids) {
	my $fname = $div . '_pgi.txt';
	$divids{$div}{'fname'} = $fname;
	$divids{$div}{'fh'} = IO::File->new($fname, "w");
	$seqprct{$div} = 0;
}

my $allseqct = 0;
open BLDB, "blastdbcmd -db $blastdb -entry all -outfmt \"%T %g\" |";
while (<BLDB>) {
	chomp;
	$allseqct++;
	my ($taxid,$gi) = split " ";
	if (defined $taxids{$taxid}) {
		my $div = $taxids{$taxid};
		my $fh = $divids{$div}{'fh'};
		print $fh "$gi\n";
		$seqprct{$div}++;
	}
}

foreach my $div (keys %divids) {
	$divids{$div}{'fh'}->close;
}


printf STDERR "Building Blast subdatabases\n";
foreach my $div (keys %divids) {
	my $gilist = $divids{$div}{'fname'};
	my ($dbname) = $gilist =~ /^(\w+)_pgi.txt$/;
	$dbname .= "_nr";
	(my $gibin = $gilist) =~ s/.txt$/.bin/;
	`blastdb_aliastool -gi_file_in $gilist -gi_file_out $gibin`;
	`blastdb_aliastool -db nr -dbtype prot -out $dbname -gilist $gibin -title "$dbname NR Database"`;
}

printf STDERR "\n";
printf STDERR "Taxonomies in DB: %d\n", scalar keys %lineage;
printf STDERR "Defined root taxonomy IDs: %d\n", scalar keys %deftax;
printf STDERR "Total taxonomy IDs: %d\n", scalar keys %taxids;
printf STDERR "Maximum depth: %d\n", $maxdepth;
print STDERR "\n";
printf STDERR "Sequences processed: %d\n", $allseqct;
print STDERR  "   by division:\n";
foreach my $div (sort {$seqprct{$b} <=> $seqprct{$a}} keys %divids) {
	printf STDERR "%10s %8d\n", "$div:", $seqprct{$div};
}
print STDERR "\n";
printf STDERR "Time: %d seconds\n", time - $stime;

exit;
