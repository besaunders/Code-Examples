#!/usr/bin/perl -w
use strict;

my $t0 = time;
my $dbdir = '/cygdrive/c/Users/bsaunders/db';

my %arpref;
%arpref = map {$_,1} qw (19AJ 21BQ 21BP 19DJ 19DG 21EL 19K 19AV);
my $arlim = scalar keys %arpref;

my %deftax;
my $taxfile = shift @ARGV;
die "Invalid taxonomy file\n" unless (-f $taxfile);

open TF, $taxfile;
while (<TF>) {
	chomp;
	my ($tax,$blastdb,$class) = split /\t/;
	$deftax{$tax} = $class;
}

my %lineage;
open LIN, "$dbdir/nodes.dmp" or die "Cannot open taxonomy tree file\n";
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

my %taxids;
my $maxdepth = 0;
foreach my $rootid (keys %deftax) {
	my $type = $deftax{$rootid};
	$taxids{$rootid} = $type;
	my @parentids = ($rootid);
	my $depth = 0;
	while ((my @childids = &getchildren(\%lineage,@parentids)) && ($depth < 100)) {
		foreach my $taxid (@childids) {
			$taxids{$taxid} = $type;
		}
	@parentids = @childids;
	$depth++;
	}
	($depth > $maxdepth) && ($maxdepth = $depth)
}

print STDERR scalar keys %deftax, " root taxonomies defined\n";
print STDERR scalar keys %taxids, " taxonomy IDs found\n";
print STDERR "\n";

my (%arinfo,%upinfo,%acc2fgi,%fgi2up);
my $arenainfo = shift @ARGV;
open ARI, $arenainfo;
my $amhead = <ARI>;
while (<ARI>) {
	s/[\r\n]+$//;
	my @arinfo = split /\t/;
	my $arenaid = $arinfo[0];
	my $exmatch = $arinfo[1];
	my $title = $arinfo[11];
	my $seqstr = $arinfo[12];
	my $upstr = $arinfo[13];

	my ($fgi) = $seqstr =~ m/^gi\|(\d+)\|/;
	$arinfo{$arenaid}{'fgi'} = $fgi;
	$arinfo{$arenaid}{'match'} = $exmatch;
	$arinfo{$arenaid}{'name'} = $title;

	foreach my $upid (split /;/, $upstr) {
		$upinfo{$upid}{'fgi'}{$fgi} = 1;
		$fgi2up{$fgi}{$upid} = 2;
	}

	my @seqids = split /;/, $seqstr;
	my @accs;
	for (my $i = 0; $i <= $#seqids; $i++) {
		my ($idtype,$pgi,$db,$acc,$entry) = split /\|/, $seqids[$i];
		if ($acc) {
			$acc2fgi{$acc} = $fgi;
			push @accs, $acc;
		}
		elsif ($entry) {
			$acc2fgi{$entry} = $fgi;
			push @accs, $entry;
		}
	}
	push @{$arinfo{$arenaid}{'accs'}}, @accs;

}
close ARI;

print STDERR "$arlim target Arena IDs\n";
print STDERR scalar keys %arinfo, " Arena proteins\n";
print STDERR scalar keys %acc2fgi, " Arena accessions found\n";
print STDERR scalar keys %fgi2up, " Arena sequences with SP\n";
print STDERR scalar keys %upinfo, " UP matches in Arena\n";
print STDERR "\n";

my (%blphits,%blpinfo);
my $arenaid = 'INITID';
my $hitnum = 0;

# my $protbl = shift @ARGV;
# open PLB, $protbl;

while (<>) {
	s/[\r\n]+$//;
	my @blval = split /\t/;
	if ($arenaid ne $blval[0]) {
		$arenaid = $blval[0];
		$hitnum = 0;
	}
	$hitnum++;
	my $hitid = sprintf("%05d",$hitnum);
	my $ablid = $arenaid . ':' . $hitid;
	my ($argroup) = $arenaid =~ m/(\d+[A-Z]+)$/;
	(! $argroup) && (($argroup) = $arenaid =~ m/Arena_(.*$)/);
	my ($fgi) = $blval[13] =~ m/^gi\|(\d+)\|/;
	$blphits{$ablid}{'fgi'} = $fgi;
	$blphits{$ablid}{'evalue'} = $blval[11];
	$blphits{$ablid}{'score'} = $blval[12];
	$blphits{$ablid}{'ident'} = $blval[1];
	$blphits{$ablid}{'qident'} = sprintf "%.2f", 100 * $blval[2] / $blval[3];
				
	if (! exists $blpinfo{$fgi}) {
		$blpinfo{$fgi}{'taxid'} = $blval[14];
		$blpinfo{$fgi}{'species'} = $blval[15];
		$blpinfo{$fgi}{'name'} = $blval[16];
		my @seqids = split /;/, $blval[13];
		my @accs;
		for (my $i = 0; $i <= $#seqids; $i++) {
			my ($idtype,$pgi,$db,$acc,$entry) = split /\|/, $seqids[$i];
			if ($db eq 'sp') {
				(my $nvacc = $acc) =~ s/\..*//;
				$upinfo{$nvacc}{'fgi'}{$fgi} = 1;
				$fgi2up{$fgi}{$nvacc} = 2;
			}
			if (defined $acc) {
				$acc2fgi{$acc} = $fgi;
				push @accs, $acc;
			}
			elsif (defined $entry) {
				$acc2fgi{$entry} = $fgi;
				push @accs, $entry;
			}
		}

		push @{$blpinfo{$fgi}{'accs'}}, @accs;
		$blpinfo{$fgi}{'levalue'} = $blval[11];
		$blpinfo{$fgi}{'lscore'} = $blval[12];
		$blpinfo{$fgi}{'lident'} = $blval[1];
		$blpinfo{$fgi}{'llen'} = $blval[5] - $blval[4] + 1;
		$blpinfo{$fgi}{'topar'} = $arenaid;
		$blpinfo{$fgi}{'qident'} = sprintf "%.2f", 100 * $blval[2] / $blval[3];
	}
	push @{$blpinfo{$fgi}{'ablid'}}, $ablid;
	$blpinfo{$fgi}{'arenaid'}{$arenaid}++;
	if (exists $arpref{$argroup}) {
		$blpinfo{$fgi}{'arpref'}{$arenaid}++;
	}
	$blpinfo{$fgi}{'arenaid'}{$arenaid}++;
	if ($blval[12] > $blpinfo{$fgi}{'lscore'}) {
		$blpinfo{$fgi}{'lscore'} = $blval[12];
		$blpinfo{$fgi}{'topar'} = $arenaid;
		$blpinfo{$fgi}{'levalue'} = $blval[11];
		$blpinfo{$fgi}{'lident'} = $blval[1];
		$blpinfo{$fgi}{'llen'} = $blval[5] - $blval[4] + 1;
		$blpinfo{$fgi}{'qident'} = sprintf "%.2f", 100 * $blval[2] / $blval[3];
	}
}

my $t1 = time;

print STDERR scalar keys %acc2fgi, " total accessions found\n";
print STDERR scalar keys %fgi2up, " NR sequences with SP\n";
print STDERR scalar keys %upinfo, " SP matches in NR\n";
print STDERR scalar keys %blphits, " total Blast hits\n";
print STDERR scalar keys %blpinfo, " unique proteins matched with Blast\n";
print STDERR $t1 - $t0, " seconds to read in Blast data\n";
print STDERR "\n";

my $fgiupfile = 'up2fgi.dat';
if (-e $fgiupfile) {
	print STDERR "Reading UP to GI info from cache\n";
	open UTA, "$fgiupfile";
	while (<UTA>) {
		chomp;
		my ($fgi,$upid) = split /\t/;
		$upinfo{$upid}{'fgi'}{$fgi}++;
		(! exists $fgi2up{$fgi}{$upid}) && ($fgi2up{$fgi}{$upid} = 1);
	}
	close UTA;
}

else {
	print STDERR "Reading UP to ACC file\n";
	open UPTOID, "zcat $dbdir/idmapping_selected.tab.gz|";
	while (<UPTOID>) {
		my @ids = split /\t/;
		my $taxid = $ids[12];
		next unless defined($taxid) and defined($taxids{$taxid});
		my $upid = $ids[0];
		if ($ids[3] =~ /\w/) {
			foreach my $acc (split /; /, $ids[3]) {
				if (defined $acc2fgi{$acc}) {
					my $fgi = $acc2fgi{$acc};
					$upinfo{$upid}{'fgi'}{$fgi}++;
					(! exists $fgi2up{$fgi}{$upid}) && ($fgi2up{$fgi}{$upid} = 1);
				}
			}
		}
		if ($ids[11] =~ /\w/) {
			foreach my $acc (split /; /, $ids[11]) {
				if (defined $acc2fgi{$acc}) {
					my $fgi = $acc2fgi{$acc};
					$upinfo{$upid}{'fgi'}{$fgi}++;
					(! exists $fgi2up{$fgi}{$upid}) && ($fgi2up{$fgi}{$upid} = 1);
				}
			}
		}
		if ($ids[17] =~ /\w/) {
			foreach my $acc (split /; /, $ids[17]) {
				if (defined $acc2fgi{$acc}) {
					my $fgi = $acc2fgi{$acc};
					$upinfo{$upid}{'fgi'}{$fgi}++;
					(! exists $fgi2up{$fgi}{$upid}) && ($fgi2up{$fgi}{$upid} = 1);
				}
			}
		}
	}
	close UPTOID;
	open UTA, ">$fgiupfile";
	foreach my $fgi (keys %fgi2up) {
		foreach my $upid (keys %{$fgi2up{$fgi}}) {
			print UTA "$fgi\t$upid\n";
		}
	}
	close UTA;
}

my $t2 = time;
print STDERR scalar keys %fgi2up, " NR sequences with Uniprot\n";
print STDERR scalar keys %upinfo, " Uniprot matches in NR\n";
print STDERR $t2 - $t1, " seconds to map Uniprot to GI\n";
print STDERR "\n";

my (%ip2up, %ipname);
my $totmatch = 0;
my $ip2upfile = 'up2ip.dat';
my $ipnfile = 'iprnames.dat';
if ((-e $ip2upfile) && (-e $ipnfile)) {
	print STDERR "Reading UP to Interpro info from cache\n";
	open UTI, $ip2upfile;
	while (<UTI>) {
		chomp;
		my ($iprid,$upid) = split /\t/;
		$upinfo{$upid}{'ipr'}{$iprid}++;
		$ip2up{$iprid}{$upid}++;
		$totmatch++;
	}
	close UTI;
	open UTN, $ipnfile;
	while (<UTN>) {
		chomp;
		my ($iprid,$name) = split /\t/;
		$ipname{$iprid} = $name;
	}
	close UTN;
}

else {
	print STDERR "Reading UP to Interpro file\n";
	open UPTOIP, "zcat $dbdir/protein2ipr.dat.gz|";
	while (<UPTOIP>) {
		my ($upid,$iprid,$fname,$match) = split /\t/;
		next unless exists $upinfo{$upid};
		(! defined $ipname{$iprid}) && ($ipname{$iprid} = $fname);
		$upinfo{$upid}{'ipr'}{$iprid}++;
		$ip2up{$iprid}{$upid}++;
		$totmatch++;
	}
	close UPTOIP;
	open UTI, ">$ip2upfile";
	foreach my $iprid (keys %ip2up) {
		foreach my $upid (keys %{$ip2up{$iprid}}) {
			print UTI "$iprid\t$upid\n";
		}
	}
	close UTI;
	open IPRN, ">$ipnfile";
	foreach my $iprid (keys %ipname) {
		print IPRN "$iprid\t$ipname{$iprid}\n";
	}
}


my $t3 = time;
print STDERR scalar keys %ip2up, " Interpro IDs found\n";
print STDERR "$totmatch total matches\n";
print STDERR $t3 - $t2, " seconds to map Uniprot to Interpro\n";
print STDERR "\n";

open PHITS, ">prefhits.txt";

my $phheader = join "\t", 'Class', 'Arena IDs', 'Match Name',
	'Evalue', 'F match', 'GPCR', "Match Rep Interpro",
	'Match Uniprot', 'Match GB Acc', 'Pct ID (query)', 'Blast Score';
print PHITS "$phheader\n";

foreach my $fgi (keys %blpinfo) {
	next unless exists $blpinfo{$fgi}{'arpref'};

	my (%aipr,%mipr);
  my ($iprlink,$upstr) = ('','');
	foreach my $upid (sort {
		$fgi2up{$fgi}{$b} <=> $fgi2up{$fgi}{$a}
			||
		$a cmp $b
		} keys %{$fgi2up{$fgi}}) {
		$upstr .= $upid . '; ';
		foreach my $iprid (keys %{$upinfo{$upid}{'ipr'}}) {
			if ($iprlink eq '') {
				$iprlink = '=HYPERLINK("http://www.ebi.ac.uk/interpro/protein/' .
					$upid . '","' . $upid . '")';
			}
			$mipr{$iprid}++;
		}
	}
	$upstr =~ s/; $//;

	foreach my $arenaid (keys %{$blpinfo{$fgi}{'arpref'}}) {
		foreach my $upid (keys %{$fgi2up{$arinfo{$arenaid}{'fgi'}}}) {
			foreach my $iprid (keys %{$upinfo{$upid}{'ipr'}}) {
				$aipr{$iprid}++;
			}
		}
	}

	my ($dmatch,$gpcr) = (0,0);
	my ($matchstr,$nmstr) = ('','');
	foreach my $iprid (keys %mipr) {
	  if ($ipname{$iprid} =~ /(GPCR)|(coupled)/) {
			$gpcr = 1;
		}
		if (exists $aipr{$iprid}) {
			$dmatch = 1;
			$matchstr .= $iprid . ': ' . $ipname{$iprid} . '; ';
		}
		else {
			$nmstr .= $iprid . ': ' . $ipname{$iprid} . '; ';
		}
	}
	$nmstr =~ s/; $//;
	$matchstr =~ s/; $//;

	my $accstr = join '; ', @{$blpinfo{$fgi}{'accs'}};
	my $parstr = join '; ', sort keys %{$blpinfo{$fgi}{'arpref'}};
	my $taxids =  $blpinfo{$fgi}{'taxid'};
	(my $taxid = $taxids) =~ s/;.*//;
	my $oclass = $taxids{$taxid};
	if ($taxids =~ /;/) {
		foreach my $itaxid (split /;/, $taxids) {
			if ($taxids{$itaxid} ne $oclass) {
				print STDERR "Warning: multiple classes for $taxids\n";
			}
		}
	}
	my $pname = $blpinfo{$fgi}{'name'};
	if (length($pname) > 32767) {
		$pname = substr $pname, 0, 32767;
		$pname =~ s/<>.*?$//;
	}
	$pname =~ s/<>/ >/g;
	
	my $prline = join "\t", $oclass, $parstr, $pname,
		$blpinfo{$fgi}{'levalue'}, $dmatch, $gpcr, $iprlink, 
		$upstr, $accstr, $blpinfo{$fgi}{'qident'}, $blpinfo{$fgi}{'lscore'};
	print PHITS "$prline\n";
}


close PHITS;

open HITS, ">allblhits.txt";

my $abheader = join "\t", 'Class', '#Arena', 'Top Arena', 'Pref',
	'Match Name', 'Evalue', 'F match', 'GPCR', "Match Rep Interpro",
	'Match Uniprot', 'Match GB Acc', 'Pct ID (query)', 'Blast Score';

print HITS "$abheader\n";

foreach my $fgi (keys %blpinfo) {
	my (%aipr,%mipr);
  my ($iprlink,$upstr) = ('','');
	foreach my $upid (sort {
		$fgi2up{$fgi}{$b} <=> $fgi2up{$fgi}{$a}
			||
		$a cmp $b
		} keys %{$fgi2up{$fgi}}) {
		$upstr .= $upid . '; ';
		foreach my $iprid (keys %{$upinfo{$upid}{'ipr'}}) {
			if ($iprlink eq '') {
				$iprlink = '=HYPERLINK("http://www.ebi.ac.uk/interpro/protein/' .
					$upid . '","' . $upid . '")';
			}
			$mipr{$iprid}++;
		}
	}
	$upstr =~ s/; $//;

	foreach my $arenaid (keys %{$blpinfo{$fgi}{'arenaid'}}) {
		foreach my $upid (keys %{$fgi2up{$arinfo{$arenaid}{'fgi'}}}) {
			foreach my $iprid (keys %{$upinfo{$upid}{'ipr'}}) {
				$aipr{$iprid}++;
			}
		}
	}

	my ($dmatch,$gpcr) = (0,0);
	my ($matchstr,$nmstr) = ('','');
	foreach my $iprid (keys %mipr) {
	  if ($ipname{$iprid} =~ /(GPCR)|(coupled)/) {
			$gpcr = 1;
		}
		if (exists $aipr{$iprid}) {
			$dmatch = 1;
			$matchstr .= $iprid . ': ' . $ipname{$iprid} . '; ';
		}
		else {
			$nmstr .= $iprid . ': ' . $ipname{$iprid} . '; ';
		}
	}
	$nmstr =~ s/; $//;
	$matchstr =~ s/; $//;

	my $arcount = scalar keys %{$blpinfo{$fgi}{'arenaid'}};
	my $accstr = join '; ', @{$blpinfo{$fgi}{'accs'}};
	my $taxids =  $blpinfo{$fgi}{'taxid'};
	(my $taxid = $taxids) =~ s/;.*//;
	my $oclass = $taxids{$taxid};
	if ($taxids =~ /;/) {
		foreach my $itaxid (split /;/, $taxids) {
			if ($taxids{$itaxid} ne $oclass) {
				print STDERR "Warning: multiple classes for $taxids\n";
			}
		}
	}
	my $pname = $blpinfo{$fgi}{'name'};
	if (length($pname) > 32767) {
		$pname = substr $pname, 0, 32767;
		$pname =~ s/<>.*?$//;
	}
	$pname =~ s/<>/ >/g;

	my $pref = 0;
	(exists $blpinfo{$fgi}{'arpref'}) && ($pref = 1);

	my $prline = join "\t", $oclass, $arcount, $blpinfo{$fgi}{'topar'},
		$pref, $pname, $blpinfo{$fgi}{'levalue'}, $dmatch, $gpcr, $iprlink, 
		$upstr, $accstr, $blpinfo{$fgi}{'qident'}, $blpinfo{$fgi}{'lscore'};

	print HITS "$prline\n";
}


close HITS;

open ARENAID, ">allmatches.txt";

my $aheader = join "\t", 'Arena ID', 'ArGr', 'Pref', 'Rep Arena Name',
	'Arena Rep Interpro', 'Arena Uniprot', 'ArIDSeq', 'Arena GB Acc',
	'Hit', 'Class', 'Match Name', 'Evalue', 'F match', 'GPCR',
	'Match Rep Interpro', 'Match Uniprot', 'Match GB Acc',
	'Pct ID (query)', 'Blast Alignment', 'CDD Alignment';
print ARENAID "$aheader\n";

foreach my $ablid (sort {$a cmp $b} keys %blphits) {
	my ($arenaid,$hitct) = split /:/, $ablid;
	$hitct =~ s/^0+//;
	my $fgi = $blphits{$ablid}{'fgi'};
	my $taxids =  $blpinfo{$fgi}{'taxid'};
	(my $taxid = $taxids) =~ s/;.*//;
	my $oclass = $taxids{$taxid};
	my ($argroup) = $arenaid =~ m/(\d+[A-Z]+)$/;
	(! $argroup) && (($argroup) = $blpinfo{$fgi}{'topar'} =~ m/Arena_(.*$)/);
	my $pref = 0;
	(exists $arpref{$argroup}) && ($pref = 1);

	my $afgi = $arinfo{$arenaid}{'fgi'};
	my $aname = $arinfo{$arenaid}{'name'};
	my $armatch = $arinfo{$arenaid}{'match'};
	$aname =~ s/<>.*//;
	
	if ($taxids =~ /;/) {
		foreach my $itaxid (split /;/, $taxids) {
			if ($taxids{$itaxid} ne $oclass) {
				print STDERR "Warning: multiple classes for $taxids\n";
			}
		}
	}
	my $pname = $blpinfo{$fgi}{'name'};
	my $accstr = join '; ', @{$blpinfo{$fgi}{'accs'}};
	my $araccstr = join '; ', @{$arinfo{$arenaid}{'accs'}};

	my (%aipr,%mipr);
	my ($iprlink,$upstr,$aiprlink,$aupstr) = ('') x 4;

	foreach my $upid (sort {
		$fgi2up{$fgi}{$b} <=> $fgi2up{$fgi}{$a}
			||
		$a cmp $b
		} keys %{$fgi2up{$fgi}}) {
		$upstr .= $upid . '; ';
		foreach my $iprid (keys %{$upinfo{$upid}{'ipr'}}) {
			if ($iprlink eq '') {
				$iprlink = '=HYPERLINK("http://www.ebi.ac.uk/interpro/protein/' .
					$upid . '","' . $upid . '")';
			}
			$mipr{$iprid}++;
		}
	}
	$upstr =~ s/; $//;

	foreach my $upid (sort {
		$fgi2up{$afgi}{$b} <=> $fgi2up{$afgi}{$a}
			||
		$a cmp $b
		} keys %{$fgi2up{$afgi}}) {
		$aupstr .= $upid . '; ';
		foreach my $iprid (keys %{$upinfo{$upid}{'ipr'}}) {
			if ($aiprlink eq '') {
				$aiprlink = '=HYPERLINK("http://www.ebi.ac.uk/interpro/protein/' .
					$upid . '","' . $upid . '")';
			}
			$aipr{$iprid}++;
		}
	}
	$aupstr =~ s/; $//;

	my ($dmatch,$gpcr) = (0,0);
	my ($matchstr,$nmstr,$amstr) = ('') x 3;
	foreach my $iprid (keys %mipr) {
	  if ($ipname{$iprid} =~ /(GPCR)|(coupled)/) {
			$gpcr = 1;
		}
		if (exists $aipr{$iprid}) {
			$dmatch = 1;
			$matchstr .= $iprid . ': ' . $ipname{$iprid} . '; ';
		}
		else {
			$nmstr .= $iprid . ': ' . $ipname{$iprid} . '; ';
		}
	}
	foreach my $iprid (keys %aipr) {
		if (! exists $mipr{$iprid}) {
			$amstr .= $iprid . ': ' . $ipname{$iprid} . '; ';
		}
	}

	$amstr =~ s/; $//;
	$nmstr =~ s/; $//;
	$matchstr =~ s/; $//;

	my $bl2link = '=HYPERLINK("http://blast.ncbi.nlm.nih.gov/' .
		'BlastAlign.cgi?PROGRAM=blastp&CMD=request&QUERY=' . $afgi .
		'&SUBJECTS=' . $fgi . '","' . $blphits{$ablid}{'score'} . '")';

	my $cddlink = '=HYPERLINK("http://www.ncbi.nlm.nih.gov/Structure/' .
		'bwrpsb/bwrpsb.cgi?queries=' . $afgi . '%0A' . $fgi .
		'","CDD Link")';

	my $prline = join "\t", $arenaid, $argroup, $pref, $aname, $aiprlink,
		$aupstr, $armatch, $araccstr, $hitct, $oclass, $pname,
		$blphits{$ablid}{'evalue'}, $dmatch, $gpcr, $iprlink, $upstr, $accstr,
		$blphits{$ablid}{'qident'}, $bl2link, $cddlink;
	print ARENAID "$prline\n";
}

my $t4 = time;
print STDERR $t4 - $t3, " seconds to process Blast data\n";
print STDERR "\n";
