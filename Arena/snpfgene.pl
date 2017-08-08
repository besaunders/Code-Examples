#!/usr/bin/perl -w

# This program uses the NCBI E-utilities to obtain a list of dbSNP
# records that pertain to a defined NCBI Gene ID, then it returns in
# XML format the detailed information for those dbSNP records. LWP is
# used to do the HTML connection

# A tab-delimited report is produced for each SNP matching the gene,
# containing the information of most interest to the project

use strict;
use LWP::UserAgent;

my $stime = time;

# my $geneid = 1269;

unless (@ARGV) {
	die "Need to provide Gene ID\n";
}

my $geneid = shift @ARGV;
unless ($geneid =~ /^\d+$/) {
	die "Invalid Gene ID $geneid\n";
}

# Find Gene ID in Arena annotation file; store all Arena and protein IDs
# that match the Gene ID

my %mprots;
my $arstr = '';
my $smatch = 1;
my $armatchf = shift @ARGV;
if ($armatchf && -e $armatchf) {
	print STDERR "Arena match file $armatchf will be used\n\n";
	# Read in Arena Match file
	open ARMATCH, $armatchf;
	my $header = <ARMATCH>;
	while (<ARMATCH>) {
		chomp;
		my @armatch = split /\t/;
		foreach my $ageneid (split /;/, $armatch[8]) {
			if ($ageneid == $geneid) {
				$arstr .= $armatch[0] . '; ';
				$smatch = $armatch[1];
				foreach my $seqid (split /;/, $armatch[12]) {
					my ($mprotid) = (split /\|/, $seqid)[3];
					$mprots{$mprotid}++;
				}
			}
		}
	}
	$arstr =~ s/; $//;
}
else {
	print STDERR "Invalid or no Arena match file defined\n";
}

# Form command to retrieve all dbSNP records pertaining to Gene id
my $gslink = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi';
my $gsargs = 'dbfrom=gene&db=snp&id=' . $geneid . '&cmd=neighbor_history';

# create HTTP user agent, request object, and send request
my $ua = new LWP::UserAgent;
$ua->agent("elink/2.0 " . $ua->agent);
my $gsreq = new HTTP::Request POST => "$gslink";
$gsreq->content_type('application/x-www-form-urlencoded');
$gsreq->content("$gsargs");
my $gsresponse = $ua->request($gsreq);
my $gskeywe= $gsresponse->content;

# Extract keys allowing access to dbSNP ID list
my ($gskey) = $gskeywe =~ m#<QueryKey>(\d+)</QueryKey>#;
my ($gswebenv) = $gskeywe =~ m#<WebEnv>([^<]+)</WebEnv>#;

if (! defined $gskey) {
	die "Gene ID $geneid has no related dbSNP records or does not exist\n";
}

# Form command to retrieve dbSNP info for IDs in previous list
my $sxlink = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
my $sxargs = '&db=snp&retmode=xml&query_key=' . $gskey .
	'&WebEnv=' . $gswebenv;

my $sxreq = new HTTP::Request POST => "$sxlink";
$sxreq->content_type('application/x-www-form-urlencoded');
$sxreq->content("$sxargs");
my $sxresponse = $ua->request($sxreq);
my $sxop= $sxresponse->content;

# Open anotation and SNP type distribution files, print annotation header,
# and loop through dbSNP records

# Loop through all SNP information. NCBI returns one record for each
# line, so a change in format would require a change in the parsing
# mechanism

# The parsing script below depends on the elements of the dbSNP XML file
# (i.e. the DTD) remaining relative stable. The user of this program should
# verify from time to time that the elements are being obtained correctly,
# and incorporate any new information of interest

my $sfile = $geneid . '_snps.txt';
my $ifile = $geneid . '_infosnp.txt';
my (%ifclass,%allfclass,%f1pclass);

open SNPF, ">$sfile";

print SNPF "SNP_ID\tdbSNP\tArena_ID\tType\tCoding\tClinical\tPubmed\t",
	"Verification\tFrequency\tFreq_Info\tAA_Change\tAA_Match\t",
	"Contig_Loc\tNT_Change\n";
my $dbscount = 0;

foreach (split /[\n\r]+/, $sxop) {
	next unless /<Rs /;
	unless (m#</Rs>#) {
		die "Error: incomplete XML record\n";
	}
	$dbscount++;

# The rsID field, i.e. the ID of the dbSNP record, is the central identifying
# element. Without it existing or being unique, this script would probably
# not work correctly
	my @rstag = /<Rs (.*?)>/;
	my ($rsid) = $rstag[0] =~ /rsId="(\d+)"/;
	if (! defined $rsid) {
		die "Warning: cannot obtain rsId in $_\n";
	}
	elsif (@rstag > 1) {
		die "Warning: more than one Rs tag in $_\n";
	}

# The FxnSet entries contain each class of the SNP, as well as the database
# entities that correspond to areas where the SNP acts. Some SNPs can
# correspond to multiple GeneIDs, where there are collisions between
# upstream and downstream entries, or even non-coding regions. We will
# only consider records correpsonding to the specified Gene ID

# Entries can have multiple classes for the same GeneID, especially in
# cases where there are multiple splice variants speicified for the gene.
# We are specificly considering all cases where amino acid information is
# provided.

	my $symbol;
	my %fclasses;
	my %protst;
	my $matchct = 0;
	while (m#<FxnSet (.*?)>#g) {
		my $fxnset = $1;
		my ($fgeneid) = $fxnset =~ /geneId="([^"]+)/;
		if (! $fgeneid) {
			print STDERR "Warning: Fxnset w/o Gene ID for $rsid\n";
			next;
		}
		elsif ($fgeneid != $geneid) {
			next;
		}
		my ($fxnclass) = $fxnset =~ /fxnClass="([^"]+)/;
		if (! defined $fxnclass) {
			print STDERR "Warning: fxnClass not defined in $rsid\n";
		}
		my ($fsymbol) = $fxnset =~ /symbol="([^"]+)/;
		if (! defined $fsymbol) {
			print STDERR "Warning: symbol not defined in $rsid\n";
		}
		elsif ((defined $symbol) && ($fsymbol ne $symbol)) {
			print STDERR "Warning: multiple gene symbols in $rsid\n";
		}
		elsif (! defined $symbol) {
			$symbol = $fsymbol;
		}
		if ($fxnclass ne 'reference') {
			$fclasses{$fxnclass}++;
			}
		if (($fxnclass eq 'reference') || ($fxnclass eq 'missense') ||
			($fxnclass eq 'stop-gained') || ($fxnclass eq 'stop-lost') ||
			($fxnclass eq 'frameshift-variant')) {
			my ($protacc) = $fxnset =~ /protAcc="([^"]+)/;
			my ($protv) = $fxnset =~ /protVer="([^"]+)/;
			my ($aapos) = $fxnset =~ /aaPosition="([^"]+)/;
			my ($residue) = $fxnset =~ /residue="([^"]+)/;
			my $protid = $protacc . '.' . $protv;
			$protst{$protid}{$fxnclass}{'aapos'} = $aapos;
			$protst{$protid}{$fxnclass}{'res'} = $residue;
			if (exists $mprots{$protid}) {
				$matchct++;
			}
		}
	}

	my $idstr = "$geneid";
	(defined $symbol) && ($idstr .= " ($symbol)");
	(length $arstr > 0) && ($idstr .= ": $arstr");
	$idstr =~ s/Arena_//g;

# Get all validation information
	my @vmatches = m#<Validation(.*?)>#g;
	if (@vmatches > 1) {
		print STDERR "Warning: more than one validation tag for $rsid\n";
	}
	my $vmatch = $vmatches[0];
	$vmatch =~ s/^\s+//;
	$vmatch =~ s/\s+$//;
	$vmatch =~ s#\s*/$##;
	$vmatch =~ s/="true"//g;
	$vmatch =~ s/^by//;
	$vmatch =~ s/ by/ /g;
	$vmatch =~ s/ /; /g;

# The component tag contains the genomic location for the SNP, and
# should be unique. We are using the start location.
	my @compon = m#(<Component.*?</Component>)#g;
	my $cloc = 0;
	if (@compon > 1) {
		print STDERR "Warning: more than Component tag for $rsid\n";
	}
	if (@compon > 0) {
		($cloc) = $compon[0] =~ /asnFrom="(\d+)/;
	}

# Parse any frequency tag that is defined (sometimes multiple tags
# exist)
	my @frqmatches = m#<Frequency(.*?)>#g;
	my $freqstr = '';
	my %freqval;
	if (@frqmatches > 0) {
		my ($samplesz) = $frqmatches[0] =~ m/sampleSize="(\d+)"/;
		my ($fallele) = $frqmatches[0] =~ m/allele="([^"]+)"/;
		my ($frequency) = $frqmatches[0] =~ m/freq="([^"]+)"/;
		$freqstr = $fallele . ': ' . $frequency . " ($samplesz)";
		$freqval{$frequency}++;
	}
	for (my $i = 1; $i <= $#frqmatches; $i++) {
		my ($samplesz) = $frqmatches[$i] =~ m/sampleSize="(\d+)"/;
		my ($fallele) = $frqmatches[$i] =~ m/allele="([^"]+)"/;
		my ($frequency) = $frqmatches[$i] =~ m/freq="([^"]+)"/;
		$freqstr .= '; ' . $fallele . ': ' . $frequency . " ($samplesz)";
		$freqval{$frequency}++;
	}
	if ((scalar keys %freqval) > 1) {
		print STDERR "Warning: more than one frequency values for $rsid\n";
	}
	my $maxfreq = '';
	if ((scalar keys %freqval) > 0) {
		$maxfreq = (sort {$freqval{$b} <=> $freqval{$a}} keys %freqval)[0];
	}

# The observed tag contains nucleic variation information
	my %observed;
	while (m#<Observed>(.*?)</Observed>#g) {
		$observed{$1}++;
	}
	my $observed = join('; ', sort keys %observed);

# Loop through all proteins that have been defined, and derive information
# that pertains to amino acid substitution
	my %aasub;
	my $trank = 0;
	foreach my $protid (keys %protst) {
		if (exists $protst{$protid}{'missense'}) {
			my $resch = $protst{$protid}{'reference'}{'res'} . '->' .
				$protst{$protid}{'missense'}{'res'};
			$aasub{$resch}{$protst{$protid}{'reference'}{'aapos'}}++;
			$aasub{$resch}{$protst{$protid}{'missense'}{'aapos'}}++;
			$trank = 1;
		}
		if (exists $protst{$protid}{'stop-gained'}) {
			my $resch = $protst{$protid}{'reference'}{'res'} . '->' .
				$protst{$protid}{'stop-gained'}{'res'};
			$aasub{$resch}{$protst{$protid}{'reference'}{'aapos'}}++;
			$aasub{$resch}{$protst{$protid}{'stop-gained'}{'aapos'}}++;
			$trank = 1;
		}
		if (exists $protst{$protid}{'stop-lost'}) {
			my $resch = $protst{$protid}{'reference'}{'res'} . '->' .
				$protst{$protid}{'stop-lost'}{'res'};
			$aasub{$resch}{$protst{$protid}{'reference'}{'aapos'}}++;
			$aasub{$resch}{$protst{$protid}{'stop-lost'}{'aapos'}}++;
			$trank = 1;
		}
		if (exists $protst{$protid}{'frameshift-variant'}) {
			my $resch = $protst{$protid}{'reference'}{'res'} . '->' .
				$protst{$protid}{'frameshift-variant'}{'res'};
			$aasub{$resch}{$protst{$protid}{'reference'}{'aapos'}}++;
			$aasub{$resch}{$protst{$protid}{'frameshift-variant'}{'aapos'}}++;
			$trank = 1;
		}
	}

	my $aastr = '';
	my $rfct = 0;
	my $aamatch = '';
	foreach my $iaasub (sort keys %aasub) {
		$aamatch = 'Yes';
		$rfct++;
		$aastr .= $iaasub;
		my $aaposstr = join '; ', sort keys %{$aasub{$iaasub}};
		$aastr .= ' (' . $aaposstr . '); ';
	}
	$aastr =~ s/; $//;

	if ((scalar keys %fclasses) == 0) {
		$fclasses{'unknown'}++;
	}
	my $fclstr = join '; ', sort keys %fclasses;

# The clinical significance field contains annotated information about
# SNPs that have been defined to cause disease effects
	my %clinsigs;
	my @clinsig = m#<ClinicalSignificance>(.*?)</ClinicalSignificance>#g;
	foreach my $clinv (@clinsig) {
		$clinsigs{$clinv}++;
	}
	my $clinstr = join '; ', sort keys %clinsigs;
	if ($clinstr ne '') {
		$clinstr = '=HYPERLINK("http://www.ncbi.nlm.nih.gov/clinvar/?' .
			'term=rs' . $rsid . '","' . $clinstr . '")';
	}

	if ($aamatch eq 'Yes') {
		if ((! $smatch) || ($rfct > 1) || (($rfct > 0) && ($matchct == 0)) ||
			($aastr =~ /;/) || ($idstr =~ /;/)) {
			$aamatch = 'No';
		}
	}

# Store all frequency information, and specify those SNPs with over a
# 1% occurence
	$allfclass{$fclstr}++;
	foreach my $fclass (keys %fclasses) {
		$ifclass{$fclass}++;
		if ($maxfreq && ($maxfreq > 0.01)) {
			$f1pclass{$fclass}++;
		}
	}

	my $snpweb = '=HYPERLINK("http://www.ncbi.nlm.nih.gov/projects/SNP/' .
		'snp_ref.cgi?rs=' . $rsid . '","' . 'dbSNP")';

	my $pmlink = '';
	if (m#<RsLinkout [^>]*?resourceId="5"#) {
		$pmlink = '=HYPERLINK("http://www.ncbi.nlm.nih.gov/pubmed/?' .
			'Db=pubmed&DbFrom=snp&Cmd=Link&LinkName=snp_pubmed_cited&' .
			'IdsFromResult=' . $rsid . '","' . 'PubMed")';
	}

 	print SNPF "$rsid\t$snpweb\t$idstr\t$fclstr\t$trank\t$clinstr\t$pmlink\t",
		"$vmatch\t$maxfreq\t$freqstr\t$aastr\t$aamatch\t$cloc\t$observed\n";

}
close SNPF;

open INFO, ">$ifile";

print STDERR "\n";
printf STDERR "%-25s %6s %6s\n", 'FXN Types Represented', 'Total', '>1%';
printf INFO "%s\t%s\t%s\n", 'FXN Types Represented', 'Total', '>1%';
print STDERR ('-') x 39, "\n";
foreach my $fclass (sort {$a cmp $b} keys %ifclass) {
	my $f1pcount = 0;
	if (exists $f1pclass{$fclass}) {
		$f1pcount =  $f1pclass{$fclass};
	}
	my $fcpr = $fclass;
	if (($fclass eq 'missense') || ($fclass eq 'stop-gained') ||
		($fclass eq 'stop-lost') || ($fclass eq 'frameshift-variant')) {
		$fcpr .= '*';
	}
	printf STDERR "%-25s %6d %6d\n", $fcpr, $ifclass{$fclass}, $f1pcount;
	printf INFO "%s\t%d\t%d\n", $fcpr, $ifclass{$fclass}, $f1pcount;
}
print STDERR "\n";
print INFO "\n";

print STDERR "Multiple FXN Types\n";
print STDERR ('-') x 57, "\n";
foreach my $fclass (sort {$a cmp $b} keys %allfclass) {
	next unless ($fclass =~ /;/);
	printf STDERR "%-60s %6d\n", $fclass, $allfclass{$fclass};
	printf INFO "%s\t%d\n", $fclass, $allfclass{$fclass};
}

close INFO;
print STDERR "\n";
print STDERR "$dbscount dbSNP records total\n";
print STDERR "\n";
print STDERR time - $stime, " seconds\n";
