#!/usr/bin/perl -w
use strict;
use LWP::UserAgent;

my $starttime = time;

my $upmap = '/cygdrive/c/Users/bsaunders/db/idmapping_selected.tab.gz';

# qseqid pident nident qlen qstart qend slen sstart send gapopen
# gaps evalue score sallseqid staxids sscinames salltitles

my %arprots;
my $hitct = 0;
while (<>) {
	chomp;
	my ($arid,$pctid,$numid,$qlen,$qst,$qend,$slen,$sst,$send,
		$gapopens,$gaps,$evalue,$score,$seqid,$staxids,$sscinames,
		$title) = split /\t/;
	$title =~ s/\s+$//;

	my %dbs;
	while ($seqid =~ m/gi\|(\d+)\|(\w+)\|/g) {
		my $pgi = $1;
		my $db = $2;
		$dbs{$db}++;
	}

	my $dbstr = join '', sort keys %dbs;
	next if (($dbstr eq 'pir') || ($dbstr eq 'prf') ||
		($dbstr eq 'pirprf'));
	
	if (! exists $arprots{$arid}) {
		$hitct = 1;
		$arprots{$arid}{'len'} = $qlen;
		$arprots{$arid}{'hit'} = $hitct;
		$arprots{$arid}{'numid'} = $numid;
		$arprots{$arid}{'slen'} = $slen;
		$arprots{$arid}{'seqids'} = $seqid;
		$arprots{$arid}{'titles'} = $title;
	}
	else {
		$hitct++;
		next if ($numid < $arprots{$arid}{'numid'});
		next if (abs($slen - $qlen) >= abs($arprots{$arid}{'slen'} - $qlen));
		$arprots{$arid}{'hit'} = $hitct;
		$arprots{$arid}{'numid'} = $numid;
		$arprots{$arid}{'slen'} = $slen;
		$arprots{$arid}{'seqids'} = $seqid;
		$arprots{$arid}{'titles'} = $title;
	}
}

my (%prot2ar,%prot2up);
my $arct = 0;
foreach my $arid (sort {$a cmp $b} keys %arprots) {
	$arct++;
	my $seqid = $arprots{$arid}{'seqids'};
	while ($seqid =~ m/gi\|(\d+)\|(\w+)\|(.*?)\|/g) {
		my $pgi = $1;
		my $db = $2;
		my $pacc = $3;
		$prot2ar{$pgi}{$arid}++;
		$prot2up{$pgi}{'db'} = $db;
		if ($db eq 'sp') {
			$pacc =~ s/\..*//;
			$prot2up{$pgi}{'up'}{$pacc}++;
		}
	}
}

my $pidstr = join '&id=', sort {$a <=> $b} keys %prot2ar;

printf STDERR "Arena proteins: %d\n", $arct;
printf STDERR "Length GI string to NCBI: %d\n", length($pidstr);

my $elink = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi';
my $eargs = 'dbfrom=protein&db=gene&id=';

#create HTTP user agent
my $ua = new LWP::UserAgent;
$ua->agent("elink/2.0 " . $ua->agent);

my $epost = $eargs . $pidstr;

#create HTTP request object
my $req = new HTTP::Request POST => "$elink";
$req->content_type('application/x-www-form-urlencoded');
$req->content("$epost");

#post the HTTP request
my $response = $ua->request($req); 
my $eres = $response->content;

printf STDERR "Length content from NCBI: %d\n", length($eres);

my $readpid = 0;
my $readlink = 0;
my $protid;
my %gene2ar;

foreach my $line (split /\n/, $eres) {
	if ($line =~ m#<IdList>#) {
		$readpid = 1;
	}
	elsif ($line =~ m#</IdList>#) {
		$readpid = 0;
	}
	elsif ($line =~ m#<Link>#) {
		$readlink = 1;
	}
	elsif ($line =~ m#</Link>#) {
		$readlink = 0;
	}
	elsif ($line =~ m#<Id>(\d+)</Id>#) {
		if ($readpid) {
			$protid = $1;
		}
		elsif ($readlink) {
			foreach my $arid (keys $prot2ar{$protid}) {
				$arprots{$arid}{'genes'}{$1}++;
				$gene2ar{$1}{$arid}++;
			}
		}
	}
	elsif ($line =~ m#</LinkSet>#) {
		undef $protid;
	}
}

my %argroups;
foreach my $protid (keys %prot2ar) {
	if (scalar keys %{$prot2ar{$protid}} > 1) {
		my $argroup = join '; ', sort {$a cmp $b} keys %{$prot2ar{$protid}};
		$argroups{$argroup}++;
		my $matchar = scalar keys %{$prot2ar{$protid}};
		foreach my $arid (keys %{$prot2ar{$protid}}) {
			$arprots{$arid}{'matchar'} = $argroup;
		}
	}
}

open UPTOID, "zcat $upmap|cut -f 1,5|";
while (<UPTOID>) {
	chomp;
	my ($upid,$pgis) = split /\t/;
	if ($pgis =~ /\w/) {
		foreach my $pgi (split /; /, $pgis) {
			if (exists $prot2up{$pgi}) {
				($prot2up{$pgi}{'up'}{$upid} = 1);
			}
		}
	}
}
close UPTOID;


print "Arena_ID\tExact\tPct_Ident\tMismatch\tQuery_Len\tTarget_Len\t",
	"HitNum\tGenes\tGene_IDs\tNCBI_Gene\tMatching_Arena\t",
	"Sequence_Titles\tSequence_ID\tRep_Uniprot\tRep_NCBI_Acc\n";
my ($idct,$ngct,$mgct) = (0) x 3;
foreach my $arid (sort {$a cmp $b} keys %arprots) {
	my $idmatch = 0;
	my $qlen = $arprots{$arid}{'len'};
	my $hitct = $arprots{$arid}{'hit'};
	my $numid = $arprots{$arid}{'numid'};
	my $slen = $arprots{$arid}{'slen'};
	my $seqid = $arprots{$arid}{'seqids'};
	my $argp = '';
	if (defined $arprots{$arid}{'matchar'}) {
		$argp = $arprots{$arid}{'matchar'};
	}

	my $genect = scalar keys %{$arprots{$arid}{'genes'}};
	my ($genestr,$genelink) = ('','');
	if ($genect == 0) {
		$ngct++;
	}
	else {
		my @genes = sort {$a <=> $b} keys %{$arprots{$arid}{'genes'}};
		$genestr = join '; ', @genes;
		$genelink = '=HYPERLINK("http://www.ncbi.nlm.nih.gov/gene/' .
			$genes[0] . '","' . $genes[0] . '")';
		($genect > 1) && $mgct++;
	}
	my $title = $arprots{$arid}{'titles'};
	my $mism = $qlen - $numid;
	if (($qlen == $slen) && ($numid == $qlen)) {
		$idmatch = 1;
		$idct++;
	}
	my $qpct = sprintf "%.2f", 100 * $numid / $qlen;

	my ($refacc,$allacc,$refup,$allup) = ('') x 4;
	my %upids;
	while ($seqid =~ m/gi\|(\d+)\|(\w+)\|(.*?)\|/g) {
		my $pgi = $1;
		my $db = $2;
		my $pacc = $3;
		if (($db eq 'ref') && ($pacc =~ /^NP_/)) {
			$refacc .= $pacc . ';';
		}
		elsif ($db eq 'sp') {
			my $spid = $pacc;
			$spid =~ s/\..*//;
			$refup .= $spid . ';';
		}
		if ($pacc) {
			$allacc .= $pacc . ';';
		}
		if (exists $prot2up{$pgi}{'up'}) {
			foreach my $upid (keys %{$prot2up{$pgi}{'up'}}) {
				$upids{$upid}++;
			}
		}
	}
	(! $refacc) && ($refacc = $allacc);
	$refacc =~ s/;$//;
	$refup =~ s/;$//;
	(! $refup) && ($refup = join ';', sort keys %upids);

 	print "$arid\t$idmatch\t$qpct\t$mism\t$qlen\t$slen\t$hitct\t",
 		"$genect\t$genestr\t$genelink\t$argp\t$title\t$seqid\t",
		"$refup\t$refacc\n";
}

print STDERR "\n";

my $genemult = 0;
foreach my $geneid (keys %gene2ar) {
	if (scalar keys %{$gene2ar{$geneid}} > 1) {
		$genemult++;
	}
}

printf STDERR "Genes with multiple Arena proteins: %d\n", $genemult;
printf STDERR "Proteins with multiple Arena matches: %d\n",
	scalar keys %argroups;
printf STDERR "Arena proteins with exact matches: %d\n", $idct;
printf STDERR "Arena proteins with multiple genes: %d\n", $mgct;
printf STDERR "Arena proteins with no genes: %d\n", $ngct;
print STDERR "\n";
printf STDERR "Time: %d seconds\n", time - $starttime;
