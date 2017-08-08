#!/usr/bin/perl -w
use strict;

my ($bestseqid,$lastarid,$besttitle);
my ($bestslen,$bestpctid,$bestnumid,$lastalen) = (0,0,0,0);
while (<>) {
	chomp;
	my ($arid,$alen,$numid,$pctid,$slen,$seqid,$title) = split /\t/;

	if ((defined $lastarid) && ($lastarid ne $arid)) {
	  my $idmatch = 0;
	  if (($lastalen == $bestslen) && ($lastalen == $bestnumid)) {
			$idmatch = 1;
		}
	  print "$lastarid\t$bestpctid\t$idmatch\t$bestseqid\t$besttitle\n";
		$lastarid = $arid;
		$lastalen = $alen;
		$bestpctid = $pctid;
		$bestseqid = $seqid;
		$besttitle = $title;
		$bestnumid = $numid;
		$bestslen = $slen;
	}

	elsif (($bestnumid < $numid) || (($bestnumid == $numid) &&
		(abs($bestslen - $alen) > abs($slen - $alen)))) {
		$bestpctid = $pctid;
		$bestseqid = $seqid;
		$besttitle = $title;
		$bestnumid = $numid;
		$bestslen = $slen;
		$lastarid = $arid;
		$lastalen = $alen;
	}

}
my $idmatch = 0;
if (($lastalen == $bestslen) && ($lastalen == $bestnumid)) {
	$idmatch = 1;
}
print "$lastarid\t$bestpctid\t$idmatch\t$bestseqid\t$besttitle\n";
