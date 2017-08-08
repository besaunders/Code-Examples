#!/usr/bin/perl -w
use strict;
use LWP::UserAgent;

my $starttime = time;

my $elink = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi';
my $eargs = 'dbfrom=protein&db=nuccore&id=';

my @testset = @ARGV;

#create HTTP user agent
my $ua = new LWP::UserAgent;
$ua->agent("elink/2.0 " . $ua->agent);

my $pidstr = join ('&id=', @testset);
my $epost = $eargs . $pidstr;

#create HTTP request object
my $req = new HTTP::Request POST => "$elink";
$req->content_type('application/x-www-form-urlencoded');
$req->content("$epost");

#post the HTTP request
my $response = $ua->request($req); 
my $eres = $response->content;

my $readpid = 0;
my $readlink = 0;
my $protid;
print $eres;

foreach my $line (split /\n/, $eres)
  {
   if ($line =~ m#<IdList>#)
     {
      $readpid = 1;
     }
   elsif ($line =~ m#</IdList>#)
     {
      $readpid = 0;
     }
   elsif ($line =~ m#<Link>#)
     {
      $readlink = 1;
     }
   elsif ($line =~ m#</Link>#)
     {
      $readlink = 0;
     }
   elsif ($line =~ m#<Id>(\d+)</Id>#)
     {
      if ($readpid)
        {
         $protid = $1;
        }
      elsif ($readlink)
        {
         my $nucid = $1;
         print STDERR "$protid\t$nucid\n";
        }
     }
   elsif ($line =~ m#</LinkSet>#)
     {
      undef $protid;
     }
  }

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
