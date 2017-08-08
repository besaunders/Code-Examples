#!/usr/bin/perl -w
use strict;

my $seqid;
my %seqhash;
while (<>)
	{ 
		chomp;
		if (/^>Arena_(.*)/)
			{
				$seqid = $1;
			}
		else
		  {
				my $seq = $_;
				my $org = 'misc';
				if ($seqid =~ /^\d+\w+/)
					{
					 $org = 'human';
					}
				elsif ($seqid =~ /^m\d+\w+/)
					{
					 $org = 'mouse';
					}
				elsif ($seqid =~ /Celegans/)
					{
					 $org = 'celegans';
					}
				elsif ($seqid =~ /drosophila/)
					{
					 $org = 'drosophila';
					}
				elsif (($seqid =~ /dog/) || ($seqid =~ /canine/) || 
					($seqid =~ /^d\d+\w+/))
					{
					 $org = 'dog';
					}
				elsif ($seqid =~ /guineapig/)
					{
					 $org = 'guineapig';
					}
				elsif ($seqid =~ /hamster/)
					{
					 $org = 'hamster';
					}
				elsif (($seqid =~ /monkey/) || ($seqid =~ /^mk\d+\w+/) ||
					($seqid =~ /^mnk\d+\w+/))
					{
					 $org = 'monkey';
					}
				elsif ($seqid =~ /^p\d+\w+/)
					{
					 $org = 'pig';
					}
				elsif ($seqid =~ /^r\d+\w+/)
					{
					 $org = 'rat';
					}
				elsif ($seqid =~ /rabbit/)
					{
					 $org = 'rabbit';
					}
				elsif ($seqid =~ /xenopus/)
					{
					 $org = 'xenopus';
					}
				elsif ($seqid =~ /zebrafish/)
					{
					 $org = 'zebrafish';
					}
				$seqhash{$org}{$seqid} = $seq;
			}
	}

print STDERR "Sequence count\n";
foreach my $org (sort {$a cmp $b} keys %seqhash)
	{
	 print STDERR "$org: ", scalar keys %{$seqhash{$org}}, "\n";
	 my $fname = 'Arena_' . $org . '.aa';
	 open FH, ">$fname";
	 foreach my $seqid (sort {$a cmp $b} keys %{$seqhash{$org}})
		{
		 print FH ">Arena_", $seqid, "\n", $seqhash{$org}{$seqid}, "\n";
		}
	 close FH;
	}

