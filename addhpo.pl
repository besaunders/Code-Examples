#!/usr/bin/perl -w
use strict;

my $hpofile = shift @ARGV;
my $idfile = shift @ARGV;

my %hpoterms;
open HPO, $hpofile;
my $hpohead = <HPO>;
while (<HPO>) {
    chomp;
    my ($code,$disorder,$hpoid,$hpdesc) = split /\t/;
    $hpoterms{$code}{$hpoid} = $hpdesc;
}

print "ExternalID\tGleesonID\tFamily\tStatus\tGender\tHPO_id\tHPO_desc\n";
open IDF, $idfile;
my $idhead = <IDF>;
while (<IDF>) {
    chomp;
    my ($extid,$glid,$family,$diseases,$affected,$sex) = split /\t/;
    my %hpoids;
    foreach my $code (split /-/, $diseases) {
        if (exists $hpoterms{$code}) {
            foreach my $hpoid (keys %{$hpoterms{$code}}) {
                $hpoids{$hpoid} = $hpoterms{$code}{$hpoid};
            }
        }
    }
    my $status = 'Unknown';
    if ($affected == 1) {
        $status = 'Index';
    } elsif ($sex eq 'F') {
        $status = 'Mother';
    } elsif ($sex eq 'M') {
        $status = 'Father';
    }
    my $hpoidp = '';
    my $hpodep = '';
    foreach my $hpoid (sort keys %hpoids) {
        $hpoidp .= $hpoid . ', ';
        $hpodep .= $hpoids{$hpoid} . ', ';
    }
    $hpoidp =~ s/, $//;
    $hpodep =~ s/, $//;
    print "$extid\t$glid\t$family\t$status\t$sex\t$hpoidp\t$hpodep\n";
}
