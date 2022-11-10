use v5.10;

use Data::Dumper;
use feature signatures;
use strict;
use warnings;

use IntervalTree;
use Text::CSV qw( csv );
use Getopt::ArgParse;
use List::Util qw(shuffle);

my $ap = Getopt::ArgParse->new_parser(
    prog => 'Design gRNAs for target cCREs using Guidescan2.',
);

$ap->add_args(
    ['--ccres', '-c', required => 1, help => 'BED file containing target cCREs.'],
    ['--grna-db', required => 1, help => 'Guidescan2 gRNA database.'],
);

my $ns = $ap->parse_args();

my @ccres = ();
my $grna_db = $ns->grna_db;

open(my $fh, '<', $ns->ccres) or die "Cannot open $ns->ccres!";
while(<$fh>) {
    chomp;
    my @fields = split;
    push(@ccres, \@fields);
}

@ccres = shuffle @ccres;
my %ccre_interval_trees = (); # chromosome indexed interval trees
for my $ccre (@ccres) {
    my $chrm = @$ccre[0];
    my $start = int(@$ccre[1]);
    my $end = int(@$ccre[2]);

    if (!(exists $ccre_interval_trees{$chrm})) {
        $ccre_interval_trees{$chrm} = IntervalTree->new();
    }

    my $it = $ccre_interval_trees{$chrm};
    $it->insert($start, $end, $ccre);
}

sub process_guide($ccre, $grna) {
    my @fields = split /\s+/, $grna;

    # parse grna features
    my $sequence = $fields[9];
    my ($grna_chrm, $grna_start, $grna_strand) = $fields[0] =~ /(\w+):(\d+):([+\-])/;
    my ($distance1_ots) = $fields[12] =~ /k\d+:i:(\d+)/;
    my ($distance2_ots) = $fields[13] =~ /k\d+:i:(\d+)/;
    my ($distance3_ots) = $fields[14] =~ /k\d+:i:(\d+)/;
    my ($distance4_ots) = $fields[15] =~ /k\d+:i:(\d+)/;
    my ($cutting_efficiency) = $fields[17] =~ /ds:f:(.*)/;
    my ($specificity) = $fields[18] =~ /cs:f:(.*)/;

    # parse ccre features
    my ($ccre_chrm, $ccre_start, $ccre_end, $ccre_name_1, $ccre_name_2, $ccre_type) = (@$ccre[0], @$ccre[1], @$ccre[2], @$ccre[3], @$ccre[4], @$ccre[5]);

    # pretty print
    say("$ccre_name_1,$ccre_type,$ccre_chrm,$ccre_start,$ccre_end,$sequence,"
        . "$grna_chrm,$grna_start,$grna_strand,$distance1_ots,$distance2_ots,"
        . "$distance3_ots,$distance4_ots,$cutting_efficiency,$specificity");
}

say("ccre_id,ccre_type,ccre_chrm,ccre_start,ccre_end,grna_sequence,grna_chrm,grna_start,grna_strand,"
    . "distance_0_matches,distance_1_matches,distance_2_matches,distance_3_matches,distance_4_matches,"
    . "cutting_efficiency,specificity");

my $pid = open(GRNA_DB, "samtools view $grna_db |") or die "Couldn't fork: $!\n";
my $counter = 0;
while (<GRNA_DB>) {
    my @fields = split /\s+/, $_;
    my ($grna_chrm, $grna_start, $grna_strand) = $fields[0] =~ /(\w+):(\d+):([+\-])/;
    $grna_start = int($grna_start);
    my $matching_ccres = ($ccre_interval_trees{$grna_chrm})->find($grna_start, $grna_start + 1);

    for my $ccre (@$matching_ccres) {
        process_guide $ccre, $_;
    }
}
