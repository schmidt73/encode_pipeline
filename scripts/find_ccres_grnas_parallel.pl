use v5.10;

use File::Basename;
use Data::Dumper;
use feature signatures;
use strict;
use warnings;

use Text::CSV qw( csv );
use Getopt::ArgParse;
use List::Util qw(shuffle);

my $ap = Getopt::ArgParse->new_parser(
    prog => 'Design gRNAs for target cCREs using Guidescan2 in parallel.',
);

$ap->add_args(
    ['--ccres', '-c', required => 1, help => 'BED file containing target cCREs.'],
);

my $ns = $ap->parse_args();
my $ccre_file = $ns->ccres;
my $output_dir = "results/candidates/" . basename($ccre_file, ".V4.bed");
`mkdir -p $output_dir`;

for my $grna_db (glob('../experiment14/organisms/k4/hg38_no_alt/results/split_dbs_scored/*.bam')) {
    my $number = basename($grna_db, ".bam");
    my $r = `bsub -R "rusage[mem=8GB]" "perl scripts/find_ccre_grnas.pl --ccres $ccre_file --grna-db $grna_db > $output_dir/$number.csv"`;
    print $r;
}
