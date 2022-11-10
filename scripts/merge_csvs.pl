use strict; use warnings; use autodie;
my $files = $#ARGV; # get number of files - 1
while (my $file = shift @ARGV) {
  open my $fh, "<", $file;
  <$fh> unless $files == @ARGV; # discard header unless first file
  print while <$fh>; # output the rest
}
