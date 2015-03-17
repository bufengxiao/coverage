#! /usr/bin/perl -W

use Getopt::Long;
#use Data::Dumper;

my @min_reads = ();
my $regions_file = '';
my @bedgraph_files;

GetOptions('regions=s' => \$regions_file,
           'cutoff=i'  => \@min_reads,
           'verbose'   => \$verbose);

if (!@min_reads) { @min_reads = (30,10); }
@min_reads = sort {$b <=> $a} @min_reads;

if (!$regions_file)	{ die "--regions option for BED specifying regions of interest is required!"; }
open REGIONFILE, "$regions_file" or die $!;
while ($line = <REGIONFILE>) {
  chomp($line);
  $x = {};
  ($x->{'chr'},$x->{'start'},$x->{'end'},$x->{'comment'}) = split ' ',$line,4;
  if (!$x->{'comment'}) { $x->{'comment'} = ''; }
  push(@regions, $x);
}
close REGIONFILE;

#$r_index = 0;
foreach $sample (@ARGV) {
  open SAMPLEFILE, "$sample" or die $!;
  while ($line = <SAMPLEFILE>) {
    chomp($line);
    $tmp = {};
    ($tmp->{'chr'},$tmp->{'start'},$tmp->{'end'},$dp,undef) = split ' ',$line,5;
    if ($dp >= $min_reads[0] ) { next; }
    foreach $region (@regions) {
      if (overlap($tmp,$region)) {
        if ($tmp->{'start'} < $region->{'start'}) { $tmp->{'start'} = $region->{'start'}; }
        if ($tmp->{'end'} > $region->{'end'}) { $tmp->{'end'} = $region->{'end'}; }
        if (!exists $region->{$sample}) {
          $region->{$sample} = ();
          for ($i=0; $i<(@min_reads); $i++) { @{$region->{$sample}[$i]} = (); }
          push(@{$region->{$sample}[0]},{%{$tmp}});
        }
        elsif (adjacent($tmp,${$region->{$sample}[0]}[-1])) {
          ${$region->{$sample}[0]}[-1]->{'end'} = $tmp->{'end'};
        }
        else {
          push(@{$region->{$sample}[0]},{%{$tmp}});
        }
#print Dumper($region);

        for ($i=1; $i<(@min_reads); $i++) {
          $threshold = $min_reads[$i];
          if ($dp < $threshold) {
            if (!@{$region->{$sample}[$i]}) {
              push(@{$region->{$sample}[$i]},{%{$tmp}});
            }
            elsif (adjacent($tmp,${$region->{$sample}[$i]}[-1])) {
              ${$region->{$sample}[$i]}[-1]->{'end'} = $tmp->{'end'};
            }
            else {
              push(@{$region->{$sample}[$i]},{%{$tmp}});
            }
          }
        }
#print Dumper($region);
      }
    }
  }
  close SAMPLEFILE;
}

print "CHR\tstart\tend\tregion size\t";
foreach $sample (@ARGV) {
  foreach $threshold (@min_reads) {
    print "$sample\@$threshold\t";
  }
}
foreach $threshold (@min_reads) {
  print "Combined low cov regions \@$threshold\tlow cov bases \@$threshold\tfraction low coverage \@$threshold\t";
}
print "\n";
@allregions = ();
foreach $region (@regions) {
  $placeholder = "$region->{'chr'}\t$region->{'start'}\t$region->{'end'}\t".($region->{'end'}-$region->{'start'})."\t";
  for ($i=0; $i<(@min_reads); $i++) { @{$allregions[$i]} = (); }
  foreach $sample (@ARGV) {
    for ($i=0; $i<(@min_reads); $i++) {
      if (exists $region->{$sample}) {
        foreach $subregion (@{$region->{$sample}[$i]}) {
          push(@{$allregions[$i]},$subregion);
          $placeholder .= "$subregion->{'chr'}:$subregion->{'start'}-$subregion->{'end'};";
        }
      }
      $placeholder .= "\t";
    }
  }

  if (@{$allregions[0]}) {
    print $placeholder;

    for ($i=0; $i<(@min_reads); $i++) {
      if (@{$allregions[$i]}) {
        @overall = mergeregions(sort { $a->{'start'} <=> $b->{'start'} } @{$allregions[$i]});
        $lowcovbases = 0;
        foreach $subregion (@overall) {
          print "$subregion->{'chr'}:$subregion->{'start'}-$subregion->{'end'};";
          $lowcovbases += $subregion->{'end'}-$subregion->{'start'};
        }
        print "\t$lowcovbases";
        $frac_low_cov = $lowcovbases / ($region->{'end'}-$region->{'start'});
        print "\t$frac_low_cov\t";
      }
      else {
        print "\t0\t0\t";
      }
    }

    print "$region->{'comment'}\n";
  }
  elsif ($verbose) {
    print $placeholder;
    for ($i=0; $i<(@min_reads); $i++) { print "\t0\t0\t"; }
    print "$region->{'comment'}\n";
  }
}

exit;

sub overlap {
  ($a, $b) = @_;
  $retval = 0;

  # You'd think the comparisons would be >= to include the case where a.end == b.start
  # But we don't actually want those since bed uses half-open coordinates
  if ($a->{'chr'} eq $b->{'chr'}) {
    if ($a->{'end'} > $b->{'start'}) {
      if ($a->{'start'} < $b->{'end'} ) {
        $retval = 1;
      }
    }
  }

  return $retval;
}

sub adjacent {
  ($a, $b) = @_;
  $retval = 0;

  # Only a.start == b.end are technically adjacent
  # a.start == b.end+1 has one base that is not included in either interval
  # However, it is inconvenient to have tiny 1 base breaks in the low-coverage intervals
  if ($a->{'chr'} eq $b->{'chr'}) {
    if ($a->{'start'} == $b->{'end'}) { $retval = 1; }
    if ($a->{'start'} == $b->{'end'}+1) { $retval = 1; }
  }

  return $retval;
}

sub mergeregions {
  @regionlist = @_;
  @retarray = ();

  push(@retarray,shift(@regionlist));

  foreach $region (@regionlist) {
    if (overlap($region,$retarray[-1]) || adjacent($region,$retarray[-1])) {
      if ($region->{'end'} > $retarray[-1]->{'end'}) {
        $retarray[-1]->{'end'} = $region->{'end'};
      }
    }
    else {
      push(@retarray,$region);
    }
  }

  return @retarray;
}
