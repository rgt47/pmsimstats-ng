#!/usr/bin/perl
# Notation linter for the pmsimstats-ng compendium.
# Checks the NOTATION.md rules against each manuscript source.
# Prose only: fenced R chunks are skipped entirely.
use strict;
use warnings;

my %count;

for my $file (@ARGV) {
  open my $fh, '<', $file or die "$file: $!";
  my ($in_chunk, $ln) = (0, 0);
  my @hits;
  while (my $line = <$fh>) {
    $ln++;
    if ($line =~ /^```\{/)              { $in_chunk = 1; next }
    elsif ($line =~ /^```\s*$/)         { $in_chunk = 0; next }
    elsif ($line =~ /\\begin\{verbatim\}/) { $in_chunk = 1; next }
    elsif ($line =~ /\\end\{verbatim\}/)   { $in_chunk = 0; next }
    next if $in_chunk;

    # Strip protected spans: \texttt{...}, `...`, and \item labels of
    # code listings. What remains is reader-facing prose and math.
    my $p = $line;
    # one level of brace nesting, so \texttt{a \textasciitilde{} b} works
    $p =~ s{\\texttt\{(?:[^{}]|\{[^{}]*\})*\}}{}g;
    $p =~ s{\\verb\{(?:[^{}]|\{[^{}]*\})*\}}{}g;
    $p =~ s/`[^`]*`//g;
    $p =~ s/\\verb\|[^|]*\|//g;

    # N1: bare Dbc outside code context (should be $D_{bc}$)
    push @hits, "$ln bare-Dbc" if $p =~ /(?<![\w\\{])Dbc(?![\w}])/;
    # N1: Sx used as a symbol outside code context (should be Y)
    push @hits, "$ln bare-Sx"  if $p =~ /(?<![\w\\{])Sx(?![\w}])/;
    # N6: legacy specification labels
    push @hits, "$ln legacy-spec-label" if $p =~ /(?<![\w])A[123](?![\w])/;
    # N3: half-life quoted in days without a units statement is
    # checked separately; here just record day-denominated half-lives
    push @hits, "$ln halflife-days" if $p =~ /t_\{1\/2\}[^.]{0,40}\bdays\b/;
    # N7: N is the total across paths, so a value of N may never be
    # given directly as a per-path count. Naming a per-path number
    # alongside the total ("$N = 140$ (70 per path)") is compliant and
    # must not fire, so only an N whose own value is immediately
    # qualified as per-path counts as a violation.
    push @hits, "$ln N-per-path"
      if $p =~ /\$?N\$?\s*(?:=|\\in)\s*\$?[\d,.{}\\\s]{1,12}\$?\s*
                (?:participants?\s+|subjects?\s+)?
                per\s+(?:randomization\s+)?path\b/x;
  }
  close $fh;
  # A day-denominated half-life is compliant if the file carries an
  # explicit units paragraph giving the day/week conversion.
  my $has_units = do {
    open my $g, '<', $file or die $!;
    local $/; my $all = <$g>; close $g;
    $all =~ /\*\*Units\.\*\*/ ? 1 : 0;
  };
  @hits = grep { !/halflife-days/ } @hits if $has_units;

  my %by;
  for (@hits) { my ($l, $k) = split ' '; push @{$by{$k}}, $l }
  my @out;
  for my $k (sort keys %by) {
    my @l = @{$by{$k}};
    push @out, sprintf('%s x%d (line%s %s)', $k, scalar @l,
      (@l > 1 ? 's' : ''), join(',', @l[0 .. ($#l > 2 ? 2 : $#l)]) . ($#l > 2 ? ',...' : ''));
    $count{$k} += scalar @l;
  }
  printf "%-54s %s\n", $file, (@out ? join('; ', @out) : 'clean');
}

print "\nTOTALS: ", join(', ', map { "$_=$count{$_}" } sort keys %count), "\n"
  if %count;
