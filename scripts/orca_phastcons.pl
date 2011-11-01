#!/usr/local/bin/perl -w

=head1 NAME

orca_phastcons.pl

=head1 SYNOPSIS

 orca_phastcons.pl -s seq_file -udb ucsc_db -track ucsc_track
        [-pos position]
        [-ex exons] [-f features]
        [-cc cons_cutoff] [-cl min_cr_len]
        [-jdb db_name] [-c collection] [-tax tax_groups] [-species species]
        [-ic min_ic]
        ([-tn tf_name] | [-pfmf pfm_file])
        [-th tfbs_threshold]
        [-start start_pos] [-end end_pos] [-nfos]
        [-hs spacer]
        [-tf tf_file] [-crf cr_file] [-csf cs_file]
        [-g out_graph] [-flip]
        [-utf ucsc_track_file]

=head1 ARGUMENTS
 
  All argument switches may be abbreviated if unique.
  -s seq_file       = FastA file containing the sequence to be analyzed.
  -udb ucsc_db      = UCSC database name, e.g. "hg18"
  -track ucsc_track = UCSC track name, e.g. "phastCons46wayPlacental"
  -pos position     = Chromosomal coordinates of input sequence in
                      format: chrN:start-end
  -ex exons         = Optional GFF file containing the exons positions
                      on the first sequence. If provided these will be
                      filtered out and excluded from the conserved
                      regions searched for TFBSs
  -f features       = Optional GFF file containing other features to plot
                      (i.e. CpG islands)
                      on the first sequence. If provided these will be
  -th threshold     = Optional. Minimum TFBS PWM score to report.
                      Specify as a relative score (i.e. 80%) or as
                      an absolute score (i.e. 11.2).
                      (default 80%)
  -cc cons_cutoff   = Optional. Min. conservation of regions in which
                      to search for TF sites. Specify as a percentage
                      in the range 0% - 100% or 0 - 1.
                      (default 0.7)
  -cl min_cr_len    = Optional. Min. conserved region length to
                      report (default 20)
  -jdb db_name      = Optional. Name of JASPAR database to use.
                      (default = JASPAR_2010)
  -c collection     = Optional. Name of JASPAR collection to use.
                      (default = CORE)
  -tax tax_groups   = Optional. Limit matrix profiles to those with these
                      taxonomic supergroups.  To specify more than a single
                      tax group, you can use multiple instances of the -tax
                      option or a comma separated string of values, e.g:
                        -tax vertebrates -tax insects
                      OR
                        -tax vertebrates,insects
  -species species  = Optional. Limit matrix profiles to those associated
                      with these species. To specify more than a single
                      species, you can use multiple instances of the
                      -species option or a comma separated string of values,
                      e.g:
                        -species human -tax mouse
                      OR
                        -species human,mouse
  -ic min_ic        = Optional. Minimum information content
                      (specificity) of TFBS PWMs to use (default 10)
  -tn tf_name       = Optional. Name of a specific TF to use in
                      the analysis.
  -pfmf file_name   = Optional. Name of a file defining one or more
                      PFM matrix to use in the analysis.
  -start start      = Optional. Starting position on sequence 1 to
                      scan for TFBSs
  -end end          = Optional. Ending position on sequence 1 to
                      scan for TFBSs
  -nfos             = Optional. Do not filter overlapping sites. Sites of
                      the same TF which overlap are filtered so that only
                      the highest scoring is reported. If specified, report
                      all sites.
  -hs spacer        = Optional. Half-site spacer. Single profile only.
                      The profile is considered to be a half site. The
                      script attempts to find pairs of sites within (<=)
                      this specified numbers of nucleotides.
  -g graph_file     = Optional name of output graphics file in PNG
                      format
  -flip             = Optional. Flip graph left to right - useful for
                      visualizing regions on the -ve strand as 5' to 3'.
  -tf tf_file       = If specified, write the putative TF sites to
                      this file, else write to STDOUT
  -crf cr_file      = If specified, write the conserved regions to
                      this file
  -csf cs_file      = If specified, write the conserved sub-sequences
                      to this file
  -utf track_file   = If specified, create a custom UCSC browser track of
                      the analysis

=head1 DESCRIPTION

Perform ConSite (http://mordor.cgb.ki.se/cgi-bin/CONSITE/consite) style
analysis. Given two sequence files, align them and find the conserved regions.
Search for TFBSs within the conserved regions (optionally restricted to the
start and end position provided). Output the putative TF sites. Optionally
output a graph consisting of the conservation profile, exons (if provided),
conserved regions and TFBS positions in PNG format.

If both masked and unmasked sequences are provided, the masked sequences are
aligned and used to compute conserved regions and the unmasked sequences are
scanned for putative TFBSs. If only masked or unmasked sequences are provided,
they are used for both conserved region computation and TFBS analysis.

=head1 REQUIRMENTS

This script has the following requirements:

  Bioperl   - at least version 1.2 or better.
  TFBS      - (http://forkhead.cgb.ki.se/TFBS/) a set of perl
              modules for the detection of TF sites (also
              provides an interface to the JASPAR database).

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia
  
  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use lib '/devel/ORCAtk/lib';

use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use TFBS::DB::JASPAR5;
use Bio::SeqIO;
use Bio::Seq;
use ORCA::Analysis::PhastCons;
use ORCA::Graphics::PhastCons;

use constant DEBUG => 0;

#
# Default minimum TFBS PWM score to report for either sequence.
# Overidden by -th argument.
#
use constant TFBS_THRESHOLD => "80%";

#
# Default absolute minimum percent identity of conserved regions to search for
# TFBSs.
# Overidden by the -cc argument.
#
use constant MIN_CONSERVATION => 0.7;

#
# Default JASPAR database
#
use constant JASPAR_DB => 'JASPAR_2010';

use constant COLLECTION => 'CORE';
use constant TAX_GROUPS => ('vertebrates');

#
# Default minimum information content (specificity) of TFBS profiles from the
# JASPAR database to search against the sequences.
# Overidden by -ic argument.
#
use constant MIN_IC     => 8;

#
# Default minimum conserved region size.
# Overidden by -cl argument.
#
use constant MIN_CR_LENGTH => 20;

#
# Mininimum amount in bp, by which a TFBS needs to overlap a conserved region
# in order to be reported by the analysis.
#
use constant MIN_TFBS_CR_OVERLAP => 1;

#
# JASPAR database connection parameters.
#
use constant JASPAR_DB_HOST => "vm5.cmmt.ubc.ca";
use constant JASPAR_DB_USER => "jaspar_r";
use constant JASPAR_DB_PASS => "";

my $seq_file;
my $ucsc_db;
my $track;
my $position;
#my $phc_file;
my $exon_file;
my $feat_file;
my $threshold;
my $min_conservation;
my $min_cr_len;
my $jaspar_db;
my $min_ic;
my $collection;
my @tax_groups;
my @species;
my $tf_name;
my $pfm_file;
my $start;
my $end;
my $nfos;
my $hs_spacer;
my $tfbs_file;
my $cr_file;
my $cs_file;
my $graph_file;
my $flip;
my $ucsc_track_file;
GetOptions(
    's=s'       => \$seq_file,
    'udb=s'     => \$ucsc_db,
    'track=s'   => \$track,
    'pos=s'     => \$position,
    #'phcf=s'   => \$phc_file,
    'ex=s'      => \$exon_file,
    'f=s'       => \$feat_file,
    'th=s'      => \$threshold,
    'cc=s'      => \$min_conservation,
    'cl=i'      => \$min_cr_len,
    'jdb=s'     => \$jaspar_db,
    'c=s'       => \$collection,
    'tax=s'     => \@tax_groups,
    'species=s' => \@species,
    'ic=i'      => \$min_ic,
    'tn=s'      => \$tf_name,
    'pfmf=s'    => \$pfm_file,
    'start=i'   => \$start,
    'end=i'     => \$end,
    'tf=s'      => \$tfbs_file,
    'crf=s'     => \$cr_file,
    'csf=s'     => \$cs_file,
    'g=s'       => \$graph_file,
    'flip'      => \$flip,
    'utf=s'     => \$ucsc_track_file,
    'nfos'      => \$nfos,
    'hs=i'      => \$hs_spacer
);

if (!$seq_file) {
    pod2usage(-verbose => 1);
}

unless ($jaspar_db) {
    $jaspar_db = JASPAR_DB;
}

unless ($min_conservation) {
    $min_conservation = MIN_CONSERVATION;
}

if ($min_conservation =~ /(.+)%/) {
    $min_conservation = $1 / 100;
}

unless ($min_ic) {
    $min_ic = MIN_IC;
}

unless ($threshold) {
    $threshold = TFBS_THRESHOLD;
}

unless ($collection) {
    $collection = COLLECTION;
}

if (@tax_groups) {
    @tax_groups = split(/\s*,\s*/, join(',', @tax_groups));
} else {
    @tax_groups = TAX_GROUPS;
}

if (@species) {
    @species = split(/\s*,\s*/, join(',', @species));
}

unless ($min_cr_len) {
    $min_cr_len = MIN_CR_LENGTH;
}

if (!$ucsc_db) {
    pod2usage(-verbose => 1, -msg => "No UCSC database name specified!\n");
}

if (!$track) {
    pod2usage(-verbose => 1, -msg => "No UCSC track name specified!\n");
}

#if (!$position) {
#    pod2usage(-verbose => 1, -msg => "No chromosomal position specified!\n");
#}

my $fos = 1;
if ($nfos) {
    $fos = 0;
}

#
# Check any files specfied exist first
#
if (!-f $seq_file) {
    die "Specified sequence file $seq_file not found!\n";
}

#if (!-f $phc_file) {
#    die "Specified phastCons score file $phc_file not found!\n";
#}

if ($exon_file && !-f $exon_file) {
    die "Specified exon file $exon_file not found!\n";
}

if ($feat_file && !-f $feat_file) {
    die "Specified feature file $feat_file not found!\n";
}

if (defined $hs_spacer && !($tf_name || $pfm_file)) {
    die "Half-site analysis requires a specific TF specified by either a"
        . " name (in JASPAR) or a PFM file\n";
}

my ($seq_chr, $seq_start, $seq_end);
if ($position) {
    $position =~ s/,//g;
    if ($position =~ /chr(\S+):(\d+)-(\d+)/) {
        $seq_chr   = $1;
        $seq_start = $2;
        $seq_end   = $3;
    } else {
        pod2usage(
            -verbose => 1,
            -msg     => "Bad format specified for position ($position)!\n"
        );
    }
}

#
# Read sequences
#
print "Reading sequence from $seq_file\n";
my $seq = read_seq($seq_file);
die "Error reading sequence\n" if !$seq;

my $display_id = $seq->display_id;
if ($display_id && $display_id =~ /chr(\S+):(\d+)-(\d+)/) {
    if ($position) {
        if ($1 ne $seq_chr || $2 != $seq_start || $3 != $seq_end) {
            die "Chromosomal position in sequence fasta header $display_id"
                . " does not match user entered position $position\n";
        }
    } else {
        $seq_chr   = $1;
        $seq_start = $2;
        $seq_end   = $3;
    }
} elsif (!$seq_chr) {
    die "Could not determine chromosomal position of sequence - please"
        . " explicitly enter a position for this sequence\n";
}

my $lseq = Bio::LocatableSeq->new(
        -primary_id     => $seq->primary_id,
        -display_id     => $seq->display_id,
        -seq            => $seq->seq,
        -start          => $seq_start,
        -end            => $seq_end,
        -strand         => 1
);

#
# Read phastCons scores
#
#print "Reading phastCons scores from $phc_file\n";
#my ($profile, $scores) = read_phastCons_scores($phc_file, $flip);
#die "Error reading phastCons scores\n" if !$scores;

#
# Optionally get exons
#
my $exons;
if ($exon_file) {
    print "Reading exons from $exon_file\n";
    $exons = read_GFF($exon_file);
    die "Could not read exons\n" if !$exons;
}

#
# Optionally get other features (e.g. CpGs)
#
my $features;
if ($feat_file) {
    print "Reading features from $feat_file\n";
    $features = read_GFF($feat_file);
    die "Could not read features\n" if !$features;
}

#
# Connect to JASPAR database
#
my $tfbs_db =
    TFBS::DB::JASPAR5->connect("dbi:mysql:$jaspar_db:" . JASPAR_DB_HOST,
    JASPAR_DB_USER, JASPAR_DB_PASS);
die "Error connecting to $jaspar_db\n" if !defined $tfbs_db;

my $matrix_set;
if ($tfbs_file) {
    if ($tf_name) {
        my $pfm = $tfbs_db->get_Matrix_by_name($tf_name);
        die "Error retrieving TFBS profile from $jaspar_db\n" if !$pfm;
        my $pwm = $pfm->to_PWM;
        die "Error converting PFM to PWM\n" if !$pwm;
        $matrix_set = TFBS::MatrixSet->new();
        die "Error creating TFBS matrix set\n" if !$matrix_set;
        $matrix_set->add_Matrix($pwm);
    } elsif ($pfm_file) {
        my $pfms = read_PFMs($pfm_file);
        die "Error reading TFBS profile(s) from PFM file $pfm_file\n"
            if !$pfms;
        $matrix_set = TFBS::MatrixSet->new();
        foreach my $pfm (@$pfms) {
            my $pwm = $pfm->to_PWM;
            die "Error converting PFM to PWM\n" if !$pwm;
            $matrix_set->add_Matrix($pwm);
        }

        if (defined $hs_spacer && $matrix_set->size > 1) {
            die "Half-site analysis requires a specific TF\n";
        }
    } else {
        my %matrix_args = (-matrixtype => 'PWM');

        $matrix_args{-min_ic}       = $min_ic if defined $min_ic;
        $matrix_args{-collection}   = $collection if defined $collection;
        $matrix_args{-tax_group}    = \@tax_groups if @tax_groups;
        $matrix_args{-species}      = \@species if @species;

        $matrix_set = $tfbs_db->get_MatrixSet(%matrix_args);
        die "Error retrieving TFBS matrix set from $jaspar_db\n"
            if !$matrix_set;
    }
}

#
# Perform phastCons analysis
#
my $phca = ORCA::Analysis::PhastCons->new(
    -seq   => $lseq,     # Use the LocatableSeq
    -db    => $ucsc_db,
    -track => $track,
    -exons => $exons,
    -chr   => $seq_chr,
    -start => $seq_start,
    -end   => $seq_end
);

my $phastCons_profile = $phca->compute_conservation_profile();

#
# Compute conserved regions
#
my %ccr_params;
$ccr_params{-min_conservation} = $min_conservation if $min_conservation;
$ccr_params{-filter_exons}     = 1                 if $exons;
$ccr_params{-min_length}       = $min_cr_len       if $min_cr_len;

my $conserved_regions = $phca->compute_conserved_regions(%ccr_params);

if ($conserved_regions) {
    if ($cr_file) {
        print "Writing conserved regions file $cr_file\n";
        write_conserved_regions_report(
            $cr_file, $min_conservation, $phca->conserved_regions
        );
    }

    if ($cs_file) {
        print "Writing conserved sub-sequences file $cs_file\n";
        write_conserved_subsequences(
            $cs_file, $phca->compute_conserved_subsequences
        );
    }
} else {
    print "No conserved regions found\n";
}

#
# Compute conserved TFBSs
#
my $tfbss;
if ($conserved_regions && $tfbs_file && $matrix_set && $matrix_set->size > 0) {
    print "Searching for conserved TFBSs...\n";

    $tfbss = $phca->compute_conserved_tfbss(
        -matrix_set               => $matrix_set,
        -min_tfbs_score           => $threshold,
        -min_tfbs_cr_overlap      => MIN_TFBS_CR_OVERLAP,
        -filter_overlapping_tfbss => $fos,
        -start                    => $start,
        -end                      => $end
    );
    #print "Conserved seq TFBS sites: " . $tfbss->size . "\n";

    if ($tfbss) {
        if ($hs_spacer) {
            my $hs_pairs = pair_half_sites($tfbss, $hs_spacer);
            write_half_sites($tfbs_file, $hs_pairs);
        } else {
            print "Writing putative TFBS sites...\n";
            write_conserved_tfbss($tfbs_file, $tfbss);
        }
    } else {
        print "No conserved TFBSs found\n";
    }
}

#
# Optionally plot results
#
if ($graph_file) {
    print "Plotting results...\n";

    #
    # Graph results.
    #
    my $graph = ORCA::Graphics::PhastCons->new(
        -analysis       => $phca,
        -other_features => $features,
        -flip           => $flip
    );

    if ($graph && $graph->{-gd_image}) {
        open PNG, ">$graph_file"
            || die "Error opening $graph_file for writing\n";
        print PNG $graph->{-gd_image}->png;
        close PNG;
    } else {
        carp "ERROR: creating phastCons graph\n";
    }
}

#   
# Optionally create a UCSC custom browser track file
#   
if ($ucsc_track_file) {
    my $track = $phca->ucsc_track;
    
    if ($track) {
        open(UTF, ">$ucsc_track_file")
            || die "Error opening $ucsc_track_file for writing\n";
        print UTF "$track\n";
        close UTF;  
    } else {
        carp "ERROR: creating UCSC track file\n";
    }
}

exit;

sub read_seq
{
    my ($seq_file) = @_;

    my $seq_in;
    $seq_in = Bio::SeqIO->new(
        -file   => $seq_file,
        -format => 'fasta'
    );
    die "SeqIO error; file: $seq_file\n" if !defined $seq_in;

    my $seq_obj = $seq_in->next_seq;
    die "Error reading sequence\n" if !$seq_obj;

    $seq_in->close;

    return $seq_obj;
}

sub read_GFF
{
    my ($file) = @_;

    my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 2);

    my @feats;
    if ($gffio) {
        while (my $feature = $gffio->next_feature) {
            push @feats, $feature;
        }
    } else {
        die "Error opening GFF file $file\n";
    }

    return @feats ? \@feats : undef;
}

sub pair_half_sites
{
    my ($tfbss, $spacer) = @_;

    print "Pairing half-sites...\n";

    my @hs_pairs;

    my $num_sites = scalar @$tfbss;

    # build hash of unique sites (based on strand and start position)
    my %unique_flag;
    my @unique_sites;
    foreach my $site (@$tfbss) {
        my $site_key = $site->strand . $site->start;
        if (!$unique_flag{$site_key}) {
            $unique_flag{$site_key} = 1;
            push @unique_sites, $site;
        }
    }

    foreach my $site1 (@unique_sites) {
        foreach my $site2 (@unique_sites) {
            next if $site2->start <= $site1->end;
            last if $site2->start > $site1->end + $spacer + 1;

            push @hs_pairs, {
                -site1 => $site1,
                -site2 => $site2
            };
        }
    }

    return @hs_pairs ? \@hs_pairs : undef;
}

sub write_conserved_regions_report
{
    my ($file, $min_conservation, $crs) = @_;

    if (!$crs || !@$crs) {
        printf STDERR "No conserved regions in report\n";
        return;
    }

    if ($file) {
        if (!open(CRS, ">$file")) {
            warn "Could not create conserved regions report output file $file"
                . " - $!\n";
            return;
        }
    } else {
        open(CRS, ">-");
    }

    printf CRS "Minimum conservation: %0.2f\n", $min_conservation
        if defined $min_conservation;
    print CRS "\n";
    foreach my $cr (@$crs) {
        printf CRS "%7d %7d %7d %.3f\n",
            $cr->start,
            $cr->end,
            $cr->end - $cr->start + 1,
            $cr->score;
    }

    close CRS;
}

sub write_conserved_subsequences
{
    my ($file, $subseqs) = @_;

    if (!$subseqs) {
        printf STDERR "No conserved sub-sequences to report\n";
        return;
    }

    if ($file) {
        if (!open(CSS, ">$file")) {
            warn "Could not create conserved regions report output file $file"
                . " - $!\n";
            return;
        }
    } else {
        open(CSS, ">-");
    }

    foreach my $css (@$subseqs) {
        printf CSS ">%s\n%s\n", $css->display_id, $css->seq;
    }

    close CSS;
}

sub write_conserved_tfbss
{
    my ($file, $tfbss) = @_;

    return if !$tfbss;

    if ($file) {
        if (!open(TFBS, ">$file")) {
            warn "Could not create TFBS output file $file - $!\n";
            return;
        }
    } else {
        open(TFBS, ">-");
    }

    print TFBS "TF Name\tStart\tEnd\tStrand\tScore\t% Score\tConservation"
        . "\tBound seq.\n";
    my $start = $phca->start;
    foreach my $site (@$tfbss) {
        my $strand;
        if ($site->strand eq '+' || $site->strand == 1) {
            $strand = '+';
        } elsif ($site->strand eq '-' || $site->strand == -1) {
            $strand = '-';
        }

        printf TFBS "%s\t%d\t%d\t%s\t%.3f\t%.1f%%\t%.1f%%\t%s\n",
            $site->pattern->name,
            # TFBSs always reported in 1-based coords even if the 
            # PhastCons analysis was in chromosmal coords
            $site->start + $start - 1,
            $site->end + $start - 1,
            $strand,
            $site->score,
            $site->rel_score * 100,
            ($site->get_tag_values('conservation'))[0] * 100,
            $site->seq->seq;
    }

    close TFBS;
}

sub write_half_sites
{
    my ($file, $hs_pairs) = @_;

    if ($file) {
        if (!open(TFBS, ">$file")) {
            warn "Could not create half-site output file $file - $!\n";
            return;
        }
    } else {
        open(TFBS, ">-");
    }

    foreach my $hs_pair (@$hs_pairs) {
        my $hs1 = $hs_pair->{'-site1'};
        my $hs2 = $hs_pair->{'-site2'};

        my $strand1;
        my $strand2;
        if ($hs1->strand eq '+' || $hs1->strand == 1) {
            $strand1 = '+';
        } elsif ($hs1->strand eq '-' || $hs1->strand == -1) {
            $strand1 = '-';
        }
        if ($hs2->strand eq '+' || $hs2->strand == 1) {
            $strand2 = '+';
        } elsif ($hs2->strand eq '-' || $hs2->strand == -1) {
            $strand2 = '-';
        }

        printf TFBS "%-16s %6d %6d %s %5.1f %5.1f %s\n",
            $hs1->pattern->name,
            $hs1->start,
            $hs1->end,
            $strand1,
            $hs1->rel_score * 100,
            ($hs1->get_tag_values('conservation'))[0] * 100,
            $hs1->seq->seq;

        printf TFBS "%-16s %6d %6d %s %5.1f %5.1f %s\n\n",
            $hs2->pattern->name,
            $hs2->start,
            $hs2->end,
            $strand2,
            $hs2->rel_score * 100,
            ($hs2->get_tag_values('conservation'))[0] * 100,
            $hs2->seq->seq;
    }

    close TFBS;
}

sub read_PFMs
{
    my ($file) = @_;

    open(FH, $file) || die "Error opening PFM file $file\n";

    my @pfms;
    my $name          = '';
    my $matrix_string = '';
    my $line_count    = 0;
    while (my $line = <FH>) {
        chomp $line;
        if ($line =~ /^>\s*(\S+)/) {
            $name = $1;
        } else {
            $line_count++;
            $matrix_string .= "$line\n";
            if ($line_count == 4) {
                push @pfms,
                    TFBS::Matrix::PFM->new(
                    -matrixstring => $matrix_string,
                    -name         => $name,
                    -ID           => $name
                    );
                $line_count    = 0;
                $name          = '';
                $matrix_string = '';
            }
        }
    }
    close(FH);

    return @pfms ? \@pfms : undef;
}
