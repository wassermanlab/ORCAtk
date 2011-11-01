#!/usr/local/bin/perl -w

=head1 NAME

orca_pairwise.pl (was consite.pl)

=head1 SYNOPSIS

 orca_pairwise.pl -s1 seq1 -s2 seq2 [-m1 masked_seq1] [-m2 masked_seq2]
        [-pos position]
        [-ex exons] [-ex2 exons2]
        [-f features] [-f2 features2]
        [-tp top_pct] [-cc cons_cutoff] [-w win_size] [-cl min_cr_len]
        [-jdb db_name] [-c collection]
        [-tax tax_groups] [-species species] [-ic min_ic]
        ([-tn tf_name] | [-pfmf pfm_file])
        [-th tfbs_threshold]
        [-s start_pos] [-e end_pos] [-nfos] [-hs spacer]
        [-a alignment] [-tf tf_file] [-crf cr_file] [-csf cs_file]
        [-g out_graph] [-flip] [-utf ucsc_track_file]

=head1 ARGUMENTS
 
  All argument switches may be abbreviated if unique.
  -s1 seq1          = FastA file containing the first sequence to be
                      analyzed.
  -s2 seq2          = FastA file containing the second sequence to be
                      analyzed.
  -m1 masked_seq1   = Optional. FastA file containing masked version
                      of sequence1.
  -m2 masked_seq2   = Optional. FastA file containing masked version
                      of sequence2.
  -pos position     = Chromosomal coordinates of input sequence in
                      format: chrN:start-end
  -ex exons         = Optional GFF file containing the exons positions
                      on the first sequence. If provided these will be
                      filtered out and excluded from the conserved
                      regions searched for TFBSs
  -ex2 exons2       = Optional GFF file containing the exons positions
                      on the second sequence.
  -f features       = Optional GFF file containing other features to plot
                      (i.e. CpG islands) on the first sequence.
  -f2 features2     = Optional GFF file containing other features to plot
                      (i.e. CpG islands) on the second sequence.
  -a alignment      = Optionally write the alignment to this file
                      (multi-FastA)
  -th threshold     = Optional. Minimum TFBS PWM score to report.
                      Specify as a relative score (i.e. 80%) or as
                      an absolute score (i.e. 11.2).
                      (default 80%)
  -tp top_pct       = Optional. Use this top percentile of conservation
                      windows to compute conservation cutoff. Specify
                      as a percentage in the range 0% - 100% or 0 - 1.
                      (default 10%)
  -cc cons_cutoff   = Optional. Min. conservation of regions in which
                      to search for TF sites. Specify as a percentage
                      in the range 0% - 100% or 0 - 1.
                      (default 70%)
  -pm pct_id_method = Optional. Method for computing the percent
                      identity. 'o' = overall, 's' = standard.
                      (default = 's')
  -cl min_cr_len    = Optional. Min. conserved region length to
                      consider for searching TFBSs (default 20)
  -jdb db_name      = Optional. Name of JASPAR database to use.
                      (default = JASPAR_2010)
  -c collection     = Optional. Name of JASPAR collection to use.
                      (default = CORE)
  -tax tax_groups   = Optional. Limit matrix profiles to those with these
                      taxonomic supergroups. To specify more than a single
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
  -w win_size       = Optional. Size of window used to compute
                      conserved regions.
  -s start          = Optional. Starting position on sequence 1 to
                      scan for TFBSs
  -e end            = Optional. Ending position on sequence 1 to
                      scan for TFBSs
  -nfos             = Optional. Do not filter overlapping sites. Overlapping
                      sites of the same TF are filtered so that only the
                      highest scoring is reported by default. With this
                      option specified all sites are reported.
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
  -utf ucsc_track_file
                    = If specified, write a file in the format of a
                      UCSC browser custom track

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
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::Tools::Run::RepeatMasker;
use TFBS::DB::JASPAR5;
use ORCA::Analysis::Pairwise;
use ORCA::Graphics::Pairwise;

use constant DEBUG => 0;

#
# Default minimum TFBS PWM score to report for either sequence.
# Overidden by -th argument.
#
use constant TFBS_THRESHOLD => "80%";

#
# Default top percentile of conserved windows to use for dynamically computing
# conservation cutoff.
# Overidden by the -tp argument.
#
use constant TOP_PERCENTILE => "10%";

#
# Default absolute minimum percent identity of conserved regions to search for
# TFBSs.
# Overidden by the -cc argument.
#
use constant MIN_CONSERVATION => "70%";

#
# Default method for computing percent identity.
#
use constant PCT_ID_METHOD => "s";

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
use constant MIN_IC => 10;

#
# Default window size used to compute conserved regions.
# Overidden by -w argument.
#
use constant WIN_SIZE => 100;

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

my $seq_file1;
my $seq_file2;
my $masked_seq_file1;
my $masked_seq_file2;
my $position;
my $exon_file;
my $exon_file2;
my $feat_file;
my $feat_file2;
my $align_file;
my $threshold;
my $top_pct;
my $min_conservation;
my $pct_id_method;
my $min_cr_len;
my $jaspar_db;
my $collection;
my @tax_groups;
my @species;
my $min_ic;
my $tf_name;
my $pfm_file;
my $win_size;
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
    's1=s'   => \$seq_file1,
    's2=s'   => \$seq_file2,
    'm1=s'   => \$masked_seq_file1,
    'm2=s'   => \$masked_seq_file2,
    'pos=s'  => \$position,
    'ex=s'   => \$exon_file,
    'ex2=s'  => \$exon_file2,
    'f=s'    => \$feat_file,
    'f2=s'   => \$feat_file2,
    'a=s'    => \$align_file,
    'th=s'   => \$threshold,
    'tp=s'   => \$top_pct,
    'cc=s'   => \$min_conservation,
    'pm=s'   => \$pct_id_method,
    'cl=i'   => \$min_cr_len,
    'jdb=s'  => \$jaspar_db,
    'c=s'    => \$collection,
    'tax=s'  => \@tax_groups,
    'species=s'  => \@species,
    'ic=i'   => \$min_ic,
    'tn=s'   => \$tf_name,
    'pfmf=s' => \$pfm_file,
    'w=i'    => \$win_size,
    's=i'    => \$start,
    'e=i'    => \$end,
    'tf=s'   => \$tfbs_file,
    'crf=s'  => \$cr_file,
    'csf=s'  => \$cs_file,
    'g=s'    => \$graph_file,
    'flip'   => \$flip,
    'utf=s'  => \$ucsc_track_file,
    'nfos'   => \$nfos,
    'hs=i'   => \$hs_spacer
);

if (!$seq_file1 || !$seq_file2) {
    pod2usage(-verbose => 1);
}

if (!$top_pct && !$min_conservation) {
    $top_pct = TOP_PERCENTILE;
}

unless (defined $threshold) {
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

unless (defined $pct_id_method) {
    $pct_id_method = PCT_ID_METHOD;
}

unless (defined $min_cr_len) {
    $min_cr_len = MIN_CR_LENGTH;
}

unless (defined $jaspar_db) {
    $jaspar_db = JASPAR_DB;
}

unless (defined $min_ic) {
    $min_ic = MIN_IC;
}

unless (defined $win_size) {
    $win_size = WIN_SIZE;
}

my $fos = 1;
if ($nfos) {
    $fos = 0;
}

#
# Check any files specfied exist first
#
if ($seq_file1 && !-f $seq_file1) {
    die "Specified sequence file 1 $seq_file1 not found!\n";
}
if ($seq_file2 && !-f $seq_file2) {
    die "Specified sequence file 2 $seq_file2 not found!\n";
}
if ($masked_seq_file1 && !-f $masked_seq_file1) {
    die "Specified masked sequence file 1 $masked_seq_file1 not found!\n";
}
if ($masked_seq_file2 && !-f $masked_seq_file2) {
    die "Specified masked sequence file 2 $masked_seq_file2 not found!\n";
}
if ($exon_file && !-f $exon_file) {
    die "Specified exon file $exon_file not found!\n";
}
if ($exon_file2 && !-f $exon_file2) {
    die "Specified exon file $exon_file2 not found!\n";
}
if ($feat_file && !-f $feat_file) {
    die "Specified feature file $feat_file not found!\n";
}
if ($feat_file2 && !-f $feat_file2) {
    die "Specified feature file $feat_file2 not found!\n";
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
print "Reading first sequence from $seq_file1\n";
my $seq1 = read_seq($seq_file1);
die "Error reading sequence 1\n" if !$seq1;

my $display_id = $seq1->display_id;
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
#} elsif (!$seq_chr) {
#    die "Could not determine chromosomal position of sequence - please"
#        . " explicitly enter a position for this sequence\n";
}

print "Reading second sequence from $seq_file2\n";
my $seq2 = read_seq($seq_file2);
die "Error reading sequence 2\n" if !$seq2;

#
# Optionally read masked sequences
#
# If masked sequence files were not specified, try to find files with the same
# names as the unmasked sequence files with a '.masked' extension added. If not
# found try the same in a subdir named 'masked'.
#
my $masked_seq1;
if ($masked_seq_file1) {
    print "Reading first masked sequence from $masked_seq_file1\n";
    $masked_seq1 = read_seq($masked_seq_file1);
    die "Error reading sequence 1\n" if !$masked_seq1;
} else {
    my ($name1, $path1, $suffix1) = fileparse($seq_file1, qr{\..*});
    $path1 = '.' if !$path1 || $path1 eq "";
    my $masked_dir1      = $path1 . "masked";
    my $masked_seq_file1 = "$path1$name1$suffix1.masked";
    if (-f $masked_seq_file1) {
        print "Found masked sequence file $masked_seq_file1\n";
        $masked_seq1 = read_seq($masked_seq_file1);
    } else {
        $masked_seq_file1 = "${path1}masked$name1$suffix1.masked";
        if (-f $masked_seq_file1) {
            print "Found masked sequence file $masked_seq_file1\n";
            $masked_seq1 = read_seq($masked_seq_file1);
        }
    }
}

my $masked_seq2;
if ($masked_seq_file2) {
    print "Reading second masked sequence from $masked_seq_file2\n";
    $masked_seq2 = read_seq($masked_seq_file2);
    die "Error reading sequence 2\n" if !$masked_seq2;
} else {
    my ($name2, $path2, $suffix2) = fileparse($seq_file2, qr{\..*});
    $path2 = '.' if !$path2 || $path2 eq "";
    my $masked_dir2      = $path2 . "masked";
    my $masked_seq_file2 = "$path2$name2$suffix2.masked";
    if (-f $masked_seq_file2) {
        print "Found masked sequence file $masked_seq_file2\n";
        $masked_seq2 = read_seq($masked_seq_file2);
    } else {
        $masked_seq_file2 = "${path2}masked$name2$suffix2.masked";
        if (-f $masked_seq_file2) {
            print "Found masked sequence file $masked_seq_file2\n";
            $masked_seq2 = read_seq($masked_seq_file2);
        }
    }
}

#
# Optionally get sequence 1 exons
#
my $seq_exons1;
if ($exon_file) {
    print "Reading base sequence exons from $exon_file\n";
    $seq_exons1 = read_exons($exon_file);
    die "Could not read sequence 1 exons\n" if !$seq_exons1;
}

#
# Optionally get sequence 2 exons
#
my $seq_exons2;
if ($exon_file2) {
    print "Reading comparison sequence exons from $exon_file2\n";
    $seq_exons2 = read_exons($exon_file2);
    die "Could not read sequence 2 exons\n" if !$seq_exons2;
}

#
# Optionally get sequence 1 features
#
my $seq_feats1;
if ($feat_file) {
    print "Reading base sequence features from $feat_file\n";
    $seq_feats1 = read_features($feat_file);
    die "Could not read features\n" if !$seq_feats1;
}

#
# Optionally get sequence 2 features
#
my $seq_feats2;
if ($feat_file2) {
    print "Reading comparison sequence features from $feat_file2\n";
    $seq_feats2 = read_features($feat_file2);
    die "Could not read features\n" if !$seq_feats2;
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
        die "Error reading TFBS profile(s) from PFM file $pfm_file\n" if !$pfms;
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

        $matrix_args{-min_ic}       = $min_ic if $min_ic;
        $matrix_args{-collection}   = $collection if defined $collection;
        $matrix_args{-tax_group}    = \@tax_groups if @tax_groups;
        $matrix_args{-species}      = \@species if @species;

        $matrix_set = $tfbs_db->get_MatrixSet(%matrix_args);
        die "Error retrieving TFBS matrix set from $jaspar_db\n"
            if !$matrix_set;
    }
}

my $pwa = ORCA::Analysis::Pairwise->new(
    -base_seq              => $seq1,
    -comparison_seq        => $seq2,
    -masked_base_seq       => $masked_seq1,
    -masked_comparison_seq => $masked_seq2,
    -base_seq_exons        => $seq_exons1,
    -chr                   => $seq_chr,
    -start                 => $seq_start,
    -end                   => $seq_end
);

#
# Create alignment
#
print "Aligning sequences...\n";
my $align = $pwa->compute_alignment();

#
# Optionally create new alignment file
#
if ($align_file) {
    print "Creating alignment file $align_file\n";
    write_align($align_file, $align);
}

if ($seq_exons2) {
    $seq_exons2 = convert_coords($pwa, $seq_exons2);
}

if ($seq_feats2) {
    $seq_feats2 = convert_coords($pwa, $seq_feats2);
}

#
# Perform conservation analysis
#
my %ccp_param;
$ccp_param{-position_type}   = 'c';
$ccp_param{-window_size}     = $win_size;
$ccp_param{-window_inc}      = 1;

my $conservation_profile = $pwa->compute_conservation_profile(%ccp_param);

my %ccr_param;
$ccr_param{-top_pct}             = $top_pct          if $top_pct;
$ccr_param{-min_conservation}    = $min_conservation if $min_conservation;
$ccr_param{-filter_exons}        = 1                 if $seq_exons1;
$ccr_param{-min_filtered_cr_len} = $min_cr_len       if $min_cr_len;
$ccr_param{-pct_id_method}       = $pct_id_method    if $pct_id_method;

my $conserved_regions = $pwa->compute_conserved_regions(%ccr_param);

if ($cr_file) {
    print "Writing conserved regions file $cr_file\n";
    write_conserved_regions_report($cr_file, $pwa->conserved_regions_report);
}

if ($cs_file) {
    print "Writing conserved sub-sequences file $cs_file\n";
    write_conserved_subsequences(
        $cs_file,
        $pwa->compute_conserved_subsequences
    );
}

#
# Perform TFBS analysis within conserved regions
#
my $tfbss;
if ($matrix_set && $matrix_set->size > 0) {
    print "Searching for conserved TFBSs...\n";

    $tfbss = $pwa->compute_conserved_tfbss(
        -matrix_set               => $matrix_set,
        -tfbs_threshold           => $threshold,
        -min_tfbs_cr_overlap      => MIN_TFBS_CR_OVERLAP,
        -filter_overlapping_sites => $fos,
        -start                    => $start,
        -end                      => $end
    );

    if ($tfbss) {
        if ($hs_spacer) {
            my $hs_pairs = pair_half_sites($tfbss, $hs_spacer);
            write_half_sites($tfbs_file, $hs_pairs);
        } else {
            print "Writing putative TFBS sites...\n";
            write_conserved_tfbss($tfbs_file, $tfbss);
        }
    } else {
        print "No conserved TFBSs found.\n";
    }
}

#
# Optionally plot TF sites
#
if ($graph_file) {
    #
    # Graph results.
    #
    my $graph = ORCA::Graphics::Pairwise->new(
        -analysis             => $pwa,
        -other_features       => $seq_feats1,
        -other_features2      => $seq_feats2,
        -flip                 => $flip
    );

    if ($graph && $graph->{-gd_image}) {
        open PNG, ">$graph_file"
            || die "Error opening $graph_file for writing\n";
        print PNG $graph->{-gd_image}->png;
        close PNG;
    } else {
        carp "ERROR: creating pairwise graph\n";
    }
}

#   
# Optionally create a UCSC custom browser track file
#   
if ($ucsc_track_file) {
    my $track = $pwa->ucsc_track;

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

sub read_exons
{
    my ($file) = @_;

    my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 2);

    my @exons;
    if ($gffio) {
        while (my $feature = $gffio->next_feature) {
            push @exons, $feature;
        }
    } else {
        die "Error opening exon file $file\n";
    }

    return @exons ? \@exons : undef;
}

sub read_features
{
    my ($file) = @_;

    my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 2);

    my @feats;
    if ($gffio) {
        while (my $feature = $gffio->next_feature) {
            push @feats, $feature;
        }
    } else {
        die "Error opening features file $file\n";
    }

    return @feats ? \@feats : undef;
}

sub convert_coords
{
    my ($pwa, $feats) = @_;

    my @cfeats;
    foreach my $f (@$feats) {
        my $cstart = $pwa->convert_seq_pos(2, 1, $f->start, 'ge');
        $cstart = $pwa->convert_seq_pos(2, 1, $f->start, 'le') if !$cstart;
        my $cend = $pwa->convert_seq_pos(2, 1, $f->end, 'ge');
        $cend = $pwa->convert_seq_pos(2, 1, $f->end, 'le') if !$cend;

        ($cstart, $cend) = ($cend, $cstart) if $cend < $cstart;

        push @cfeats, Bio::SeqFeature::Generic->new(
            -primary_id   => $f->primary_id,
            -display_name => $f->display_name,
            -source_tag   => $f->source_tag,
            -strand       => $f->strand,
            -score        => $f->score,
            -start        => $cstart,
            -end          => $cend,
        );
    }

    return @cfeats ? \@cfeats : undef;
}

sub pair_half_sites
{
    my ($tfbss, $spacer) = @_;

    print "Pairing half-sites...\n";

    my @hs_pairs;

    my $num_sites = scalar @$tfbss;

    # build hash of unique sites (based on strand and start position)
    my %unique_flag;
    my @unique_site_pairs;
    foreach my $sp (@$tfbss) {
        my $site_key = $sp->site1->strand . $sp->site1->start;
        if (!$unique_flag{$site_key}) {
            $unique_flag{$site_key} = 1;
            push @unique_site_pairs, $sp;
        }
    }

    foreach my $sp1 (@unique_site_pairs) {
        foreach my $sp2 (@unique_site_pairs) {
            next if $sp2->site1->start <= $sp1->site1->end;
            last if $sp2->site1->start > $sp1->site1->end + $spacer + 1;

            push @hs_pairs, {
                -site_pair1 => $sp1,
                -site_pair2 => $sp2
            };
        }
    }

    return @hs_pairs ? \@hs_pairs : undef;
}

sub write_align
{
    my ($align_file, $align) = @_;

    my $alignIO = Bio::AlignIO->new(
        -file   => ">$align_file",
        -format => 'fasta'
    );

    if ($alignIO) {
        if (!$alignIO->write_aln($align)) {
            warn "Could not create alignment file $align_file\n";
        }
        $alignIO->close;
    } else {
        warn "Could not create alignment file $align_file\n";
    }
}

sub write_conserved_regions_report
{
    my ($file, $report) = @_;

    return if !$report;

    my $crs = $report->conserved_regions;
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

    printf CRS "Window Size         : %d bp\n", $report->param('window_size');
    printf CRS "Window Increment    : %d bp\n", $report->param('window_inc');
    printf CRS "Top X%% identities   : %.1f\n", $report->param('top_pct')
        if defined $report->param('top_pct');
    printf CRS "Minimum %% identity  : %.1f\n", $report->param('min_pct_id')
        if defined $report->param('min_pct_id');
    printf CRS "Computed %% identity : %.1f\n", $report->param('comp_pct_id')
        if defined $report->param('comp_pct_id');
    printf CRS "Effective %% identity: %.1f\n", $report->param('cutoff')
        if defined $report->param('cutoff');
    print CRS "\n";

    foreach my $cr (@$crs) {
        printf CRS "%7d %7d %7d %7d %7d %7d %7d %7d %.4f\n",
            $cr->seq1_start,
            $cr->seq1_end,
            $cr->seq1_end - $cr->seq1_start + 1,
            $cr->seq2_start,
            $cr->seq2_end,
            $cr->seq2_end - $cr->seq2_start + 1,
            $cr->align_start,
            $cr->align_end,
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
    my ($file, $tfbs) = @_;

    return if !$tfbs;

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
    my $start = $pwa->start;
    foreach my $site_pair (@$tfbss) {
        my $site1 = $site_pair->site1;
        #my $site2 = $site_pair->site2;

        my $strand1;
        if ($site1->strand eq '+' || $site1->strand == 1) {
            $strand1 = '+';
        } elsif ($site1->strand eq '-' || $site1->strand == -1) {
            $strand1 = '-';
        }

        #my $strand2;
        #if ($site2->strand eq '+' || $site2->strand == 1) {
        #    $strand2 = '+';
        #} elsif ($site2->strand eq '-' || $site2->strand == -1) {
        #    $strand2 = '-';
        #}

        #printf TFBS "%-16s %6d %6d %s %5.1f %-25s\n",
        #    $site1->pattern->name,
        #    $site1->start,
        #    $site1->end,
        #    $strand1,
        #    $site1->rel_score * 100,
        #    $site1->seq->seq;
        #printf TFBS "\t\t %6d %6d %s %5.1f %-25s %5.1f\n",
        #    $site2->start,
        #    $site2->end,
        #    $strand2,
        #    $site2->rel_score * 100,
        #    $site2->seq->seq,
        #    ($site1->get_tag_values('conservation'))[0] * 100;

        printf TFBS "%s\t%d\t%d\t%s\t%.3f\t%.1f%%\t%.1f%%\t%s\n",
            $site1->pattern->name,
            $site1->start + $start - 1,
            $site1->end + $start - 1,
            $strand1,
            $site1->score,
            $site1->rel_score * 100,
            ($site1->get_tag_values('conservation'))[0] * 100,
            $site1->seq->seq;
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
        my $hs1 = $hs_pair->{'-site_pair1'}->site1;
        my $hs2 = $hs_pair->{'-site_pair2'}->site1;

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
