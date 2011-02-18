#!/usr/local/bin/perl -w

=head1 NAME

cons_graph.pl

=head1 SYNOPSIS

  cons_graph.pl -s1 seq_file1 (fasta)
  		-s2 seq_file2 (fasta)
		[-ex exon_file (GFF)]
  		[-a alignment_file (multi-fasta)]
  		[-mus muscle_file]
  		[-liv liver_file]
		[-mc muscle_cons_file]
		[-lc liver_cons_file]
  		[-tf tfbs_file]
  		[-th tfbs_threshold (%)]
  		[-cc cons_cutoff (%)]
		[-tp top_pct (%)]
  		[-cl min_cr_length]
  		[-ic min_ic]
  		[-w win_size]
		[-s start_pos]
		[-e end_pos]
		-o out_png_file

=head1 ARGUMENTS
 
  All argument switches may be abbreviated if unique.
  -s1 seq_file1		= FastA file containing the first sequence to be
  			  analyzed.
  -s2 seq_file2		= FastA file containing the second sequence to be
  			  analyzed.
  -ex exon_file		= Optional. GFF file containing the exons positions
  			  on the first sequence. If provided these will be
			  filtered out excluded from the conserved regions
			  and reported as a track on the graph.
  -a align_file		= Optional. Write the alignment to this file
  			  (multi-FastA)
  -mus muscle_file	= Optional. File containing position and scores
  			  for muscle model analysis as output by lra_step
  -liv liver_file	= Optional. File containing position and scores
  			  for liver model analysis as output by lra_step
  -mc muscle_cons_file	= Optional. File containing position and scores
  			  of the product of the muscle model and conservation
			  profile (as output by lra_cons).
  -lc liver_cons_file	= Optional. File containing position and scores
  			  of the product of the liver model and conservation
			  profile (as output by lra_cons).
  -tf tfbs_file		= Optional. Write the TFBS sites to this file
  -th threshold		= Optional. Minimum TFBS score to report.
  			  (default 70; range 0 - 100)
  -tp top_pct		= Optional. Use this top % of conservation windows
			  to compute conservation cutoff
			  (default 10; range 0 - 100)
  -cc cons_cutoff	= Optional. Min. conservation of regions in which
  			  to search for TF sites
  			  (default 70; range 0 - 100)
  -cl			= Optional. Min. conserved region length to
  			  consider for searching TFBSs (default 20)
  -ic			= Optional. Min. information content of TFBS PWM
  			  matrices to use (default 8)
  -s			= Optional. Starting position on sequence 1 to plot
  			  (and perform TFBS analysis)
  -e			= Optional. Ending position on sequence 1 to plot
  			  (and perform TFBS analysis)
  -o out_png_file	= Name of output graphics file in PNG format

=head1 DESCRIPTION

Given two sequence files, align them and find the conservation. Search
for TFBSs in the conserved regions. Output a graph consisting of the
conservation profile, exons (if provided), conserved regions and TFBS
positions in PNG format.

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia
  
  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use Getopt::Long;
use Pod::Usage;
use Orca;
use File::Basename;
use TFBS::DB::JASPAR2;
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Align::AlignI;
use Bio::Tools::Run::RepeatMasker;
use ConservationAnalysis::ConservationAnalysis;
use ConservationAnalysis::TFBSConservation;
use ConservationAnalysis::Graphics::Graph;

use constant DEBUG => 0;

#
# Default min. score to report for either sequence.
# Overidden by -threshold argument.
#
use constant TFBS_THRESHOLD => "70%";

#
# Mininimum amount in bp, by which a TFBS needs to overlap a conserved region
# in order to be included in analysis
#
use constant MIN_TFBS_CONSERVATION_OVERLAP => 1;

#
# JASPAR database connection parameters.
#
#use constant JASPAR_DB_HOST => "forkhead.cgb.ki.se";
#use constant JASPAR_DB_USER => "krivan";
#use constant JASPAR_DB_PASS => "wk3003";
#use constant JASPAR_DB_HOST => "mordor.cgb.ki.se";
#use constant JASPAR_DB_USER => "cmmt";
#use constant JASPAR_DB_PASS => "vancouver";
use constant JASPAR_DB_HOST => "napa.cmmt.ubc.ca";
use constant JASPAR_DB_USER => "jaspar_r";
use constant JASPAR_DB_PASS => "";
use constant JASPAR_DB_NAME => "JASPAR2";

my $PROG = 'cons_graph.pl';

my $seq_file1;
my $seq_file2;
my $exon_file;
my $align_file;
my $muscle_file;
my $liver_file;
my $muscons_file;
my $livcons_file;
my $tfbs_file;
my $threshold = TFBS_THRESHOLD;
my $top_pct;
my $min_conservation;
my $min_cr_len;
my $min_ic;
my $win_size;
my $start;
my $end;
my $out_file;
GetOptions(
    's1=s'		=> \$seq_file1,
    's2=s'		=> \$seq_file2,
    'ex=s'		=> \$exon_file,
    'a=s'		=> \$align_file,
    'mus=s'		=> \$muscle_file,
    'liv=s'		=> \$liver_file,
    'mc=s'		=> \$muscons_file,
    'lc=s'		=> \$livcons_file,
    'tf=s'		=> \$tfbs_file,
    'th=f'		=> \$threshold,
    'tp=f'		=> \$top_pct,
    'cc=f'		=> \$min_conservation,
    'cl=i'		=> \$min_cr_len,
    'ic=i'		=> \$min_ic,
    'w=i'		=> \$win_size,
    's=i'		=> \$start,
    'e=i'		=> \$end,
    'o=s'		=> \$out_file
);

if (!$seq_file1 || !$seq_file2 || !$out_file) {
    pod2usage(-verbose => 1);
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
if ($exon_file && !-f $exon_file) {
    die "Specified exon file $exon_file not found!\n";
}
if ($muscle_file && !-f $muscle_file) {
    die "Specified muscle model file $muscle_file not found!\n";
}
if ($liver_file && !-f $liver_file) {
    die "Specified liver model file $liver_file not found!\n";
}
if ($muscons_file && !-f $muscons_file) {
    die "Specified muscle model conservation file $muscons_file not found!\n";
}
if ($livcons_file && !-f $livcons_file) {
    die "Specified liver model conservation file $livcons_file not found!\n";
}

#
# Read sequences
#
print "Reading first sequence from $seq_file1\n";
my $seq1 = read_seq($seq_file1);
die "Error reading sequence 1\n" if !$seq1;

print "Reading second sequence from $seq_file2\n";
my $seq2 = read_seq($seq_file2);
die "Error reading sequence 2\n" if !$seq2;

#
# Optionally get first sequence exons
#
my $seq1_exons;
if ($exon_file) {
    print "Reading base sequence exons from $exon_file...\n";
    $seq1_exons = read_exons($exon_file);
    die "Could not read exons\n" if !$seq1_exons;
}

#
# Optionally read muscle model profile
#
my $muscle_profile;
if ($muscle_file) {
    $muscle_profile = read_profile($muscle_file);
    die "Error reading muscle model profile\n" if !$muscle_profile;
    #
    # Currently lra_step outputs positions as the end point of the window
    # whereas the default for align_cons is to output the centre of the
    # window. Therefore shift the lra_step profile to the centre of the
    # window for consistency.
    #
    $muscle_profile = shift_profile($muscle_profile);
}

#
# Optionally read liver model profile
#
my $liver_profile;
if ($liver_file) {
    $liver_profile = read_profile($liver_file);
    die "Error reading liver model profile\n" if !$liver_profile;
    #
    # Currently lra_step outputs positions as the end point of the window
    # whereas the default for align_cons is to output the centre of the
    # window. Therefore shift the lra_step profile to the centre of the
    # window for consistency.
    #
    $liver_profile = shift_profile($liver_profile);
}

#
# Optionally read muscle model conservation profile
#
my $muscons_profile;
if ($muscons_file) {
    $muscons_profile = read_profile($muscons_file);
    die "Error reading muscle model conservation profile\n"
    	if !$muscons_profile;
}

#
# Optionally read liver model conservation profile
#
my $livcons_profile;
if ($livcons_file) {
    $livcons_profile = read_profile($livcons_file);
    die "Error reading liver model conservation profile\n"
    	if !$livcons_profile;
}

#
# Connect to JASPAR2 database
#
my $tfbs_db = TFBS::DB::JASPAR2->connect(
			"dbi:mysql:" . JASPAR_DB_NAME . ":" . JASPAR_DB_HOST,
		    	JASPAR_DB_USER, JASPAR_DB_PASS);
die "Error connecting to TFBS database\n" if !defined $tfbs_db;

my %matrix_args = (-matrixtype	=> 'PWM');
$matrix_args{-min_ic} = $min_ic if $min_ic;

my $matrix_set = $tfbs_db->get_MatrixSet(%matrix_args);
die "Error retrieving matrix set from JASPAR2\n" if !$matrix_set;


#
# Try to read masked sequence files. Try to find files with the same names as
# the two unmasked sequence files with a '.masked' extension added. If not
# found try the same in a subdir named 'masked'.
#
my ($masked_seq1, $masked_seq2, $seq1_repeats);
my ($name1, $path1, $suffix1) = fileparse($seq_file1, qr{\..*});
my ($name2, $path2, $suffix2) = fileparse($seq_file2, qr{\..*});
$path1 = '.' if !$path1 || $path1 eq "";
$path2 = '.' if !$path2 || $path2 eq "";
my $masked_dir1 = $path1 . "masked";
my $masked_dir2 = $path2 . "masked";
my $masked_seq_file1 = "$path1$name1$suffix1.masked";
my $masked_seq_file2 = "$path2$name2$suffix2.masked";

print "Checking for masked sequences...\n";
if (-f $masked_seq_file1) {
    print "Found masked sequence file $masked_seq_file1\n";
    $masked_seq1 = read_seq($masked_seq_file1);
} else {
    $masked_seq_file1 = "${path1}masked$name1$suffix1.masked";
    if (-f $masked_seq_file1) {
	print "Found masked sequence file $masked_seq_file1\n";
	$masked_seq1 = read_seq($masked_seq_file1);
    } else {
	print "Masked sequence 1 not found.\n";
    }
}
if (-f $masked_seq_file2) {
    print "Found masked sequence file $masked_seq_file2\n";
    $masked_seq2 = read_seq($masked_seq_file2);
} else {
    $masked_seq_file2 = "${path2}masked$name2$suffix2.masked";
    if (-f $masked_seq_file2) {
	print "Found masked sequence file $masked_seq_file2\n";
	$masked_seq2 = read_seq($masked_seq_file2);
    } else {
	print "Masked sequence 2 not found.\n";
    }
}

#
# If repeat masked sequences were not found...
#
if (!$masked_seq1) {
    #
    # ...create masked sequence
    #
    print "Repeat masking sequence 1\n";
    ($masked_seq1, $seq1_repeats) = repeat_mask($seq1);
    if ($masked_seq1) {
	#
	# Write repeat masked sequence to a file
	#
	print "Creating masked sequence file 1 $path1$name1$suffix1.masked\n";
	my $seqIO = Bio::SeqIO->new(-file => ">$path1$name1$suffix1.masked",
				    -format => 'fasta');
	$seqIO->write_seq($masked_seq1);
    } else {
	warn "Could not repeat mask sequence 1\n" if !$masked_seq1;
    }
}

if (!$masked_seq2) {
    #
    # ...create masked sequence
    #
    print "Repeat masking sequence 2\n";
    ($masked_seq2) = repeat_mask($seq2, "mus");
    if ($masked_seq2) {
	#
	# Write repeat masked sequence to a file
	#
	print "Creating masked sequence file 2 $path2$name2$suffix2.masked\n";
	my $seqIO = Bio::SeqIO->new(-file => ">$path2$name2$suffix2.masked",
				    -format => 'fasta');
	$seqIO->write_seq($masked_seq2);
    } else {
	warn "Could not repeat mask sequence 2\n" if !$masked_seq2;
    }
}

#
# Create alignment
#
print "Aligning masked sequences...\n";
my $align;
my $orca = Orca->new();
die "Error initiating ORCA\n" if !$orca;
$align = $orca->Align(
	    -seq1	=> $masked_seq1 || $seq1,
	    -seq2	=> $masked_seq2 || $seq2);

#
# Create new alignment file
#
if ($align_file) {
    print "Creating alignment file $align_file\n";
    write_align($align_file, $align);
}

my $orient;
my $seq2_strand = $align->get_seq_by_pos(2)->strand;
if ($seq2_strand eq '-' || $seq2_strand eq "-1") {
    $orient = 0;
} else {
    $orient = 1;
}

#
# Perform conservation analysis
#
my $ca = ConservationAnalysis::ConservationAnalysis->new(
			base_seq                => $seq1,
			comparison_seq          => $seq2,
			alignment               => $align,
			masked_base_seq         => $masked_seq1,
			masked_comparison_seq   => $masked_seq2,
			base_seq_exons          => $seq1_exons);
$ca->param('top_pct' => $top_pct) if $top_pct;
$ca->param('min_conservation' => $min_conservation) if $min_conservation;
$ca->param('window_size' => $win_size) if $win_size;

my %cr_params = (filter_exons => 1);
$cr_params{'min_filtered_cr_len'} = $min_cr_len if $min_cr_len;
my $conserved_regions = $ca->compute_conserved_regions(%cr_params);

my $conservation_profile = $ca->compute_conservation_profile(
    						position_type => 'c');
#
# Perform TFBS analysis within conserved regions
#
print "Computing conserved TFBSs...\n";
my $tfbs_cons = ConservationAnalysis::TFBSConservation->new(
	matrix_set                  	=> $matrix_set,
	tfbs_threshold              	=> $threshold,
	conservation_cutoff         	=> $min_conservation,
	min_tfbs_conservation_overlap	=> MIN_TFBS_CONSERVATION_OVERLAP);

my $tf_sites = $tfbs_cons->find_conserved_tf_sites(
		conservation_analysis       => $ca,
		filter_overlapping_sites    => 1,
		start                       => $start,
		end                         => $end);

#
# Convert TF site pair set to feature list
#
my @tf_feats;
my $iter = $tf_sites->Iterator(-sort_by => 'start');
while (my $site_pair = $iter->next) {
    push @tf_feats, $site_pair->feature1;
}

if ($tfbs_file) {
    print "Creating TFBS file $tfbs_file\n";
    write_tfbs($tfbs_file, $tf_sites, $conserved_regions);
}

#
# Graph results.
#
my $cutoff;
if ($min_conservation) {
    if ($min_conservation =~ /(.+)%/) {
	$cutoff = $1 / 100;
    } else {
	$cutoff = $min_conservation;
    }
}
my $cons_graph = ConservationAnalysis::Graphics::Graph->new(
    			base_seq		=> $seq1,
			conservation_profile	=> $conservation_profile,
    			conserved_regions	=> $conserved_regions,
			cutoff			=> $cutoff,
    			tf_sites		=> \@tf_feats,
			exons			=> $seq1_exons,
		    	repeat_regions		=> $seq1_repeats,
		    	muscle_profile		=> $muscle_profile,
		    	liver_profile		=> $liver_profile,
		    	mus_cons_profile	=> $muscons_profile,
		    	liv_cons_profile	=> $livcons_profile,
			start			=> $start,
			end			=> $end
		    );

open PNG, ">$out_file" || die "Error opening $out_file for writing\n";
print PNG $cons_graph->{gd_image}->png;
close PNG;

exit;

sub read_seq
{
    my ($seq_file) = @_;

    my $seq_in;
    $seq_in = Bio::SeqIO->new(
			    -file	=> $seq_file,
			    -format	=> 'fasta');
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

sub read_profile
{
    my ($file) = @_;

    if (!open(FH, $file)) {
	die "Error opening profile file $file = $!\n";
    }

    my @profile;
    while (<FH>) {
	if (/^\s*(\d+)\s+(\d+\.\d+)\s*$/) {
	    push @profile, {position => $1, score => $2};
	}
    }

    close FH;

    return @profile ? \@profile : undef;
}

#
# Some profiles may be reported as the start or end position of the window
# (currently lra_step output the end point of the window). This routine
# shifts the position of a profile by the specified amount. If no amount
# is specified, try to guess whether profile is reported as start or end
# and shift accordingly
#
sub shift_profile
{
    my ($profile, $shift) = @_;

    if (!$shift) {
	my $pos1 = $profile->[0]->{position};
	my $pos2 = $profile->[1]->{position};
	if ($pos1 == 0 || $pos1 == 1) {
	    # positions are start of window
	    $shift = int($pos2 / 2);
	} else {
	    # positions are end of window
	    $shift = -(int($pos1 / 2));
	}
    }
    foreach my $point (@$profile) {
	$point->{position} += $shift;
    }

    return $profile;
}

sub repeat_mask
{
    my ($seq, $spec_mask) = @_;

    my %params = (
		    w		=> 1,
		    s		=> 1,
		    no_is	=> 1,
		    frag	=> 1000000);
    if ($spec_mask) {
	$params{$spec_mask} = 1;
    }

    my $rm = Bio::Tools::Run::RepeatMasker->new(%params);
    die "Error running RepeatMasker\n" if !$rm;

    my $masked_seq;
    my @feats;
    if ($rm->run($seq)) {
	$masked_seq = $rm->masked_seq;
	@feats = $rm->repeat_features;
    }

    return ($masked_seq, @feats ? \@feats : undef);
}

sub write_align
{
    my ($align_file, $align) = @_;

    my $alignIO = Bio::AlignIO->new(
				    -file	=> ">$align_file",
				    -format	=> 'fasta');
    if ($alignIO) {
	if (!$alignIO->write_aln($align)) {
	    warn "Could not create alignment file $align_file\n"
	}
	$alignIO->close;
    } else {
	warn "Could not create alignment file $align_file\n"
    }
}

sub write_tfbs
{
    my ($file, $tfbs, $conserved_regions) = @_;

    if (!open(TFBS, ">$file")) {
	warn "Could not create TFBS output file - $!\n";
	return;
    }

    my $iter = $tfbs->Iterator(-sort_by => 'start');
    while (my $site_pair = $iter->next()) {
	my $site1 = $site_pair->site1;
	my $site2 = $site_pair->site2;

	my $strand1;
	my $strand2;
	if ($site1->strand eq '+' || $site1->strand == 1) {
	    $strand1 = '+';
	} elsif ($site1->strand eq '-' || $site1->strand == -1) {
	    $strand1 = '-';
	}
	if ($site2->strand eq '+' || $site2->strand == 1) {
	    $strand2 = '+';
	} elsif ($site2->strand eq '-' || $site2->strand == -1) {
	    $strand2 = '-';
	}

	printf(TFBS "%s %6d %6d %s %s %.1f %6d %6d %s %s %.1f %.1f\n", 
	    $site1->pattern->name,
	    $site1->start,
	    $site1->end,
	    $strand1,
	    $site1->seq->seq,
	    $site1->rel_score * 100,
	    $site2->start,
	    $site2->end,
	    $strand2,
	    $site2->seq->seq,
	    $site2->rel_score * 100,
	    ($site1->get_tag_values('conservation'))[0] * 100
	);
    }

    close TFBS;
}
