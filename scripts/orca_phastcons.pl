#!/usr/local/bin/perl -w

=head1 NAME

orca_phastcons.pl

=head1 SYNOPSIS

 orca_phastcons.pl -sf seq_file -udb ucsc_db -track ucsc_track
 		-pos position
		[-ex exons] [-f features] [-rc]
		[-th tfbs_threshold]
		[-cc cons_cutoff] [-cl min_cr_len]
		[-db db_name] [-ic min_ic]
		([-tn tf_name] | [-pfmf pfm_file]) [-w win_size]
		[-s start_pos] [-e end_pos] [-fos] [-hs spacer]
		[-g out_graph] [-tf tf_file] [-crf cr_file] [-csf cs_file]

=head1 ARGUMENTS
 
  All argument switches may be abbreviated if unique.
  -sf seq_file		= FastA file containing the sequence to be analyzed.
  -udb ucsc_db		= UCSC database name, e.g. "hg18"
  -track ucsc_track	= UCSC track name, e.g. "phastCons28way"
  -pos position		= Chromosomal coordinates of input sequence in
  			  format: chrN:start-end
  -ex exons		= Optional GFF file containing the exons positions
  			  on the first sequence. If provided these will be
			  filtered out and excluded from the conserved
			  regions searched for TFBSs
  -f features		= Optional GFF file containing other features to plot
			  (i.e. CpG islands)
  			  on the first sequence. If provided these will be
  -rc			= Reverse complement phastCons profile
  -th threshold		= Optional. Minimum TFBS PWM score to report.
  			  Specify as a relative score (i.e. 80%) or as
			  an absolute score (i.e. 11.2).
  			  (default 80%)
  -cc cons_cutoff	= Optional. Min. conservation of regions in which
  			  to search for TF sites. Specify as a percentage
			  in the range 0% - 100% or 0 - 1.
  			  (default 70%)
  -cl min_cr_len	= Optional. Min. conserved region length to
  			  report (default 20)
  -db db_name		= Optional. Name of JASPAR database to use.
  			  (default = JASPAR_CORE_2008)
  -ic min_ic		= Optional. Minimum information content
  			  (specificity) of TFBS PWMs to use (default 10)
  -tn tf_name		= Optional. Name of a specific TF to use in
  			  the analysis.
  -pfmf file_name	= Optional. Name of a file defining one or more
  			  PFM matrix to use in the analysis.
  -w win_size		= Optional. Size of window used to compute
  			  conserved regions.
  -s start		= Optional. Starting position on sequence 1 to
  			  scan for TFBSs
  -e end		= Optional. Ending position on sequence 1 to
  			  scan for TFBSs
  -fos			= Optional. Graph indicating that overlapping sites
			  of the same TF should be filtered so that only
			  the highest scoring is reported
  -hs spacer		= Optional. Half-site spacer. Single profile only.
			  The profile is considered to be a half site. The
			  script attempts to find pairs of sites within (<=)
			  this specified numbers of nucleotides.
  -g graph_file		= Optional name of output graphics file in PNG
  			  format
  -tf tf_file		= If specified, write the putative TF sites to
  			  this file, else write to STDOUT
  -crf cr_file		= If specified, write the conserved regions to
  			  this file
  -csf cs_file		= If specified, write the conserved sub-sequences
			  to this file

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

  Bioperl		- at least version 1.2 or better.
  JASPAR4		- (http://jaspar.cgb.ki.se/cgi-bin/jaspar_db.pl)
			  a MySQL database of TF binding site profiles.
  TFBS			- (http://forkhead.cgb.ki.se/TFBS/) a set of perl
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

#use lib '/space/devel/ORCAtk/lib';

use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use TFBS::DB::JASPAR4;
use Bio::SeqIO;
use Bio::Seq;
use ORCA::ConservationAnalysis::PhastCons;
use ORCA::TFBSConservation::PhastCons;
use ORCA::Graphics::PhastConsGraph;

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
use constant MIN_CONSERVATION => "70%";

#
# Default JASPAR database
#
use constant JASPAR_DB => 'JASPAR_CORE_2008';

#
# Default minimum information content (specificity) of TFBS profiles from the
# JASPAR database to search against the sequences.
# Overidden by -ic argument.
#
use constant MIN_IC => 10;

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
use constant JASPAR_DB_HOST => "napa.cmmt.ubc.ca";
use constant JASPAR_DB_USER => "jaspar_r";
use constant JASPAR_DB_PASS => "";

my $seq_file;
my $ucsc_db;
my $track;
my $position;
#my $phc_file;
my $exon_file;
my $feat_file;
my $revcom = 0;
my $threshold = TFBS_THRESHOLD;
my $min_conservation;
my $min_cr_len = MIN_CR_LENGTH;
my $jaspar_db = JASPAR_DB;
my $min_ic = MIN_IC;
my $tf_name;
my $pfm_file;
my $start;
my $end;
my $fos = 0;
my $hs_spacer;
my $tfbs_file;
my $cr_file;
my $cs_file;
my $graph_file;
GetOptions(
    'sf=s'		=> \$seq_file,
    'udb=s'		=> \$ucsc_db,
    'track=s'		=> \$track,
    'pos=s'		=> \$position,
    #'phcf=s'		=> \$phc_file,
    'ex=s'		=> \$exon_file,
    'f=s'		=> \$feat_file,
    'rc'		=> \$revcom,
    'th=s'		=> \$threshold,
    'cc=s'		=> \$min_conservation,
    'cl=i'		=> \$min_cr_len,
    'db=s'		=> \$jaspar_db,
    'ic=i'		=> \$min_ic,
    'tn=s'		=> \$tf_name,
    'pfmf=s'		=> \$pfm_file,
    's=i'		=> \$start,
    'e=i'		=> \$end,
    'tf=s'		=> \$tfbs_file,
    'crf=s'		=> \$cr_file,
    'csf=s'		=> \$cs_file,
    'g=s'		=> \$graph_file,
    'fos'		=> \$fos,
    'hs=i'		=> \$hs_spacer
);

if (!$seq_file) {
    pod2usage(-verbose => 1);
}

if (!$min_conservation) {
    $min_conservation = MIN_CONSERVATION;
}

if (!$ucsc_db) {
    pod2usage(-verbose => 1, -msg => "No UCSC database name specified!\n");
}

if (!$track) {
    pod2usage(-verbose => 1, -msg => "No UCSC track name specified!\n");
}

if (!$position) {
    pod2usage(-verbose => 1, -msg => "No chromosomal position specified!\n");
}

my ($seq_chr, $seq_start, $seq_end);
$position =~ s/,//g;
if ($position =~ /chr(\S+):(\d+)-(\d+)/) {
    $seq_chr = $1;
    $seq_start = $2;
    $seq_end = $3;
} else {
    pod2usage(-verbose => 1,
    		-msg => "Bad format specified for position ($position)!\n");
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

#
# Read sequences
#
print "Reading sequence from $seq_file\n";
my $seq = read_seq($seq_file);
die "Error reading sequence\n" if !$seq;

#
# Read phastCons scores
#
#print "Reading phastCons scores from $phc_file\n";
#my ($profile, $scores) = read_phastCons_scores($phc_file, $revcom);
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
my $tfbs_db = TFBS::DB::JASPAR4->connect(
			"dbi:mysql:$jaspar_db:" . JASPAR_DB_HOST,
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
	my %matrix_args = (-matrixtype	=> 'PWM');
	$matrix_args{-min_ic} = $min_ic if $min_ic;
	$matrix_set = $tfbs_db->get_MatrixSet(%matrix_args);
	die "Error retrieving TFBS matrix set from $jaspar_db\n"
		if !$matrix_set;
    }
}

#
# Perform phastCons analysis
#
my $ca = ORCA::ConservationAnalysis::PhastCons->new(
					-seq		=> $seq,
					-db		=> $ucsc_db,
					-track		=> $track,
					-exons		=> $exons,
					-chr_name	=> $seq_chr,
					-start		=> $seq_start,
					-end		=> $seq_end);

my $phastCons_profile = $ca->fetch_conservation_profile();

#
# Compute conserved regions
#
my %ccr_params;
$ccr_params{-min_conservation} = $min_conservation if $min_conservation;
$ccr_params{-filter_exons} = 1 if $exons;
$ccr_params{-min_cr_len} = $min_cr_len if $min_cr_len;
my $conserved_regions = $ca->compute_conserved_regions(%ccr_params);

if ($conserved_regions && $cr_file) {
    print "Writing conserved regions file $cr_file\n";
    write_conserved_regions_report($cr_file, $min_conservation,
				    $ca->conserved_regions);
}

if ($cs_file) {
    print "Writing conserved sub-sequences file $cs_file\n";
    write_conserved_subsequences($cs_file, $ca->extract_conserved_subsequences);
}

my $tf_sites;
if ($tfbs_file && $matrix_set && $matrix_set->size > 0) {
    print "Searching for conserved TFBSs...\n";

    my $tfbs_phc = ORCA::TFBSConservation::PhastCons->new(
						-conservation_analysis => $ca);

    $tf_sites = $tfbs_phc->find_conserved_tf_sites(
			    -matrix_set                 => $matrix_set,
			    -tfbs_threshold             => $threshold,
			    -min_tfbs_cr_overlap	=> MIN_TFBS_CR_OVERLAP,
			    -filter_overlapping_sites	=> $fos,
			    -start			=> $start,
			    -end			=> $end
			);
    #print "Conserved seq TFBS sites: " . $tf_sites->size . "\n";

    if ($hs_spacer) {
	my $hs_pairs = pair_half_sites($tf_sites, $hs_spacer);
	write_half_sites($tfbs_file, $hs_pairs);
    } else {
	print "Writing putative TFBS sites...\n";
	write_tf_sites($tfbs_file, $tf_sites);
    }
}

#
# Optionally plot results
#
if ($graph_file) {
    print "Plotting results...\n";

    #
    # Convert TF site pair set to feature list
    #
    my @tf_feats;
    if ($tf_sites) {
	my $iter = $tf_sites->Iterator(-sort_by => 'start');
	while (my $site = $iter->next) {
	    push @tf_feats, $site;
	}
    }

    #
    # Graph results.
    #
    my $cutoff = 0;
    if ($min_conservation) {
	if ($min_conservation =~ /(.+)%/) {
	    $cutoff = $1 / 100;
	} else {
	    $cutoff = $min_conservation;
	}
    }

    my $cons_graph = ORCA::Graphics::PhastConsGraph->new(
			-seq			=> $seq,
			-conservation_profile	=> $phastCons_profile,
			-conserved_regions	=> $conserved_regions,
			-cutoff			=> $cutoff,
			-tf_sites		=> @tf_feats
							? \@tf_feats : undef,
			-exons			=> $exons,
			-features		=> $features,
			-start			=> $start,
			-end			=> $end
		    );

    if ($cons_graph && $cons_graph->{-gd_image}) {
	open PNG, ">$graph_file"
			    || die "Error opening $graph_file for writing\n";
	print PNG $cons_graph->{-gd_image}->png;
	close PNG;
    } else {
    	carp "ERROR: creating phastCons graph\n";
    }
}

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
    my ($tf_sites, $spacer) = @_;

    print "Pairing half-sites...\n";

    my @hs_pairs;

    my $num_sites = $tf_sites->size;

    # build hash of unique sites (based on strand and start position)
    my %unique_flag;
    my @unique_sites;
    my $iter = $tf_sites->Iterator(-sort_by => 'start');
    while (my $site = $iter->next()) {
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
				-site1	=> $site1,
				-site2	=> $site2};
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

    printf CRS "Percentage identity: %s\n", $min_conservation
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

sub write_tf_sites
{
    my ($file, $tfbs) = @_;

    if ($file) {
	if (!open(TFBS, ">$file")) {
	    warn "Could not create TFBS output file $file - $!\n";
	    return;
	}
    } else {
	open(TFBS, ">-");
    }

    my $iter = $tfbs->Iterator(-sort_by => 'start');
    while (my $site = $iter->next()) {
	my $strand;
	if ($site->strand eq '+' || $site->strand == 1) {
	    $strand = '+';
	} elsif ($site->strand eq '-' || $site->strand == -1) {
	    $strand = '-';
	}

	printf TFBS "%-16s %6d %6d %s %5.1f %5.1f %s\n", 
	    $site->pattern->name,
	    $site->start,
	    $site->end,
	    $strand,
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
    my $name = '';
    my $matrix_string = '';
    my $line_count = 0;
    while (my $line = <FH>) {
    	chomp $line;
	if ($line =~ /^>\s*(\S+)/) {
	    $name = $1;
	} else {
	    $line_count++;
	    $matrix_string .= "$line\n";
	    if ($line_count == 4) {
		push @pfms, TFBS::Matrix::PFM->new(
                                    -matrixstring       => $matrix_string,
                                    -name               => $name,
                                    -ID                 => $name);
	    	$line_count = 0;
	    	$name = '';
		$matrix_string = '';
	    }
	}
    }
    close(FH);

    return @pfms ? \@pfms : undef;
}
