=head1 NAME

ORCA::TFBSConservation::Pairwise - Object for facilitating
analysis of TFBSs in conserved regions of a DNA sequence.

=head1 SYNOPSIS

    use ORCA::TFBSConservation::Pairwise;

    # Create a new ORCA::TFBSConservation::Pairwise object.
    $tcpw = ORCA::TFBSConservation::Pairwise->new(
    			-conservation_analysis	=> $ca);

    # Find TFBSs which meet the conservation criteria.
    my $sites = $tcpw->find_conserved_tf_sites(
			-matrix_set			=> $matrix_set,
			-tfbs_threshold			=> $tfbs_threshold,
			-min_tfbs_cr_overlap		=> 1,
			-filter_overlapping_sites	=> 1,
			-start				=> $start,
			-end				=> $end);

    # Once computed, the conserved TFBSs can be obtained later...
    my $sites = $tcpw->conserved_tf_sites();
    
=head1 DESCRIPTION

ORCA::TFBSConservation::Pairwise is an object used for the purpose of
finding TFBS sites within conserved regions of two aligned orthologous
DNA sequences.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::TFBSConservation::Pairwise;

use strict;

use constant DEBUG => 0;

use constant DFLT_MIN_TFBS_CA_OVERLAP	=> 1;

use Carp;
use ORCA::ConservationAnalysis::Pairwise;
use TFBS::SitePairSet;

=head2 new

 Title    : new
 Usage    : $tcpw = ORCA::TFBSConservation::Pairwise->new(
					    -conservation_analysis => $ca);
 Function : Create a new ORCA::TFBSConservation::Pairwise object.
 Returns  : A new ORCA::TFBSConservation::Pairwise object.
 Args     : -conservation_analysis = reference to an
 				     ORCA::ConservationAnalysis::Pairwise
				     object

=cut

sub new
{
    my ($class, %args) = @_;
    my $self = bless {
			%args
		    }, ref $class || $class;

    my $ca = $self->conservation_analysis;

    if (!defined $ca) {
	carp "No ORCA::ConservationAnalysis::Pairwise object provided\n";
	return;
    }

    if (!$ca->isa("ORCA::ConservationAnalysis::Pairwise")) {
	carp "-conservation_analysis argument is not an"
		. " ORCA::ConservationAnalysis::Pairwise object\n";
	return;
    }


    return $self;
}

################################################################################
#
# Analysis methods follow.
#
################################################################################

=head2 find_conserved_tf_sites

 Title    : find_conserved_tf_sites
 Usage    : $sites = $tcpw->find_conserved_tf_sites(
			-matrix_set			=> $matrix_set,
			-tfbs_threshold			=> $tfbs_thresh,
			-min_tfbs_cr_overlap		=> $overlap,
			-filter_overlapping_sites	=> 1,
			-start				=> $start,
			-end				=> $end);
 Function : Compute and return the list of conserved TF sites.
 Returns  : A reference to a TFBS::SitePairSet object.
 Args     :
	    -matrix_set			=> A TFBS::MatrixSet object
	    -tfbs_threshold		=> Report only TFBSs with at least
	    				   this score on both sequences
	    -min_tfbs_cr_overlap	=> Report only TFBSs which overlap
					   a conserved region by at least
					   this number of nucleotides.
	    -filter_overlapping_sites	=> boolean; if true, for any given
	    				   matrix, report only the best
					   scoring of any two overlapping
					   site pairs
	    -start			=> restrict TFBS search region to
					   start at this position
	    -end			=> restrict TFBS search region to
					   end at this position

=cut

sub find_conserved_tf_sites
{
    my ($self, %args) = @_;

    my $matrix_set = $self->matrix_set($args{-matrix_set});
    if (!$matrix_set || !$matrix_set->isa("TFBS::MatrixSet")) {
	carp "-matrix_set argument not provided or value is not a"
		. " TFBS::MatrixSet object\n";
	return;
    }

    if ($matrix_set->size == 0) {
	carp "matrix set is empty\n";
	return;
    }

    my $threshold = $self->tfbs_threshold($args{-tfbs_threshold});
    if (!defined $threshold) {
	carp "-tfbs_threshold argument not provided\n";
    }
    
    my $overlap = $self->min_tfbs_cr_overlap(
					$args{-min_tfbs_cr_overlap}
					|| DFLT_MIN_TFBS_CA_OVERLAP);

    my $filter_overlaps = $args{-filter_overlapping_sites} || 0;

    my $start1 = $args{-start};
    my $end1 = $args{-end};

    my $ca = $self->conservation_analysis;
    if (!defined $ca) {
	carp "ORCA::ConservationAnalysis::Pairwise object has not been set\n";
	return;
    }

    my $seq1 = $ca->base_seq;
    if (!defined $seq1) {
	carp "ORCA::ConservationAnalysis::Pairwise object has no base"
		. " sequence\n";
	return;
    }

    my $seq2 = $ca->comparison_seq;
    if (!defined $seq2) {
	carp "ORCA::ConservationAnalysis::Pairwise object has no comparison"
		. " sequence\n";
	return;
    }

    my $align = $ca->alignment;
    if (!defined $align) {
	carp "ORCA::ConservationAnalysis::Pairwise object has no aligment\n";
	return;
    }

    my $regions = $ca->conserved_regions;
    if (!defined $regions) {
	carp "ORCA::ConservationAnalysis::Pairwise object has no conserved"
		. " regions\n";
	return;
    }

    my $cutoff = $ca->conserved_regions_report->param('cutoff');
    if (defined $cutoff) {
	$cutoff /= 100;
	if ($cutoff < 0 || $cutoff > 1) {
	    carp "ERROR: conservation cutoff is out of range\n";
	    return;
	}
    }

    my $aln_orient = $ca->alignment_orientation;

    $start1 = 1 if !$start1;
    $end1 = $seq1->length if !$end1;
    my $start2;
    my $end2;
    if ($aln_orient == -1) {
	$start2 = $ca->convert_seq_pos(1, 2, $end1, 'le');
	$end2 = $ca->convert_seq_pos(1, 2, $start1, 'ge');
    } else {
	$start2 = $ca->convert_seq_pos(1, 2, $start1, 'le');
	$end2 = $ca->convert_seq_pos(1, 2, $end1, 'ge');
    }
    $start2 = 1 if !$start2;
    $end2 = $seq2->length if !$end2 || $end2 > $seq2->length;

    my $cons_tf_count1 = 0;
    my $cons_tf_count2 = 0;

    my $site_pairs = TFBS::SitePairSet->new();

    #
    # Changed the logic below so that matrix loop is within the region loop
    # as opposed to the other way around. This assumes that a TFBS may only
    # be contained within a single region. This assumption makes certain
    # computations simpler and thus much more efficient.
    # Previously the case was considered where a TFBS
    # could potentially overlap more than one conserved region. However this
    # is extremely unlikely (impossible?) as it would require two conserved
    # regions to be within a number of base pairs of each other which is less
    # than the width of the TFBS without them having been combined into a single
    # region. This would seem to imply that the intervening short region is of
    # such poor conservation that there would be no chance of finding a TFBS
    # encompassing it anyway. DJA 050329
    #
#    foreach my $region (@$regions) {
#	if ($cutoff) {
#	    next if $region->feature1->score < $cutoff;
#	}
#
#	my $reg_start1 = $region->feature1->start;
#	my $reg_end1 = $region->feature1->end;
#
#	#
#	# If conserved region falls outside the sequence region of
#	# interest exclude it from analysis
#	#
#	next if $reg_start1 > $end1;
#	next if $reg_end1 < $start1;
#
#	my $reg_start2 = $region->feature2->start;
#	if (!$reg_start2) {
#	    if ($aln_orient == -1) {
#		$reg_start2 = $ca->convert_seq_pos(
#					1, 2, $reg_end1, 'le');
#	    } else {
#		$reg_start2 = $ca->convert_seq_pos(
#					1, 2, $reg_start1, 'le');
#	    }
#	}
#	$reg_start2 = $start2 if !$reg_start2
#					|| $reg_start2 < $start2;
#
#	my $reg_end2 = $region->feature2->end;
#	if (!$reg_end2) {
#	    if ($aln_orient == -1) {
#		$reg_end2 = $ca->convert_seq_pos(
#					1, 2, $reg_start1, 'ge');
#	    } else {
#		$reg_end2 = $ca->convert_seq_pos(
#					1, 2, $reg_end1, 'ge');
#	    }
#	}
#	$reg_end2 = $end2 if !$reg_end2 || $reg_end2 > $end2;
#
#	if (DEBUG) {
#	    printf "\nRegion %d  %d  %d  %d  %.2f\n",
#		$reg_start1,
#		$reg_end1,
#		$reg_start2,
#		$reg_end2,
#		$region->feature1->score;
#	}
#
#	my $ms_iter = $matrix_set->Iterator();
#	while (my $matrix = $ms_iter->next) {
#	    my $matrix_width = $matrix->length;
#
#	    my $m_start1 = $reg_start1 - $matrix_width + $overlap;
#	    $m_start1 = $start1 if $m_start1 < $start1;
#
#	    my $m_end1 = $reg_end1 + $matrix_width - $overlap;
#	    $m_end1 = $end1 if $m_end1 > $end1;
#
#	    next if $m_end1 - $m_start1 + 1 < $matrix_width;
#
#	    my $m_start2 = $reg_start2 - $matrix_width + $overlap;
#	    $m_start2 = $start2 if $m_start2 < $start2;
#
#	    my $m_end2 = $reg_end2 + $matrix_width - $overlap;
#	    $m_end2 = $end2 if $m_end2 > $end2;
#
#	    next if $m_end2 - $m_start2 + 1 < $matrix_width; 
#
#	    if (DEBUG) {
#		printf "\nTranscription Factor %s  %s\n",
#			    $matrix->name, $matrix->class;
#	    }
#
#	    my $m_sites1 = _find_seq_tf_sites($matrix, $seq1,
#					$m_start1, $m_end1,
#					$threshold);
#
#	    next if !$m_sites1 || $m_sites1->size == 0;
#	    if (DEBUG) {
#		_print_site_set($ca, 1, $m_sites1) if defined $m_sites1;
#	    }
#
#	    my $m_sites2 = _find_seq_tf_sites($matrix, $seq2,
#					$m_start2, $m_end2,
#					$threshold);
#
#	    next if !$m_sites2 || $m_sites2->size == 0;
#	    if (DEBUG) {
#		_print_site_set($ca, 2, $m_sites2) if defined $m_sites2;
#	    }
#
#	    my $m_site_pairs = _tf_sites_to_site_pairs($ca, $m_sites1,
#							$m_sites2, 1);
#
#	    next if !$m_site_pairs || $m_site_pairs->size == 0;
#	    if (DEBUG) {
#		print "\nUnfiltered matrix site pair set\n";
#		_print_site_pair_set($m_site_pairs)
#	    }
#
#	    if ($filter_overlaps) {
#		my $filt_m_site_pairs = _tf_filter_overlapping_site_pairs(
#								$m_site_pairs);
#		if (DEBUG) {
#		    print "\nFiltered matrix site pair set\n";
#		    _print_site_pair_set($filt_m_site_pairs)
#		}
#		#
#		# We are now assuming that a TFBS can only be contained within
#		# the current region, so set the conservation score to that
#		# of the region (as opposed to the more complex computation
#		# below). DJA 050329
#		#
#		$self->_set_site_pair_conservation($filt_m_site_pairs,
#						    $region->feature1->score);
#		$site_pairs->add_site_pair_set($filt_m_site_pairs);
#	    } else {
#		#
#		# We are now assuming that a TFBS can only be contained within
#		# the current region, so set the conservation score to that
#		# of the region (as opposed to the more complex computation
#		# below). DJA 050329
#		#
#		$self->_set_site_pair_conservation($m_site_pairs,
#						    $region->feature1->score);
#		$site_pairs->add_site_pair_set($m_site_pairs);
#	    }
#	}
#    }
#
#    #
#    # We are no longer considering that a TFBS can be contained within more
#    # than one search region so remove this more complex computation in favour
#    # of the simpler one above. DJA 050329
#    #
#    #$self->_set_max_conservation($site_pairs, $ca);


    #
    # The logic above did not hold as there were cases (albeit rare) where
    # TFBSs did overlap two adjacent conserved regions resulting in that
    # TFBS being reported twice. Thus go back to having the regions loop
    # within the matrix loop. However put the site pairs together for each
    # region separately (much more efficient) but then modify the overlap
    # filtering routine so that it also removes these duplicates based on
    # conservation (more conserved one is preserved).
    #
    my $ms_iter = $matrix_set->Iterator();
    while (my $matrix = $ms_iter->next) {
	if (DEBUG) {
	    printf "\nTranscription Factor %s  %s\n",
			$matrix->name, $matrix->class;
	}

	my $matrix_width = $matrix->length;

	my $m_site_pairs = TFBS::SitePairSet->new();
	foreach my $region (@$regions) {
	    if ($cutoff) {
		next if $region->feature1->score < $cutoff;
	    }

	    my $reg_start1 = $region->feature1->start;
	    my $reg_end1 = $region->feature1->end;

	    #
	    # If conserved region falls outside the sequence region of
	    # interest exclude it from analysis
	    #
	    next if $reg_start1 > $end1;
	    next if $reg_end1 < $start1;

	    my $m_start1 = $reg_start1 - $matrix_width + $overlap;
	    $m_start1 = $start1 if $m_start1 < $start1;

	    my $m_end1 = $reg_end1 + $matrix_width - $overlap;
	    $m_end1 = $end1 if $m_end1 > $end1;

	    next if $m_end1 - $m_start1 + 1 < $matrix_width;

	    my $reg_start2 = $region->feature2->start;
	    if (!$reg_start2) {
		if ($aln_orient == -1) {
		    $reg_start2 = $ca->convert_seq_pos(
					    1, 2, $reg_end1, 'le');
		} else {
		    $reg_start2 = $ca->convert_seq_pos(
					    1, 2, $reg_start1, 'le');
		}
	    }
	    $reg_start2 = $start2 if !$reg_start2
					    || $reg_start2 < $start2;

	    my $reg_end2 = $region->feature2->end;
	    if (!$reg_end2) {
		if ($aln_orient == -1) {
		    $reg_end2 = $ca->convert_seq_pos(
					    1, 2, $reg_start1, 'ge');
		} else {
		    $reg_end2 = $ca->convert_seq_pos(
					    1, 2, $reg_end1, 'ge');
		}
	    }
	    $reg_end2 = $end2 if !$reg_end2 || $reg_end2 > $end2;

	    my $m_start2 = $reg_start2 - $matrix_width + $overlap;
	    $m_start2 = $start2 if $m_start2 < $start2;

	    my $m_end2 = $reg_end2 + $matrix_width - $overlap;
	    $m_end2 = $end2 if $m_end2 > $end2;

	    next if $m_end2 - $m_start2 + 1 < $matrix_width; 

	    if (DEBUG) {
		printf "\nRegion %d  %d  %d  %d  %.2f\n",
		    $reg_start1,
		    $reg_end1,
		    $reg_start2,
		    $reg_end2,
		    $region->feature1->score;
	    }

	    my $reg_m_sites1 = _find_seq_tf_sites($matrix, $seq1,
					$m_start1, $m_end1,
					$threshold);

	    next if !$reg_m_sites1 || $reg_m_sites1->size == 0;
	    if (DEBUG) {
		_print_site_set($ca, 1, $reg_m_sites1) if defined $reg_m_sites1;
	    }

	    my $reg_m_sites2 = _find_seq_tf_sites($matrix, $seq2,
					$m_start2, $m_end2,
					$threshold);

	    next if !$reg_m_sites2 || $reg_m_sites2->size == 0;
	    if (DEBUG) {
		_print_site_set($ca, 2, $reg_m_sites2) if defined $reg_m_sites2;
	    }

	    my $reg_m_site_pairs = _tf_sites_to_site_pairs($ca, $reg_m_sites1,
							    $reg_m_sites2, 1);

	    next if !$reg_m_site_pairs || $reg_m_site_pairs->size == 0;
	    if (DEBUG) {
		print "\nRegion matrix site pair set\n";
		_print_site_pair_set($reg_m_site_pairs)
	    }

	    $self->_set_site_pair_conservation($reg_m_site_pairs,
						$region->feature1->score);

	    $m_site_pairs->add_site_pair_set($reg_m_site_pairs);
	}

	next if !$m_site_pairs || $m_site_pairs->size == 0;
	if ($filter_overlaps) {
	    my $filt_m_site_pairs = _tf_filter_overlapping_site_pairs(
							$m_site_pairs);
	    if (DEBUG) {
		print "\nFiltered site pair set\n";
		_print_site_pair_set($filt_m_site_pairs)
	    }
	    $site_pairs->add_site_pair_set($filt_m_site_pairs);
	} else {
	    $site_pairs->add_site_pair_set($m_site_pairs);
	}
    }

    #print "seq1 conserved TF sites: $cons_tf_count1\n";
    #print "seq2 conserved TF sites: $cons_tf_count2\n";

    return $site_pairs;
}

################################################################################
#
# Getter/Setter methods follow.
#
################################################################################

=head2 matrix_set

 Title    : matrix_set
 Usage    : $matrix_set = $tcpw->matrix_set($matrix_set);
 Function : Get/Set the TFBS::MatrixSet.
 Returns  : A TFBS::MatrixSet object.
 Args     : Optionally a new TFBS::MatrisSet object.

=cut

sub matrix_set
{
    my ($self, $matrix_set) = @_;

    if (defined ($matrix_set)) {
	if (!$matrix_set->isa("TFBS::MatrixSet")) {
	    carp "matrix_set is not a TFBS::MatrixSet object\n";
	    return undef;
	}
	$self->{-matrix_set} = $matrix_set;
    }

    return $self->{-matrix_set};
}

=head2 min_base_tfbs_score

 Title    : min_base_tfbs_score
 Usage    : $score = $tcpw->min_base_tfbs_score($score);
 Function : Get/Set the minimum base sequence TFBS score.
 Returns  : An integer.
 Args     : Optionally a new minimum base sequence TFBS score.

=cut

sub min_base_tfbs_score
{
    my ($self, $score) = @_;

    if (defined ($score)) {
	$self->{-min_base_tfbs_score} = $score;
    }

    return $self->{-min_base_tfbs_score};
}

=head2 min_comparison_tfbs_score

 Title    : min_comparison_tfbs_score
 Usage    : $score = $tcpw->min_comparison_tfbs_score($score);
 Function : Get/Set the minimum comparison sequence TFBS score.
 Returns  : An integer.
 Args     : Optionally a new minimum comparison sequence TFBS score.

=cut

sub min_comparison_tfbs_score
{
    my ($self, $score) = @_;

    if (defined ($score)) {
	$self->{-min_comparison_tfbs_score} = $score;
    }

    return $self->{-min_comparison_tfbs_score};
}

=head2 tfbs_threshold

 Title    : tfbs_threshold
 Usage    : $threshold = $tcpw->tfbs_threshold($threshold);
 Function : Get/Set the TFBS threshold.
 Returns  : An integer.
 Args     : Optionally a new TFBS threshold.

=cut

sub tfbs_threshold
{
    my ($self, $threshold) = @_;

    if (defined ($threshold)) {
	$self->{-tfbs_threshold} = $threshold;
    }

    return $self->{-tfbs_threshold};
}

=head2 conservation_cutoff

 Title    : conservation_cutoff
 Usage    : $cutoff = $tcpw->conservation_cutoff($cutoff);
 Function : Get/Set the conservation cutoff.
 Returns  : An real number between 0 and 100.
 Args     : Optionally a new conservation cutoff.

=cut

sub conservation_cutoff
{
    my ($self, $cutoff) = @_;

    if (defined ($cutoff)) {
	$self->{-conservation_cutoff} = $cutoff;
    }

    return $self->{-conservation_cutoff};
}

=head2 min_tfbs_cr_overlap

 Title    : min_tfbs_cr_overlap
 Usage    : $overlap = $tcpw->min_tfbs_cr_overlap($overlap);
 Function : Get/Set the minimum TFBS conservation overlap.
 Returns  : An integer.
 Args     : Optionally a new minimum TFBS conservation overlap.

=cut

sub min_tfbs_cr_overlap
{
    my ($self, $overlap) = @_;

    if (defined ($overlap)) {
	$self->{-min_tfbs_cr_overlap} = $overlap;
    }

    return $self->{-min_tfbs_cr_overlap};
}

=head2 conservation_analysis

 Title    : conservation_analysis
 Usage    : $ca = $tcpw->conservation_analysis();
 Function : Get the conservation analysis object.
 Returns  : A ORCA::ConservationAnalysis::Pairwise object.
 Args     : None.

=cut

sub conservation_analysis
{
    my ($self, $ca) = @_;

    if (defined $ca) {
	$self->{-conservation_analysis} = $ca;
    }

    return $self->{-conservation_analysis};
}

=head2 conserved_tf_sites

 Title    : conserved_tf_sites
 Usage    : $sites = $tcpw->conserved_tf_sites();
 Function : Get the conserved TF sites.
 Returns  : A TFBS::SitePairSet reference.
 Args     : None.

=cut

sub conserved_tf_sites
{
    $_[0]->{-conserved_tf_sites};
}

################################################################################
#
# Internal methods follow.
#
################################################################################

sub _find_seq_tf_sites
{
    my ($matrix, $seqobj, $start, $end, $threshold) = @_;

    my %search_args = (
	    -seqobj	=> $seqobj,
	    -threshold	=> $threshold,
	    -subpart	=> {
			    -start		=> $start,
			    -end		=> $end
	    },
    );

    return $matrix->search_seq(%search_args);
}

#
# XXX not used
#
# XXX This may have to be revisited for more sophisticated filtering.
# Take a TFBS::SiteSet and filter overlapping sites such that only
# the highest scoring site of any mutually overlapping sites is kept.
# In the event that sites score equally, the first site is kept, i.e.
# bias is towards the site with the lowest starting position.
#
sub _filter_overlapping_sites
{
    my ($sites) = @_;

    return if !defined $sites || $sites->size == 0;

    my $filtered_sites = TFBS::SiteSet->new();

    my $iter = $sites->Iterator(-sort_by => 'start');
    my $prev_site;
    while (my $site = $iter->next) {
	if (defined $prev_site) {
	    if ($site->overlaps($prev_site)) {
		if ($site->score > $prev_site->score) {
		    $prev_site = $site;
		}
	    } else {
		$filtered_sites->add_site($prev_site);
	    }
	} else {
	    $prev_site = $site;
	}
    }
    $filtered_sites->add_site($prev_site);

    return $filtered_sites;
}

#
# This may have to be revisited for more sophisticated filtering.
# Take a TFBS::SitePairSet where each site pair in the set corresponds to the
# same transcription factor and filter overlapping site pairs such that only
# the highest scoring site pair of any mutually overlapping site pairs is kept.
# In the event that site pairs score equally, the first site pair is kept, i.e.
# bias is towards the site pair with the lowest starting position.
#
sub _tf_filter_overlapping_site_pairs
{
    my ($site_pairs) = @_;

    return if !defined $site_pairs || $site_pairs->size == 0;

    my $filtered_pairs = TFBS::SitePairSet->new();

    my $iter = $site_pairs->Iterator(-sort_by => 'start');
    my $prev_pair = $iter->next;
    if ($prev_pair) {
	while (my $pair = $iter->next) {
	    if ($pair->overlaps($prev_pair)) {
		#
		# Bias is toward the site pair with the lower start
		# site (i.e. if the scores are equal).
		# 
		if (_cmp_site_pair_score($pair, $prev_pair) > 0) {
		    $prev_pair = $pair;
		}
	    } else {
		$filtered_pairs->add_site_pair($prev_pair);
		$prev_pair = $pair;
	    }
	}
	$filtered_pairs->add_site_pair($prev_pair);
    }

    return $filtered_pairs;
}

#
# XXX - no longer used
#
# This may have to be revisited for more sophisticated filtering.
# Take a TFBS::SitePairSet and filter overlapping site pairs such that only
# the highest scoring site pair of any mutually overlapping site pairs is kept.
# In the event that site pairs score equally, the first site pair is kept, i.e.
# bias is towards the site pair with the lowest starting position.
#
sub _filter_overlapping_site_pairs
{
    my ($site_pairs) = @_;

    return if !defined $site_pairs || $site_pairs->size == 0;

    #
    # Create a temp. list of all site pairs. Process this list such that for
    # any mutually overlapping site pairs corresponding to the same pattern,
    # the one with the highest score is kept.
    #
    my @temp_pairs;
    my $iter = $site_pairs->Iterator();
    while (my $pair = $iter->next) {
	push @temp_pairs, $pair;
    }
    for (my $i = 0; $i < (scalar @temp_pairs) - 1; $i++) {
	for (my $j = $i + 1; $j < scalar @temp_pairs; $j++) {
	    if (defined $temp_pairs[$i] && defined $temp_pairs[$j]) {
		if ($temp_pairs[$i]->site1->pattern->name
			eq $temp_pairs[$j]->site1->pattern->name)
		{
		    if ($temp_pairs[$i]->site1->overlaps(
			    $temp_pairs[$j]->site1))
		    {
			#
			# Bias is toward the site pair with the lower start
			# site (i.e. if the scores are equal). Note that
			# strand is not taken into account so this is not
			# necessarily the most upstream site.
			# 
			if (_cmp_site_pair_score(
				$temp_pairs[$i], $temp_pairs[$j]) < 0)
			{
			    $temp_pairs[$i] = undef;
			} else {
			    $temp_pairs[$j] = undef;
			}
		    }
		}
	    }
	}
    }

    my $filtered_pairs = TFBS::SitePairSet->new();
    my $i = 0;
    while ($i < scalar @temp_pairs) {
	if (defined $temp_pairs[$i]) {
	    $filtered_pairs->add_site_pair($temp_pairs[$i]);
	}
	$i++;
    }

    return $filtered_pairs;
}

#
# Take two site sets of a given transcription factor and combine them into
# a site pair set by reconciling their positions within the alignment.
#
sub _tf_sites_to_site_pairs
{
    my ($ca, $site_set1, $site_set2, $exact) = @_;

    return if !defined $ca || !defined $site_set1 || !defined $site_set2
    		|| $site_set1->size == 0 || $site_set2->size == 0;

    my $aln_orient = $ca->alignment_orientation;

    my $site_pair_set = TFBS::SitePairSet->new();

    if ($exact) {
	my $iter1 = $site_set1->Iterator();
	while (my $site1 = $iter1->next) {
	    my $site1_aln_start = $ca->seq_to_align_pos(1, $site1->start);
	    my $iter2 = $site_set2->Iterator();
	    while (my $site2 = $iter2->next) {
		my $site2_aln_start;
		if ($aln_orient == -1) {
		    $site2_aln_start = $ca->seq_to_align_pos(2, $site2->end);
		} else {
		    $site2_aln_start = $ca->seq_to_align_pos(2, $site2->start);
		}
		if ($site1_aln_start == $site2_aln_start) {
		    $site_pair_set->add_site_pair(
					TFBS::SitePair->new(
							    $site1,
							    $site2));
		}
	    }
	}

#	my $iter1 = $site_set1->Iterator(-sort_by => 'start');
#	my $iter2;
#	if ($cmp_strand == -1 || $cmp_strand eq '-') {
#	    $iter2 = $site_set2->Iterator(-sort_by => 'end');
#	} else {
#	    $iter2 = $site_set2->Iterator(-sort_by => 'start');
#	}
#	my $site1 = $iter1->next;
#	my $site2 = $iter2->next;
#	while (defined $site1 && defined $site2) {
#	    my $site1_aln_start = $ca->seq_to_align_pos(1, $site1->start);
#	    my $site2_aln_start;
#	    if ($cmp_strand == -1 || $cmp_strand eq '-') {
#		$site2_aln_start = $ca->seq_to_align_pos(2, $site2->end);
#	    } else {
#		$site2_aln_start = $ca->seq_to_align_pos(2, $site2->start);
#	    }
#	    #
#	    # Note: added $site1->strand <=> $site2->strand check as we
#	    # no longer care about +/- and -/+ matches.
#	    #
#	    my $cmp = $site1_aln_start <=> $site2_aln_start;
#	    #|| $site1->strand <=> $site2->strand;
#	    if ($cmp == 0) {
#		$site_pair_set->add_site_pair(
#				    TFBS::SitePair->new($site1, $site2));
#		
#		#
#		# We have to take into account different strand combos
#		# with same start position, i.e. +/+, +/-, -/+, -/+
#		# OK we no longer want to deal with +/- or -/+ and the
#		# addition of $site1->strand <=> $site2->strand to the $cmp
#		# should take care of everything so we no longer have to do
#		# the following.
#		#
#		my $prev_site1 = $site1;
#		my $prev_site2 = $site2;
#		my $prev_site1_aln_start = $site1_aln_start;
#		my $prev_site2_aln_start = $site2_aln_start;
#		
#		$site1 = $iter1->next;
#		$site2 = $iter2->next;
#		#
#		# Note: the case where both of the new site start positions
#		# correspond to both the previous site start positions
#		# will be handled in the next pass through the loop.
#		#
#		if (defined $site1) {
#		    $site1_aln_start = $ca->seq_to_align_pos(1, $site1->start);
#		    if ($site1_aln_start == $prev_site2_aln_start) {
#			$site_pair_set->add_site_pair(
#			    TFBS::SitePair->new($site1, $prev_site2));
#		    }
#		}
#		if (defined $site2) {
#		    if ($cmp_strand == -1 || $cmp_strand eq '-') {
#			$site2_aln_start
#				= $ca->seq_to_align_pos(2, $site2->end);
#		    } else {
#			$site2_aln_start
#				= $ca->seq_to_align_pos(2, $site2->start);
#		    }
#		    if ($site2_aln_start == $prev_site1_aln_start) {
#			$site_pair_set->add_site_pair(
#			    TFBS::SitePair->new($prev_site1, $site2));
#		    }
#		}
#	    } elsif ($cmp < 0) {
#		$site1 = $iter1->next;
#	    } elsif ($cmp > 0) {
#		$site2 = $iter2->next;
#	    }
#	}
    } else {
	#
	# XXX - not implemented
	#
    }

    return $site_pair_set;
}

#
# For each site pair set the conservation score as a tag value to site1.
#
sub _set_site_pair_conservation
{
    my ($self, $sps, $score) = @_;

    return if !$sps || $sps->size == 0;

    my $iter = $sps->Iterator();
    return if !$iter;
    while (my $sp = $iter->next) {
	$sp->site1->add_tag_value('conservation', $score);
    }
}

#
# For each site pair add the maximum conservation score as a tag value to
# site1. The maximum convervation score is defined as the maximum of each of
# the conservation scores of any conserved regions which the site pair
# overlaps by the min_tfbs_cr_overlap amount.
#
sub _set_max_conservation
{
    my ($self, $sites, $ca) = @_;

    return if !$sites || !$ca || $sites->size == 0;

    my $regions = $ca->conserved_regions;
    return if !$regions;

    my $min_ol = $self->min_tfbs_cr_overlap || 0;

    my $iter = $sites->Iterator();
    return if !$iter;
    while (my $sp = $iter->next) {
	my $site1 = $sp->site1;
	my $conservation = 0;
	foreach my $reg (@$regions) {
	    if (($site1->end - $reg->feature1->start + 1 >= $min_ol)
		    && ($reg->feature1->end - $site1->start + 1 >= $min_ol))
	    { 
		$conservation = $reg->feature1->score
				    if $reg->feature1->score > $conservation;
	    }
	}
	$site1->add_tag_value('conservation', $conservation);
    }
}

#
# XXX - no longer used
#
# Take two site sets and combine them into a site pair set by reconciling
# their names and positions within the alignment.
#
sub _sites_to_site_pairs
{
    my ($ca, $site_set1, $site_set2, $exact) = @_;

    return if !defined $ca || !defined $site_set1 || !defined $site_set2
    		|| $site_set1->size == 0 || $site_set2->size == 0;

    my $site_pair_set = TFBS::SitePairSet->new();

    if ($exact) {
	my $iter1 = $site_set1->Iterator();
	while (defined (my $site1 = $iter1->next)) {
	    my $iter2 = $site_set2->Iterator();
	    while (defined (my $site2 = $iter2->next)) {
		if ($site1->pattern->name eq $site2->pattern->name) {
		    my $site1_aln_start = $ca->seq_to_align_pos(
						    1, $site1->start);
		    my $site2_aln_start = $ca->seq_to_align_pos(
						    2, $site2->start);

		    #
		    # XXX Should we also look at the end positions???
		    #
		    if ($site1_aln_start == $site2_aln_start) {
			$site_pair_set->add_site_pair(
					    TFBS::SitePair->new(
								$site1,
								$site2));
		    }
		}
	    }
	}
    } else {
	#
	# XXX
	# This will create every possible pair of overlapping sites. We really
	# want to only create pairs where the comparison sequence site is the
	# highest scoring of all overlapping sites.
	#
	# Not currently using this anyway.
	#
	my $iter1 = $site_set1->Iterator();
	while (defined (my $site1 = $iter1->next)) {
	    my $iter2 = $site_set2->Iterator();
	    while (defined (my $site2 = $iter2->next)) {
		if ($site1->pattern->name eq $site2->pattern->name) {
		    my $site1_aln_start = $ca->seq_to_align_pos(
						    1, $site1->start);
		    my $site1_aln_end = $ca->seq_to_align_pos(
						    1, $site1->end);
		    my $site2_aln_start = $ca->seq_to_align_pos(
						    2, $site2->start);
		    my $site2_aln_end = $ca->seq_to_align_pos(
						    2, $site2->end);

		    if ($site1_aln_start <= $site2_aln_end
			    && $site1_aln_end >= $site2_aln_start)
		    {
			$site_pair_set->add_site_pair(
					 TFBS::SitePair->new(
							    $site1,
							    $site2));
		    }
		}
	    }
	}
    }

    return $site_pair_set;
}

#
# This may be changed in the future.
# Compares two site pairs. If the first pair 'scores higher' than the
# second pair, returns 1, if the scores are the same, returns 0 and if the
# first pair scores lower, returns -1.
#
# 'Scoring higher' is defined to mean the site pair with the highest TFBS score
# on the base sequence. Then, if the base sequence scores are equal, the
# highest TFBS score on the comparison sequence.
#
# Changed so that we now check for the case in which the overlapping pairs are
# actually the same pair falling into two adjacent conserved regions and thus
# counted twice. In this case, keep the one with the higher conservation score.
# Otherwise the behaviour is as before. If for example, two non-identical
# TFBS pairs overlap one another but fall into different conserved regions,
# count the higher scoring one. (Maybe we should count the more highly
# conserved in this case also).
# DJA 050331
#
sub _cmp_site_pair_score
{
    my ($pair1, $pair2) = @_;

    my $feat11 = $pair1->feature1;
    my $feat12 = $pair1->feature2;
    my $feat21 = $pair2->feature1;
    my $feat22 = $pair2->feature2;
    my $cons1 = ($feat11->get_tag_values('conservation'))[0];
    my $cons2 = ($feat21->get_tag_values('conservation'))[0];

    if (defined $cons1 && defined $cons2
	    && $cons1 != $cons2
	    && $feat11->start == $feat21->start
	    && $feat11->strand == $feat21->strand
	    && $feat12->strand == $feat22->strand)
    {
	#
	# Case of the same TFBS pair overlapping two different conserved
	# regions and being counted twice.
	#
	return ($cons1 <=> $cons2);
    }

    #
    # All other cases
    #
    return (($feat11->score <=> $feat21->score)
		|| ($feat12->score <=> $feat22->score));
}

#
# For debugging purposes
#
sub _print_site_set
{
    my ($ca, $site_seq, $site_set) = @_;

    my $seq_type = 'Unknown';
    if ($site_seq == 1) {
	$seq_type = 'Base';
    } elsif ($site_seq == 2) {
	$seq_type = 'Comparison';
    }

    my $iter = $site_set->Iterator(-sort_by => 'start');

    while (my $site = $iter->next) {
	printf "%-10s: %-12s %5d %5d %2d %5d %5d %0.2f\n",
		$seq_type,
		$site->pattern->name,
		$site->start,
		$site->end,
		$site->strand,
		$ca->seq_to_align_pos($site_seq, $site->start),
		$ca->seq_to_align_pos($site_seq, $site->end),
		$site->rel_score * 100;
    }
}

#
# For debugging purposes
#
sub _print_site_pair_set
{
    my ($site_pair_set) = @_;

    print "\n";
    my $iter = $site_pair_set->Iterator(-sort_by => 'start');
    while (my $pair = $iter->next) {
	my $site1 = $pair->site1;
	my $site2 = $pair->site2;
	printf "%-12s %5d %5d %2d %0.2f %-12s %5d %5d %2d %.2f\n",
		$site1->pattern->name,
		$site1->start,
		$site1->end,
		$site1->strand,
		$site1->rel_score * 100,
		$site2->pattern->name,
		$site2->start,
		$site2->end,
		$site2->strand,
		$site2->rel_score * 100;
    }
    print "\n";
}

1;
