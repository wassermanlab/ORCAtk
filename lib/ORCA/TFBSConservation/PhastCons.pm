=head1 NAME

ORCA::TFBSConservation::PhastCons - Object for facilitating analysis of TFBSs
in (phastCons) conserved regions of a DNA sequence.

=head1 SYNOPSIS

    use ORCA::TFBSConservation::PhastCons;

    # Create a new ORCA::TFBSConservation::PhastCons object.
    $tcphc = ORCA::TFBSConservation::PhastCons->new(
					    -conservation_analysis => $ca);

    # Find TFBSs which meet the conservation criteria.
    my $sites = $tcphc->find_conserved_tf_sites(
			    -matrix_set			=> $matrix_set,
			    -tfbs_threshold		=> $tfbs_threshold,
			    -min_tfbs_cr_overlap	=> $overlap);
			    -filter_overlapping_sites	=> 1,
			    -start			=> $start,
			    -end			=> $end);

    # Once computed, the conserved TFBSs can be obtained later...
    my $sites = $tcphc->conserved_tf_sites();
    
=head1 DESCRIPTION

ORCA::TFBSConservation::PhastCons is an object used for the purpose of finding
TFBS sites within (phastCons) conserved regions of a DNA sequence.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::TFBSConservation::PhastCons;

use strict;

use constant DEBUG => 0;

use constant DFLT_MIN_TFBS_CA_OVERLAP	=> 1;

use Carp;
use ORCA::ConservationAnalysis::PhastCons;
use TFBS::SitePairSet;

=head2 new

 Title    : new
 Usage    : $tcphc = ORCA::TFBSConservation::PhastCons->new(
					    -conservation_analysis => $ca);
 Function : Create a new ORCA::TFBSConservation::PhastCons object.
 Returns  : A new ORCA::TFBSConservation::PhastCons object.
 Args     : -conservation_analysis => An ORCA::ConservationAnalysis::PhastCons
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
    	carp "No ORCA::ConservationAnalysis::PhastCons object provided\n";
	return;
    }

    if (!$ca->isa("ORCA::ConservationAnalysis::PhastCons")) {
    	carp "-conservation_analysis argument is not an"
		. " ORCA::ConservationAnalysis::PhastCons object\n";
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
 Usage    : $sites = $tcphc->find_conserved_tf_sites(
			    -matrix_set			=> $matrix_set,
			    -tfbs_threshold		=> $tfbs_threshold,
			    -min_tfbs_cr_overlap	=> 1
			    -filter_overlapping_sites	=> 1,
			    -start			=> $start,
			    -end			=> $end);
 Function : Compute and return the list of conserved TF sites.
 Returns  : A reference to a TFBS::SiteSet object.
 Args     : -matrix_set			=> a TFBS::MatrixSet object 
	    -tfbs_threshold		=> Report only TFBSs with at least
	    				   this score
	    -min_tfbs_cr_overlap	=> Report only TFBSs which overlap a
					   conserved region by at least this
					   number of nucleotides.
	    -filter_overlapping_sites	=> boolean; if true, for any given
	    				   matrix, report only the best
					   scoring of any two overlapping
					   site pairs
	    -start			=> restrict TFBS search region to
	    				   start at position $start
	    -end			=> restrict TFBS search region to
	    				   end at position $end

=cut

sub find_conserved_tf_sites
{
    my ($self, %args) = @_;

    my $matrix_set = $self->matrix_set($args{-matrix_set});
    if (!defined $matrix_set || !$matrix_set->isa("TFBS::MatrixSet")) {
	carp "matrix_set not provided or is not a TFBS::MatrixSet\n";
	return;
    }

    if ($matrix_set->size == 0) {
        carp "matrix set is empty\n";
        return;
    }

    my $filter_overlaps = $args{-filter_overlapping_sites} || 0;

    my $ca = $self->conservation_analysis;
    if (!defined $ca) {
        carp "ORCA::ConservationAnalysis::PhastCons object not set\n";
        return;
    }

    if (!$ca->isa("ORCA::ConservationAnalysis::PhastCons")) {
	carp "first argument is not a"
		. " ORCA::ConservationAnalysis::PhastCons object\n";
	return;
    }

    my $seq = $ca->seq;
    if (!defined $seq) {
        carp "PhastConsAnalysis object has no sequence\n";
        return;
    }

    my $start = $args{-start};
    my $end = $args{-end};
    $start = 1 if !$start;
    $end = $seq->length if !$end;

    my $regions = $ca->conserved_regions;
    if (!defined $regions) {
        carp "PhastConsAnalysis object has no conserved regions\n";
        return;
    }
    
    my $overlap = $self->min_tfbs_cr_overlap($args{-min_tfbs_cr_overlap}
						|| DFLT_MIN_TFBS_CA_OVERLAP);

    my $threshold = $self->tfbs_threshold($args{-tfbs_threshold});
    if (!defined $threshold) {
    	carp "-tfbs_threshold argument not provided\n";
        return;
    }

    my $cutoff = $self->conservation_cutoff;
    if ($cutoff) {
    	if ($cutoff =~ /(.+)%/) {
	    $cutoff = $1 / 100;
	}
	if ($cutoff < 0 || $cutoff > 1) {
	    carp "conservation cutoff is out of range; please specify as a
		    number between 0 and 1 or as a string between 0% and 100%";
	    return;
	}
    }

    my $cons_tf_count = 0;

    my $sites = TFBS::SiteSet->new();

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

	my $m_sites = TFBS::SiteSet->new();
	foreach my $region (@$regions) {
	    if ($cutoff) {
		next if $region->score < $cutoff;
	    }

	    my $reg_start = $region->start;
	    my $reg_end = $region->end;

	    #
	    # If conserved region falls outside the sequence region of
	    # interest exclude it from analysis
	    #
	    next if $reg_start > $end;
	    next if $reg_end < $start;

	    my $m_start = $reg_start - $matrix_width + $overlap;
	    $m_start = $start if $m_start < $start;

	    my $m_end = $reg_end + $matrix_width - $overlap;
	    $m_end = $end if $m_end > $end;

	    next if $m_end - $m_start + 1 < $matrix_width;

	    if (DEBUG) {
		printf "\nRegion %d  %d  %.2f\n",
		    $reg_start,
		    $reg_end,
		    $region->score;
	    }

	    my $reg_m_sites = _find_seq_tf_sites($matrix, $seq,
					$m_start, $m_end,
					$threshold);

	    next if !$reg_m_sites || $reg_m_sites->size == 0;
	    if (DEBUG) {
		_print_site_set($reg_m_sites) if defined $reg_m_sites;
	    }

	    $self->_set_site_conservation($reg_m_sites, $region->score);

	    $m_sites->add_siteset($reg_m_sites);
	}

	next if !$m_sites || $m_sites->size == 0;
	if ($filter_overlaps) {
	    my $filt_m_sites = _tf_filter_overlapping_sites($m_sites);
	    if (DEBUG) {
		print "\nFiltered site set\n";
		_print_site_set($filt_m_sites)
	    }
	    $sites->add_siteset($filt_m_sites);
	} else {
	    $sites->add_siteset($m_sites);
	}
    }

    return $sites;
}

################################################################################
#
# Getter/Setter methods follow.
#
################################################################################

=head2 matrix_set

 Title    : matrix_set
 Usage    : $matrix_set = $tcphc->matrix_set($matrix_set);
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
 Usage    : $score = $tcphc->min_base_tfbs_score($score);
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
 Usage    : $score = $tcphc->min_comparison_tfbs_score($score);
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
 Usage    : $threshold = $tcphc->tfbs_threshold($threshold);
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
 Usage    : $cutoff = $tcphc->conservation_cutoff($cutoff);
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
 Usage    : $overlap = $tcphc->min_tfbs_cr_overlap($overlap);
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
 Usage    : $phca = $tcphc->conservation_analysis();
 Function : Get the phastcons conservation analysis object.
 Returns  : A ORCA::ConservationAnalyis::PhastCons object.
 Args     : None.

=cut

sub conservation_analysis
{
    $_[0]->{-conservation_analysis};
}

=head2 conserved_tf_sites

 Title    : conserved_tf_sites
 Usage    : $sites = $tcphc->conserved_tf_sites();
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
# XXX This may have to be revisited for more sophisticated filtering.
# Take a TFBS::SiteSet and filter overlapping sites such that only
# the highest scoring site of any mutually overlapping sites is kept.
# In the event that sites score equally, the first site is kept, i.e.
# bias is towards the site with the lowest starting position.
#
sub _tf_filter_overlapping_sites
{
    my ($sites) = @_;

    return if !defined $sites || $sites->size == 0;

    my $filtered_sites = TFBS::SiteSet->new();

    my $iter = $sites->Iterator(-sort_by => 'start');
    my $prev_site = $iter->next;
    if ($prev_site) {
	while (my $site = $iter->next) {
	    if ($site->overlaps($prev_site)) {
		if ($site->score > $prev_site->score) {
		    $prev_site = $site;
		}
	    } else {
		$filtered_sites->add_site($prev_site);
		$prev_site = $site;
	    }
	}
	$filtered_sites->add_site($prev_site);
    }

    return $filtered_sites;
}

#
# For each site in the given site set, set the conservation score as a
# tag value for the site.
#
sub _set_site_conservation
{
    my ($self, $site_set, $score) = @_;

    return if !$site_set || $site_set->size == 0;

    my $iter = $site_set->Iterator();
    return if !$iter;
    while (my $site = $iter->next) {
	$site->add_tag_value('conservation', $score);
    }
}

#
# For debugging purposes
#
sub _print_site_set
{
    my ($site_set) = @_;

    my $iter = $site_set->Iterator(-sort_by => 'start');

    while (my $site = $iter->next) {
	printf "%-12s %5d %5d %2d %0.2f\n",
		$site->pattern->name,
		$site->start,
		$site->end,
		$site->strand,
		$site->rel_score * 100;
    }
}

1;
