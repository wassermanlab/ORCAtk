
=head1 NAME

ORCA::Analysis::Pairwise - Object for facilitating pairwise conservation
analysis of two aligned DNA sequences.

=head1 SYNOPSIS

    use ORCA::Analysis::Pairwise;

    # Minimally, create a new ORCA::Analysis::Pairwise object with
    # the given sequences.
    $pwa = ORCA::Analysis::Pairwise->new(
			-base_seq		    => $seq1,
			-comparison_seq		=> $seq2
    );

    # Optionally also provide positional information, masked sequences,
    # alignment and sequence exons.
    $pwa = ORCA::Analysis::Pairwise->new(
            -chr                    => $chr,
            -start                  => $start,
            -end                    => $end,
			-base_seq		        => $seq1,
			-comparison_seq		    => $seq2,
			-masked_base_seq	    => $masked_seq1,
			-masked_comparison_seq	=> $masked_seq2,
			-alignment		        => $align,
			-base_seq_exons	     	=> $exons1,
			-comparison_seq_exons	=> $exons2
    );

    # If not explicitly provided, compute the alignment.
    $pwa->compute_alignment();

    # Explicitly compute and retrieve the conservation profile.
    my $profile = $pwa->compute_conservation_profile(
        -position_type  => 'c',
        -window_size    => 100,
        -window_inc     => 1
    );

    # where
        position_type	= scoring window positions are reported as
			  nucleotide position on base_seq as either:
				c = center of scoring window
				s = start of scoring window
				e = end of scoring window
    	window_size	= width of scoring window in nucleotides on the
			  base_seq
    	window_inc	= amount by which to slide scoring window along
			  the base_seq in nucleotides for each iteration
    

    # The computed conservation profile may be retrieved again later
    # (if not already explicitly computed, automatically compute it) 
    my $profile = $pwa->conservation_profile();

    # Explicitly compute and retrieve the conserved regions.
    my $regions = $pwa->compute_conserved_regions(
        -top_pct            => '10%',
        -min_conservation   => '70%',
        -filter_exons       => 1,
        -min_cr_len         => 20,
        -pct_id_method      => 's'
    );

    # where
    	top_pct			    = report the top X percentile of conserved
                              regions (specify as a number between
                              0 and 1 or a string such as '10%')
    	min_conservation    = the absolute minimum percent identity of
                              regions to report (specify as a number
                              between 0 and 1 or a string such as
                              '70%')
        filter_exons        = (bool) exons are to be filtered out of
                              the conserved regions (if base_seq_exons
                              has been set).
        min_cr_len          = (int) min. length of a conserved region
                              to report                         
        pct_id_method       = method used to compute the percent
                              identity ('s' = standard, 'o' = overall)

    # The conserved regions can be obtained later
    # (if not already explicitly computed, automatically compute them) 
    my $regions = $pwa->conserved_regions();

    # Get the conserved sub-sequences corresponding to the conserved
    # regions.
    my $subseqs = $pwa->compute_conserved_subsequences();
    
    # The conserved sub-sequences can be obtained later
    # (compute if not already explicitly computed)
    my $subseqs = $pwa->conserved_subsequences();

    # Get the conserved sub-alignments corresponding to the conserved
    # regions.
    my $subaligns = $pwa->compute_conserved_subalignments();
    
    # The conserved sub-alignments can be obtained later
    # (compute if not already explicitly computed)
    my $subaligns = $pwa->conserved_subalignments();

    # Find TFBSs which meet the conservation criteria.
    my $sites = $pwa->compute_conserved_tfbss(
        -matrix_set                 => $matrix_set,
        -tfbs_threshold             => $tfbs_threshold,
        -min_tfbs_cr_overlap        => 1,
        -filter_overlapping_sites   => 1,
        -start                      => $start,
        -end                        => $end
    );

    # Once computed, the conserved TFBSs can be obtained later...
    my $sites = $pwa->conserved_tfbss();

=head1 DESCRIPTION

ORCA::Analysis::Pairwise is an object used for the purpose
of performing various conservation analysis tasks on a pair of aligned
orthologous DNA sequences. These tasks include finding the conserved regions,
filtering exon features out of the conserved regions, generating a conservation
profile, and extracting the conserved sub-alignments/sub-sequences.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::Analysis::Pairwise;

use strict;

use Carp;
use Bio::SeqFeature::Generic;
use ORCA::Analysis::Run::AlignCons;
use ORCA::Analysis::ConservationReport;
use TFBS::SitePairSet;

use constant DEBUG => 0;
use constant DFLT_MIN_TFBS_CA_OVERLAP => 1;

=head2 new

 Title    : new
 Usage    : $pwa = ORCA::Analysis::Pairwise->new();
 Function : Create a new ORCA::Analysis::Pairwise object.
 Returns  : A new ORCA::Analysis::Pairwise object.
 Args     : Optional named parameters:
             -base_seq          => Bio::Seq object defining the base
                                   sequence
             -comparison_seq	=> Bio::Seq object defining the
                                   comparison sequence
             -masked_base_seq   => Bio::Seq object defining repeat
                                   masked base sequence
             -masked_comparison_seq
                                => Bio::Seq object defining repeat
                                   masked comparison sequence
             -alignment         => A Bio::SimpleAlign object
                                   defining the alignment between
                                   base_seq and comparison_seq
             -base_seq_exons    => A listref of Bio::SeqFeature
                                   objects defining exon positions
                                   on base_seq
             -comparison_seq_exons
                                => A listref of Bio::SeqFeature
                                   objects defining exon positions
                                   on comparison_seq

=cut

sub new
{
    my ($class, %args) = @_;
    my $self = bless {
        -chr                   => undef,
        -start                 => undef,
        -end                   => undef,
        -base_seq              => undef,
        -comparison_seq        => undef,
        -masked_base_seq       => undef,
        -masked_comparison_seq => undef,
        -alignment             => undef,
        -base_seq_exons        => undef,
        -comparison_seq_exons  => undef,
        -conserved_regions     => undef,
        -conservation_profile  => undef,
        -conserved_subseqs     => undef,
        -conserved_subaligns   => undef,
        %args
    }, ref $class || $class;

    my $base_seq = $self->base_seq;
    if ($base_seq) {
        if (   !$base_seq->isa("Bio::SeqI")
            && !$base_seq->isa("Bio::PrimarySeqI"))
        {
            carp "Base seq is not a Bio::SeqI or Bio::PrimarySeqI compliant"
                . " object\n";
            return;
        }

        if ($base_seq->isa('Bio::LocatableSeq')) {
            if (defined $self->start) {
                if ($self->start != $base_seq->start) {
                    carp "Base sequence start does not match start argument"
                        . " provided\n";
                    return;
                }
            } else {
                $self->start($base_seq->start);
            }

            if (defined $self->end) {
                if ($self->end != $base_seq->end) {
                    carp "Base sequence end does not match end argument"
                        . " provided\n";
                    return;
                }
            } else {
                $self->end($base_seq->end);
            }
        } else {
            # using relative sequence coords
            $self->start(1) if !$self->start;
            $self->end($base_seq->length()) if !$self->end;
        }
    } else {
        carp "Must provide a base sequence\n";
    }

    my $comparison_seq = $self->comparison_seq;
    if (defined $comparison_seq) {
        if (   !$comparison_seq->isa("Bio::SeqI")
            && !$comparison_seq->isa("Bio::PrimarySeqI"))
        {
            carp "Comparison_seq is not a Bio::SeqI or Bio::PrimarySeqI"
                . " compliant object\n";
            return;
        }
    } else {
        carp "must provide a comparison sequence\n";
    }

    if (defined $self->{-masked_base_seq}) {
        if (   !$self->{-masked_base_seq}->isa("Bio::SeqI")
            && !$self->{-masked_base_seq}->isa("Bio::PrimarySeqI"))
        {
            carp "masked_base_seq is not a Bio::SeqI or Bio::PrimarySeqI"
                . " compliant object\n";
            return;
        }
    }

    if (defined $self->{-masked_comparison_seq}) {
        if (   !$self->{-masked_comparison_seq}->isa("Bio::SeqI")
            && !$self->{-masked_comparison_seq}->isa("Bio::PrimarySeqI"))
        {
            carp
                "masked_comparison_seq is not a Bio::SeqI or Bio::PrimarySeqI"
                . " compliant object\n";
            return;
        }
    }

    if (defined $self->{-alignment}) {
        if (!$self->{-alignment}->isa("Bio::SimpleAlign")) {
            carp "alignment is not a Bio::SimpleAlign compliant object\n";
            return;
        }
    }

    if (defined $self->{-base_seq_exons}) {
        if (   ref $self->{-base_seq_exons} ne "ARRAY"
            || !$self->{-base_seq_exons}->[0]->can("start")
            || !$self->{-base_seq_exons}->[0]->can("end"))
        {
            carp "base_seq_exons is not a list ref of objects with"
                . " start/end methods\n";
            return;
        }
    }

    if (defined $self->{-comparison_seq_exons}) {
        if (   ref $self->{-comparison_seq_exons} ne "ARRAY"
            || !$self->{-comparison_seq_exons}->[0]->can("start")
            || !$self->{-comparison_seq_exons}->[0]->can("end"))
        {
            carp "comparison_seq_exons is not a list ref of objects with"
                . " start/end methods\n";
            return;
        }
    }

    if ($self->{-alignment}) {
        $self->_create_pos_conversion_arrays;
    }

    return $self;
}

################################################################################
#
# Getter/Setter methods follow.
#
################################################################################

=head2 chr

 Title    : chr
 Usage    : $chr = $pwa->chr($chr);
 Function : Get the chromosom name.
 Returns  : The chromosom name.
 Args     : A new chromosome name.

=cut

sub chr
{
    my ($self, $chr) = @_;

    if ($chr) {
        $self->{-chr} = $chr;
    }

    return $self->{-chr};
}

=head2 start

 Title    : start
 Usage    : $start = $pwa->start();
 Function : Get/set the chromosomal start coordinate of the base sequence.
 Returns  : The chromosomal start coordinate of the base sequence.
 Args     : A new chromosomal start coordinate of the base sequence.

=cut

sub start
{
    my ($self, $start) = @_;

    if (defined($start)) {
        $self->{-start} = $start;
    }

    return $self->{-start};
}

=head2 end

 Title    : end
 Usage    : $end = $pwa->end();
 Function : Get/set the chromosomal end coordinate of the base sequence.
 Returns  : The chromosomal end coordinate of the base sequence.
 Args     : A new chromosomal end coordinate of the base sequence.

=cut

sub end
{
    my ($self, $end) = @_;

    if (defined($end)) {
        $self->{-end} = $end;
    }

    return $self->{-end};
}

=head2 base_seq

 Title    : base_seq
 Usage    : $base_seq = $pwa->base_seq($seq);
 Function : Get/Set the base sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new base sequence (Bio::Seq object)

=cut

sub base_seq
{
    my ($self, $base_seq) = @_;

    if (defined($base_seq)) {
        if (   !$base_seq->isa("Bio::SeqI")
            && !$base_seq->isa("Bio::PrimarySeqI"))
        {
            carp "base_seq is not a Bio::SeqI or Bio::PrimarySeqI compliant"
                . " object\n";
            return undef;
        }
        $self->{-base_seq} = $base_seq;
    }

    return $self->{-base_seq};
}

=head2 comparison_seq

 Title    : comparison_seq
 Usage    : $comparison_seq = $pwa->comparison_seq($seq);
 Function : Get/Set the comparison sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new comparison sequence (Bio::Seq object)

=cut

sub comparison_seq
{
    my ($self, $comparison_seq) = @_;

    if (defined($comparison_seq)) {
        if (   !$comparison_seq->isa("Bio::SeqI")
            && !$comparison_seq->isa("Bio::PrimarySeqI"))
        {
            carp "comparison_seq is not a Bio::SeqI or Bio::PrimarySeqI"
                . " compliant object\n";
            return undef;
        }
        $self->{-comparison_seq} = $comparison_seq;
    }

    return $self->{-comparison_seq};
}

=head2 masked_base_seq

 Title    : masked_base_seq
 Usage    : $masked_base_seq = $pwa->masked_base_seq($seq);
 Function : Get/Set the masked base sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new masked base sequence (Bio::Seq object)

=cut

sub masked_base_seq
{
    my ($self, $masked_base_seq) = @_;

    if (defined($masked_base_seq)) {
        if (   !$masked_base_seq->isa("Bio::SeqI")
            && !$masked_base_seq->isa("Bio::PrimarySeqI"))
        {
            carp "masked_base_seq is not a Bio::SeqI or Bio::PrimarySeqI"
                . " compliant object\n";
            return undef;
        }
        $self->{-masked_base_seq} = $masked_base_seq;
    }

    return $self->{-masked_base_seq};
}

=head2 masked_comparison_seq

 Title    : masked_comparison_seq
 Usage    : $masked_comparison_seq = $pwa->masked_comparison_seq($seq);
 Function : Get/Set the masked comparison sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new masked comparison sequence (Bio::Seq object)

=cut

sub masked_comparison_seq
{
    my ($self, $masked_comparison_seq) = @_;

    if (defined($masked_comparison_seq)) {
        if (   !$masked_comparison_seq->isa("Bio::SeqI")
            && !$masked_comparison_seq->isa("Bio::PrimarySeqI"))
        {
            carp
                "masked_comparison_seq is not a Bio::SeqI or Bio::PrimarySeqI"
                . " compliant object\n";
            return undef;
        }
        $self->{-masked_comparison_seq} = $masked_comparison_seq;
    }

    return $self->{-masked_comparison_seq};
}

=head2 alignment

 Title    : alignment
 Usage    : $align = $pwa->alignment($align);
 Function : Get/Set the alignment.
 Returns  : A Bio::SimpleAlign object or undef.
 Args     : Optionally a new alignment (Bio::SimpleAlign) object.

=cut

sub alignment
{
    my ($self, $align) = @_;

    if (defined($align)) {
        if (!$align->isa("Bio::SimpleAlign")) {
            carp "alignment is not a Bio::SimpleAlign compliant object\n";
            return undef;
        }

        $self->{-alignment} = $align;
        $self->_destroy_position_conversion_arrays;
    } elsif (!$self->{-alignment}) {
        $self->compute_alignment();
    }

    return $self->{-alignment};
}

=head2 base_seq_exons

 Title    : base_seq_exons
 Usage    : $exons = $pwa->base_seq_exons();
 Function : Get/Set the base sequence exons.
 Returns  : A reference to a list of Bio::SeqFeatureI objects.
 Args     : Optionally new base sequence exons (reference to a
            list of Bio::SeqFeatureI objects)

=cut

sub base_seq_exons
{
    my ($self, $exons) = @_;

    if (defined $exons) {
        if (    ref $exons ne "ARRAY"
            || !$exons->[0]->can("start")
            || !$exons->[0]->can("end"))
        {
            carp "Base sequence exons is not a list ref of objects with"
                . " start/end methods\n";
            return undef;
        }
        $self->{-base_seq_exons} = $exons;
    }

    return $self->{-base_seq_exons};
}

sub comparison_seq_exons
{
    my ($self, $exons) = @_;

    if (defined $exons) {
        if (   ref $exons ne "ARRAY"
            || !$exons->[0]->can("start")
            || !$exons->[0]->can("end"))
        {
            carp "Comparison sequence exons is not a list ref of objects"
                . " with start/end methods\n";
            return undef;
        }
        $self->{-comparison_seq_exons} = $exons;
    }

    return $self->{-comparison_seq_exons};
}

=head2 alignment_orientation

 Title    : alignment_orientation
 Usage    : $align_orient = $pwa->alignment_orientation();
 Function : Get the alignment orientation.
 Returns  : +1 or -1 depending on whether the comparison sequence is
	    reverse complemented in the alignment.
 Args     : None.

=cut

sub alignment_orientation
{
    my $self = shift;

    return $self->{_alignment_orientation};
}

=head2 position_type

 Title    : position_type
 Usage    : $position_type = $pwa->position_type($position_type);
 Function : Get/set the position type. Scoring window positions are
            reported as nucleotide position on base_seq as either:
				c = center of scoring window
				s = start of scoring window
				e = end of scoring window
 Returns  : The position type.
 Args     : Optional new position type.

=cut

sub position_type
{
    my ($self, $position_type) = @_;

    if (defined $position_type) {
        unless (   $position_type eq 'c'
                || $position_type eq 's'
                || $position_type eq 'e'
        ) {
            carp "Invalid position type $position_type"
                . " - should be 'c', 's' or 'e'\n";
        } else {
            $self->{-position_type} = $position_type;
        }
    }

    return $self->{-position_type};
}

=head2 window_size

 Title    : window_size
 Usage    : $window_size = $pwa->window_size($window_size);
 Function : Get/set the size of the alignment scoring window. Generally will
            be somewhere in the range of 50 to 200.
 Returns  : The size of the alignment scoring window.
 Args     : Optional new alignment scoring window size.

=cut

sub window_size
{
    my ($self, $window_size) = @_;

    if (defined $window_size) {
        $self->{-window_size} = $window_size;
    }

    return $self->{-window_size};
}

=head2 window_inc

 Title    : window_inc
 Usage    : $window_inc = $pwa->window_inc($window_inc);
 Function : Get/set the alignment scoring window increment amount.
 Returns  : The size of the alignment scoring window increment.
 Args     : Optional new alignment scoring window  increment.

=cut

sub window_inc
{
    my ($self, $window_inc) = @_;

    if (defined $window_inc) {
        $self->{-window_inc} = $window_inc;
    }

    return $self->{-window_inc};
}

=head2 pct_id_method

 Title    : pct_id_method
 Usage    : $pct_id_method = $pwa->pct_id_method($pct_id_method);
 Function : Get/set the method of computing the percent identity of
            conserved regions ('s' = standard, 'o' = overall).
 Returns  : The percent ID method.
 Args     : A percent ID method; either 's' or 'o'.

=cut

sub pct_id_method
{
    my ($self, $method) = @_;

    if (defined $method) {
        unless ($method eq 's' || $method eq 'o') {
            carp "Incorrect percent ID method specified - '$method'.\n" 
                . " Should be specified as either 's' (standard)"
                . " or 'o' (overall)\n";
            return;
        }

        $self->{-pct_id_method} = $method;
    }

    return $self->{-pct_id_method};
}

=head2 min_conservation

 Title    : min_conservation
 Usage    : $min_conservation = $pwa->min_conservation($min_conservation);
 Function : Get/set the minimum conservation.
 Returns  : 
 Args     : A minimum conservation specified as a float between 0 and 1 or
            as a percentage string, e.g. '70%'.

=cut

sub min_conservation
{
    my ($self, $min_conservation) = @_;

    if (defined $min_conservation) {
        if ($min_conservation =~ /(.+)%/) {
            $min_conservation = $1 / 100;
        }

        if ($min_conservation < 0 || $min_conservation > 1) {
            carp "Min. conservation is out of range; please specify as a number"
                . " between 0 and 1 or as a string between '0%' and '100%'";
            return;
        }

        $self->{-min_conservation} = $min_conservation;
    }

    return $self->{-min_conservation};
}

=head2 top_pct

 Title    : top_pct
 Usage    : $top_pct = $pwa->top_pct($top_pct);
 Function : Get/set the top percentile of scoring windows that is used to
            compute the conservation cutoff dynamically.
 Returns  : The top percentile.
 Args     : A top percentile specified as a float between 0 and 1 or
            as a percentage string, e.g. '70%'.

=cut

sub top_pct
{
    my ($self, $top_pct) = @_;

    if (defined $top_pct) {
        if ($top_pct =~ /(.+)%/) {
            $top_pct = $1 / 100;
        }

        if ($top_pct < 0 || $top_pct > 1) {
            carp "Top percentile is out of range; please specify as a number"
                . " between 0 and 1 or as a string between '0%' and '100%'";
            return;
        }

        $self->{-top_pct} = $top_pct;
    }

    return $self->{-top_pct};
}

=head2 conservation_cutoff

 Title    : conservation_cutoff
 Usage    : $conservation_cutoff = $pwa->conservation_cutoff();
 Function : Get the conservation cutoff. This is the actual conservation
            cutoff value computed from the min. conservation / top
            percentile combination.
 Returns  : The conservation cutoff.
 Args     : None

=cut

sub conservation_cutoff
{
    my ($self) = @_;

    my $crr = $self->conserved_regions_report();

    if ($crr) {
        # Conserved regions report cutoff is in rainge 0-100, convert to 0-1
        return $crr->param('cutoff') / 100;
    }

    return undef;
}

=head2 filter_exons

 Title    : filter_exons
 Usage    : $filter_exons = $pwa->filter_exons($filter_exons);
 Function : Get/set whether to filter exons from the conserved regions
 Returns  : Boolean indicating whether exons are filtered from the
            conserved regions.
 Args     : Optional boolean value indicating whether exons should be
            filtered from conserved regions.

=cut

sub filter_exons
{
    my ($self, $filter_exons) = @_;

    if (defined $filter_exons) {
        $self->{-filter_exons} = $filter_exons;
    }

    return $self->{-filter_exons};
}

=head2 min_cr_len

 Title    : min_cr_len
 Usage    : $min_cr_len = $pwa->min_cr_len($min_cr_len);
 Function : Get/set the minimum length of conserved regions.
 Returns  : The minimum length of conserved regions.
 Args     : Optional new minimum conserved region length.

=cut

sub min_cr_len
{
    my ($self, $min_cr_len) = @_;

    if (defined $min_cr_len) {
        $self->{-min_cr_len} = $min_cr_len;
    }

    return $self->{-min_cr_len};
}

=head2 conservation_profile

 Title    : conservation_profile
 Usage    : Synonym for conservation_profile_array() method

=cut

sub conservation_profile
{
    my $self = shift;

    return $self->conservation_profile_array;
}

=head2 conservation_profile_array

 Title    : conservation_profile_array
 Usage    : $profile = $pwa->conservation_profile_array();
 Function : Get the conservation profile.
 Returns  : An array of phastCons scores where the value at each array
            index corresponds to the score at the location in the range of
            start..end
 Args     : None.

=cut

sub conservation_profile_array
{
    my $self = shift;

    if (!$self->{_profile_array}) {
        $self->compute_conservation_profile();
    }

    return $self->{_profile_array};
}


=head2 conservation_profile_hash

 Title    : conservation_profile_hash
 Usage    : $profile = $pwa->conservation_profile_hash();
 Function : Get the conservation profile as a hash
 Returns  : A hashref of conservation scores where the hash key is the
            position (in the range start..end) and the hash value is the
            conservation score at that position.
 Args     : None.

=cut

sub conservation_profile_hash
{
    my $self = shift;

    if (!$self->{_profile_hash}) {
        $self->compute_conservation_profile();
    }

    return $self->{_profile_hash};
}

=head2 conserved_regions

 Title    : conserved_regions
 Usage    : $regions = $pwa->conserved_regions()
 	    OR $pwa->conserved_regions($regions);
 Function : Get/set the conserved regions. If getting conserved regions
 	    and they have not already been computed, compute them.
 Returns  : A reference to a list of Bio::SeqFeature::FeaturePair objects
	    or undef.
 Args     : Optionally, a listref of Bio::SeqFeature::FeaturePair objects.

=cut

sub conserved_regions
{
    my ($self, $regions) = @_;

    if ($regions) {
        if (ref $regions ne "ARRAY"
            || !$regions->[0]->isa("Bio::SeqFeature::FeaturePair"))
        {
            carp "regions argument is not a listref of"
                . " Bio::SeqFeature::FeaturePair objects\n";
            return undef;
        }
        $self->{-conserved_regions} = $regions;
    } elsif (!$self->{-conserved_regions}) {
        $self->compute_conserved_regions();
    }

    return $self->{-conserved_regions};
}

=head2 conservation_profile_report

 Title    : conservation_profile_report
 Usage    : $cp_report = $pwa->conservation_profile_report();
 Function : Get the conservation profile report.
 Returns  : A ORCA::Analysis::ConservationReport object or
	    undef.
 Args     : None.

=cut

sub conservation_profile_report
{
    $_[0]->{-conservation_profile_report};
}

=head2 conserved_regions_report

 Title    : conserved_regions_report
 Usage    : $cr_report = $pwa->conserved_regions_report();
 	    OR $pwa->conserved_regions_report($report);
 Function : Get/set the conserved regions report.
 Returns  : A ORCA::Analysis::ConservationReport object or
	    undef.
 Args     : Optionally a ORCA::Analysis::ConservationReport.

=cut

sub conserved_regions_report
{
    my ($self, $report) = @_;

    if ($report) {
        if (!$report->isa("ORCA::Analysis::ConservationReport")) {
            carp "report argument is not an "
                . " ORCA::Analysis::ConservationReport"
                . " object\n";
            return undef;
        }
        $self->{-conserved_regions_report} = $report;
    }

    return $self->{-conserved_regions_report};
}

=head2 conserved_subsequences

 Title    : conserved_subsequences
 Usage    : $subseqs = $pwa->conserved_subsequences();
 Function : Get the list of conserved sub-sequences.
 Returns  : A reference to a list of Bio::Seq objects.
 Args     : None.

=cut

sub conserved_subsequences
{
    my $self = shift;

    if (!$self->{-conserved_subseqs}) {
        $self->compute_conserved_subsequences();
    }

    return $self->{-conserved_subseqs};
}

=head2 conserved_subalignments

 Title    : conserved_subalignments
 Usage    : $subaligns = $pwa->conserved_subalignments();
 Function : Get the list of conserved sub-alignments.
 Returns  : A reference to a list of Bio::SimpleAlign objects.
 Args     : None.

=cut

sub conserved_subalignments
{
    my $self = shift;

    if (!$self->{-conserved_subaligns}) {
        $self->compute_conserved_subalignments();
    }

    return $self->{-conserved_subaligns};
}

################################################################################
#
# Analysis methods follow.
#
################################################################################

=head2 compute_alignment

 Title    : compute_alignment
 Usage    : $alignment = $pwa->compute_alignment();
 Function : Compute and return the pairwise alignment of the base and
 	    comparison sequences.
 Returns  : A Bio::SimpleAlign object.
 Args     : None.

=cut

sub compute_alignment
{
    my ($self) = @_;

    require ORCA::Aligner;

    my $aligner = ORCA::Aligner->new();
    if (!$aligner) {
        carp "could not initialize ORCA::Aligner\n" if !$aligner;
        return undef;
    }

    my $aln = $aligner->align(
        -seq1 => $self->masked_base_seq       || $self->base_seq,
        -seq2 => $self->masked_comparison_seq || $self->comparison_seq
    );
    if (!$aln) {
        carp "could not align sequences with ORCA::Aligner\n" if !$aln;
        return undef;
    }

    $self->{-alignment} = $aln;

    $self->_create_pos_conversion_arrays;

    return $self->{-alignment};
}

=head2 compute_conservation_profile

 Title    : compute_conservation_profile
 Usage    : $profile = $pwa->compute_conservation_profile();
 Function : Compute and return the conservation profile.
 Returns  : A listref of hashes containing '-position' and '-score'
            key/value pairs.
 Args     : Optionally
            -window_size    => Size of the conservation scoring window 
            -window_inc     => Number of positions to slide the conservation
                               scoring window 
            -position_type	=> ('c', 's', or 'e') indicating positions
				   should be reported as the center, start
				   or end position of the scoring window.

=cut

sub compute_conservation_profile
{
    my ($self, %args) = @_;

    my $start = $self->start;
    if (!$start) {
        carp "Error: chromosome start coordinate not specified\n";
        return;
    }

    my $end = $self->end;
    if (!$end) {
        carp "Error: chromosome end coordinate not specified\n";
        return;
    }

    my $align = $self->alignment;
    if (!$align) {
        $align = $self->compute_alignment;

        if (!$align) {
            carp "error computing alignment\n";
            return undef;
        }
    }

    my $win_size = $self->window_size($args{-window_size});
    my $win_inc  = $self->window_inc($args{-window_inc});
    my $pos_type = $self->position_type($args{-position_type});

    my $align_cons = ORCA::Analysis::Run::AlignCons->new();
    if (!$align_cons) {
        carp "error initializing ORCA::Analysis::Run::AlignCons\n";
        return undef;
    }

    my %ac_run_params = (-alignment => $align, -r => 'p');
    $ac_run_params{-w} = $win_size if $win_size;
    $ac_run_params{-n} = $win_inc  if $win_inc;
    $ac_run_params{-f} = $pos_type if $pos_type;

    my $report = $align_cons->run(%ac_run_params);

    $self->{-conservation_profile_report} = $report;

    #
    # Conservation profile returned by report is an array of hashes
    # with -position and -score tags.
    #
    my $report_cp = $report->conservation_profile;

    #
    # Store conservation scores as position/value pairs. Positions will
    # be stored in the start..end range, i.e. chromosomal coordinates.
    # May be non-continuous data. DJA 11/10/21
    #
    my %profile_hash;
    foreach my $cp (@$report_cp) {
        $profile_hash{$cp->{-position} + $start - 1} = $cp->{-score};
    }
    $self->{_profile_hash} = \%profile_hash;

    #
    # Generate a profile as continuous array data (0-based coords). Fill
    # gaps with 0 values. DJA 11/10/21
    #
    my @profile_array;
    foreach my $idx (0 .. $end - $start) {
        my $pos = $idx + $start;
        $profile_array[$idx] = $profile_hash{$pos} || 0;
    }

    $self->{_profile_array} = @profile_array ? \@profile_array : undef;
}

=head2 compute_conserved_regions

 Title    : compute_conserved_regions
 Usage    : $regions = $pwa->compute_conserved_regions();
 Function : Compute and return the list of conserved regions.
 Returns  : A reference to a list of Bio::SeqFeature::Generic objects.
 Args     : Optionally:
            -top_pct        => report the top X percentile of conserved
                               regions (specify as a number between
                               0 and 1 or a string such as '10%')
            -min_conservation
                            => the absolute minimum percent identity of
                               regions to report (specify as a number
                               between 0 and 1 or a string such as
                               '70%')
            -filter_exons   => indicates exons should be
                               filtered out of conserved
                               regions.
            -min_cr_len     => minimum length of conserved
                               region to keep when exons are
                               filtered.
            -pct_id_method  => method used to compute the percent
                               identity ('s' = standard,
                               'o' = overall)

=cut

sub compute_conserved_regions
{
    my ($self, %args) = @_;

    my $aln = $self->alignment;
    if (!$aln) {
        $aln = $self->compute_alignment;

        if (!$aln) {
            carp "error computing alignment\n";
            return undef;
        }
    }

    my $top_pct = $self->top_pct($args{-top_pct});
    my $min_conservation = $self->min_conservation($args{-min_conservation});

    unless (defined $top_pct || defined $min_conservation) {
        carp "No top percent or minimum conservation defined\n";
        return;
    }

    my $win_size = $self->window_size;
    my $win_inc  = $self->window_inc;
    my $pos_type = $self->position_type;
    my $filter_exons = $self->filter_exons($args{-filter_exons});

    #
    # Now set minimum conserved region length regardless of whether
    # filter_exons is set. DJA 2008/01/10
    #
    my $min_cr_len = $self->min_cr_len($args{-min_cr_len});

    my $pct_id_method = $self->pct_id_method($args{-pct_id_method});

    my $align_cons = ORCA::Analysis::Run::AlignCons->new();
    if (!$align_cons) {
        carp "error initializing ORCA::Analysis::Run::AlignCons\n";
        return undef;
    }

    my %ac_run_params = (-alignment => $aln);

    if ($filter_exons && $self->base_seq_exons) {
        my $start = $self->start;
        if ($start > 1) {
            #
            # Assume if start > 1 that we are using a chromosomal coordinate
            # system and that exons are passed in with chromosomal coordinates.
            # If so we have to convert the coords to 1-based for align_cons.
            #
            my @exons;
            foreach my $exon (@{$self->base_seq_exons()}) {
                push @exons, Bio::SeqFeature::Generic->new(
                    -primary_tag    => $exon->primary_tag,
                    -source_tag     => $exon->source_tag,
                    -display_name   => $exon->display_name,
                    -strand         => $exon->strand,
                    -start          => $exon->start - $start + 1,
                    -end            => $exon->end   - $start + 1
                );
            }
            $ac_run_params{-features} = \@exons;
        } else {
            $ac_run_params{-features} = $self->base_seq_exons;
        }
    }

    $ac_run_params{-r} = $pos_type         if $pos_type;
    $ac_run_params{-w} = $win_size         if $win_size;
    $ac_run_params{-n} = $win_inc          if $win_inc;
    $ac_run_params{-s} = $top_pct          if $top_pct;
    $ac_run_params{-t} = $min_conservation if $min_conservation;
    $ac_run_params{-l} = $min_cr_len       if $min_cr_len;
    $ac_run_params{-m} = $pct_id_method    if $pct_id_method;

    #
    # Now set the aligment reverse complemented flag as an ac_params value
    # rather than letting the ORCA::Analysis::Run::AlignCons module
    # computed based on the raw alignment. DJA 2008/01/10
    #
    if ($self->alignment_orientation == -1) {
        $ac_run_params{-c} = 1;
    }

    my $report = $align_cons->run(%ac_run_params);

    if ($report) {
        $self->conserved_regions_report($report);

        # check if undef before assigning to prevent infinite recursion
        my $crfp = $report->conserved_regions_as_feature_pairs;

        if ($crfp) {
            $self->conserved_regions($crfp);

            return $crfp;
        }
    }

    return undef;
}

=head2 compute_conserved_subsequences

 Title    : compute_conserved_subsequences
 Usage    : $subseqs = $pwa->compute_conserved_subsequences($seq_num);
 Function : Return the sub-sequences corresponding to the conserved regions.
 Returns  : A reference to a list of Bio::Seq objects.
 Args     : Optionally specify a sequence number; default = 1.

=cut

sub compute_conserved_subsequences
{
    my ($self, $seq_num) = @_;

    $seq_num = 1 if !$seq_num;

    my $seq;
    if ($seq_num == 1) {
        $seq = $self->base_seq;
    } elsif ($seq_num == 2) {
        $seq = $self->comparison_seq;
    } else {
        carp "invalid sequence number $seq_num";
        return;
    }
    return if !defined $seq;

    my $seq_display_id = $seq->display_id;
    #
    # Not used and will not work for Bio::Seq objects (only Bio::LocatableSeq)
    # DJA 090820
    #
    #my $seq_start = $seq->start;
    #my $seq_end = $seq->end;

    my $crs = $self->conserved_regions;
    if (!defined $crs) {
        carp "no conserved regions defining sub-sequences to extract\n";
        return undef;
    }

    my @subseqs;
    foreach my $cr (@$crs) {
        my $cr_start;
        my $cr_end;
        if ($seq_num == 1) {
            $cr_start = $cr->feature1->start;
            $cr_end   = $cr->feature1->end;
        } elsif ($seq_num == 2) {
            $cr_start = $cr->feature2->start;
            $cr_end   = $cr->feature2->end;
        }

        my $subseq = $seq->trunc($cr_start, $cr_end);

        $seq_display_id = "SUBSEQ" if !$seq_display_id;
        my $subseq_display_id = $seq_display_id . '/' . "$cr_start-$cr_end";
        $subseq->display_id($subseq_display_id);
        push @subseqs, $subseq;
    }

    $self->{-conserved_subseqs} = @subseqs ? \@subseqs : undef;
}

=head2 compute_conserved_subalignments

 Title    : compute_conserved_subalignments
 Usage    : $subaligns = $pwa->compute_conserved_subalignments();
 Function : Return the sub-alignments corresponding to the conserved
 	    regions.
 Returns  : A reference to a list of Bio::SimpleAlign objects.
 Args     : None.

=cut

sub compute_conserved_subalignments
{
    my ($self) = @_;

    my $align = $self->alignment;
    return undef if !defined $align;

    my $regions = $self->conserved_regions;
    if (!defined $regions) {
        carp "no conserved regions defining sub-alignments to extract\n";
        return undef;
    }

    my @subaligns;
    foreach my $reg (@$regions) {
        #
        # Get alignment coordinates directly DJA 051021
        #
        my $align_start = ($reg->feature1->get_tag_values('align_start'))[0];
        my $align_end   = ($reg->feature1->get_tag_values('align_end'))[0];

        #
        # Just in case
        #
        if (!$align_start) {
            $align_start = $self->seq_to_align_pos(1, $reg->feature1->start);
        }
        if (!$align_end) {
            $align_end = $self->seq_to_align_pos(1, $reg->feature1->end);
        }

        my $subalign;
        #
        # Just in case start and end coords are invalid
        #
        eval {
            $subalign = $align->slice($align_start, $align_end);

            my $subseq1;
            my $subseq2;
            eval {
                #
                # Bioperl complains if one of the sequences contains only
                # gap chars and does not include it in the alignment which
                # will cause an exception to be thrown by the get_seq_by_pos
                # method. This sub-alignment is not useful, so just don't keep
                # it at all.
                # DJA 051021
                #
                $subseq1 = $subalign->get_seq_by_pos(1);
                $subseq2 = $subalign->get_seq_by_pos(2);
            };

            if (!$subseq1 || !$subseq2) {
                undef $subalign;
            } else {
                # because the slice method doesn't seem to preserve this
                # information for some reason
                $subseq1->alphabet('dna');
                $subseq2->alphabet('dna');
            }
        };

        if ($subalign) {
            push @subaligns, $subalign;
        }
    }

    $self->{-conserved_subaligns} = @subaligns ? \@subaligns : undef;
}

################################################################################
#
# Utility methods follow.
#
################################################################################

=head2 seq_to_align_pos

 Title    : seq_to_align_pos
 Usage    : $align_pos = $pwa->seq_to_align_pos(
				    $which_seq, $seq_pos, $coord_type);
 Function : Return the position within the alignment corresponding to the
 	    position within the specified sequence.
 Returns  : The alignment position (int) on success, otherwise 0.
 Args     : $which_seq  - either a 1 or 2 or 'base' or 'comparison'
 			  indicating which sequence the sequence positon
			  refers to
	    $seq_pos    - the nucleotide position within the sequence
	    $coord_type - either true (non-zero) indicating seq_pos is
	    		  given in absolute coordinates, or 0 indicating
			  relative coordinates

=cut

sub seq_to_align_pos
{
    my ($self, $which_seq, $seq_pos, $coord_type) = @_;

    return 0 if !$seq_pos || $seq_pos < 1;

    #
    # New way - faster???
    #
    if (!$self->{_pos_conversion_arrays_created}) {
        $self->_create_pos_conversion_arrays;
    }

    if ($coord_type) {
        # convert back to relative position
        $seq_pos = $seq_pos 
            - $self->alignment->get_seq_by_pos($which_seq)->start + 1;
    }

    my $align_pos = 0;
    if ($which_seq == 1) {
        if ($seq_pos > 0 && $seq_pos < scalar @{$self->{_base_to_align_pos}})
        {
            $align_pos = $self->{_base_to_align_pos}->[$seq_pos];
        }
    } elsif ($which_seq == 2) {
        if (   $seq_pos > 0
            && $seq_pos < scalar @{$self->{_comparison_to_align_pos}})
        {
            $align_pos = $self->{_comparison_to_align_pos}->[$seq_pos];
        }
    }

    return $align_pos;

    #
    #     Old way
    #
    #    my $seq_nucs;
    #    if ($which_seq eq '1' || lc $which_seq eq 'base') {
    #	$seq_nucs = $self->_aligned_base_nucs;
    #    } elsif ($which_seq eq '2' || lc $which_seq eq 'comparison') {
    #	$seq_nucs = $self->_aligned_comparison_nucs;
    #    } else {
    #	return 0;
    #    }
    #
    #    return 0 if !$seq_nucs;
    #
    #    my $seq_idx = 0;
    #    my $nuc_count = 0;
    #    my $aln_pos = 0;
    #    while ($seq_idx < scalar @$seq_nucs) {
    #	if ($seq_nucs->[$seq_idx] ne '-') {
    #	    $nuc_count++;
    #	    if ($nuc_count == $seq_pos) {
    #		$aln_pos = $seq_idx + 1;
    #		last;
    #	    }
    #	}
    #	$seq_idx++;
    #    }
    #
    #    return $aln_pos;
}

=head2 align_to_seq_pos

 Title    : align_to_seq_pos
 Usage    : $seq_pos = $pwa->align_to_seq_pos($which_seq, $align_pos,
						$match_type, $coord_type);
 Function : Return the position within the specified sequence corresponding
 	    to the position within the alignment.
 Returns  : The sequence position (int) on success, otherwise 0. See
 	    the description of the match_type parameter below.
 Args     : $which_seq  - either a 1 or 2 or 'base' or 'comparison'
 			  indicating which sequence the sequence positon
			  referse to
	    $align_pos  - the position within the alignment
	    $match_type - either 'eq', 'le' or 'ge' indicating the exactness
	    		  of the position match. If the alignment position
			  corresponds to a gap then a match_type of 'eq'
			  results in a return value of 0, a match_type of
			  'le' returns the position of the nearest
			  nucleotide less than or equal to the gap position
			  and a match_type of 'ge' returns the position
			  of the nearest nucleotide greater than or equal
			  to the gap position.
	    $coord_type - either true (non-zero) indicating absolute
			  coordinates, or 0 indicating relative coordinates

=cut

sub align_to_seq_pos
{
    my ($self, $which_seq, $aln_pos, $match_type, $coord_type) = @_;

    return 0 if !$aln_pos || $aln_pos < 1;

    $match_type = 'eq' if !defined $match_type;
    $match_type = lc $match_type;
    $match_type = 'le'
        if $match_type eq 'lt'
            || $match_type eq '<'
            || $match_type eq '<=';
    $match_type = 'ge'
        if $match_type eq 'gt'
            || $match_type eq '>'
            || $match_type eq '>=';

    $coord_type = 0 if !defined $coord_type;

    #
    # New way - faster???
    #
    if (!$self->{_pos_conversion_arrays_created}) {
        $self->_create_pos_conversion_arrays;
    }

    my $seq_pos = 0;
    if ($which_seq == 1) {
        if ($aln_pos < scalar @{$self->{_align_to_base_pos}}) {
            if ($match_type eq 'eq') {
                $seq_pos = $self->{_align_to_base_pos}->[$aln_pos];
            } elsif ($match_type eq 'le') {
                while (!$seq_pos && $aln_pos > 0) {
                    $seq_pos = $self->{_align_to_base_pos}->[$aln_pos];
                    $aln_pos--;
                }
            } elsif ($match_type eq 'ge') {
                while (!$seq_pos
                    && $aln_pos < scalar @{$self->{_align_to_base_pos}})
                {
                    $seq_pos = $self->{_align_to_base_pos}->[$aln_pos];
                    $aln_pos++;
                }
            }
        }
    } elsif ($which_seq == 2) {
        if ($aln_pos < scalar @{$self->{_align_to_comparison_pos}}) {
            if ($match_type eq 'eq') {
                $seq_pos = $self->{_align_to_comparison_pos}->[$aln_pos];
            } elsif ($match_type eq 'le') {
                while (!$seq_pos && $aln_pos > 0) {
                    $seq_pos = $self->{_align_to_comparison_pos}->[$aln_pos];
                    $aln_pos--;
                }
            } elsif ($match_type eq 'ge') {
                while (!$seq_pos
                    && $aln_pos < scalar @{$self->{_align_to_comparison_pos}})
                {
                    $seq_pos = $self->{_align_to_comparison_pos}->[$aln_pos];
                    $aln_pos++;
                }
            }
        }
    }

    if ($coord_type && $seq_pos) {
        # convert to absolute coordinates
        $seq_pos = $seq_pos 
            + $self->alignment->get_seq_by_pos($which_seq)->start - 1;
    }

    return $seq_pos;

   #
   #    Old way
   #
   #    my $seq_nucs;
   #    if ($which_seq eq '1' || lc $which_seq eq 'base') {
   #	$seq_nucs = $self->_aligned_base_nucs;
   #    } elsif ($which_seq eq '2' || lc $which_seq eq 'comparison') {
   #	$seq_nucs = $self->_aligned_comparison_nucs;
   #    } else {
   #	return 0;
   #    }
   #
   #    return 0 if !$seq_nucs;
   #
   #    return 0 if ($aln_pos > scalar @$seq_nucs);
   #
   #    return 0 if ($match_type eq 'eq' && $seq_nucs->[$aln_pos - 1] eq '-');
   #
   #    my $nuc_count = 0;
   #    my $aln_idx = 0;
   #    if ($match_type eq 'le' || $match_type eq 'eq') {
   #	while ($aln_idx < $aln_pos) {
   #	    if ($seq_nucs->[$aln_idx] ne '-') {
   #		$nuc_count++;
   #	    }
   #	    $aln_idx++;
   #	}
   #    } elsif ($match_type eq 'ge') {
   #	while ($aln_idx < scalar @$seq_nucs) {
   #	    if ($seq_nucs->[$aln_idx] ne '-') {
   #		$nuc_count++;
   #		last if $aln_idx >= $aln_pos - 1;
   #	    }
   #	    $aln_idx++;
   #	}
   #	$nuc_count = 0 if $aln_idx >= scalar @$seq_nucs;
   #    }
   #
   #    return $nuc_count;
}

=head2 convert_seq_pos

 Title    : convert_seq_pos
 Usage    : $to_pos = $pwa->convert_seq_pos($from, $to, $from_pos,
					    $match_type, $coord_type);
 Function : Convert nucleotide position on one sequence to the aligned
 	    sequence nucleotide position on the other sequence.
 Returns  : The aligned sequence nucleotide position (int) on success,
            otherwise 0. See discussion of match_type parameter below.
 Args     : $from       - either a 1 or 2 or 'base' or 'comparison'
 			  indicating from which sequence to convert the
			  nucleotide positon
	    $to         - either a 1 or 2 or 'base' or 'comparison'
	    		  indicating to which sequence to convert the
			  nucleotide positon
	    $from_pos   - the nucleotide position on the 'from' sequence
	    $match_type - either 'eq', 'le' or 'ge' indicating the exactness
	    		  of the position match. If the 'from' sequence
			  position corresponds to a gap in the aligment on
			  the 'to' sequence, then a match_type of 'eq'
			  results in a return value of 0, a match_type of
			  'le' returns the position of the nearest
			  nucleotide less than or equal to the gap position
			  and a match_type of 'ge' returns the position
			  of the nearest nucleotide greater than or equal
			  to the gap position.
	    $coord_type - either true (non-zero) indicating from_pos and
	    		  returned position are in absolute coordinates
			  or 0 indicating they are in relative coordinates

=cut

sub convert_seq_pos
{
    my ($self, $from, $to, $from_pos, $match_type, $coord_type) = @_;

    return $self->align_to_seq_pos($to,
        $self->seq_to_align_pos($from, $from_pos, $coord_type),
        $match_type, $coord_type);
}

=head2 exon1_offset

 Title    : exon1_offset
 Usage    : $offset = $pwa->exon1_offset();
 Function : Return a the difference in alignment coordinates between the 3'
 	    end of the first exon of the comparison sequence and the 3' end
	    of the first exons of the base sequence, or undef if this could
	    not be computed.
	    This can be used to check the quality of the alignment in terms
	    of mis-annotation of first exons or alternate promoters
	    etc.
 Returns  : Difference in position of 3' end of first exons of two aligned
 	    sequences in alignment coordinates. Value is computed as
	    comparison sequence exon1 3' position minus base_sequence exon1
	    3' position and converted according to strand of first sequence.
	    Thus the return value is signed. A positive value indicates that
	    comparison sequence first exon is more downstream. Return undef
	    on error.
 Args     : None.

=cut

sub exon1_offset
{
    my $self = shift;

    my $align   = $self->alignment;
    my $bexons  = $self->base_seq_exons;
    my $bstrand = $self->base_seq->strand;
    my $cexons  = $self->comparison_seq_exons;
    my $cstrand = $self->comparison_seq->strand;

    return undef if !$align || !$bexons || !$cexons;

    my $bexon1;
    my $bexon1_pos;
    if ($bstrand == 1) {
        $bexon1 = $bexons->[0];
        return undef if !$bexon1;
        $bexon1_pos = $self->seq_to_align_pos(1, $bexon1->end);
    } elsif ($bstrand == -1) {
        $bexon1 = $bexons->[(scalar @$bexons) - 1];
        return undef if !$bexon1;
        $bexon1_pos = $self->seq_to_align_pos(1, $bexon1->start);
    } else {
        carp "unknown strand $bstrand\n";
        return undef;
    }

    my $cexon1;
    my $cexon1_pos;
    if ($cstrand == 1) {
        $cexon1 = $cexons->[0];
        return undef if !$cexon1;
        $cexon1_pos = $self->seq_to_align_pos(2, $cexon1->end);
    } elsif ($cstrand == -1) {
        $cexon1 = $cexons->[(scalar @$cexons) - 1];
        return undef if !$cexon1;
        $cexon1_pos = $self->seq_to_align_pos(2, $cexon1->start);
    } else {
        carp "unknown strand $bstrand\n";
        return undef;
    }

    return undef if !$bexon1_pos || !$cexon1_pos;

    return ($cexon1_pos - $bexon1_pos) * $bstrand;
}

################################################################################
#
# Internal (private) methods follow. These should not be called directly.
#
################################################################################

sub _create_pos_conversion_arrays
{
    my ($self) = @_;

    return if $self->{_pos_conversion_arrays_created};

    my @base_seq_to_aln;
    my @base_aln_to_seq;
    $base_seq_to_aln[0] = 0;
    $base_aln_to_seq[0] = 0;
    my $base_nucs      = $self->_aligned_base_nucs;
    my $base_nuc_count = 0;
    my $base_aln_idx   = 0;
    foreach my $base_nuc (@$base_nucs) {
        if ($base_nuc ne '-') {
            $base_nuc_count++;
            $base_seq_to_aln[$base_nuc_count] = $base_aln_idx + 1;
            $base_aln_to_seq[$base_aln_idx + 1] = $base_nuc_count;
        } else {
            $base_aln_to_seq[$base_aln_idx + 1] = 0;
        }
        $base_aln_idx++;
    }
    $self->{_base_to_align_pos} = \@base_seq_to_aln;
    $self->{_align_to_base_pos} = \@base_aln_to_seq;

    #my $cmp_strand = $self->alignment->get_seq_by_pos(2)->strand;

   #
   # Compute the alignment orientation. If the comparison sequence was reverse
   # complemented during the alignment process, then the orientation is -1,
   # otherwise it's +1. We are assuming that the orientation of the base
   # sequence is always the same as the aligned base sequence.
   #
    my $aln_orient;
    my $seq2     = $self->comparison_seq;
    my $aln_seq1 = $self->alignment->get_seq_by_pos(1);
    my $aln_seq2 = $self->alignment->get_seq_by_pos(2);
    if ($seq2->isa('Bio::LocatableSeq')) {
        #
        # Check if the second sequence is reverse complemented in the
        # alignment.
        #
        $aln_orient = $seq2->strand * $aln_seq2->strand;
    } else {
        #
        # We don't have strand info for the comparison sequence. Assume then
        # that if the two aligned sequences are on opposite strands that the
        # second sequence was reverse complemented in the alignment.
        #
        $aln_orient = $aln_seq1->strand * $aln_seq2->strand;
    }

    my @cmp_seq_to_aln;
    my @cmp_aln_to_seq;
    $cmp_seq_to_aln[0] = 0;
    $cmp_aln_to_seq[0] = 0;
    my $cmp_nucs      = $self->_aligned_comparison_nucs;
    my $cmp_nuc_count = 0;
    my $cmp_aln_idx   = 0;
    #if ($cmp_strand && ($cmp_strand eq '-' || $cmp_strand == -1)) {
    if ($aln_orient == -1) {
        my $ncmp_nucs   = $self->comparison_seq->length;
        my $cmp_seq_pos = 0;
        foreach my $cmp_nuc (@$cmp_nucs) {
            if ($cmp_nuc ne '-') {
                $cmp_nuc_count++;
                $cmp_seq_pos = $ncmp_nucs - $cmp_nuc_count + 1;
                $cmp_seq_to_aln[$cmp_seq_pos] = $cmp_aln_idx + 1;
                $cmp_aln_to_seq[$cmp_aln_idx + 1] = $cmp_seq_pos;
            } else {
                $cmp_aln_to_seq[$cmp_aln_idx + 1] = 0;
            }
            $cmp_aln_idx++;
        }
    } else {    # plus strand
        foreach my $cmp_nuc (@$cmp_nucs) {
            if ($cmp_nuc ne '-') {
                $cmp_nuc_count++;
                $cmp_seq_to_aln[$cmp_nuc_count] = $cmp_aln_idx + 1;
                $cmp_aln_to_seq[$cmp_aln_idx + 1] = $cmp_nuc_count;
            } else {
                $cmp_aln_to_seq[$cmp_aln_idx + 1] = 0;
            }
            $cmp_aln_idx++;
        }
    }
    $self->{_comparison_to_align_pos} = \@cmp_seq_to_aln;
    $self->{_align_to_comparison_pos} = \@cmp_aln_to_seq;

    $self->{_alignment_orientation} = $aln_orient;

    $self->{_pos_conversion_arrays_created} = 1;
}

sub _destroy_position_conversion_arrays
{
    my $self = shift;

    undef $self->{_alignment_orientation};
    undef $self->{_aligned_base_nucs};
    undef $self->{_aligned_comparison_nucs};
    undef $self->{_comparison_to_align_pos};
    undef $self->{_align_to_comparison_pos};
    undef $self->{_base_to_align_pos};
    undef $self->{_align_to_base_pos};

    $self->{_pos_conversion_arrays_created} = 0;
}

#
# Return a reference to an array of the aligned base sequence nucleotide
# characters. Create the array if it has not already been created.
#
sub _aligned_base_nucs
{
    my ($self) = @_;

    if (!defined $self->{_aligned_base_nucs}) {
        my $aligned_base_seq = $self->alignment->get_seq_by_pos(1)->seq;
        if (defined $aligned_base_seq) {
            $self->{_aligned_base_nucs} = [split //, $aligned_base_seq];
        }
    }

    return $self->{_aligned_base_nucs};
}

#
# Return a reference to an array of the comparison sequence nucleotide
# characters.  Create the array if it has not already been created.
#
sub _aligned_comparison_nucs
{
    my ($self) = @_;

    if (!defined $self->{_aligned_comparison_nucs}) {
        my $aligned_comparison_seq = $self->alignment->get_seq_by_pos(2)->seq;
        if (defined $aligned_comparison_seq) {
            $self->{_aligned_comparison_nucs} =
                [split //, $aligned_comparison_seq];
        }
    }

    return $self->{_aligned_comparison_nucs};
}

=head2 compute_conserved_tfbss

 Title    : compute_conserved_tfbss
 Usage    : $sites = $pwa->compute_conserved_tfbss(
                -matrix_set			        => $matrix_set,
                -tfbs_threshold			    => $tfbs_thresh,
                -min_tfbs_cr_overlap		=> $overlap,
                -filter_overlapping_sites	=> 1,
                -start				        => $start,
                -end				        => $end
            );
 Function : Compute and return the list of conserved TF sites.
 Returns  : A reference to a TFBS::SitePairSet object.
 Args     : -matrix_set                 => A TFBS::MatrixSet object
            -tfbs_threshold		        => Report only TFBSs with at least
                                           this score on both sequences
            -min_tfbs_cr_overlap	    => Report only TFBSs which overlap
                                           a conserved region by at least
                                           this number of nucleotides.
            -filter_overlapping_sites   => boolean; if true, for any given
                                           matrix, report only the best
                                           scoring of any two overlapping
                                           site pairs
            -start			            => restrict TFBS search region to
                                           start at this position
            -end			            => restrict TFBS search region to
                                           end at this position

=cut

sub compute_conserved_tfbss
{
    my ($self, %args) = @_;

    my $start = $self->start();
    if (!defined $start) {
        carp "Pairwise analysis object has no start defined\n";
        return;
    }

    my $end = $self->end();
    if (!defined $start) {
        carp "Pairwise analysis object has no end defined\n";
        return;
    }

    my $seq1 = $self->base_seq;
    if (!defined $seq1) {
        carp "ORCA::Analysis::Pairwise object has no base sequence defined\n";
        return;
    }

    my $seq2 = $self->comparison_seq;
    if (!defined $seq2) {
        carp "ORCA::Analysis::Pairwise object has no comparison sequence"
            . " defined\n";
        return;
    }

    my $align = $self->alignment;
    if (!defined $align) {
        carp "ORCA::Analysis::Pairwise object has no aligment defined\n";
        return;
    }

    my $conserved_regions = $self->conserved_regions;
    if (!defined $conserved_regions) {
        carp "ORCA::Analysis::Pairwise object has no conserved regions\n";
        return;
    }

    my $cutoff = $self->conservation_cutoff;
    if (defined $cutoff) {
        $cutoff /= 100;
        if ($cutoff < 0 || $cutoff > 1) {
            carp "ERROR: conservation cutoff is out of range\n";
            return;
        }
    }

    my $threshold = $self->tfbs_threshold($args{-tfbs_threshold});
    if (!defined $threshold) {
        carp "-tfbs_threshold argument not provided\n";
    }

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

    my $overlap =
        $self->min_tfbs_cr_overlap($args{-min_tfbs_cr_overlap}
            || DFLT_MIN_TFBS_CA_OVERLAP);

    my $filter_overlaps = $args{-filter_overlapping_sites} || 0;

    my $search_start1 = $args{-start};
    my $search_end1   = $args{-end};

    #
    # XXX if using relative coords
    #
    $search_start1  = 1 if !defined $search_start1;
    $search_end1    = $seq1->length() if !defined $search_end1;

    if ($search_end1 < 1 || $search_start1 > $seq1->length()) {
        carp "TF site search start/end outside of base sequence start/end"
            . " - make sure start/end is specified in relative coordinates\n";
        return;
    }
    #
    # Otherwise if using chromosomal coords
    #
    #$search_start1  = $start if !defined $search_start1;
    #$search_end1    = $end if !defined $search_end1;
    #
    #if ($search_end1 < $start || $search_start1 > $end) {
    #    carp "TF site search start/end outside of base sequence start/end"
    #        . " - make sure coordinate systems are consistent\n";
    #    return;
    #}

    my $aln_orient = $self->alignment_orientation;

    my $search_start2;
    my $search_end2;
    if ($aln_orient == -1) {
        $search_start2 = $self->convert_seq_pos(1, 2, $search_end1,   'le');
        $search_end2   = $self->convert_seq_pos(1, 2, $search_start1, 'ge');
    } else {
        $search_start2 = $self->convert_seq_pos(1, 2, $search_start1, 'le');
        $search_end2   = $self->convert_seq_pos(1, 2, $search_end1,   'ge');
    }
    $search_start2 = 1 if !$search_start2;
    $search_end2 = $seq2->length
        if !$search_end2 || $search_end2 > $seq2->length;

    my $cons_tf_count1 = 0;
    my $cons_tf_count2 = 0;

    my $site_pairs = TFBS::SitePairSet->new();

    #
    # The previous logic did not hold as there were cases (albeit rare) where
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
        foreach my $cr (@$conserved_regions) {
            if ($cutoff) {
                next if $cr->feature1->score < $cutoff;
            }

            my $cr_start1 = $cr->feature1->start;
            my $cr_end1   = $cr->feature1->end;

            #
            # If conserved region falls outside the sequence region of
            # interest exclude it from analysis
            #
            next if $cr_start1 > $search_end1;
            next if $cr_end1 < $search_start1;

            my $m_start1 = $cr_start1 - $matrix_width + $overlap;
            $m_start1 = $search_start1 if $m_start1 < $search_start1;

            my $m_end1 = $cr_end1 + $matrix_width - $overlap;
            $m_end1 = $search_end1 if $m_end1 > $search_end1;

            next if $m_end1 - $m_start1 + 1 < $matrix_width;

            my $cr_start2 = $cr->feature2->start;
            if (!$cr_start2) {
                if ($aln_orient == -1) {
                    $cr_start2 = $self->convert_seq_pos(1, 2, $cr_end1, 'le');
                } else {
                    $cr_start2 = $self->convert_seq_pos(1, 2, $cr_start1, 'le');
                }
            }
            $cr_start2 = $search_start2
                if !$cr_start2 || $cr_start2 < $search_start2;

            my $cr_end2 = $cr->feature2->end;
            if (!$cr_end2) {
                if ($aln_orient == -1) {
                    $cr_end2 = $self->convert_seq_pos(1, 2, $cr_start1, 'ge');
                } else {
                    $cr_end2 = $self->convert_seq_pos(1, 2, $cr_end1, 'ge');
                }
            }
            $cr_end2 = $search_end2 if !$cr_end2 || $cr_end2 > $search_end2;

            my $m_start2 = $cr_start2 - $matrix_width + $overlap;
            $m_start2 = $search_start2 if $m_start2 < $search_start2;

            my $m_end2 = $cr_end2 + $matrix_width - $overlap;
            $m_end2 = $search_end2 if $m_end2 > $search_end2;

            next if $m_end2 - $m_start2 + 1 < $matrix_width;

            if (DEBUG) {
                printf "\nRegion %d  %d  %d  %d  %.2f\n",
                    $cr_start1,
                    $cr_end1,
                    $cr_start2,
                    $cr_end2,
                    $cr->feature1->score;
            }

            my $cr_m_sites1 = _find_seq_tf_sites(
                $matrix, $seq1, $m_start1, $m_end1, $threshold
            );

            next if !$cr_m_sites1 || $cr_m_sites1->size == 0;
            if (DEBUG) {
                $self->_print_site_set(1, $cr_m_sites1) if defined $cr_m_sites1;
            }

            my $cr_m_sites2 = _find_seq_tf_sites(
                $matrix, $seq2, $m_start2, $m_end2, $threshold
            );

            next if !$cr_m_sites2 || $cr_m_sites2->size == 0;
            if (DEBUG) {
                $self->_print_site_set(2, $cr_m_sites2) if defined $cr_m_sites2;
            }

            my $cr_m_site_pairs =
                $self->_tf_sites_to_site_pairs($cr_m_sites1, $cr_m_sites2, 1);

            next if !$cr_m_site_pairs || $cr_m_site_pairs->size == 0;
            if (DEBUG) {
                print "\nRegion matrix site pair set\n";
                _print_site_pair_set($cr_m_site_pairs);
            }

            $self->_set_site_pair_conservation(
                $cr_m_site_pairs, $cr->feature1->score
            );

            $m_site_pairs->add_site_pair_set($cr_m_site_pairs);
        }

        next if !$m_site_pairs || $m_site_pairs->size == 0;
        if ($filter_overlaps) {
            my $filt_m_site_pairs =
                _tf_filter_overlapping_site_pairs($m_site_pairs);
            if (DEBUG) {
                print "\nFiltered site pair set\n";
                _print_site_pair_set($filt_m_site_pairs);
            }
            $site_pairs->add_site_pair_set($filt_m_site_pairs);
        } else {
            $site_pairs->add_site_pair_set($m_site_pairs);
        }
    }

    #print "seq1 conserved TF sites: $cons_tf_count1\n";
    #print "seq2 conserved TF sites: $cons_tf_count2\n";

    my @tfbss;
    my $iter = $site_pairs->Iterator(-sort_by => 'start');
    while (my $site_pair = $iter->next()) {
        push @tfbss, $site_pair;
    }

    #
    # Changed to return a listref of TFBS::SitePair objects rather than a
    # TFBS::SitePairSet object. DJA 11/10/19
    #
    $self->{-conserved_tfbss} = @tfbss ? \@tfbss : undef;
}

=head2 ucsc_track

 Title    : ucsc_track
 Usage    : $track = $pwa->ucsc_track();
 Function : Create a UCSC custom browser track
 Returns  : A custom UCSC browser track formatted string
 Args     : None

=cut

sub ucsc_track
{
    my $self = shift;

    my $chr     = $self->chr;
    my $start   = $self->start;
    my $end     = $self->end;

    unless ($chr && $start && $end) {
        carp "Can't create UCSC browser track without positional information\n";
        return;
    }

    if (   $chr !~ /scaffold/i
        && $chr !~ /contig/i
        && $chr !~ /ultra/i
        && $chr !~ /super/i)
    {
        $chr = "chr$chr";
    }

    my $track = '';
    $track = sprintf "browser position %s:%d-%d\n", $chr, $start, $end;

    my $crs = $self->conserved_regions;
    if ($crs) {
        $track .= "track name=\"Conserved Regions\" description=\"ORCA Conserved Regions\" color=0,255,255\n";

        my $num = 0;
        foreach my $cr (@$crs) {
            $num++;
            $track .= sprintf "%s\t%d\t%d\t%s\t%.3f\n",
                $chr,
                # Conserved regions are in relative coords - convert
                # to chromosomal
                $cr->start + $start - 1,
                $cr->end + $start - 1,
                "CR$num",
                $cr->score;
        }
    }

    my $tfbss = $self->conserved_tfbss;
    if ($tfbss) {
        $track .= "track name=\"TFBSs\" description=\"TF Binding Sites\" color=0,0,255 visibility=3\n";

        foreach my $site_pair (@$tfbss) {
            my $site1 = $site_pair->site1;
            $track .= sprintf "%s\t%d\t%d\t'%s'\t%.1f\n",
                $chr,
                # TFBSs are in relative coords - convert to chromosomal
                $site1->start + $start - 1,
                $site1->end + $start - 1,
                $site1->pattern->name,
                $site1->rel_score;
        }
    }

    my $cp = $self->conservation_profile;
    if ($cp) {
        $track .= "track type=wiggle_0 name=\"Conservation\" description=\"ORCA Conservation (pairwise)\" color=255,0,0 graphType=bar yLineMark=0.0 yLineOnOff=on visibility=2\n";

        $track .= sprintf "fixedStep chrom=%s start=%d step=1\n",
            $chr, $start;
        foreach my $score (@$cp) {
            $track .= sprintf "%.3f\n", $score;
        }
    }

    return $track;
}

################################################################################
#
# Getter/Setter methods follow.
#
################################################################################

=head2 matrix_set

 Title    : matrix_set
 Usage    : $matrix_set = $pwa->matrix_set($matrix_set);
 Function : Get/Set the TFBS::MatrixSet.
 Returns  : A TFBS::MatrixSet object.
 Args     : Optionally a new TFBS::MatrisSet object.

=cut

sub matrix_set
{
    my ($self, $matrix_set) = @_;

    if (defined($matrix_set)) {
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
 Usage    : $score = $pwa->min_base_tfbs_score($score);
 Function : Get/Set the minimum base sequence TFBS score.
 Returns  : An integer.
 Args     : Optionally a new minimum base sequence TFBS score.

=cut

sub min_base_tfbs_score
{
    my ($self, $score) = @_;

    if (defined($score)) {
        $self->{-min_base_tfbs_score} = $score;
    }

    return $self->{-min_base_tfbs_score};
}

=head2 min_comparison_tfbs_score

 Title    : min_comparison_tfbs_score
 Usage    : $score = $pwa->min_comparison_tfbs_score($score);
 Function : Get/Set the minimum comparison sequence TFBS score.
 Returns  : An integer.
 Args     : Optionally a new minimum comparison sequence TFBS score.

=cut

sub min_comparison_tfbs_score
{
    my ($self, $score) = @_;

    if (defined($score)) {
        $self->{-min_comparison_tfbs_score} = $score;
    }

    return $self->{-min_comparison_tfbs_score};
}

=head2 tfbs_threshold

 Title    : tfbs_threshold
 Usage    : $threshold = $pwa->tfbs_threshold($threshold);
 Function : Get/Set the TFBS threshold.
 Returns  : An integer.
 Args     : Optionally a new TFBS threshold.

=cut

sub tfbs_threshold
{
    my ($self, $threshold) = @_;

    if (defined($threshold)) {
        $self->{-tfbs_threshold} = $threshold;
    }

    return $self->{-tfbs_threshold};
}

=head2 min_tfbs_cr_overlap

 Title    : min_tfbs_cr_overlap
 Usage    : $overlap = $pwa->min_tfbs_cr_overlap($overlap);
 Function : Get/Set the minimum TFBS conservation overlap.
 Returns  : An integer.
 Args     : Optionally a new minimum TFBS conservation overlap.

=cut

sub min_tfbs_cr_overlap
{
    my ($self, $overlap) = @_;

    if (defined($overlap)) {
        $self->{-min_tfbs_cr_overlap} = $overlap;
    }

    return $self->{-min_tfbs_cr_overlap};
}

=head2 conserved_tfbss

 Title    : conserved_tfbss
 Usage    : $site_pairs = $pwa->conserved_tfbss();
 Function : Get the conserved TFBSs.
 Returns  : A list ref of TFBS::SitePair objects.
 Args     : None.

=cut

sub conserved_tfbss
{
    $_[0]->{-conserved_tfbss};
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
        -seqobj    => $seqobj,
        -threshold => $threshold,
        -subpart   => {
            -start => $start,
            -end   => $end
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
                if ($temp_pairs[$i]->site1->pattern->name eq
                    $temp_pairs[$j]->site1->pattern->name)
                {
                    if ($temp_pairs[$i]
                        ->site1->overlaps($temp_pairs[$j]->site1))
                    {
                        #
                        # Bias is toward the site pair with the lower start
                        # site (i.e. if the scores are equal). Note that
                        # strand is not taken into account so this is not
                        # necessarily the most upstream site.
                        #
                        if (
                            _cmp_site_pair_score($temp_pairs[$i],
                                $temp_pairs[$j]) < 0
                            )
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
    my $i              = 0;
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
    my ($self, $site_set1, $site_set2, $exact) = @_;

    return if !defined $site_set1
        || !defined $site_set2
        || $site_set1->size == 0
        || $site_set2->size == 0;

    my $aln_orient = $self->alignment_orientation;

    my $site_pair_set = TFBS::SitePairSet->new();

    if ($exact) {
        my $iter1 = $site_set1->Iterator();
        while (my $site1 = $iter1->next) {
            my $site1_aln_start = $self->seq_to_align_pos(1, $site1->start);
            my $iter2 = $site_set2->Iterator();
            while (my $site2 = $iter2->next) {
                my $site2_aln_start;
                if ($aln_orient == -1) {
                    $site2_aln_start = $self->seq_to_align_pos(2, $site2->end);
                } else {
                    $site2_aln_start =
                        $self->seq_to_align_pos(2, $site2->start);
                }
                if ($site1_aln_start == $site2_aln_start) {
                    $site_pair_set->add_site_pair(
                        TFBS::SitePair->new($site1, $site2));
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
        #	    my $site1_aln_start = $self->seq_to_align_pos(1, $site1->start);
        #	    my $site2_aln_start;
        #	    if ($cmp_strand == -1 || $cmp_strand eq '-') {
        #		$site2_aln_start = $self->seq_to_align_pos(2, $site2->end);
        #	    } else {
        #		$site2_aln_start = $self->seq_to_align_pos(2, $site2->start);
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
        #		    $site1_aln_start = $self->seq_to_align_pos(1, $site1->start);
        #		    if ($site1_aln_start == $prev_site2_aln_start) {
        #			$site_pair_set->add_site_pair(
        #			    TFBS::SitePair->new($site1, $prev_site2));
        #		    }
        #		}
        #		if (defined $site2) {
        #		    if ($cmp_strand == -1 || $cmp_strand eq '-') {
        #			$site2_aln_start
        #				= $self->seq_to_align_pos(2, $site2->end);
        #		    } else {
        #			$site2_aln_start
        #				= $self->seq_to_align_pos(2, $site2->start);
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
    my ($self, $sites) = @_;

    return if !$sites || $sites->size == 0;

    my $regions = $self->conserved_regions;
    return if !$regions;

    my $min_ol = $self->min_tfbs_cr_overlap || 0;

    my $iter = $sites->Iterator();
    return if !$iter;
    while (my $sp = $iter->next) {
        my $site1        = $sp->site1;
        my $conservation = 0;
        foreach my $reg (@$regions) {
            if (   ($site1->end - $reg->feature1->start + 1 >= $min_ol)
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
    my ($self, $site_set1, $site_set2, $exact) = @_;

    return if !defined $site_set1
            || !defined $site_set2
            || $site_set1->size == 0
            || $site_set2->size == 0;

    my $site_pair_set = TFBS::SitePairSet->new();

    if ($exact) {
        my $iter1 = $site_set1->Iterator();
        while (defined(my $site1 = $iter1->next)) {
            my $iter2 = $site_set2->Iterator();
            while (defined(my $site2 = $iter2->next)) {
                if ($site1->pattern->name eq $site2->pattern->name) {
                    my $site1_aln_start =
                        $self->seq_to_align_pos(1, $site1->start);
                    my $site2_aln_start =
                        $self->seq_to_align_pos(2, $site2->start);

                    #
                    # XXX Should we also look at the end positions???
                    #
                    if ($site1_aln_start == $site2_aln_start) {
                        $site_pair_set->add_site_pair(
                            TFBS::SitePair->new($site1, $site2));
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
        while (defined(my $site1 = $iter1->next)) {
            my $iter2 = $site_set2->Iterator();
            while (defined(my $site2 = $iter2->next)) {
                if ($site1->pattern->name eq $site2->pattern->name) {
                    my $site1_aln_start =
                        $self->seq_to_align_pos(1, $site1->start);
                    my $site1_aln_end = $self->seq_to_align_pos(1, $site1->end);
                    my $site2_aln_start =
                        $self->seq_to_align_pos(2, $site2->start);
                    my $site2_aln_end = $self->seq_to_align_pos(2, $site2->end);

                    if (   $site1_aln_start <= $site2_aln_end
                        && $site1_aln_end >= $site2_aln_start)
                    {
                        $site_pair_set->add_site_pair(
                            TFBS::SitePair->new($site1, $site2));
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
    my $cons1  = ($feat11->get_tag_values('conservation'))[0];
    my $cons2  = ($feat21->get_tag_values('conservation'))[0];

    if (   defined $cons1
        && defined $cons2
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
    return (   ($feat11->score <=> $feat21->score)
            || ($feat12->score <=> $feat22->score));
}

#
# For debugging purposes
#
sub _print_site_set
{
    my ($self, $site_seq, $site_set) = @_;

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
            $self->seq_to_align_pos($site_seq, $site->start),
            $self->seq_to_align_pos($site_seq, $site->end),
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
