=head1 NAME

ORCA::ConservationAnalysis::Pairwise - Object for facilitating pairwise
conservation analysis of two aligned DNA sequences.

=head1 SYNOPSIS

    use ORCA::ConservationAnalysis::Pairwise;

    # Minimally, create a new ORCA::ConservationAnalysis::Pairwise object with
    # the given sequences.
    $ca = ORCA::ConservationAnalysis::Pairwise->new(
			-base_seq		=> $seq1,
			-comparison_seq		=> $seq2);

    # Optionally also provide masked sequences, alignment and sequence exons.
    $ca = ORCA::ConservationAnalysis::Pairwise->new(
			-base_seq		=> $seq1,
			-comparison_seq		=> $seq2,
			-masked_base_seq	=> $masked_seq1,
			-masked_comparison_seq	=> $masked_seq2,
			-alignment		=> $align,
			-base_seq_exons		=> $exons1,
			-comparison_seq_exons	=> $exons2);

    # If not explicitly provided, compute the alignment.
    $ca->align();

    # Set conservation profile parameters
    $ca->param('position_type', 'c');
    $ca->param('window_size', 100);
    $ca->param('window_inc', 1);
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
    
    # Explicitly compute and retrieve the conservation profile.
    my $profile = $ca->compute_conservation_profile();

    # Optionally provide the parameters directly to the method (see above).
    my $profile = $ca->compute_conservation_profile(
					    -position_type	=> 'c',
					    -window_size	=> 100,
					    -window_inc		=> 1
					);

    # The computed conservation profile may be retrieved again later
    # (if not already explicitly computed, automatically compute it) 
    my $profile = $ca->conservation_profile();

    # Set conserved region parameters
    $ca->param('top_pct', '10%');
    $ca->param('min_conservation', '70%');
    $ca->param('filter_exons', 1);
    $ca->param('min_cr_len', 20);

    # where
    	top_pct			= report the top X percentile of conserved
				  regions (specify as a number between
				  0 and 1 or a string such as '10%')
    	min_conservation	= the absolute minimum percent identity of
				  regions to report (specify as a number
				  between 0 and 1 or a string such as
				  '70%')
	filter_exons		= (bool) exons are to be filtered out of
				  the conserved regions (if base_seq_exons
				  has been set).
	min_cr_len		= (int) min. length of a conserved region
				  to report.
	pct_id_method		= method used to compute the percent
				  identity ('s' = standard, 'o' = overall)

    # Explicitly compute and retrieve the conserved regions.
    my $regions = $ca->compute_conserved_regions();

    # Optionally provide the conserved regions parameters directly to the
    # method (see above).
    my $regions = $ca->compute_conserved_regions(
					    -top_pct		=> '10%');
					    -min_conservation	=> '70%');
					    -filter_exons	=> 1,
					    -min_cr_len		=> 20,
					    -pct_id_method	=> 's');

    # The conserved regions can be obtained later
    # (if not already explicitly computed, automatically compute them) 
    my $regions = $ca->conserved_regions();

    # Get the conserved sub-sequences corresponding to the conserved
    # regions.
    my $subseqs = $ca->extract_conserved_subsequences();
    
    # The conserved sub-sequences can be obtained later
    # (compute if not already explicitly computed)
    my $subseqs = $ca->conserved_subsequences();

    # Get the conserved sub-alignments corresponding to the conserved
    # regions.
    my $subaligns = $ca->extract_conserved_subalignments();
    
    # The conserved sub-alignments can be obtained later
    # (compute if not already explicitly computed)
    my $subaligns = $ca->conserved_subalignments();

=head1 DESCRIPTION

ORCA::ConservationAnalysis::Pairwise is an object used for the purpose
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

package ORCA::ConservationAnalysis::Pairwise;

use strict;

use Carp;
use ORCA::ConservationAnalysis::Run::AlignCons;
use ORCA::ConservationAnalysis::ConservationReport;

=head2 new

 Title    : new
 Usage    : $ca = ORCA::ConservationAnalysis::Pairwise->new();
 Function : Create a new ORCA::ConservationAnalysis::Pairwise object.
 Returns  : A new ORCA::ConservationAnalysis::Pairwise object.
 Args     : Optional named parameters:
	     -base_seq			=> Bio::Seq object defining the base
	      			           sequence
	     -comparison_seq		=> Bio::Seq object defining the
	      			           comparison sequence
	     -masked_base_seq		=> Bio::Seq object defining repeat
	     				   masked base sequence
	     -masked_comparison_seq	=> Bio::Seq object defining repeat
	     				   masked comparison sequence
	     -alignment			=> A Bio::SimpleAlign object
	     				   defining the alignment between
					   base_seq and comparison_seq
	     -base_seq_exons		=> A listref of Bio::SeqFeature
	     				   objects defining exon positions
					   on base_seq
	     -comparison_seq_exons	=> A listref of Bio::SeqFeature
	     				   objects defining exon positions
					   on comparison_seq

=cut

sub new
{
    my ($class, %args) = @_;
    my $self = bless {
		    -base_seq			=> undef,
		    -comparison_seq		=> undef,
		    -masked_base_seq		=> undef,
		    -masked_comparison_seq	=> undef,
		    -alignment			=> undef,
		    -base_seq_exons		=> undef,
		    -comparison_seq_exons	=> undef,
		    -conserved_regions		=> undef,
		    -conservation_profile	=> undef,
		    -conserved_subseqs		=> undef,
		    -conserved_subaligns	=> undef,
		    -parameters			=> {
			window_size		=> undef,
			window_inc		=> undef,
			position_type		=> undef,
			top_pct			=> undef,
			min_conservation	=> undef,
			filter_exons		=> undef,
			min_cr_len		=> undef
		    },
		    %args
		}, ref $class || $class;

    if (defined $self->{-base_seq}) {
	if (!$self->{-base_seq}->isa("Bio::SeqI")
		    && !$self->{-base_seq}->isa("Bio::PrimarySeqI"))
	{
	    carp "base_seq is not a Bio::SeqI or Bio::PrimarySeqI compliant"
		    . " object\n";
	    return;
	}
    } else {
    	carp "must provide a base sequence\n";
    }

    if (defined $self->{-comparison_seq}) {
	if (!$self->{-comparison_seq}->isa("Bio::SeqI")
		    && !$self->{-comparison_seq}->isa("Bio::PrimarySeqI"))
	{
	    carp "comparison_seq is not a Bio::SeqI or Bio::PrimarySeqI"
		    . " compliant object\n";
	    return;
	}
    } else {
    	carp "must provide a comparison sequence\n";
    }

    if (defined $self->{-masked_base_seq}) {
	if (!$self->{-masked_base_seq}->isa("Bio::SeqI")
		&& !$self->{-masked_base_seq}->isa("Bio::PrimarySeqI"))
	{
	    carp "masked_base_seq is not a Bio::SeqI or Bio::PrimarySeqI"
	    	. " compliant object\n";
	    return;
	}
    }

    if (defined $self->{-masked_comparison_seq}) {
	if (!$self->{-masked_comparison_seq}->isa("Bio::SeqI")
		&& !$self->{-masked_comparison_seq}->isa("Bio::PrimarySeqI"))
	{
	    carp "masked_comparison_seq is not a Bio::SeqI or Bio::PrimarySeqI"
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
	if (ref $self->{-base_seq_exons} ne "ARRAY"
		|| !$self->{-base_seq_exons}->[0]->isa("Bio::SeqFeatureI")) {
	    carp "base_seq_exons is not a list ref of Bio::SeqFeatureI"
	    		. " compliant objects\n";
	    return;
	}
    }

    if (defined $self->{-comparison_seq_exons}) {
	if (ref $self->{-comparison_seq_exons} ne "ARRAY"
	    || !$self->{-comparison_seq_exons}->[0]->isa("Bio::SeqFeatureI"))
	{
	    carp "comparison_seq_exons is not a list ref of Bio::SeqFeatureI"
	    		. " compliant objects\n";
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

=head2 base_seq

 Title    : base_seq
 Usage    : $base_seq = $ca->base_seq($seq);
 Function : Get/Set the base sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new base sequence (Bio::Seq object)

=cut

sub base_seq
{
    my ($self, $base_seq) = @_;

    if (defined ($base_seq)) {
	if (!$base_seq->isa("Bio::SeqI")
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
 Usage    : $comparison_seq = $ca->comparison_seq($seq);
 Function : Get/Set the comparison sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new comparison sequence (Bio::Seq object)

=cut

sub comparison_seq
{
    my ($self, $comparison_seq) = @_;

    if (defined ($comparison_seq)) {
	if (!$comparison_seq->isa("Bio::SeqI")
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
 Usage    : $masked_base_seq = $ca->masked_base_seq($seq);
 Function : Get/Set the masked base sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new masked base sequence (Bio::Seq object)

=cut

sub masked_base_seq
{
    my ($self, $masked_base_seq) = @_;

    if (defined ($masked_base_seq)) {
	if (!$masked_base_seq->isa("Bio::SeqI")
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
 Usage    : $masked_comparison_seq = $ca->masked_comparison_seq($seq);
 Function : Get/Set the masked comparison sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new masked comparison sequence (Bio::Seq object)

=cut

sub masked_comparison_seq
{
    my ($self, $masked_comparison_seq) = @_;

    if (defined ($masked_comparison_seq)) {
	if (!$masked_comparison_seq->isa("Bio::SeqI")
	    	&& !$masked_comparison_seq->isa("Bio::PrimarySeqI"))
	{
	    carp "masked_comparison_seq is not a Bio::SeqI or Bio::PrimarySeqI"
	    	. " compliant object\n";
	    return undef;
	}
	$self->{-masked_comparison_seq} = $masked_comparison_seq;
    }

    return $self->{-masked_comparison_seq};
}

=head2 alignment

 Title    : alignment
 Usage    : $align = $ca->alignment($align);
 Function : Get/Set the alignment.
 Returns  : A Bio::SimpleAlign object or undef.
 Args     : Optionally a new alignment (Bio::SimpleAlign) object.

=cut

sub alignment
{
    my ($self, $align) = @_;

    if (defined ($align)) {
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
 Usage    : $exons = $ca->base_seq_exons();
 Function : Get/Set the base sequence exons.
 Returns  : A reference to a list of Bio::SeqFeatureI objects.
 Args     : Optionally new base sequence exons (reference to a
            list of Bio::SeqFeatureI objects)

=cut

sub base_seq_exons
{
    my ($self, $exons) = @_;

    if (defined $exons) {
	if (ref $exons ne "ARRAY" || !$exons->[0]->isa("Bio::SeqFeatureI")) {
	    carp "exons is not a list ref of Bio::SeqFeatureI compliant"
		    . " objects\n";
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
	if (ref $exons ne "ARRAY" || !$exons->[0]->isa("Bio::SeqFeatureI")) {
	    carp "exons is not a list ref of Bio::SeqFeatureI compliant"
		    . " objects\n";
	    return undef;
	}
	$self->{-comparison_seq_exons} = $exons;
    }

    return $self->{-comparison_seq_exons};
}

=head2 alignment_orientation

 Title    : alignment_orientation
 Usage    : $align_orient = $ca->alignment_orientation();
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

=head2 conservation_cutoff

 Title    : conservation_cutoff
 Usage    : $cutoff = $ca->conservation_cutoff();
 Function : Get the conservation cutoff (the effective % identity
 	    threshold computed by AlignCons from the min_conservation
	    and top_pct).
 Returns  : A float in the range 0 - 1.
 Args     : None.

=cut

sub conservation_cutoff
{
    my $self = shift;

    return $self->{_conservation_cutoff};
}

=head2 param

 Title    : param
 Usage    : $value = $ca->param($param, $value); OR $params = $ca->param;
 Function : Get/set the value of a conservation parameter or get a
	    listref of all parameter key/value pairs.
 Returns  : If a parameter name is provided, the value of the parameter.
 	    Otherwise a listref of parameter key/value pairs.
 Args     : Optionally the name and new value of a parameter.

=cut

sub param
{
    my ($self, $key, $value) = @_;

    if ($key && defined $value) {
	$self->{-parameters}->{-$key} = $value;
	return $value;
    }
    
    if ($key) {
	return $self->{-parameters}->{-$key};
    }

    return $self->{-parameters};
}

=head2 conservation_profile

 Title    : conservation_profile
 Usage    : $profile = $ca->conservation_profile();
 Function : Get the conservation profile.
 Returns  : A reference to a list of hashes containing '-position' and
 	    '-score' key/value pairs.
 Args     : None.

=cut

sub conservation_profile
{
    my $self = shift;

    if (!$self->{-conservation_profile}) {
    	$self->compute_conservation_profile();
    }

    return $self->{-conservation_profile};
}

=head2 conserved_regions

 Title    : conserved_regions
 Usage    : $regions = $ca->conserved_regions()
 	    OR $ca->conserved_regions($regions);
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
 Usage    : $cp_report = $ca->conservation_profile_report();
 Function : Get the conservation profile report.
 Returns  : A ORCA::ConservationAnalysis::ConservationReport object or
	    undef.
 Args     : None.

=cut

sub conservation_profile_report
{
    $_[0]->{-conservation_profile_report};
}

=head2 conserved_regions_report

 Title    : conserved_regions_report
 Usage    : $cr_report = $ca->conserved_regions_report();
 	    OR $ca->conserved_regions_report($report);
 Function : Get/set the conserved regions report.
 Returns  : A ORCA::ConservationAnalysis::ConservationReport object or
	    undef.
 Args     : Optionally a ORCA::ConservationAnalysis::ConservationReport.

=cut

sub conserved_regions_report
{
    my ($self, $report) = @_;

    if ($report) {
	if (!$report->isa("ORCA::ConservationAnalysis::ConservationReport")) {
	    carp "report argument is not an "
		    . " ORCA::ConservationAnalysis::ConservationReport"
		    . " object\n";
	    return undef;
	}
	$self->{-conserved_regions_report} = $report;
    }

    return $self->{-conserved_regions_report};
}

=head2 conserved_subsequences

 Title    : conserved_subsequences
 Usage    : $subseqs = $ca->conserved_subsequences();
 Function : Get the list of conserved sub-sequences.
 Returns  : A reference to a list of Bio::Seq objects.
 Args     : None.

=cut

sub conserved_subsequences
{
    my $self = shift;

    if (!$self->{-conserved_subseqs}) {
	$self->extract_conserved_subsequences();
    }

    return $self->{-conserved_subseqs};
}

=head2 conserved_subalignments

 Title    : conserved_subalignments
 Usage    : $subaligns = $ca->conserved_subalignments();
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
 Usage    : $alignment = $ca->compute_alignment();
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
				-seq1	=> $self->masked_base_seq
						|| $self->base_seq,
				-seq2	=> $self->masked_comparison_seq
						|| $self->comparison_seq);
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
 Usage    : $profile = $ca->compute_conservation_profile();
 Function : Compute and return the conservation profile.
 Returns  : A listref of hashes containing '-position' and '-score'
	    key/value pairs.
 Args     : Optionally
 		-position_type	=> ('c', 's', or 'e') indicating positions
				   should be reported as the center, start
				   or end position of the scoring window.

=cut

sub compute_conservation_profile
{
    my ($self, %args) = @_;

    my $align = $self->alignment;
    if (!$align) {
    	$align = $self->compute_alignment;

	if (!$align) {
	    carp "error computing alignment\n";
	    return undef;
	}
    }
    
    my $win_size = $self->param('window_size', $args{-window_size});
    my $win_inc =  $self->param('window_inc', $args{-window_inc});
    my $pos_type = $self->param('position_type', $args{-position_type});

    my $align_cons = ORCA::ConservationAnalysis::Run::AlignCons->new();
    if (!$align_cons) {
	carp "error initializing ORCA::ConservationAnalysis::Run::AlignCons\n";
	return undef;
    }
 
    my %ac_run_params = (-alignment => $align, -r => 'p');
    $ac_run_params{-w} = $win_size if $win_size;
    $ac_run_params{-n} = $win_inc if $win_inc;
    $ac_run_params{-f} = $pos_type if $pos_type;

    my $report = $align_cons->run(%ac_run_params);

    $self->{-conservation_profile_report} = $report;
    $self->{-conservation_profile} = $report->conservation_profile;
}

=head2 compute_conserved_regions

 Title    : compute_conserved_regions
 Usage    : $regions = $ca->compute_conserved_regions();
 Function : Compute and return the list of conserved regions.
 Returns  : A reference to a list of Bio::SeqFeature::Generic objects.
 Args     : Optionally:
		-top_pct	=> report the top X percentile of conserved
				   regions (specify as a number between
				   0 and 1 or a string such as '10%')
		-min_conservation
				=> the absolute minimum percent identity of
				  regions to report (specify as a number
				  between 0 and 1 or a string such as
				  '70%')
		-filter_exons	=> indicates exons should be
				   filtered out of conserved
				   regions.
		-min_cr_len 	=> minimum length of conserved
				   region to keep when exons are
				   filtered.
		-pct_id_method	=> method used to compute the percent
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
    
    my $win_size = $self->param('window_size');
    my $win_inc = $self->param('window_inc');
    my $top_pct = $self->param('top_pct', $args{-top_pct});
    my $min_conservation = $self->param('min_conservation',
					    $args{-min_conservation});
    my $filter_exons = $self->param('filter_exons', $args{-filter_exons});

    #
    # Now set minimum conserved region length regardless of whether
    # filter_exons is set. DJA 2008/01/10
    #
    my $min_cr_len = $self->param('min_cr_len', $args{-min_cr_len});

    my $pct_id_method = $self->param('pct_id_method', $args{-pct_id_method});

    my $align_cons = ORCA::ConservationAnalysis::Run::AlignCons->new();
    if (!$align_cons) {
	carp "error initializing ORCA::ConservationAnalysis::Run::AlignCons\n";
	return undef;
    }

    my %ac_run_params = (-alignment => $aln);
    if ($filter_exons && $self->base_seq_exons) {
	$ac_run_params{-features} = $self->base_seq_exons;
    }
    $ac_run_params{-r} = 'c';
    $ac_run_params{-w} = $win_size if $win_size;
    $ac_run_params{-n} = $win_inc if $win_inc;
    $ac_run_params{-s} = $top_pct if $top_pct;
    $ac_run_params{-t} = $min_conservation if $min_conservation;
    $ac_run_params{-l} = $min_cr_len if $min_cr_len;
    $ac_run_params{-m} = $pct_id_method if $pct_id_method;

    #
    # Now set the aligment reverse complemented flag as an ac_params value
    # rather than letting the ORCA::ConservationAnalysis::Run::AlignCons module
    # computed based on the raw alignment. DJA 2008/01/10
    #
    if ($self->alignment_orientation == -1) {
    	$ac_run_params{-c} = 1;
    }

    my $report = $align_cons->run(%ac_run_params);

    $self->conserved_regions_report($report) if $report;

    my $conservation_cutoff = $min_conservation if $min_conservation;
    $conservation_cutoff = $report->param('cutoff') / 100
					    if $report->param('cutoff'); 
    $self->{_conservation_cutoff} = $conservation_cutoff;

    # check if undef before assigning to prevent infinite recursion
    my $crfp = $report->conserved_regions_as_feature_pairs;
    $self->conserved_regions($crfp) if $crfp;
}

=head2 extract_conserved_subsequences

 Title    : extract_conserved_subsequences
 Usage    : $subseqs = $ca->extract_conserved_subsequences($seq_num);
 Function : Return the sub-sequences corresponding to the conserved regions.
 Returns  : A reference to a list of Bio::Seq objects.
 Args     : Optionally specify a sequence number; default = 1.

=cut

sub extract_conserved_subsequences
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
	    $cr_end = $cr->feature1->end;
	} elsif ($seq_num == 2) {
	    $cr_start = $cr->feature2->start;
	    $cr_end = $cr->feature2->end;
	}

	my $subseq = $seq->trunc($cr_start, $cr_end);
    
	$seq_display_id = "SUBSEQ" if !$seq_display_id;
	my $subseq_display_id = $seq_display_id . '/' . "$cr_start-$cr_end";
	$subseq->display_id($subseq_display_id);
	push @subseqs, $subseq;
    }

    $self->{-conserved_subseqs} = @subseqs ? \@subseqs : undef;
}

=head2 extract_conserved_subalignments

 Title    : extract_conserved_subalignments
 Usage    : $subaligns = $ca->extract_conserved_subalignments();
 Function : Return the sub-alignments corresponding to the conserved
 	    regions.
 Returns  : A reference to a list of Bio::SimpleAlign objects.
 Args     : None.

=cut

sub extract_conserved_subalignments
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
	my $align_end = ($reg->feature1->get_tag_values('align_end'))[0];

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
 Usage    : $align_pos = $ca->seq_to_align_pos(
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
	if ($seq_pos > 0 && $seq_pos < scalar @{$self->{_base_to_align_pos}}) {
	    $align_pos = $self->{_base_to_align_pos}->[$seq_pos];
	}
    } elsif ($which_seq == 2) {
	if ($seq_pos > 0
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
 Usage    : $seq_pos = $ca->align_to_seq_pos($which_seq, $align_pos,
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
    $match_type = 'le' if  $match_type eq 'lt' || $match_type eq '<'
    			|| $match_type eq '<=';
    $match_type = 'ge' if  $match_type eq 'gt' || $match_type eq '>'
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
		while (!$seq_pos &&
		    	$aln_pos < scalar @{$self->{_align_to_base_pos}})
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
		while (!$seq_pos &&
		    	$aln_pos < scalar @{$self->{_align_to_comparison_pos}})
		{
		    $seq_pos = $self->{_align_to_comparison_pos}->[$aln_pos];
		    $aln_pos++;
		}
	    }
	}
    }

    if ($coord_type && $seq_pos) {
    	# convert to absolute coordinates
    	$seq_pos = $seq_pos + $self->alignment->get_seq_by_pos(
					$which_seq)->start - 1;
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
 Usage    : $to_pos = $ca->convert_seq_pos($from, $to, $from_pos,
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

    return $self->align_to_seq_pos(
    			$to,
    			$self->seq_to_align_pos($from, $from_pos, $coord_type),
						    $match_type, $coord_type);
}

=head2 exon1_offset

 Title    : exon1_offset
 Usage    : $offset = $ca->exon1_offset();
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

    my $align = $self->alignment;
    my $bexons = $self->base_seq_exons;
    my $bstrand = $self->base_seq->strand;
    my $cexons = $self->comparison_seq_exons;
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
    my $base_nucs = $self->_aligned_base_nucs;
    my $base_nuc_count = 0;
    my $base_aln_idx = 0;
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
    my $seq2 = $self->comparison_seq;
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
    my $cmp_nucs = $self->_aligned_comparison_nucs;
    my $cmp_nuc_count = 0;
    my $cmp_aln_idx = 0;
    #if ($cmp_strand && ($cmp_strand eq '-' || $cmp_strand == -1)) {
    if ($aln_orient == -1) {
	my $ncmp_nucs = $self->comparison_seq->length;
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
    } else {	# plus strand
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
	my $aligned_base_seq
		= $self->alignment->get_seq_by_pos(1)->seq;
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
	my $aligned_comparison_seq
		= $self->alignment->get_seq_by_pos(2)->seq;
	if (defined $aligned_comparison_seq) {
	    $self->{_aligned_comparison_nucs}
		    = [split //, $aligned_comparison_seq];
	}
    }

    return $self->{_aligned_comparison_nucs};
}

1;
