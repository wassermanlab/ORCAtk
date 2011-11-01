=head1 NAME

ORCA::Analysis::PhastCons - Object for facilitating conserved TFBS analysis
of a DNA sequence using phastCons.

=head1 SYNOPSIS

 use ORCA::Analysis::PhastCons;

 #
 # Create a new ORCA::Analysis::PhastCons object with
 # the given sequence and (UCSC) database and track names.
 # Optionally also provide exons positions to mask.
 #
 $phca = ORCA::Analysis::PhastCons->new(
     -db         => 'hg18',
     -track      => 'phastCons28way',
     -chr        => '18',
     -start      => 55033380,
     -end        => 55048980,
     -seq        => $seq,
     -exons      => $exons
 );

 #
 # Compute and return the conservation profile.
 #
 my $profile = $phca->compute_conservation_profile();
 
 #
 # The conservation profile can be retrieved again later.
 # (Automatically computed if not already)
 #
 my $profile = $phca->conservation_profile();

 #
 # Compute and return the conserved regions.
 #
 my $regions = $phca->compute_conserved_regions(
     -min_conservation   => '70%',
     -min_length         => 20,
     -filter_exons       => 1,
     -flank_size         => 100
 );

 #
 # The conserved regions can be retrieved again later.
 # (Automatically computed if not already)
 #
 my $regions = $phca->conserved_regions();

 #
 # Compute and return the sub-sequences corresponding to the conserved
 # regions. If the regions have been filtered, returns the sub-sequences
 # corresponding to the filtered conserved regions, otherwise returns
 # the sub-sequence corresponding to the unfiltered conserved regions.
 #
 my $subseqs = $phca->compute_conserved_subsequences();
 
 #
 # The conserved sub-sequences can be retrieved again later
 #
 my $subseqs = $phca->conserved_subsequences();

 #
 # Search TFBSs which meet the conservation criteria.
 #
 my $tfbss = $phca->compute_conserved_tfbss(
     -matrix_set                => $matrix_set,
     -min_tfbs_score            => $min_tfbs_score,
     -min_tfbs_cr_overlap       => $overlap,
     -filter_overlapping_tfbss  => 1,
     -start                     => $start,
     -end                       => $end
 );

 # Once computed, the conserved TFBSs can be obtained later...
 my $tfbss = $phca->conserved_tfbss();

=head1 DESCRIPTION

ORCA::Analysis::PhastCons is an object used for the purpose
of performing various analysis tasks on a DNA sequence by using phastCons
scores as a conservation filter. These tasks include extracting the
phastCons score profile from a UCSC database, finding the conserved
regions, filtering exon features out of the conserved regions and
extracting the conserved sub-sequences.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::Analysis::PhastCons;

use strict;

use Carp;
use Bio::SeqFeature::Generic;
use ORCA::Analysis::Run::hgWiggle;

use constant DEBUG                    => 0;
use constant DFLT_MIN_TFBS_CA_OVERLAP => 1;

=head2 new

 Title    : new
 Usage    : $phca = ORCA::Analysis::PhastCons->new();
 Function : Create a new ORCA::Analysis::PhastCons object.
 Returns  : A new ORCA::Analysis::PhastCons object.
 Args     : Optional named parameters:
            -db         => Name of the UCSC database containing phastCons
                           scores
            -track      => Name of the UCSC track (table name)  containing
                           the phastCons scores
            -chr        => Name of chromosome for the region
            -start      => Optional chromosomal start position of the
                           region. If not provided and sequence is provided
                           and is a Bio::LocatableSeq object, set to start
                           of sequence, otherwise set to 1 (relative coord).
            -end        => Optional chromosomal end position of the region.
                           If not provided and sequence is provided and is
                           a Bio::LocatableSeq object, set to end of
                           sequence, otherwise set to length of sequence
                           (relative coord).
            -seq        => Optional Bio::Seq object defining the sequence
            -exons      => Optional reference to a list of exon objects
                           (takes any objects which implement start/end
                           methods)

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        -db         => undef,
        -track      => undef,
        -chr        => undef,
        -start      => undef,
        -end        => undef,
        -seq        => undef,
        -exons      => undef,
        %args
    }, ref $class || $class;

    #
    # Do not require phastCons DB and track as these are only required for
    # computing conservation. For some analysis we may want to us this
    # object for searching for TFBSs using pre-computed conserved regions.
    #
    #if (!defined $self->db()) {
    #    carp "Must provide a UCSC database name\n";
    #    return;
    #}
    #
    #if (!defined $self->track()) {
    #    carp "Must provide a UCSC track (table) name\n";
    #    return;
    #}

    if (!defined $self->chr()) {
        if (defined $args{-chr_name}) {
            carp "Deprecated tag -chr_name - please use -chr instead\n";
            $self->chr($args{-chr_name});
        } else {
            carp "Must provide a chromosome name\n";
            return;
        }
    }

    my $seq = $self->seq;
    if ($seq) {
        if (   !$seq->isa("Bio::SeqI")
            && !$seq->isa("Bio::PrimarySeqI"))
        {
            carp "Seq is not a Bio::SeqI or Bio::PrimarySeqI compliant"
                . " object\n";
            return;
        }

        if ($seq->isa('Bio::LocatableSeq')) {
            if (defined $self->start) {
                if ($self->start != $seq->start) {
                    carp "Sequence start does not match start argument"
                        . " provided\n";
                    return;
                }
            } else {
                $self->start($seq->start);
            }

            if (defined $self->end) {
                if ($self->end != $seq->end) {
                    carp "Sequence end does not match end argument"
                        . " provided\n";
                    return;
                }
            } else {
                $self->end($seq->end);
            }
        } else {
            # using relative sequence coords
            $self->start(1) if !$self->start;
            $self->end($seq->length()) if !$self->end;
        }
    }

    if (!defined $self->start) {
        carp "Must provide a start argument\n";
        return;
    }

    if (!defined $self->end) {
        carp "Must provide an end argument\n";
        return;
    }

    if (defined $self->exons) {
        if (   ref $self->exons ne "ARRAY"
            || !$self->exons->[0]->can("start")
            || !$self->exons->[0]->can("end"))
        {
            carp
                "Exons is not a list ref of objects with start/end methods\n";
            return;
        }
    }

    return $self;
}

################################################################################
#
# Getter/Setter methods follow.
#
################################################################################

=head2 seq

 Title    : seq
 Usage    : $seq = $phca->seq($seq);
 Function : Get/Set the base sequence.
 Returns  : A Bio::Seq object.
 Args     : Optionally a new base sequence (Bio::Seq object)

=cut

sub seq
{
    my ($self, $seq) = @_;

    if (defined($seq)) {
        if (   !$seq->isa("Bio::SeqI")
            && !$seq->isa("Bio::PrimarySeqI"))
        {
            carp "Seq is not a Bio::SeqI or Bio::PrimarySeqI compliant"
                . " object\n";
            return undef;
        }
        $self->{-seq} = $seq;
    }

    return $self->{-seq};
}

=head2 db

 Title    : db
 Usage    : $db = $phca->db($db);
 Function : Get/Set the UCSC database name.
 Returns  : The UCSC database name.
 Args     : Optionally a new UCSC database name.

=cut

sub db
{
    my ($self, $db) = @_;

    if (defined($db)) {
        $self->{-db} = $db;
    }

    return $self->{-db};
}

=head2 track

 Title    : track
 Usage    : $track = $phca->track($track);
 Function : Get/Set the UCSC track name.
 Returns  : The UCSC track name.
 Args     : Optionally a new UCSC track name.

=cut

sub track
{
    my ($self, $track) = @_;

    if (defined($track)) {
        $self->{-track} = $track;
    }

    return $self->{-track};
}

=head2 chr

 Title    : chr
 Usage    : $chr = $phca->chr($chr);
 Function : Get/Set the UCSC chromosome name.
 Returns  : The chromosome name.
 Args     : Optionally a new chromosome name.

=cut

sub chr
{
    my ($self, $chr) = @_;

    if (defined($chr)) {
        $self->{-chr} = $chr;
    }

    return $self->{-chr};
}

=head2 chr_name

 Title    : chr_name
 Usage    : Deprecated method; use chr() instead.

=cut

sub chr_name
{
    my ($self, $chr) = @_;

    carp "Deprecated method chr_name() - please use chr() method instead\n";

    return $self->chr($chr);
}

=head2 start

 Title    : start
 Usage    : $start = $phca->start($start);
 Function : Get/Set the start coordinate on the chromosome.
 Returns  : The start coordinate on the chromosome.
 Args     : Optionally a new start coordinate on the chromosome.

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
 Usage    : $end = $phca->end($end);
 Function : Get/Set the end coordinate on the chromosome.
 Returns  : The end coordinate on the chromosome.
 Args     : Optionally a new end coordinate on the chromosome.

=cut

sub end
{
    my ($self, $end) = @_;

    if (defined($end)) {
        $self->{-end} = $end;
    }

    return $self->{-end};
}

=head2 exons

 Title    : exons
 Usage    : $exons = $phca->exons();
 Function : Get/Set the sequence exons.
 Returns  : A reference to a list of Bio::SeqFeatureI objects.
 Args     : Optionally new sequence exons (reference to a
            list of Bio::SeqFeatureI objects)

=cut

sub exons
{
    my ($self, $exons) = @_;

    if (defined $exons) {
        if (ref $exons ne "ARRAY" || !$exons->[0]->isa("Bio::SeqFeatureI")) {
            carp "Exons is not a list ref of Bio::SeqFeatureI compliant"
                . " objects\n";
            return undef;
        }
        $self->{-exons} = $exons;
    }

    return $self->{-exons};
}

=head2 min_conservation

 Title    : min_conservation
 Usage    : $min_conservation = $phca->min_conservation($min_conservation);
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

=head2 min_conserved_region_length

 Title    : min_conserved_region_length
 Usage    : $min_len = $phca->min_conserved_region_length($min_len);
 Function : Get/set the minimum conserved region length.
 Returns  : 
 Args     : The minimum conserved region length to report.

=cut

sub min_conserved_region_length
{
    my ($self, $len) = @_;

    if (defined $len) {
        $self->{-min_length} = $len;
    }

    return $self->{-min_length};
}

=head2 filter_exons

 Title    : filter_exons
 Usage    : $filter_exons = $phca->filter_exons($filter_exons);
 Function : Get/set whether to filter exons from the conserved regions.
 Returns  : Boolean value
 Args     : Boolean value

=cut

sub filter_exons
{
    my ($self, $bool) = @_;

    if (defined $bool) {
        $self->{-filter_exons} = $bool;
    }

    return $self->{-filter_exons};
}

=head2 flank_size

 Title    : flank_size
 Usage    : $size = $phca->flank_size($size);
 Function : Get/set the amount of flanking sequence added to the conserved
            regions.
 Returns  : 
 Args     : The flank size used in computing the conserved regions.

=cut

sub flank_size
{
    my ($self, $len) = @_;

    if (defined $len) {
        $self->{-flank_size} = $len;
    }

    return $self->{-flank_size};
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
 Usage    : $profile = $phca->conservation_profile_array();
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
 Usage    : $profile = $phca->conservation_profile_hash();
 Function : Get the conservation profile as a hash.
 Returns  : A hashref of phastCons scores where the hash key is the position
            (in the range start..end) and the hash value is the phastCons
            score at that position.
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
 Usage    : $regions = $phca->conserved_regions()
            OR $phca->conserved_regions($regions);
 Function : Get/set the conserved regions.
 Returns  : A reference to a list of Bio::SeqFeature::FeaturePair objects
            or undef.
 Args     : Optionally, a listref of Bio::SeqFeature::FeaturePair objects.

=cut

sub conserved_regions
{
    my ($self, $regions) = @_;

    if ($regions) {
        if (ref $regions ne "ARRAY"
            || !$regions->[0]->isa("Bio::SeqFeature::Generic"))
        {
            carp "Conserved regions argument is not a listref of"
                . " Bio::SeqFeature::Generic objects\n";
            return undef;
        }
        $self->{-conserved_regions} = $regions;
    # XXX If not computed, do not compute. Causes problems if no conservation
    # parameters set
    #} else {
    #    if (!$self->{-conserved_regions} && $self->min_conservation()) {
    #        $self->compute_conserved_regions();
    #    }
    }

    return $self->{-conserved_regions};
}

=head2 conserved_subsequences

 Title    : conserved_subsequences
 Usage    : $subseqs = $phca->conserved_subsequences();
 Function : Get the list of conserved sub-sequences.
 Returns  : A reference to a list of Bio::Seq objects.
 Args     : None.

=cut

sub conserved_subsequences
{
    my $self = shift;

    if (!$self->{-conserved_subseqs} && $self->{-conserved_regions}) {
        $self->compute_conserved_subsequences();
    }

    return $self->{-conserved_subseqs};
}

=head2 matrix_set

 Title    : matrix_set
 Usage    : $matrix_set = $phca->matrix_set($matrix_set);
 Function : Get/Set the TFBS::MatrixSet.
 Returns  : A TFBS::MatrixSet object.
 Args     : Optionally a new TFBS::MatrisSet object.

=cut

sub matrix_set
{
    my ($self, $matrix_set) = @_;

    if (defined($matrix_set)) {
        if (!$matrix_set->isa("TFBS::MatrixSet")) {
            carp "Matrix_set is not a TFBS::MatrixSet object\n";
            return undef;
        }
        $self->{-matrix_set} = $matrix_set;
    }

    return $self->{-matrix_set};
}

=head2 tfbs_threshold

 Title    : tfbs_threshold
 Usage    : Deprecated method - use min_tfbs_score instead.

=cut

sub tfbs_threshold
{
    my ($self, $score) = @_;

    carp "Deprecated method tfbs_threshold() - please use min_tfbs_score()"
        . " instead\n";

    return $self->min_tfbs_score($score);
}

=head2 min_tfbs_score

 Title    : min_tfbs_score
 Usage    : $score = $phca->min_tfbs_score($score);
 Function : Get/Set the minimum TFBS score.
 Returns  : The minimum TFBS score.
 Args     : An absolute score as a float or a percentage value string in
            the range '0%' to '100%'.

=cut

sub min_tfbs_score
{
    my ($self, $score) = @_;

    if (defined($score)) {
        $self->{-min_tfbs_score} = $score;
    }

    return $self->{-min_tfbs_score};
}

=head2 min_tfbs_cr_overlap

 Title    : min_tfbs_cr_overlap
 Usage    : $overlap = $phca->min_tfbs_cr_overlap($overlap);
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

=head2 filter_overlapping_tfbss

 Title    : filter_overlapping_tfbss
 Usage    : $filter = $phca->filter_overlapping_tfbss($bool);
 Function : Get/set whether to filter overlapping TFBSs so that only the
            highest scoring binding site of a group of sites for a given
            TF is reported by compute_conserved_tfbss().
 Returns  : Boolean value
 Args     : Boolean value

=cut

sub filter_overlapping_tfbss
{
    my ($self, $bool) = @_;

    if (defined $bool) {
        $self->{-filter_overlapping_tfbss} = $bool;
    }

    return $self->{-filter_overlapping_tfbss};
}

=head2 tfbs_search_start

 Title    : tfbs_search_start
 Usage    : $start = $phca->tfbs_search_start($start);
 Function : Get/Set the start coordinate of the TFBS search region.
 Returns  : The start coordinate of the TFBS search region.
 Args     : Optionally a new start coordinate of the TFBS search region.

=cut

sub tfbs_search_start
{
    my ($self, $start) = @_;

    if (defined($start)) {
        $self->{-tfbs_search_start} = $start;
    }

    return $self->{-tfbs_search_start};
}

=head2 tfbs_search_end

 Title    : tfbs_search_end
 Usage    : $end = $phca->tfbs_search_end($end);
 Function : Get/Set the end coordinate of the TFBS search region.
 Returns  : The end coordinate of the TFBS search region.
 Args     : Optionally a new end coordinate of the TFBS search region.

=cut

sub tfbs_search_end
{
    my ($self, $end) = @_;

    if (defined($end)) {
        $self->{-tfbs_search_end} = $end;
    }

    return $self->{-tfbs_search_end};
}

=head2 conserved_tf_sites

 Title    : conserved_tf_sites
 Usage    : Deprecated - use conserved_tfbss() method.

=cut

sub conserved_tf_sites
{
    my $self = shift;

    carp "Deprecated method conserved_tf_sites() - please use"
        . " conserved_tfbss() method instead\n";

    return $self->conserved_tfbss;
}

=head2 conserved_tfbss

 Title    : conserved_tfbss
 Usage    : $sites = $phca->conserved_tfbss();
 Function : Get the conserved TFBSs.
 Returns  : A TFBS::SiteSet reference.
 Args     : None.

=cut

sub conserved_tfbss
{
    my $self = shift;

    # XXX Do not compute - causes problems if parameters are not already set
    #if (!$self->{-conserved_tfbss}) {
    #    $self->compute_conserved_tfbss();
    #}

    return $self->{-conserved_tfbss};
}

################################################################################
#
# Analysis methods follow.
#
################################################################################

=head2 fetch_conservation_profile

 Title    : fetch_conservation_profile
 Usage    : Deprecated method - use compute_conservation_profile() instead

=cut

sub fetch_conservation_profile
{
    my ($self, %args) = @_;

    carp "Deprecated method fetch_conservation_profile() - please use"
        . " compute_conservation_profile() instead\n";

    return $self->compute_conservation_profile(%args);
}

=head2 compute_conservation_profile

 Title    : compute_conservation_profile
 Usage    : $profile = $phca->compute_conservation_profile();
 Function : Fetch the phastCons profile from the appropriate UCSC
            db/flat file resources.
 Returns  : A reference to an array of hashes containing '-position' and
            '-score' key/value pairs.
 Args     : None.

=cut

sub compute_conservation_profile
{
    my ($self, %args) = @_;

    #my $seq = $self->seq;
    #if (!$seq) {
    #    carp "Sequence not set\n";
    #    return;
    #}

    my $db = $self->db;
    if (!$db) {
        carp "UCSC database name not set\n";
        return;
    }

    my $track = $self->track;
    if (!$track) {
        carp "UCSC track name not set\n";
        return;
    }

    my $chr = $self->chr;
    if (!$chr) {
        carp "chromosome name not specified\n";
        return;
    }

    if ($chr !~ /^chr/) {
        $chr = "chr$chr";
    }

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

    my $hgw = ORCA::Analysis::Run::hgWiggle->new();
    if (!$hgw) {
        carp "Error initializing ORCA::Analysis::Run::hgWiggle\n";
        return;
    }

    my $out = $hgw->run(
        {
            -db       => $db,
            -position => "$chr:$start-$end",
            # change coords from 0-based to 1-based
            -lift => 1
        }, $track
    );

    my @lines = split "\n", $out;

    #
    # Store conservation scores as position/value pairs. Positions will
    # be stored in the start..end range, i.e. chromosomal coordinates.
    # May be non-continuous data. DJA 09/12/11
    #
    my %profile_hash;
    foreach my $line (@lines) {
        if ($line =~ /^\s*(\d+)\s+(\d+\.*\d*)/) {
            $profile_hash{$1} = $2;
        }
    }
    $self->{_profile_hash} = \%profile_hash;

    #
    # Generate a profile as continuous array data (0-based coords). Fill
    # gaps with 0 values. DJA 09/12/11
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
 Usage    : $regions = $phca->compute_conserved_regions();
 Function : Compute and return the list of conserved regions.
 Returns  : A reference to a list of Bio::SeqFeature::Generic objects.
 Args     : Optionally:
            -min_conservation   => the minimum percent identity of
                                   regions to report (specify as a
                                   number between 0 and 1 or a string
                                   such as '70%')
            -filter_exons       => (bool) indicates exons should be
                                   filtered out of conserved regions.
            -min_length         => (int) minimum length of conserved
                                   region to keep
            -flank_size         => (int) if specified, add this much
                                   flanking sequence to either side of the
                                   conserved regions

=cut

sub compute_conserved_regions
{
    my ($self, %args) = @_;

    my $min_conservation = $args{-min_conservation};
    if (defined $min_conservation) {
        if ($min_conservation =~ /(.+)%/) {
            $min_conservation = $1 / 100;
        }

        if ($min_conservation <= 0 || $min_conservation > 1) {
            carp
                "min_conservation is out of range; please specify as a number"
                . " between 0 and 1 or as a string between 0% and 100%";
            return;
        }

        $self->min_conservation($min_conservation);
    }

    $min_conservation = $self->min_conservation();

    if (!defined $min_conservation) {
        carp "No minimum conservation defined\n";
        return;
    }

    $self->filter_exons($args{-filter_exons}) if defined $args{-filter_exons};
    my $filter_exons = $self->filter_exons();

    my $min_cr_len = 0;
    if (defined $args{-min_length}) {
        $min_cr_len = $self->min_conserved_region_length($args{-min_length});
    } elsif (defined $args{-min_cr_len}) {
        carp "Deprecated argument -min_cr_len; please use -min_length\n";
        $min_cr_len = $self->min_conserved_region_length($args{-min_cr_len})
    }

    my $flank_size = 0;
    if (defined $args{-flank_size}) {
        $flank_size = $self->flank_size($args{-flank_size});
    }

    my $exons = $self->exons;

    @$exons = sort {$a->start <=> $b->start} @$exons;

    my @conservation;
    if ($self->conservation_profile) {
        @conservation = @{$self->conservation_profile};
    } else {
        $self->{-conserved_regions} = undef;
        return;
    }

    my $start     = $self->start;
    my $end       = $self->end;

    my $rel_start = 1;
    my $rel_end   = $end - $start + 1;

    if ($filter_exons && $exons) {
        foreach my $exon (@$exons) {
            #
            # XXX Assume exon coords are relative to start of sequence
            #
            #my $mask_start = $exon->start + $start - 1;
            #my $mask_end   = $exon->end + $start - 1;
            #
            # Otherwise if exon coords are in same coord system as sequence
            #
            my $mask_start = $exon->start;
            my $mask_end   = $exon->end;

            next if ($mask_end < $start || $mask_start > $end);
 
            $mask_start = $start if $mask_start < $start;
            $mask_end   = $end   if $mask_end > $end;

            # Convert to 0-based array index coords
            my $idx_start = $mask_start - $start;
            my $idx_end   = $mask_end - $start;

            # Zero out exon positions in the profile
            map {$conservation[$_] = 0} ($idx_start .. $idx_end);
        }
    }

    my $cr_start_idx = undef;
    my $cr_end_idx   = undef;
    my $counter      = 0;
    my @crs;
    foreach my $i (0 .. $#conservation - 1) {
        if ($conservation[$i] >= $min_conservation) {
            if (!defined $cr_start_idx) {
                $cr_start_idx = $i;
            }
        } else {
            if (defined $cr_start_idx) {
                $cr_end_idx = $i - 1;

                my $score = _score_region(
                    \@conservation, $cr_start_idx, $cr_end_idx
                );

                my $cr = Bio::SeqFeature::Generic->new(
                    -primary_id   => sprintf("CR%d", $counter + 1),
                    -display_name => sprintf("CR%d", $counter + 1),
                    -source_tag   => "ORCA",
                    # Leave in 0-based index coords until after
                    # merging regions below
                    -start => $cr_start_idx,
                    -end   => $cr_end_idx,
                    -score => $score
                );

                push @crs, $cr;

                $cr_start_idx = undef;
                $cr_end_idx   = undef;
                $counter++;
            }
        }
    }

    if (defined $cr_start_idx) {
        # Started a region and didn't finish it (scores remained above
        # min. conservation until end of sequence)
        $cr_end_idx = $#conservation;

        my $score = _score_region(\@conservation, $cr_start_idx, $cr_end_idx);

        my $cr = Bio::SeqFeature::Generic->new(
            -primary_id   => sprintf("CR%d", $counter + 1),
            -display_name => sprintf("CR%d", $counter + 1),
            -source_tag   => "ORCA",
            # Leave in 0-based index coords until after
            # merging regions below
            -start => $cr_start_idx,
            -end   => $cr_end_idx,
            -score => $score
        );

        push @crs, $cr;
    }

    #
    # Combine regions into larger regions which still score
    # above min. conservation
    #
    if ($filter_exons && $exons) {
        while (my $combined_crs =
            _combine_conserved_regions_exluding_exons(
                \@conservation, \@crs, $exons, $min_conservation, $min_cr_len,
                $self->start
            )
        ) {
            @crs = ();
            @crs = @{$combined_crs};
        }
    } else {
        while (my $combined_crs =
            _combine_conserved_regions(
                \@conservation, \@crs, $min_conservation, $min_cr_len
            )
        ) {
            @crs = ();
            @crs = @{$combined_crs};
        }
    }

    #
    # Explicitly remove exons from conserved regions (above combining of
    # regions can result in conserved region spanning an exon position)
    # NOTE: this may result in regions which now score below min. conservation
    #
    #if ($filter_exons && $exons && @crs) {
    #    my $cut_crs = _cut_exons_from_conserved_regions(
    #        \@conservation, \@crs, $exons, $self->start
    #    );
    #
    #    if ($cut_crs) {
    #        @crs = @$cut_crs;
    #    } else {
    #        @crs = undef;
    #    }
    #}

    #
    # Keep only regions > min. region length and which still score above min.
    # conservation.
    #
    my @signif_crs;
    foreach my $cr (@crs) {
        if (
            ($cr->end - $cr->start + 1) > $min_cr_len
            && $cr->score >= $min_conservation
        )
        {
            push @signif_crs, $cr;
        }
    }

    @crs = ();
    @crs = @signif_crs;

    if ($flank_size) {
        my $flanked_crs = _add_conserved_region_flanks(
            \@conservation, \@crs, $exons, $self->start, $flank_size
        );

        if ($flanked_crs) {
            @crs = @$flanked_crs;
        } else {
            @crs = undef;
        }

        if (@crs) {
            while (my $combined_crs = _combine_flanked_conserved_regions(\@crs))            {
                @crs = ();
                @crs = @{$combined_crs};
            }
        }
    }

    #
    # XXX Change conserved region coordinates from 0-based (array index) coords
    # back relative (1-based) coords DJA 09/12/13
    #
    foreach my $cr (@crs) {
        $cr->start($cr->start + 1);
        $cr->end($cr->end + 1);
    }
    #
    # Otherwise if using chromosomal coords
    #
    #foreach my $cr (@crs) {
    #    $cr->start($cr->start + $start);
    #    $cr->end($cr->end + $start);
    #}

    $self->{-conserved_regions} = @crs ? \@crs : undef;
}

=head2 extract_conserved_subsequences

 Title    : extract_conserved_subsequences
 Usage    : Deprecated method - use compute_conserved_subsequences.

=cut

sub extract_conserved_subsequences
{
    my ($self) = @_;

    carp "Deprecated method extract_conserved_subsequences()"
        . " - use compute_conserved_subsequences() instead\n";

    $self->compute_conserved_subsequences();
}

=head2 compute_conserved_subsequences

 Title    : compute_conserved_subsequences
 Usage    : $subseqs = $phca->compute_conserved_subsequences($seq_num);
 Function : Return the sub-sequences corresponding to the conserved regions.
 Returns  : A reference to a list of Bio::Seq objects.
 Args     : Optionally specify a sequence number; default = 1.

=cut

sub compute_conserved_subsequences
{
    my ($self) = @_;

    my $seq = $self->seq;

    return if !defined $seq;

    my $crs = $self->conserved_regions;
    if (!defined $crs) {
        carp "no conserved regions defining sub-sequences to extract\n";
        return undef;
    }

    my $start = $self->start;
    my $end   = $self->end;

    my $seq_display_id = $seq->display_id;

    my @subseqs;
    foreach my $cr (@$crs) {
        # XXX If using relative coords
        my $ss_start = $cr->start;
        my $ss_end   = $cr->end;
        # Otherwise if using chromosomal coords, convert to relative coords
        #my $ss_start = $cr->start - $start + 1;
        #my $ss_end   = $cr->end - $start + 1;

        # NOTE: the trunc method expects relative coords even if $seq is a
        # Bio::LocatableSeq
        my $subseq = $seq->trunc($ss_start, $ss_end);

        $seq_display_id = "SUBSEQ" if !$seq_display_id;
        my $subseq_display_id = $seq_display_id . '/' . "$ss_start-$ss_end";
        $subseq->display_id($subseq_display_id);

        push @subseqs, $subseq;
    }

    $self->{-conserved_subseqs} = @subseqs ? \@subseqs : undef;
}

=head2 find_conserved_tf_sites

 Title    : find_conserved_tf_sites
 Usage    : Deprecated - use compute_conserved_tfbss

=cut

sub find_conserved_tf_sites
{
    my ($self, %args) = @_;

    carp "Deprecated method find_conserved_tf_sites() - please use"
        . " compute_conserved_tfbss() method instead\n";

    return $self->compute_conserved_tfbss(%args);
}

=head2 compute_conserved_tfbss

 Title    : compute_conserved_tfbss
 Usage    : $sites = $phca->compute_conserved_tfbss(
                -matrix_set                 => $matrix_set,
                -min_tfbs_score             => $tfbs_threshold,
                -min_tfbs_cr_overlap        => 1
                -filter_overlapping_tfbss   => 1,
                -start                      => $start,
                -end                        => $end
            );
 Function : Compute and return the list of conserved TF sites.
 Returns  : A reference to a TFBS::SiteSet object.
 Args     : -matrix_set                 => a TFBS::MatrixSet object 
            -min_tfbs_score             => Report only TFBSs with at least
                                           this score
            -min_tfbs_cr_overlap        => Report only TFBSs which overlap
                                           a conserved region by at least
                                           this number of nucleotides.
            -filter_overlapping_tfbss   => boolean; if true, for any given
                                           matrix, report only the best
                                           scoring of any two overlapping
                                           site pairs
            -start                      => restrict TFBS search region to
                                           start at position $start
            -end                        => restrict TFBS search region to
                                           end at position $end

=cut

sub compute_conserved_tfbss
{
    my ($self, %args) = @_;

    my $seq = $self->seq();
    if (!defined $seq) {
        carp "PhastCons analysis object has no sequence defined\n";
        return;
    }

    my $start = $self->start();
    if (!defined $start) {
        carp "PhastCons analysis object has no start defined\n";
        return;
    }

    my $end = $self->end();
    if (!defined $start) {
        carp "PhastCons analysis object has no end defined\n";
        return;
    }

    my $conserved_regions = $self->conserved_regions;
    if (!defined $conserved_regions) {
        carp "PhastCons analysis object has no conserved regions defined\n";
        return;
    }

    $self->matrix_set($args{-matrix_set}) if defined $args{-matrix_set};
    my $matrix_set = $self->matrix_set();
    if (!defined $matrix_set || !$matrix_set->isa("TFBS::MatrixSet")) {
        carp "Matrix set not provided or is not a TFBS::MatrixSet\n";
        return;
    }

    if ($matrix_set->size == 0) {
        carp "Matrix set is empty\n";
        return;
    }

    if (defined $args{-filter_overlapping_tfbss}) {
        $self->filter_overlapping_tfbss($args{-filter_overlapping_tfbss});
    } elsif (defined $args{-filter_overlapping_sites}) {
        carp "Deprecated argument -filter_overlapping_sites - please use"
            . " -filter_overlapping_tfbss instead\n";
        $self->filter_overlapping_tfbss($args{-filter_overlapping_sites});
    }
    my $filter_overlaps = $self->filter_overlapping_tfbss || 0;

    my $search_start    = $args{-start};
    my $search_end      = $args{-end};

    #
    # XXX if using relative coords
    #
    $search_start   = 1 if !defined $search_start;
    $search_end     = $seq->length() if !defined $search_end;
 
    if ($search_end < 1 || $search_start > $seq->length()) {
        carp "TF site search start/end outside of sequence start/end"
            . " - make sure start/end is specified in relative coordinates\n";
        return;
    }

    #
    # Otherwise if using chromosomal coords
    #
    #$search_start   = $start if !defined $search_start;
    #$search_end     = $end if !defined $search_end;
    #
    #if ($search_end < $start || $search_start > $end) {
    #    carp "TF site search start/end outside of sequence start/end"
    #        . " - make sure coordinate systems are consistent\n";
    #    return;
    #}

    $self->tfbs_search_start($search_start);
    $self->tfbs_search_end($search_end);

    $self->min_tfbs_cr_overlap($args{-min_tfbs_cr_overlap})
        if defined $args{-min_tfbs_cr_overlap};

    my $overlap = $self->min_tfbs_cr_overlap();
    $overlap = DFLT_MIN_TFBS_CA_OVERLAP if !defined $overlap;

    if (defined $args{-min_tfbs_score}) {
        $self->min_tfbs_score($args{-min_tfbs_score});
    } elsif (defined $args{-tfbs_threshold}) {
        carp "Deprecated argument -tfbs_threshold - please use"
            . " -min_tfbs_score instead\n";
        $self->min_tfbs_score($args{-tfbs_threshold});
    }

    my $threshold = $self->min_tfbs_score();
    if (!defined $threshold) {
        carp "Min. TFBS score not provided\n";
        return;
    }

    my $cons_tf_count = 0;

    my $site_set = TFBS::SiteSet->new();

    my $ms_iter = $matrix_set->Iterator();
    while (my $matrix = $ms_iter->next) {
        if (DEBUG) {
            printf "\nTranscription Factor %s  %s\n",
                $matrix->name, $matrix->class;
        }

        my $matrix_width = $matrix->length;

        my $m_site_set = TFBS::SiteSet->new();
        foreach my $cr (@$conserved_regions) {
            my $cr_start = $cr->start;
            my $cr_end   = $cr->end;

            #
            # If conserved region falls outside the sequence region of
            # interest exclude it from analysis
            #
            next if $cr_start > $search_end;
            next if $cr_end < $search_start;

            my $m_start = $cr_start - $matrix_width + $overlap;
            $m_start = $search_start if $m_start < $search_start;

            my $m_end = $cr_end + $matrix_width - $overlap;
            $m_end = $search_end if $m_end > $search_end;

            next if $m_end - $m_start + 1 < $matrix_width;

            if (DEBUG) {
                printf "\nTFBS search region %d  %d  %.2f\n",
                    $cr_start,
                    $cr_end,
                    $cr->score;
            }

            my $cr_m_site_set = _find_seq_tf_site_set(
                $matrix, $seq,
                # XXX if using relative coords
                $m_start, $m_end,
                # Otherwise if using chromomal coords
                #$m_start - $start + 1, $m_end - $start + 1,
                $threshold
            );

            next if !$cr_m_site_set || $cr_m_site_set->size == 0;
            if (DEBUG) {
                _print_site_set($cr_m_site_set) if defined $cr_m_site_set;
            }

            _set_site_conservation($cr_m_site_set, $cr->score);

            $m_site_set->add_siteset($cr_m_site_set);
        }

        next if !$m_site_set || $m_site_set->size == 0;
        if ($filter_overlaps) {
            my $filt_m_site_set =
                _tf_filter_overlapping_site_set($m_site_set);
            if (DEBUG) {
                print "\nFiltered site set\n";
                _print_site_set(
                    $filt_m_site_set);
            }
            $site_set->add_siteset($filt_m_site_set);
        } else {
            $site_set->add_siteset($m_site_set);
        }
    }

    my @tfbss;
    my $iter = $site_set->Iterator(-sort_by => 'start');
    while (my $site = $iter->next()) {
        #
        # XXX if using chromosomal coords but...
        # Unfortunately converting to chromosomal coords causes problems
        # because if you later call $site->seq(), Bioperl throws an exception
        # if the end coordinate of the site is greater than the total sequence
        # length - it doesn't seem to handle LocatableSeqs!!!
        #
        #$site->start($site->start + $start - 1);
        #$site->end($site->end + $start - 1);
        #

        push @tfbss, $site;
    }

    #
    # Changed to return a listref of TFBS::Site objects rather than a
    # TFBS::SiteSet object. DJA 10/01/07
    #
    $self->{-conserved_tfbss} = @tfbss ? \@tfbss : undef;
}

=head2 ucsc_track

 Title    : ucsc_track
 Usage    : $track = $phca->ucsc_track();
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

        foreach my $tf (@$tfbss) {
            $track .= sprintf "%s\t%d\t%d\t%s\t%.3f\n",
                $chr,
                # TFBSs are in relative coords - convert to chromosomal
                $tf->start + $start - 1,
                $tf->end + $start - 1,
                $tf->pattern->name,
                $tf->rel_score;
        }
    }

    my $cp = $self->conservation_profile;
    if ($cp) {
        $track .= sprintf "track type=wiggle_0 name=\"Conservation\" description=\"ORCA Conservation (%s)\" color=255,0,0 graphType=bar yLineMark=0.0 yLineOnOff=on visibility=2\n", $self->track() || 'phastCons';

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
# Internal methods follow.
#
################################################################################

sub _find_seq_tf_site_set
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
# XXX This may have to be revisited for more sophisticated filtering.
# Take a TFBS::SiteSet and filter overlapping sites such that only
# the highest scoring site of any mutually overlapping sites is kept.
# In the event that sites score equally, the first site is kept, i.e.
# bias is towards the site with the lowest starting position.
#
sub _tf_filter_overlapping_site_set
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
    my ($site_set, $score) = @_;

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

sub _combine_conserved_regions
{
    my ($conservation, $crs, $min_conservation, $min_len) = @_;

    my @crs2;
    my $combined = 0;
    my $i        = 0;
    while (defined $crs->[$i + 1]) {
        # 0-based index coords
        my $rstart = $crs->[$i]->start;
        my $rend   = $crs->[$i + 1]->end;

        my $score = _score_region($conservation, $rstart, $rend);

        if ($score >= $min_conservation) {
            #printf
            #  "combining region %s (%d-%d) and %s (%d-%d); score = %.3f\n",
            #        $crs->[$i]->display_name, $rstart, $crs->[$i]->end,
            #        $crs->[$i+1]->display_name, $crs->[$i+1]->start, $rend,
            #        $score;

            my $cr = Bio::SeqFeature::Generic->new(
                -primary_id   => $crs->[$i]->display_name,
                -display_name => $crs->[$i]->display_name,
                -source_tag   => "ORCA",
                -start        => $rstart,
                -end          => $rend,
                -score        => $score
            );
            push @crs2, $cr;
            $combined = 1;
            $i++;
        } else {
            #printf "NOT combining region %s (%d-%d) and %s (%d-%d);"
            #        . " score = %.3f\n",
            #        $crs->[$i]->display_name, $rstart, $crs->[$i]->end,
            #        $crs->[$i+1]->display_name, $crs->[$i+1]->start, $rend,
            #        $score;

            push @crs2, $crs->[$i];
        }
        $i++;
    }

    if (defined $crs->[$i]) {
        push @crs2, $crs->[$i];
    }

    return $combined ? \@crs2 : undef;
}

sub _combine_conserved_regions_exluding_exons
{
    my ($conservation, $crs, $exons, $min_conservation, $min_len, $start) = @_;

    my @crs2;
    my $combined = 0;
    my $i        = 0;
    while (defined $crs->[$i + 1]) {
        # 0-based index coords
        my $lstart = $crs->[$i]->start;
        my $lend   = $crs->[$i]->end;

        my $rstart = $crs->[$i + 1]->start;
        my $rend   = $crs->[$i + 1]->end;

        my $exon_between = 0;
        foreach my $exon (@$exons) {
            #
            # Exons are in seq region coords but conserved regions are in
            # 0-based coords
            #
            my $xstart = $exon->start - $start;
            my $xend   = $exon->end - $start;

            last if $xstart > $rstart;
            next if $xend < $lend;

            if ($xend >= $lend && $xstart <= $rstart) {
                $exon_between = 1;
            }
        }

        if ($exon_between) {
            #
            # If there is an exon between the two conserved regions,
            # don't combine them.
            #
            push @crs2, $crs->[$i];
        } else {
            my $score = _score_region($conservation, $lstart, $rend);

            if ($score >= $min_conservation) {
                #printf
                #  "combining region %s (%d-%d) and %s (%d-%d); score = %.3f\n",
                #        $crs->[$i]->display_name, $rstart, $crs->[$i]->end,
                #        $crs->[$i+1]->display_name, $crs->[$i+1]->start, $rend,
                #        $score;

                my $cr = Bio::SeqFeature::Generic->new(
                    -primary_id   => $crs->[$i]->display_name,
                    -display_name => $crs->[$i]->display_name,
                    -source_tag   => "ORCA",
                    -start        => $lstart,
                    -end          => $rend,
                    -score        => $score
                );
                push @crs2, $cr;
                $combined = 1;
                $i++;
            } else {
                #printf "NOT combining region %s (%d-%d) and %s (%d-%d);"
                #        . " score = %.3f\n",
                #        $crs->[$i]->display_name, $rstart, $crs->[$i]->end,
                #        $crs->[$i+1]->display_name, $crs->[$i+1]->start, $rend,
                #        $score;

                push @crs2, $crs->[$i];
            }
        }
        $i++;
    }

    if (defined $crs->[$i]) {
        push @crs2, $crs->[$i];
    }

    return $combined ? \@crs2 : undef;
}

sub _add_conserved_region_flanks
{
    my ($conservation, $crs, $exons, $start, $flank_size) = @_;

    my $last_idx = $#{$conservation};

    # assumes sorted in calling method
    #@$exons = sort {$a->start <=> $b->start} @$exons;

    foreach my $cr (@$crs) {
        my $cr_start = $cr->start - $flank_size;
        $cr_start = 0 if $cr_start < 0;

        my $cr_end = $cr->end + $flank_size;
        $cr_end = $last_idx if $cr_end > $last_idx;

        foreach my $ex (@$exons) {
            #
            # Exons are in seq region coords but conserved regions are in
            # 0-based coords
            #
            my $ex_start = $ex->start - $start;
            my $ex_end   = $ex->end - $start;

            next if $ex_end < $cr_start;
            last if $ex_start > $cr_end;

            if ($ex_end <= $cr->end && $ex_end > $cr_start) {
                $cr_start = $ex_end + 1;
            }

            if ($ex_start >= $cr->start && $ex_start < $cr_end) {
                $cr_end = $ex_start - 1;
            }
        }

        $cr->start($cr_start);
        $cr->end($cr_end);
    }

    return $crs;
}

sub _combine_flanked_conserved_regions
{
    my ($crs) = @_;

    my @crs2;
    my $combined = 0;
    my $i        = 0;
    while (defined $crs->[$i + 1]) {
        # 0-based index coords
        my $lstart = $crs->[$i]->start;
        my $lend   = $crs->[$i]->end;
        my $rstart = $crs->[$i + 1]->start;
        my $rend   = $crs->[$i + 1]->end;

        if ($lend >= $rstart && $lstart <= $lend) {
            #
            # Score of the combined region takes the lower score of the two
            # regions being combined
            #
            my $score = $crs->[$i]->score;
            if ($crs->[$i + 1]->score < $score) {
                $score = $crs->[$i + 1]->score;
            }

            my $cr = Bio::SeqFeature::Generic->new(
                -primary_id   => $crs->[$i]->display_name,
                -display_name => $crs->[$i]->display_name,
                -source_tag   => "ORCA",
                -start        => $lstart,
                -end          => $rend,
                -score        => $score
            );

            push @crs2, $cr;
            $combined = 1;
            $i++;
        } else {
            push @crs2, $crs->[$i];
        }
        $i++;
    }

    if (defined $crs->[$i]) {
        push @crs2, $crs->[$i];
    }

    return $combined ? \@crs2 : undef;
}

sub _cut_exons_from_conserved_regions
{
    my ($conservation, $crs, $exons, $start) = @_;

    my @cut_crs;

    foreach my $cr (@$crs) {
        my $cr_start = $cr->start;
        my $cr_end   = $cr->end;

        foreach my $ex (@$exons) {
            #
            # Exons are in seq region coords but conserved regions are in
            # 0-based coords
            #
            my $ex_start = $ex->start - $start;
            my $ex_end   = $ex->end - $start;

            last if $ex_start > $cr_end;
            next if $ex_end < $cr_start;

            if ($ex_start <= $cr_start && $ex_end >= $cr_start) { 
                $cr_start = $ex_end + 1 if $ex_end + 1 > $cr_start;
            }

            if ($ex_start <= $cr_end && $ex_end >= $cr_end) { 
                $cr_end = $ex_start - 1 if $ex_start - 1 < $cr_end;
            }
        }

        if ($cr_start <= $cr_end) {
            my $cut_cr = Bio::SeqFeature::Generic->new(
                -start  => $cr_start,
                -end    => $cr_end,
                -score  => _score_region($conservation, $cr_start, $cr_end)
            );
            push @cut_crs, $cut_cr;
        }
    }

    my $split = 1;
    while ($split) {
        $split = 0;

        foreach my $cr (@cut_crs) {
            my $cr_start = $cr->start;
            my $cr_end   = $cr->end;

            my $cr2;
            foreach my $ex (@$exons) {
                #
                # Exons are in seq region coords but conserved regions are
                # in 0-based coords
                #
                my $ex_start = $ex->start - $start;
                my $ex_end   = $ex->end - $start;

                last if $ex_start > $cr_end;
                next if $ex_end < $cr_start;

                if ($ex_start > $cr_start && $ex_end < $cr_end) {
                    $cr2 = Bio::SeqFeature::Generic->new(
                        -start  => $ex_end + 1,
                        -end    => $cr->end,
                        -score  => _score_region(
                            $conservation, $ex_end + 1, $cr->end
                        )
                    );

                    $cr->end($ex_start - 1);
                    $cr->score(
                        _score_region($conservation, $cr->start, $cr->end)
                    );

                    $split = 1;
                    last;
                }
            }

            if ($split) {
                push @cut_crs, $cr2;
                last;
            }
        }
    }

    return @cut_crs ? \@cut_crs : undef;
}

#
# Score conserved region defined by start/end on conservation array.
# NOTE: start and end should be specified in index (0-based) coords
#
sub _score_region
{
    my ($conservation, $start, $end) = @_;

    my $score = 0;
    foreach my $i ($start .. $end) {
        $score += $conservation->[$i];
    }
    $score /= ($end - $start + 1);

    return $score;
}

1;
