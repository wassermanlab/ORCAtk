
=head1 NAME

ORCA::ConservationAnalysis::PhastCons - Object for facilitating phastCons
conservation analysis of a DNA sequence.

=head1 SYNOPSIS

    use ORCA::ConservationAnalysis::PhastCons;

    # Create a new ORCA::ConservationAnalysis::PhastCons object with
    # the given sequence and (UCSC) database and track names.
    # Optionally also provide exons positions to mask.
    $ca = ORCA::ConservationAnalysis::PhastCons->new(
                    -seq        => $seq,
                    -exons      => $exons,
                    -db         => 'hg18',
                    -track      => 'phastCons28way',
                    -chr_name   => '18',
                    -start      => 55033380,
                    -end        => 55048980
    );

    # Retrieve the phastCons score profile.
    my $profile = $ca->fetch_conservation_profile();

    # Set the conserved regions parameters
    $ca->param('min_conservation', '70%');
    $ca->param('filter_exons', 1);
    $ca->param('min_cr_len', 20);

    # where
    min_conservation    = the absolute minimum percent identity of
                          regions to report (specify as a number
                          between 0 and 1 or a string such as
                          '70%')
    filter_exons        = (bool) exons are to be filtered out of
                          the conserved regions (if exons
                          has been set).
    min_cr_len          = (int) min. length of a conserved region
                          to report
    
    # Compute the conserved regions.
    my $regions = $ca->compute_conserved_regions();
    
    # Optionally provide the conserved regions parameters directly to the
    # method.
    my $regions = $ca->compute_conserved_regions(
        -min_conservation   => '70%',
        -filter_exons       => 1,
        -min_cr_len         => 20
    );
    
    
    # The phastCons score profile can be obtained later.
    # (Automatically computed if not already)
    my $profile = $ca->conservation_profile();
    
    # The conserved regions can be obtained later.
    # (Automatically computed if not already)
    my $regions = $ca->conserved_regions();

    # Get the conserved sub-sequences corresponding to the conserved
    # regions. If the regions have been filtered, returns the sub-sequences
    # corresponding to the filtered conserved regions, otherwise returns
    # the sub-sequence corresponding to the unfiltered conserved regions.
    my $subseqs = $ca->extract_conserved_subsequences();
    
    # The conserved sub-sequences can be obtained later
    my $subseqs = $ca->conserved_subsequences();

=head1 DESCRIPTION

ORCA::ConservationAnalysis::PhastCons is an object used for the purpose
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

package ORCA::ConservationAnalysis::PhastCons;

use strict;

use Carp;
use Bio::SeqFeature::Generic;
use ORCA::ConservationAnalysis::Run::hgWiggle;

=head2 new

 Title    : new
 Usage    : $ca = ORCA::ConservationAnalysis::PhastCons->new();
 Function : Create a new ORCA::ConservationAnalysis::PhastCons object.
 Returns  : A new ORCA::ConservationAnalysis::PhastCons object.
 Args     : Optional named parameters:
            -seq        => Bio::Seq object defining the sequence
            -exons      => a reference to a list of Bio::SeqFeature objects
                           defining exons on the sequence
            -db         => name of the UCSC database containing phastCons
             scores
            -track      => name of the UCSC track (table name)  containing
                           the phastCons scores
            -chr_name   => name of chromosome for the sequence provided
            -start      => chromosomal start position of the sequence
            -end        => chromosomal end position of the sequence

=cut

sub new
{
    my ($class, %args) = @_;

    my $self = bless {
        -seq               => undef,
        -exons             => undef,
        -db                => undef,
        -track             => undef,
        -conserved_regions => undef,
        -conserved_subseqs => undef,
        -parameters        => {
            min_conservation => undef,
            filter_exons     => undef,
            min_cr_len       => undef,
            min_cr_gap       => undef
        }, %args
    }, ref $class || $class;

    my $seq = $self->{-seq};
    if ($seq) {
        if (   !$seq->isa("Bio::SeqI")
            && !$seq->isa("Bio::PrimarySeqI"))
        {
            carp "seq is not a Bio::SeqI or Bio::PrimarySeqI compliant"
                . " object\n";
            return;
        }

        if ($seq->isa('Bio::LocatableSeq')) {
            if (defined $self->{-start}) {
                if ($self->{-start} != $seq->start) {
                    carp "sequence start does not match start argument"
                        . " provided\n";
                    return;
                }
            } else {
                $self->{-start} = $seq->start;
            }

            if (defined $self->{-end}) {
                if ($self->{-end} != $seq->end) {
                    carp "sequence end does not match end argument"
                        . " provided\n";
                    return;
                }
            } else {
                $self->{-end} = $seq->end;
            }
        }
    }

    if (!defined $self->{-start}) {
        carp "must provide a start argument\n";
        return;
    }

    if (!defined $self->{-end}) {
        carp "must provide an end argument\n";
        return;
    }

    if (defined $self->{-exons}) {
        if (   ref $self->{-exons} ne "ARRAY"
            || !$self->{-exons}->[0]->can("start")
            || !$self->{-exons}->[0]->can("end"))
        {
            carp
                "exons is not a list ref of objects with start/end methods\n";
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
 Usage    : $seq = $ca->seq($seq);
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
            carp "seq is not a Bio::SeqI or Bio::PrimarySeqI compliant"
                . " object\n";
            return undef;
        }
        $self->{-seq} = $seq;
    }

    return $self->{-seq};
}

=head2 db

 Title    : db
 Usage    : $db = $ca->db($db);
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
 Usage    : $track = $ca->track($track);
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

=head2 chr_name

 Title    : chr_name
 Usage    : $chr_name = $ca->chr_name($chr_name);
 Function : Get/Set the UCSC track name.
 Returns  : The chromosome name.
 Args     : Optionally a new chromosome name.

=cut

sub chr_name
{
    my ($self, $chr_name) = @_;

    if (defined($chr_name)) {
        $self->{-chr_name} = $chr_name;
    }

    return $self->{-chr_name};
}

=head2 start

 Title    : start
 Usage    : $start = $ca->start($start);
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

=head2 start

 Title    : end
 Usage    : $end = $ca->end($end);
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
 Usage    : $exons = $ca->exons();
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
            carp "exons is not a list ref of Bio::SeqFeatureI compliant"
                . " objects\n";
            return undef;
        }
        $self->{-exons} = $exons;
    }

    return $self->{-exons};
}

=head2 param

 Title    : param
 Usage    : $value = $ca->param($param, $value);
 Function : Get/Set the value of a conservation parameter or return a list
            of parameter names.
 Returns  : If a parameter name is provided, the value of the parameter if
            defined. Otherwise a list of defined parameter names.
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

    if (!$self->{_profile_array}) {
        $self->fetch_conservation_profile();
    }

    return $self->{_profile_array};
}

=head2 conserved_regions

 Title    : conserved_regions
 Usage    : $regions = $ca->conserved_regions()
            OR $ca->conserved_regions($regions);
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
            carp "conserved regions argument is not a listref of"
                . " Bio::SeqFeature::Generic objects\n";
            return undef;
        }
        $self->{-conserved_regions} = $regions;
    } else {
        if (!$self->{-conserved_regions}) {
            $self->compute_conserved_regions();
        }
    }

    return $self->{-conserved_regions};
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

################################################################################
#
# Analysis methods follow.
#
################################################################################

=head2 fetch_conservation_profile

 Title    : fetch_conservation_profile
 Usage    : $profile = $ca->fetch_conservation_profile();
 Function : Fetch the phastCons profile from the appropriate UCSC
            db/flat file resources.
 Returns  : A reference to an array of hashes containing '-position' and
            '-score' key/value pairs.
 Args     : None.

=cut

sub fetch_conservation_profile
{
    my ($self, %args) = @_;

    my $seq = $self->seq;
    carp "Error: sequence not set\n" if !$seq;

    my $db = $self->db;
    carp "Error: UCSC database name not set\n" if !$db;

    my $track = $self->track;
    carp "Error: UCSC track name not set\n" if !$track;

    my $chr = $self->chr_name;
    carp "Error: chromosome name not specified\n" if !$chr;
    if ($chr !~ /^chr/) {
        $chr = "chr$chr";
    }

    my $start = $self->start;
    carp "Error: chromosome start coordinate not specified\n" if !$start;

    my $end = $self->end;
    carp "Error: chromosome end coordinate not specified\n" if !$end;

    my $hgw = ORCA::ConservationAnalysis::Run::hgWiggle->new();
    carp "Error initializing ORCA::ConservationAnalysis::Run::hgWiggle\n"
        if !$hgw;

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
    foreach my $idx (0..$end - $start) {
        my $pos = $idx + $start;
        $profile_array[$idx] = $profile_hash{$pos} || 0;
    }

    $self->{_profile_array} = @profile_array ? \@profile_array : undef;
}

=head2 compute_conserved_regions

 Title    : compute_conserved_regions
 Usage    : $regions = $ca->compute_conserved_regions();
 Function : Compute and return the list of conserved regions.
 Returns  : A reference to a list of Bio::SeqFeature::Generic objects.
 Args     : Optionally:
            -min_conservation   => the minimum percent identity of
                                   regions to report (specify as a
                                   number between 0 and 1 or a string
                                   such as '70%')
            -filter_exons       => (bool) indicates exons should be
                                   filtered out of conserved regions.
            -min_cr_len         => (int) minimum length of conserved
                                   region to keep

=cut

sub compute_conserved_regions
{
    my ($self, %args) = @_;

    my $min_conservation =
        $self->param('min_conservation', $args{-min_conservation});
    if ($min_conservation) {
        if ($min_conservation =~ /(.+)%/) {
            $min_conservation = $1 / 100;
        }
        if ($min_conservation <= 0 || $min_conservation > 1) {
            carp
                "min_conservation is out of range; please specify as a number"
                . " between 0 and 1 or as a string between 0% and 100%";
            return;
        }
        $self->param('min_conservation', $min_conservation);
    } else {
        carp "No minimum conservation defined\n";
        return;
    }

    my $filter_exons = $self->param('filter_exons', $args{-filter_exons});
    my $min_cr_len   = $self->param('min_cr_len',   $args{-min_cr_len});

    my $exons = $self->exons;

    my @conservation = @{$self->conservation_profile};

    my $start   = $self->start;
    my $end     = $self->end;

    if ($filter_exons && $exons) {
        foreach my $exon (@$exons) {
            my ($exon_start, $exon_end) = ($exon->start, $exon->end);

            next if ($exon_end < $start || $exon_start > $end);

            $exon_start = $start if $exon_start < $start;
            $exon_end   = $end   if $exon_end > $end;

            my $idx_start = $exon_start - $start;
            my $idx_end   = $exon_end - $start;

            # Zero out exon positions in the profile
            map {$conservation[$_] = 0} (($idx_start - 1) .. ($idx_end - 1));
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

                my $score = $self->_score_region(
                                \@conservation, $cr_start_idx, $cr_end_idx);

                my $cr = Bio::SeqFeature::Generic->new(
                    -primary_id   => sprintf("CR%d", $counter + 1),
                    -display_name => sprintf("CR%d", $counter + 1),
                    -source_tag   => "ORCA",
                    # leave in 0-based index coords for now
                    -start        => $cr_start_idx,
                    -end          => $cr_end_idx,
                    -score        => $score
                );

                push @crs, $cr;

                $cr_start_idx = undef;
                $cr_end_idx   = undef;
                $counter++;
            }
        }
    }

    if (defined $cr_start_idx) {
        # started a region and didn't finish it (scores remained above
        # min. conservation until end of sequence)
        $cr_end_idx = $#conservation;

        my $score =
            $self->_score_region(\@conservation, $cr_start_idx, $cr_end_idx);

        my $cr = Bio::SeqFeature::Generic->new(
            -primary_id   => sprintf("CR%d", $counter + 1),
            -display_name => sprintf("CR%d", $counter + 1),
            -source_tag   => "ORCA",
            # leave in 0-based index coords for now
            -start        => $cr_start_idx,
            -end          => $cr_end_idx,
            -score        => $score
        );

        push @crs, $cr;
    }

    # Combine regions into larger regions which still score
    # above min. conservation
    while (my $combined_crs =
        $self->_combine_conserved_regions(\@conservation, \@crs))
    {
        @crs = ();
        @crs = @{$combined_crs};
    }

    # Keep only regions > min. region length.
    my @final_crs;
    foreach my $cr (@crs) {
        if (($cr->end - $cr->start + 1) > $min_cr_len) {
            push @final_crs, $cr;
        }
    }

    #
    # Change back to absolute (chromosomal) coords DJA 09/12/11
    #
    #foreach my $cr (@final_crs) {
    #    $cr->start($cr->start + $start);
    #    $cr->end($cr->end + $start);
    #}

    #
    # Change back to relative 1-based coords DJA 09/12/11
    #
    foreach my $cr (@final_crs) {
        $cr->start($cr->start + 1);
        $cr->end($cr->end + 1);
    }


    # check if undef before assigning to prevent infinite recursion
    $self->{-conserved_regions} = @final_crs ? \@final_crs : undef;
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
    my ($self) = @_;

    my $seq = $self->seq;

    return if !defined $seq;

    my $crs = $self->conserved_regions;
    if (!defined $crs) {
        carp "no conserved regions defining sub-sequences to extract\n";
        return undef;
    }

    my $seq_display_id = $seq->display_id;
    #my $seq_start = $seq->start;
    #my $seq_end = $seq->end;
    my $seq_start = $self->start;
    my $seq_end   = $self->end;

    my @subseqs;
    foreach my $cr (@$crs) {
        # if CRs are using relative coords
        my $cr_start = $cr->start;
        my $cr_end   = $cr->end;
        # if CRs are using chromosomal coords
        #my $cr_start = $cr->start - $seq_start + 1;
        #my $cr_end   = $cr->end - $seq_start + 1;

        # trunc expects relative coords even if $seq is a Bio::LocatableSeq
        my $subseq = $seq->trunc($cr_start, $cr_end);

        $seq_display_id = "SUBSEQ" if !$seq_display_id;
        my $subseq_display_id = $seq_display_id . '/' . "$cr_start-$cr_end";
        $subseq->display_id($subseq_display_id);
        push @subseqs, $subseq;
    }

    $self->{-conserved_subseqs} = @subseqs ? \@subseqs : undef;
}

sub _combine_conserved_regions
{
    my ($self, $conservation, $crs) = @_;

    my $min_conservation = $self->param('min_conservation');

    my @crs2;
    my $combined = 0;
    my $i        = 0;
    while (defined $crs->[$i + 1]) {
        # 0-based index coords
        my $rstart = $crs->[$i]->start;
        my $rend   = $crs->[$i + 1]->end;

        my $score = $self->_score_region($conservation, $rstart, $rend);

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

sub _score_region
{
    my ($self, $conservation, $start, $end) = @_;

    my $score = 0;
    foreach my $i ($start - 1 .. $end - 1) {    # convert coords to 0-based
        $score += $conservation->[$i];
    }
    $score /= ($end - $start + 1);

    return $score;
}

1;
