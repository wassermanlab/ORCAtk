package ORCA::ConservationAnalysis::ConservedRegion;

use strict;

sub new
{
    my ($class, %args) = @_;
    my $self = bless {
			-seq1_start	=> undef,
			-seq1_end	=> undef,
			-seq2_start	=> undef,
			-seq2_end	=> undef,
			-align_start	=> undef,
			-align_end	=> undef,
			-score		=> undef,
			%args
		    }, ref $class || $class;

    return $self;
}

#
# For backwards compatibility
#
sub seq_start
{
    my ($self, $seq_start) = @_;

    return $self->seq1_start($seq_start);
}

#
# For backwards compatibility
#
sub seq_end
{
    my ($self, $seq_end) = @_;

    return $self->seq1_end($seq_end);
}

sub seq1_start
{
    my ($self, $seq_start) = @_;

    if ($seq_start) {
	$self->{-seq1_start} = $seq_start;
    }

    return $self->{-seq1_start};
}

sub seq1_end
{
    my ($self, $seq_end) = @_;

    if ($seq_end) {
	$self->{-seq1_end} = $seq_end;
    }

    return $self->{-seq1_end};
}

sub seq2_start
{
    my ($self, $seq_start) = @_;

    if ($seq_start) {
	$self->{-seq2_start} = $seq_start;
    }

    return $self->{-seq2_start};
}

sub seq2_end
{
    my ($self, $seq_end) = @_;

    if ($seq_end) {
	$self->{-seq2_end} = $seq_end;
    }

    return $self->{-seq2_end};
}

sub align_start
{
    my ($self, $align_start) = @_;

    if ($align_start) {
	$self->{-align_start} = $align_start;
    }

    return $self->{-align_start};
}

sub align_end
{
    my ($self, $align_end) = @_;

    if ($align_end) {
	$self->{-align_end} = $align_end;
    }

    return $self->{-align_end};
}

sub score
{
    my ($self, $score) = @_;

    if ($score) {
	$self->{-score} = $score;
    }

    return $self->{-score};
}

#
# For backwards compatibility
#
sub length
{
    my ($self) = @_;

    return $self->length1;
}

sub length1
{
    my ($self) = @_;

    my $length = undef;
    my $seq_start = $self->seq1_start;
    my $seq_end = $self->seq1_end;
    if ($seq_start && $seq_end) {
	$length = $seq_end - $seq_start + 1;
    }
    return $length;
}

sub length2
{
    my ($self) = @_;

    my $length = undef;
    my $seq_start = $self->seq2_start;
    my $seq_end = $self->seq2_end;
    if ($seq_start && $seq_end) {
	$length = $seq_end - $seq_start + 1;
    }
    return $length;
}

sub to_string
{
    my ($self) = @_;

    my $str = undef;
    
    $str = sprintf("%d\t%d\t%d\t%d\t%d\t%d\t%.4f",
		    $self->seq1_start || 0, $self->seq1_end || 0,
		    $self->seq2_start || 0, $self->seq2_end || 0,
		    $self->align_start || 0, $self->align_end || 0,
		    $self->score || 0);

    return $str;
}

1;
