=head1 NAME

ORCA::Graphics::PhastCons - Object for graphing the results of an ORCA
phastCons analysis.

=head1 SYNOPSIS

 use ORCA::Graphics::PhastCons;

 my $graph = ORCA::Graphics::PhastCons->new(
     # An ORCA::Analysis::PhastCons object
     -analysis       => $analysis,

     -repeat_regions => $repeats,
     -cpg_islands    => $cpgs,
     -other_features => $features,

     # If set to true, graph is flipped left
     # to right. Useful for -ve strand genes
     # so 5' end is on the left.
     -flip           => $flip_graph
 );

 open(PLOT, ">$plot_file);
 print PLOT $graph->{-gd_image}->png;
 close PLOT;

=head1 DESCRIPTION

ORCA::Graphics::PhastCons is an object used for the purpose of performing
graphical output of the results of an ORCA phastCons analysis.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::Graphics::PhastCons;

use Bio::Graphics::Panel;
use Bio::SeqFeature::Generic;
use GD;
use Carp;

use strict;
use warnings;

use constant IMAGE_WIDTH    => 1024;
#use constant IMAGE_HEIGHTi  => 1024;
use constant LEFT_MARGIN    => 25;
use constant RIGHT_MARGIN   => 50;    # need extra room for TFBS labels
use constant TOP_MARGIN     => 5;
use constant BOTTOM_MARGIN  => 10;

=head2

 Title    : new
 Usage    : my $graph = ORCA::Graphics::PhastCons->new(
                -analysis       => $phca,
                -repeat_regions => $repeats,
                -cpg_islands    => $cpgs,
                -other_features => $features,
                -flip           => $flip_graph
            );
 Function : Create a new ORCA::Graphics::PhastCons object.
 Returns  : An ORCA::Graphics::PhastCons object.
 Args     : -analysis       => A reference to an ORCA::Analysis::PhastCons
                               object.
            -repeat_regions => Optional listref of Bio::SeqFeature::Generic
                               objects indicating repeat regions.
            -cpg_islands    => Optional listref of Bio::SeqFeature::Generic
                               objects indicating CpG islands.
            -other_features => Optional listref of Bio::SeqFeature::Generic
                               objects indicating additional features to
                               plot.
            -flip           => Optional boolean indicating wether to reverse
                               graph (useful for plotting gene sequences
                               which fall on the reverse strand so that
                               5' to 3' is oriented left to right).

=cut

sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
        -seq_no         => 1,
        -width          => IMAGE_WIDTH,
        #-height         => IMAGE_HEIGHT,
        -left_margin    => LEFT_MARGIN,
        -right_margin   => RIGHT_MARGIN,
        -top_margin     => TOP_MARGIN,
        -bottom_margin  => BOTTOM_MARGIN,

        -analysis       => undef,
        -flip           => undef,
        %args
    }, ref $caller || $caller;

    if (!$self->analysis()) {
        croak "No ORCA::Analysis::PhastCons object provided!";
    }

    $self->_draw_panel;

    return $self;
}

=head2

 Title    : analysis
 Usage    : my $analysis = ORCA::Graphics::PhastCons->analysis($phca);
 Function : Get/set the phastCons analysis.
 Returns  : An ORCA::Analysis::PhastCons object.
 Args     : A reference to an ORCA::Analysis::PhastCons object.

=cut

sub analysis
{
    my ($self, $phca) = @_;

    if ($phca) {
        if (!$phca->isa('ORCA::Analysis::PhastCons')) {
            carp "Not an ORCA::Analysis::PhastCons object\n";
            return;
        }

        $self->{-analysis} = $phca;
    }

    return $self->{-analysis};
}

sub _draw_panel
{
    my ($self, %args) = @_;

    my $phca = $self->analysis;
    if (!$phca) {
        carp "No phastCons analysis object\n";
        return;
    }

    my $seq = $phca->seq;
    if (!$seq) {
        carp "No phastCons analysis sequence\n";
        return;
    }

    my $seq_start = $phca->start || 1;
    my $seq_end   = $phca->end   || $seq->length();

    my $panel = Bio::Graphics::Panel->new(
        -width      => $self->{-width},
        -key_style  => "between",
        -pad_left   => $self->{-left_margin},
        -pad_right  => $self->{-right_margin},
        -pad_top    => $self->{-top_margin},
        -pad_bottom => $self->{-bottom_margin},
        -start      => $seq_start,
        -stop       => $seq_end,
        -grid       => 1,
        -flip       => $self->{-flip},
        #-gridcolor  => 'lightgrey'
    );

    my $seqarea = Bio::SeqFeature::Generic->new(
        -start => $seq_start,
        -end   => $seq_end
    );

    #
    # Ruler
    #
    #print STDERR "ADDING TRACK: human seq\n";
    my $display_name = $seq->display_name;
    $panel->add_track(
        arrow => [$seqarea],
        -tick => 2,
        -key  => "$display_name"
    );

    #
    # Repeats
    #
    if ($self->{-repeat_regions}) {
        #print STDERR "ADDING TRACK: repeats\n";

        my $repeats = $self->{-repeat_regions};
        foreach my $repeat (@$repeats) {
            $repeat->start($repeat->start + $seq_start - 1);
            $repeat->end($repeat->end + $seq_start - 1);
        }

        $panel->add_track(
            generic    => [$repeats],
            -key       => "Repeats",
            -fgcolor   => "black",
            -bgcolor   => "black",
            -bump      => 0,
            -connector => "none",
            #-label => sub { shift->primary_tag }
        );
    }

    #
    # Exons track
    #
    if ($phca->exons()) {
        #print STDERR "ADDING TRACK: exons\n";

        #my @exons;
        #foreach my $exon (@{$phca->exons}) {
        #    push @exons, Bio::SeqFeature::Generic->new(
        #        #-primary_tag    => $exon->primary_tag,
        #        #-source_tag     => $exon->source_tag,
        #        #-display_name   => $exon->display_name,
        #        #-strand         => $exon->strand,
        #        -start          => $exon->start + $seq_start - 1,
        #        -end            => $exon->end + $seq_start - 1
        #    );
        #}

        $panel->add_track(
            generic    => [$phca->exons],
            #generic    => [\@exons],
            -key       => "Exons",
            -bgcolor   => "black",
            -fgcolor   => "black",
            -bump      => 0,
            -connector => "none",
            -label     => sub {$_[0]->seq_id}
        );
    }

    #
    # CpG islands track
    #
    if ($self->{-cpg_islands}) {
        #print STDERR "ADDING TRACK: CpG islands\n";

        my $cpg_islands = $self->{-cpg_islands};
        my @cpgs;
        foreach my $cpg (@$cpg_islands) {
            push @cpgs, Bio::SeqFeature::Generic->new(
                -start  => $cpg->start + $seq_start - 1,
                -end    => $cpg->end + $seq_start - 1
            );
        }

        $panel->add_track(
            generic    => [\@cpgs],
            -key       => "CpG Islands",
            -bgcolor   => "green",
            -fgcolor   => "green",
            -bump      => 1,
            -connector => "none",
            -label     => sub {$_[0]->seq_id}
        );
    }

    #
    # Other features track
    #
    if ($self->{-other_features}) {
        #print STDERR "ADDING TRACK: other features\n";

        my $other_features = $self->{-other_features};
        foreach my $other_feature (@$other_features) {
            $other_feature->start($other_feature->start + $seq_start - 1);
            $other_feature->end($other_feature->end + $seq_start - 1);
        }

        $panel->add_track(
            generic    => [$other_features],
            -key       => "Other Features",
            -bgcolor   => "purple",
            -fgcolor   => "purple",
            -bump      => 1,
            -connector => "none",
            -label     => sub {$_[0]->seq_id}
        );
    }

    #
    # Conserved regions track
    #
    if ($phca->conserved_regions) {
        #print STDERR "ADDING TRACK: conserved regions\n";

        my @crs;
        foreach my $cr (@{$phca->conserved_regions}) {
            push @crs, Bio::SeqFeature::Generic->new(
                -start      => $cr->start + $seq_start - 1,
                -end        => $cr->end + $seq_start - 1
            );
        }

        $panel->add_track(
            generic    => [$phca->conserved_regions],
            generic    => [\@crs],
            -key       => "Conserved Regions",
            -bgcolor   => "turquoise",
            -fgcolor   => "turquoise",
            -bump      => 0,
            -connector => "none",
            #-label     => sub {$_[0]->primary_tag}
        );
    }

    # Conserved TFBSs
    if ($phca->conserved_tfbss) {
        #print STDERR "ADDING TRACK: Conserved TFBSs\n";

        my @tfbss;
        foreach my $tfbs (@{$phca->conserved_tfbss}) {
            my $new_tfbs = Bio::SeqFeature::Generic->new(
                -start      => $tfbs->start + $seq_start - 1,
                -end        => $tfbs->end + $seq_start - 1,
                -score      => $tfbs->score,
                -strand     => $tfbs->strand
            );

            $new_tfbs->add_tag_value(
                'TF', ($tfbs->each_tag_value("TF"))[0]);

            push @tfbss, $new_tfbs;
        }

        $panel->add_track(
            generic    => [\@tfbss],
            -key       => "Conserved TFBSs",
            -fgcolor   => 'blue',
            -bgcolor   => 'blue',
            -bump      => 1,
            -connector => "none",
            -label     => sub {
                $_[0]->has_tag("TF") && ($_[0]->each_tag_value("TF"))[0]
                    || "?";
            }
        );
    }

    #
    # PhastCons conservation profile
    #
    #print STDERR "ADDING TRACK: phastCons conservation profile\n";
    #
    # Conservation profile is now a 0-based array of scores DJA 09/12/11
    #
    my $cp  = $phca->conservation_profile;
    my $idx = 0;
    foreach my $score (@$cp) {
        my $pos = $idx + $seq_start;
        if ($pos <= $seq_end) {
            $score = 0.001 if $score < 0.001;
            $seqarea->add_sub_SeqFeature(
                Bio::SeqFeature::Generic->new(
                    -start => $pos,
                    -end   => $pos,
                    -score => $score
                )
            );
        }
        $idx++;
    }

    $panel->{_conservation_cutoff} = $phca->min_conservation();
    $panel->add_track(
        [$seqarea],
        -glyph      => "my_xyplot",
        -key        => sprintf("Conservation (%s score ranging from 0 to 1)",
                               $phca->track),
        -height     => 200,
        -min_score  => 0,
        -max_score  => 1,
        -graph_type => "line",
        -fgcolor    => "red",
        #-scale      => "left"
    );

    $self->{-gd_image} = $panel->gd;
}

1;
