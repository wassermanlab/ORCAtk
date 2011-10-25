
=head1 NAME

ORCA::Graphics::Pairwise - Object for graphing the results of an
ORCA pairwise analysis

=head1 SYNOPSIS

    use ORCA::Graphics::Pairwise;

    my $graph = ORCA::Graphics::Pairwise->new(
        -analysis               => $analysis,
        -cpg_islands            => \@cpgs,
        -cpg_islands2           => \@cpgs2,
        -other_features         => \@feats,
        -other_features2        => \@feats2,
        -flip                   => $flip_graph
    );

    open(PLOT, ">$plot_file);
    print PLOT $graph->{-gd_image}->png;
    close PLOT;

=head1 DESCRIPTION

ORCA::Graphics::Pairwise is an object used for the purpose of performing
graphical output of the results of an ORCA pairwise analysis.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::Graphics::Pairwise;

use Bio::Graphics::Panel;
use GD;
use Carp;

use strict;
use warnings;

use constant IMAGE_WIDTH => 1024;
#use constant IMAGE_HEIGHTi => 1024;
use constant LEFT_MARGIN   => 25;
use constant RIGHT_MARGIN  => 50;    # need extra room for TFBS labels
use constant TOP_MARGIN    => 5;
use constant BOTTOM_MARGIN => 10;

=head2

 Title    : new
 Usage    : my $graph = ORCA::Graphics::Pairwise->new(
                -analysis               => $analysis,
                -repeat_regions         => \@repeats, 
                -cpg_islands            => \@cpgs,
                -cpg_islands2           => \@cpgs2,
                -other_features         => \@feats,
                -other_features2        => \@feats2,
                -flip                   => $flip_graph
            );
 Function : Create a new ORCA::Graphics::Pairwise object.
 Returns  : An ORCA::Graphics::Pairwise object.
 Args     : -analysis        => An ORCA::Analysis::Pairwise object
            -cpg_islands     => An listre of Bio::SeqFeature objects
            -repeat_regions  => Optional listref of Bio::SeqFeature::Generic
                                objects indicating repeat regions.
            -cpg_islands     => Optional listref of Bio::SeqFeature::Generic
                                objects indicating CpG islands on the base
                                sequence.
            -cpg_islands2    => Optional listref of Bio::SeqFeature::Generic
                                objects indicating CpG islands on the
                                comparison sequence.
            -other_features  => Optional listref of Bio::SeqFeature::Generic
                                objects indicating additional features on the
                                base sequence to plot.
            -other_features2 => Optional listref of Bio::SeqFeature::Generic
                                objects indicating additional features on the
                                comparison sequence to plot.
            -flip            => Boolean indicating whether to reverse graph
                                (useful for plotting gene sequences which
                                fall on the reverse strand so that 5' to 3'
                                is oriented left to right.

=cut

sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
        -seq_no              => 1,
        -width               => IMAGE_WIDTH,
        #-height             => IMAGE_HEIGHT,
        -left_margin         => LEFT_MARGIN,
        -right_margin        => RIGHT_MARGIN,
        -top_margin          => TOP_MARGIN,
        -bottom_margin       => BOTTOM_MARGIN,

        -analysis            => undef,
        -flip                => undef,

        %args
    }, ref $caller || $caller;

    if (!$self->{-analysis}) {
        croak "No ORCA::Analysis::Pairwise object provided!";
    }

    $self->_draw_panel;

    return $self;
}

=head2

 Title    : analysis
 Usage    : my $analysis = ORCA::Graphics::Pairwise->analysis($pwa);
 Function : Get/set the pairwise analysis.
 Returns  : An ORCA::Analysis::Pairwise object.
 Args     : A reference to an ORCA::Analysis::Pairwise object.

=cut

sub analysis
{
    my ($self, $pwa) = @_;

    if ($pwa) {
        if (!$pwa->isa('ORCA::Analysis::Pairwise')) {
            carp "Not an ORCA::Analysis::Pairwise object\n";
            return;
        }

        $self->{-analysis} = $pwa;
    }

    return $self->{-analysis};
}

sub _draw_panel
{
    my ($self, %args) = @_;

    my $pwa = $self->analysis;
    if (!$pwa) {
        carp "No pairwise analysis object\n";
        return;
    }

    my $base_seq = $pwa->base_seq;
    if (!$base_seq) {
        carp "No pairwise analysis base sequence\n";
        return;
    }

    my $seq_start = $pwa->start || 1;
    my $seq_end   = $pwa->end   || $base_seq->length;

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
        #-gridcolor => 'lightgrey'
    );

    my $seqarea = Bio::SeqFeature::Generic->new(
        -start => $seq_start,
        -end   => $seq_end
    );

    #
    # Ruler
    #
    #print STDERR "ADDING TRACK: human seq\n";
    my $display_name = $base_seq->display_name;
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
    # Base seq exons track
    #
    if ($pwa->base_seq_exons) {
        #print STDERR "ADDING TRACK: exons\n";
        $panel->add_track(
            generic    => [$pwa->base_seq_exons],
            -key       => "Exons",
            -bgcolor   => "black",
            -fgcolor   => "black",
            -bump      => 0,
            -connector => "none",
            -label     => sub {$_[0]->seq_id}
        );
    }

    #
    # Comparison seq exons track
    #
    if ($pwa->comparison_seq_exons) {
        #print STDERR "ADDING TRACK: exons\n";
        $panel->add_track(
            generic    => [$pwa->comparison_seq_exons],
            -key       => "Exons2",
            -bgcolor   => "grey",
            -fgcolor   => "grey",
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
    # CpG islands 2 track
    #
    if ($self->{-cpg_islands2}) {
        #print STDERR "ADDING TRACK: CpG islands 2\n";
        my $cpg_islands = $self->{-cpg_islands2};
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
    # Other features 2 track
    #
    if ($self->{-other_features2}) {
        #print STDERR "ADDING TRACK: other features 2\n";
        my $other_features = $self->{-other_features2};
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
    if ($pwa->conserved_regions) {
        #print STDERR "ADDING TRACK: conserved regions\n";
        my @crs;
        foreach my $cr (@{$pwa->conserved_regions}) {
            push @crs, Bio::SeqFeature::Generic->new(
                -start      => $cr->start + $seq_start - 1,
                -end        => $cr->end + $seq_start - 1
            );
        }

        $panel->add_track(
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
    if ($pwa->conserved_tfbss) {
        #print STDERR "ADDING TRACK: Conserved TFBSs\n";
        my @tfbss;
        foreach my $tfbs (@{$pwa->conserved_tfbss}) {
            my $site1 = $tfbs->site1;
            my $new_tfbs = TFBS::Site->new(
                -start      => $site1->start + $seq_start - 1,
                -end        => $site1->end + $seq_start - 1,
                -score      => $site1->score,
                -strand     => $site1->strand
            );

            $new_tfbs->add_tag_value(
                'TF', ($site1->each_tag_value("TF"))[0]);

            push @tfbss, $new_tfbs;
        }

        $panel->add_track(
            generic    => [\@tfbss],
            -key       => "TF Sites",
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
    # Conservation profile
    #
    #print STDERR "ADDING TRACK: conservation profile\n";
    #
    # Conservation profile is now a 0-based array of scores DJA 09/12/11
    #
    my $cp  = $pwa->conservation_profile;
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

    $panel->{_conservation_cutoff} = $pwa->conservation_cutoff;
    $panel->add_track(
        [$seqarea],
        -glyph      => "my_xyplot",
        -key        => "Pairwise Conservation (score ranging from 0 to 1)",
        -height     => 200,
        -min_score  => 0,
        -max_score  => 1,
        -graph_type => "line",
        -fgcolor    => "red",
        #-scale      => "left"
    );

    #
    # Add other profile tracks
    #
#    $self->_add_profile_track(
#        $panel,
#        $self->{-muscle_profile},
#        (   -title  => 'Muscle',
#            -colour => 'blue'
#        )
#    ) if $self->{-muscle_profile};
#    $self->_add_profile_track(
#        $panel,
#        $self->{-mus_cons_profile},
#        (   -title  => 'Muscle Conservation',
#            -colour => 'magenta'
#        )
#    ) if $self->{-mus_cons_profile};
#    $self->_add_profile_track(
#        $panel,
#        $self->{-liver_profile},
#        (   -title  => 'Liver',
#            -colour => 'green'
#        )
#    ) if $self->{-liver_profile};
#    $self->_add_profile_track(
#        $panel,
#        $self->{-liv_cons_profile},
#        (   -title  => 'Liver Conservation',
#            -colour => 'cyan'
#        )
#    ) if $self->{-liv_cons_profile};

    $self->{-gd_image} = $panel->gd;
}

sub _add_profile_track
{
    my ($self, $panel, $profile, %args) = @_;

    return if !$profile;

    my $title  = $args{-title}  || "";
    my $colour = $args{-colour} || 'black';
    my $glyph  = $args{-glyph}  || 'xyplot';

    my $start = $self->{-start} || 1;
    my $end   = $self->{-end}   || $self->{-base_seq}->length;

    my $area = Bio::SeqFeature::Generic->new(
        -start => $start,
        -end   => $end
    );
    #print STDERR "ADDING TRACK: $title\n" if $title;
    foreach my $point (@$profile) {
        my $pos = $point->{-position};
        if ($pos >= $start && $pos <= $end) {
            #
            # Convert to percentage
            #
            my $score = $point->{-score} * 100;
            $score = 0.001 if $score < 0.001;
            $area->add_sub_SeqFeature(
                Bio::SeqFeature::Generic->new(
                    -start => $pos,
                    -end   => $pos,
                    -score => $score
                )
            );
        }
    }

    $panel->add_track(
        [$area],
        -glyph      => $glyph,
        -key        => $title,
        -height     => 200,
        -graph_type => "line",
        -fgcolor    => $colour,
        -scale      => "left"
    );
}

1;
