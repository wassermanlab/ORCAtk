=head1 NAME

ORCA::Graphics::PhastConsGraph - Object for graphing the results of an
ORCA phastCons analysis

=head1 SYNOPSIS

    use ORCA::Graphics::PhastConsGraph;

    my $graph = ORCA::Graphics::PhastConsGraph->new(
			-seq               	=> $seq_obj,
			-cutoff                 => $conservation_cutoff,
			-conservation_profile   => \@conservation_profile,
			-conserved_regions      => \@conserved_regions,
			-tf_sites               => \@tf_sites,
			-exons                  => \@exons,
			-cpg_islands            => \@cpgs,
			-start                  => $start,
			-end                    => $end,
			-flip                   => $flip_graph
        );

    open(PLOT, ">$plot_file);
    print PLOT $graph->{-gd_image}->png;
    close PLOT;

=head1 DESCRIPTION

ORCA::Graphics::PhastConsGraph is an object used for the purpose of performing
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

package ORCA::Graphics::PhastConsGraph;

use Bio::Graphics::Panel;
use GD;
use Carp;

use strict;
use warnings;

use constant IMAGE_WIDTH	=> 1024;
#use constant IMAGE_HEIGHTi	=> 1024;
use constant LEFT_MARGIN	=> 25;
use constant RIGHT_MARGIN	=> 50; # need extra room for TFBS labels
use constant TOP_MARGIN		=> 5;
use constant BOTTOM_MARGIN	=> 10;

=head2

 Title    : new
 Usage    : my $graph = ORCA::Graphics::PhastConsGraph->new(
			-seq               	=> $seq_obj,
			-cutoff                 => $conservation_cutoff,
			-conservation_profile   => \@conservation_profile,
			-conserved_regions      => \@conserved_regions,
			-tf_sites               => \@tf_sites,
			-exons                  => \@exons,
			-cpg_islands            => \@cpgs,
			-start                  => $start,
			-end                    => $end,
			-flip                   => $flip_graph
		    );
 Function : Create a new ORCA::Graphics::PhastConsGraph object.
 Returns  : An ORCA::Graphics::PhastConsGraph object.
 Args     : -seq               		=> Base sequence as a Bio::Seq or
					   Bio::LocatableSeq object
	    -cutoff			=> Conservation cutoff in the range
					   0-1
	    -conservation_profile	=> Conservation profile as an array of scores
                                   Changed from an array of hashes DJA 09/12/11
	    -conserved_regions		=> Conserved regions as an arrayref
					   of Bio::SeqFeature::Generic
					   objects
	    -tf_sites			=> TFBSs as an arrayref of
					   TFBS::Site objects
	    -exons			=> Exons as a listref of
					   Bio::SeqFeature::Gene::Exon
	    -cpg_islands		   objects
	    -start			=> Chromosomal start coordinate
					   of graph
	    -end			=> Chromosomal end coordinate of
					   graph
	    -flip			=> Boolean indicating wether to
					   reverse graph (useful for
					   plotting gene sequences which
					   fall on the reverse strand so
					   that 5' to 3' is oriented left
					   to right

=cut


sub new
{
    my ($caller, %args) = @_;

    my $self = bless {
        -seq_no			=> 1,
        -width			=> IMAGE_WIDTH,
	#-height    		=> IMAGE_HEIGHT,
        -left_margin		=> LEFT_MARGIN,
        -right_margin		=> RIGHT_MARGIN,
        -top_margin		=> TOP_MARGIN,
        -bottom_margin		=> BOTTOM_MARGIN,

	-seq			=> undef,
        -alignment		=> undef,
        -conservation_profile	=> undef,
        -conserved_regions	=> undef,
	-cutoff			=> undef,
	-tf_sites		=> undef,
	-exons			=> undef,
	-exons2			=> undef,
	-cpg_islands		=> undef,
	-features		=> undef,
	-repeat_regions		=> undef,
        -start			=> undef,
        -end			=> undef,
	-flip			=> undef,

        %args
    }, ref $caller || $caller;

    if (!$self->{-seq}) {
        croak "No base sequence provided!";
    }

    if (!$self->{-conservation_profile}) {
        croak "No conservation profile provided!";
    }

    $self->_draw_panel;

    return $self;
}

sub _draw_panel
{
    my ($self, %args) = @_;

    my $start = $self->{-start} || 1;
    my $end = $self->{-end} || $self->{-seq}->length;

    my $panel = Bio::Graphics::Panel->new(
        -width		=> $self->{-width},
        -key_style	=> "between",
        -pad_left	=> $self->{-left_margin},
        -pad_right	=> $self->{-right_margin},
        -pad_top	=> $self->{-top_margin},
        -pad_bottom	=> $self->{-bottom_margin},
        -start		=> $start,
        -stop		=> $end,
        -grid		=> 1,
	-flip		=> $self->{-flip},
	#-gridcolor	=> 'lightgrey'
    );

    my $seqarea = Bio::SeqFeature::Generic->new(
        -start => $start,
        -end   => $end
    );

    #
    # Ruler
    #
    #print STDERR "ADDING TRACK: human seq\n";
    my $display_name = $self->{-seq}->display_name;
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
	$panel->add_track(
	    generic	=> [$self->{-repeat_regions}],
	    -key	=> "Repeats",
	    -fgcolor	=> "black",
	    -bgcolor	=> "black",
	    -bump	=> 0,
	    -connector	=> "none",
	    #-label => sub { shift->primary_tag }
	);
    }

    #
    # Base seq exons track
    #
    if ($self->{-exons}) {
	#print STDERR "ADDING TRACK: exons\n";
	$panel->add_track (
	    generic	=> [$self->{-exons}],
	    -key	=> "Exons",
	    -bgcolor	=> "black",
	    -fgcolor	=> "black",
	    -bump	=> 0,
	    -connector	=> "none",
	    -label	=> sub {$_[0]->seq_id}
	);
    }

    #
    # Comparison seq exons track
    #
    if ($self->{-exons2}) {
	#print STDERR "ADDING TRACK: exons\n";
	$panel->add_track (
	    generic	=> [$self->{-exons2}],
	    -key	=> "Exons2",
	    -bgcolor	=> "grey",
	    -fgcolor	=> "grey",
	    -bump	=> 0,
	    -connector	=> "none",
	    -label	=> sub {$_[0]->seq_id}
	);
    }

    #
    # CpG islands track
    #
    if ($self->{-cpg_islands}) {
	#print STDERR "ADDING TRACK: CpG islands\n";
	$panel->add_track (
	    generic	=> [$self->{-cpg_islands}],
	    -key	=> "CpG Islands",
	    -bgcolor	=> "green",
	    -fgcolor	=> "green",
	    -bump	=> 1,
	    -connector	=> "none",
	    -label	=> sub {$_[0]->seq_id}
	);
    }

    #
    # Other features track
    #
    if ($self->{-features}) {
	#print STDERR "ADDING TRACK: other features\n";
	$panel->add_track (
	    generic	=> [$self->{-features}],
	    -key	=> "Other Features",
	    -bgcolor	=> "purple",
	    -fgcolor	=> "purple",
	    -bump	=> 1,
	    -connector	=> "none",
	    -label	=> sub {$_[0]->seq_id}
	);
    }

    #
    # Conserved regions track
    #
    if ($self->{-conserved_regions}) {
	#print STDERR "ADDING TRACK: conserved regions\n";
	$panel->add_track (
	    generic	=> [$self->{-conserved_regions}],
	    -key	=> "Conserved Regions",
	    -bgcolor	=> "turquoise",
	    -fgcolor	=> "turquoise",
	    -bump	=> 0,
	    -connector	=> "none",
	    #-label		=> sub {$_[0]->primary_tag}
	);
    }

    # TF sites
    if ($self->{-tf_sites}) {
	#print STDERR "ADDING TRACK: TF sites\n";
	$panel->add_track(
	    generic  	=> [$self->{-tf_sites}],
	    -key     	=> "TF Sites",
	    -fgcolor 	=> 'blue',
	    -bgcolor 	=> 'blue',
	    -bump	=> 1,
	    -connector	=> "none",
	    -label 	=> sub {$_[0]->has_tag("TF")
			&& ($_[0]->each_tag_value("TF"))[0] || "?"}
	);
    }

    #
    # PhastCons profile
    #
    #print STDERR "ADDING TRACK: phastCons profile\n";
    #
    # Conservation profile is now a 0-based array of scores DJA 09/12/11
    #
    my $cp = $self->{-conservation_profile};
    my $idx = 0;
    foreach my $score (@$cp) {
        my $pos = $idx + $start;
        if ($pos >= $start && $pos <= $end) {
            $score = 0.001 if $score < 0.001;
            $seqarea->add_sub_SeqFeature(
                Bio::SeqFeature::Generic->new(
                    -start	=> $pos,
                    -end	=> $pos,
                    -score	=> $score
                )
            );
        }
        $idx++;
    }

    $panel->{_conservation_cutoff} = $self->{-cutoff};
    $panel->add_track(
        [$seqarea],
	-glyph      => "my_xyplot",
        -key        => "Conservation Profile (PhastCons score ranging from 0 to 1)",
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
