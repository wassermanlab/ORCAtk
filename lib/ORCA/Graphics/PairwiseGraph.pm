=head1 NAME

ORCA::Graphics::PairwiseGraph - Object for graphing the results of an
ORCA pairwise analysis

=head1 SYNOPSIS

    use ORCA::Graphics::PairwiseGraph;

    my $graph = ORCA::Graphics::PairwiseGraph->new(
			-base_seq               => $seq_obj,
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

ORCA::Graphics::PairwiseGraph is an object used for the purpose of performing
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

package ORCA::Graphics::PairwiseGraph;

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
 Usage    : my $graph = ORCA::Graphics::PairwiseGraph->new(
			-base_seq               => $seq_obj,
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
 Function : Create a new ORCA::Graphics::PairwiseGraph object.
 Returns  : An ORCA::Graphics::PairwiseGraph object.
 Args     : -base_seq               	=> Base sequence as a Bio::Seq or
					   Bio::LocatableSeq object
	    -cutoff			=> Conservation cutoff in the range
					   0-100
	    -conservation_profile	=> Conservation profile as an
					   arrayref of hashes with -position
					   and -score elements
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

	-base_seq		=> undef,
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
        -muscle_profile		=> undef,
        -liver_profile		=> undef,
        -mus_cons_profile	=> undef,
        -liv_cons_profile	=> undef,
        -start			=> undef,
        -end			=> undef,
	-flip			=> undef,

        %args
    }, ref $caller || $caller;

    if (!$self->{-base_seq}) {
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
    my $end = $self->{-end} || $self->{-base_seq}->length;

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
    my $display_name = $self->{-base_seq}->display_name;
    if ($self->{-flip}) {
    	$display_name .= " (reversed)";
    }
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
    # CpG islands 2 track
    #
    if ($self->{-cpg_islands2}) {
	#print STDERR "ADDING TRACK: CpG islands 2\n";
	$panel->add_track (
	    generic	=> [$self->{-cpg_islands2}],
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
    # Other features 2 track
    #
    if ($self->{-features2}) {
	#print STDERR "ADDING TRACK: other features 2\n";
	$panel->add_track (
	    generic	=> [$self->{-features2}],
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
    # Conservation profile
    #
    #print STDERR "ADDING TRACK: conservation profile\n";
    my $cp = $self->{-conservation_profile};
    foreach my $point (@$cp) {
	my $pos = $point->{-position};
	if ($pos >= $start && $pos <= $end) {
	    my $score = $point->{-score} * 100;
	    $score = 0.001 if $score < 0.001;
	    $seqarea->add_sub_SeqFeature(
		Bio::SeqFeature::Generic->new(
		    -start	=> $pos,
		    -end	=> $pos,
		    -score	=> $score
		)
	    );
	}
    }

    $panel->{_conservation_cutoff} = $self->{-cutoff};
    $panel->add_track(
        [$seqarea],
	-glyph      => "my_xyplot",
        -key        => "Conservation Profile",
        -height     => 200,
        -max_score  => 100,
        -min_score  => 0,
        -graph_type => "line",
        -fgcolor    => "red",
        -scale      => "left"
    );

    #
    # Add other profile tracks
    #
    $self->_add_profile_track($panel, $self->{-muscle_profile},
				(
				    -title => 'Muscle',
				    -colour => 'blue'
				)
			    ) if $self->{-muscle_profile};
    $self->_add_profile_track($panel, $self->{-mus_cons_profile},
				(
				    -title => 'Muscle Conservation',
				    -colour => 'magenta'
				)
			    ) if $self->{-mus_cons_profile};
    $self->_add_profile_track($panel, $self->{-liver_profile},
				(
				    -title => 'Liver',
				    -colour => 'green'
				)
			    ) if $self->{-liver_profile};
    $self->_add_profile_track($panel, $self->{-liv_cons_profile},
				(
				    -title => 'Liver Conservation',
				    -colour => 'cyan'
				)
			    ) if $self->{-liv_cons_profile};

    $self->{-gd_image} = $panel->gd;
}

sub _add_profile_track
{
    my ($self, $panel, $profile, %args) = @_;

    return if !$profile;

    my $title = $args{-title} || "";
    my $colour = $args{-colour} || 'black';
    my $glyph = $args{-glyph} || 'xyplot';

    my $start = $self->{-start} || 1;
    my $end = $self->{-end} || $self->{-base_seq}->length;

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
		    -start	=> $pos,
		    -end	=> $pos,
		    -score	=> $score
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
