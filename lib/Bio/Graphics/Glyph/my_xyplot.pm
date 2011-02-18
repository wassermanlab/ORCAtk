package Bio::Graphics::Glyph::my_xyplot;

use strict;
use Bio::Graphics::Glyph::xyplot;
use vars '@ISA';

@ISA = 'Bio::Graphics::Glyph::xyplot';

use constant DEFAULT_POINT_RADIUS=>1;

my %SYMBOLS = (
	       triangle => \&draw_triangle,
	       square   => \&draw_square,
	       disc     => \&draw_disc,
	       point    => \&draw_point,
	      );


sub draw {
  my $self = shift;
  my ($gd,$dx,$dy) = @_;
  my ($left,$top,$right,$bottom) = $self->calculate_boundaries($dx,$dy);
  my $gray = $gd->colorAllocate(220,220,220);
  my $black = $gd->colorAllocate(0,0,0);

  my $cut = $self->panel->{_conservation_cutoff};
  $self->filled_box($gd,$left,$top+($bottom-$top)*(1-$cut),
		    $right,$bottom,$gray, $gray);
  $gd->rectangle($left,$top-2,$right,$bottom, $black);

  return $self->SUPER::draw(@_);

}

1;

