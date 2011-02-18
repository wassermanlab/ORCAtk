
=head1 NAME

ORCA::Web::State - Object for keeping track of the state information of the
ORCA web application.

=head1 SYNOPSIS

    use ORCA::Web::State;

    my $state;
    my $sid = $q->param('sid');
    if ($sid) {
        #
        # If session ID is already defined, load the state from file
        #
        my $filename = ABS_TMP_PATH . "/" . $sid;
        $state = ORCA::Web::State->load(__Fn => $filename);
    } else {
        #
        # Create a new session with a unique session ID (in this example the
        # process ID plus current date/time
        #
        $sid = $$ . time;
        my $filename = ABS_TMP_PATH . "/" . $sid;
        $state = ORCA::Web::State->new(__Fn => $filename, -sid => $sid);

        # Set some parameters
        $state->sid($sid);
        $state->species1('human');
        $state->species2('mouse');
    }

    # Save the state to a file
    $state->dumper->Purity(1);
    $state->dumper->Deepcopy(1);
    $state->commit();

=head1 DESCRIPTION

ORCA::Web::State is an object used for the purpose of keeping track of the
web application's state between pages. Inherits from
Persistence::Object::Simple to save and load state information to/from a flat
file between succesive CGI calls based on a unique session ID. It uses the
AUTOLOAD facility to define it's members/methods dynamically.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::Web::State;

use strict;

use vars qw/@ISA $AUTOLOAD/;

use Carp;
use Persistence::Object::Simple;

@ISA = qw/Persistence::Object::Simple/;


sub AUTOLOAD
{
    my $self = shift;

    my $type = ref($self) || croak "$self is not an object";

    my $name = $AUTOLOAD;
    $name =~ s/.*://;    # strip fully-qualified portion

    if (@_) {
        $self->{-$name} = shift;
    }

    return $self->{-$name};
}

1;
