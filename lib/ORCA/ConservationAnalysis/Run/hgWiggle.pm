=head1 NAME

ORCA::ConservationAnalysis::Run::hgWiggle - Wrapper for Jim Kent source tree
(UCSC) hgWiggle program

=head1 SYNOPSIS

    use ORCA::ConservationAnalysis::Run::hgWiggle;

    my $hgw = ORCA::ConservationAnalysis::Run::hgWiggle->new();

    my $db = "hg18";
    my $pos = "chr6:108593955-108616706";
    my $opts = {-db => $db, -position => $pos}

    my $track = "phastCons28way";
    my $output = $hgw->run($opts, $track);

=head1 REQUIREMENTS

The hgWiggle program requires the database parameters to be stored in a file
called ".hg.conf" in the user's home directory. This file MUST HAVE a mode of
600 (-rw-------). An example of this file is:

    db.host=sonoma.cmmt.ubc.ca
    db.user=ucsc_r
    db.password=

The hgWiggle program also requires the MySQL libs, so the path to these
libs must be set (e.g. in LD_LIBRARY_PATH).

=head1 DESCRIPTION

Call the hgWiggle program with the given options and track name and return
the output.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::ConservationAnalysis::Run::hgWiggle;

use strict;

use constant PROG		=> "hgWiggle";

use Carp;
use File::Spec::Functions;	# for catfile

=head2

 Title    : new
 Usage    : $hgw = ORCA::ConservationAnalysis::Run::hgWiggle->new();
 Function : Create a new ORCA::ConservationAnalysis::Run::hgWiggle object.
 Returns  : An ORCA::ConservationAnalysis::Run::hgWiggle object.
 Args	  : None

=cut

sub new
{
    my ($class, %args) = @_;

    if (!_is_executable(PROG)) {
	carp "Error: could not locate " . PROG . " executable\n";
	return;
    }

    my $self = bless {
			%args
		    }, ref $class || $class;

    return $self;
}

=head2

 Title    : run
 Usage    : $output = $hgw->run($opts, $track);
 Function : Run hgWiggle and return results from STDOUT.
 Returns  : A scalar containing contents of STDOUT.
 Args	  : A hashref of valid command line arguments to hgWiggle,
 	    A scalar containing the name of the track for which to retrieve
	    data.

=cut

sub run
{
    my ($self, $opts, $track) = @_;

    if (!$opts || $opts == {}) {
	carp "Must provide hgWiggle command line options\n";
	return undef;
    }

    if (!$track) {
	carp "Must provide track name to hgWiggle\n";
	return undef;
    }

    my $opt_str = "";
    while (my ($key, $val) = (each %$opts)) {
	if ($key !~ /^-/) {
	    $key = "-$key";
	}

	if ($val) {
	    $opt_str .= " $key=$val";
	} else {
	    $opt_str .= " $key";
	}
    }

    my $ld_library_path = $ENV{LD_LIBRARY_PATH};
    if (!$ld_library_path) {
    	$ld_library_path = '/usr/local/mysql/lib/mysql';
    } elsif ($ld_library_path !~ /mysql\/lib\/mysql/) {
    	$ld_library_path .= ':/usr/local/mysql/lib/mysql';
    }
    $ENV{LD_LIBRARY_PATH} = $ld_library_path;

    my $cmd = PROG . "$opt_str $track";
    my $output = `$cmd`;
    my $rval = $? >> 8;

    if ($rval) {
	carp "Error running $cmd\n";
    }

    return $output
}

#
# Check to see if the specified program exists and is executable anywhere
# within the user's path.
#
sub _is_executable
{
	my ($prog) = @_;

	my $is_executable = 0;

	return $is_executable if !defined $prog;

	if (-x $prog) {
		$is_executable = 1;
	} else {
		my $path = $ENV{PATH};
		my @pathdirs = split(':', $path);
		foreach my $dir (@pathdirs) {
			my $exe = catfile($dir, $prog);
			if (-x $exe) {
				$is_executable = 1;
				last;
			}
		}
	}

	return $is_executable;
}

1;
