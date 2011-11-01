=head1 NAME

ORCA::Analysis::Run::AlignCons - Wrapper for align_cons program

=head1 SYNOPSIS

    use ORCA::Analysis::Run::AlignCons;

    my $alignIO = Bio::AlignIO->new(
				    -file	=> "align.orca",
				    -format	=> "fasta");
    my $aln = $alignIO->next_aln;

    my $ac = ORCA::Analysis::Run::AlignCons->new();

    my %params = (
		   -r	=> 'p',
		   -w	=> 100,
		   -n	=> 1,
		   -f	=> 'c');

    my $report = $ac->run(-alignment => $aln, %params);

=head1 DESCRIPTION

The align_cons program is a program which analysis an alignment of two DNA
sequences and generates one of two types of reports. The first is a simple
report containing a list of positions and the corresponding percent identity
score for that position. The second is a list of conserved regions defined as
regions whose percent identity score at or above a certain threshold over a
certain minimum number of bases on the first sequence in the alignment
(base sequence).

This is accomplished by sliding a window of a given number of nucleotides on
the base sequence over the alignment and counting the number of nucleotide
matches on the comparison sequence to compute a percent identity score for the
window.

For the first type of report, the position and score are simply output for each
window. The -f parameter determines whether the position output is the start,
end of center of the window. For the second type of report the conserved
regions are computed by merging each of the scoring windows such that the
overall percent identity of each of the merged regions remains at or above
threshold. This threshold can be specified explicitly with the -t parameter.
Alternatively a stringency may be specified using the -s parameter. In this
case the threshold is computed dynamically as the value above which X percent
of windows score. I.e. If the -s parameter has a value of 0.1, the threshold T
is computed such that 10% of windows score greater than or equal to T. If both
the -t and -s parameters are specified, a dynamic threshold is computed, but if
it falls below the threshold specified by the -t argument, the specified
threshold argument is used.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::Analysis::Run::AlignCons;

use strict;

use constant PROG		=> "align_cons";
#use constant DFLT_REPORT_TYPE	=> 'c';	# should be 'c' or 'p'
#use constant DFLT_WIN_SIZE	=> 100;
#use constant DFLT_WIN_INC	=> 1;
#use constant DFLT_STRINGENCY	=> 0.1;
#use constant DFLT_REPORT_FORMAT	=> 'c'; # should be 's', 'e' or 'c'

use Carp;
use File::Temp qw/tempfile tempdir/;
use File::Path;			# for rmtree
use File::Spec::Functions;	# for catfile
use IO::String;
use Bio::AlignIO;
use ORCA::Analysis::ConservedRegion;
use ORCA::Analysis::ConservationReport;

=head2

 Title    : new
 Usage    : $ac = ORCA::Analysis::Run::AlignCons->new(%params);
 Function : Create a new ORCA::Analysis::Run::AlignCons object.
 Returns  : An ORCA::Analysis::Run::AlignCons object.
 Args	  : None.

=cut

sub new
{
    my ($class, %args) = @_;

    if (!_is_executable(PROG)) {
	carp "error: could not locate " . PROG . " executable\n";
	return;
    }

    my $self = bless {
			_tempdir	=> undef,
			_tempalignfile	=> undef,
			%args
		    }, ref $class || $class;

    return $self;
}

=head2

 Title    : run
 Usage    : $report = $ac->run(%args);
 Function : Run conservation analysis and report result.
 Returns  : An ORCA::Analysis::ConservationReport object.
 Args	  :   -alignment => A Bio::SimpleAlign object.
 	      -features	=> A list of Bio::SeqFeatureI compliant objects.
	    		   If provided these features will be filtered out
			   of the conserved regions.
	      -r	=> report type: 'c' for conserved regions, 'p' for
	                   conservation profile (position/score)
              -w	=> size of scoring window in bp
	      -n	=> window increment - number of bp by which to slide
	                   scoring window on each iteration
	      -t	=> threshold - min. identity score of consserved
	      		   regions to report; range: 0 <= t <= 1 (may also
			   specify as a string, e.g. -t => '70%').
			   only relevant to conserved regions analysis
	      -s	=> stringency - top X windows are used to dynamically
	                   calculate threshold; range: 0 < s < 1 (may also
			   specify as a string, e.g. -s => '10%').
			   results in a threshold calulated such that 10% of
			   all scoring windows fall above that threshold.
			   only relevant to conserved regions analysis
	      -f	=> format 's', 'e' or 'c' - for conservation profile
	                   analysis specifies whether scoring window
			   positions are returned as start, end or center
	      -m	=> method for computing percent identity,
	      		   's' = standard, 'o' = overall. 
	      -c	=> flag indicating sequence 2 was reverse complemented
	      		   in alignment

=cut

sub run
{
    my ($self, %args) = @_;

    my $alignment = $args{-alignment};
    if (!defined $alignment) {
	carp "Must provide alignment\n";
	return undef;
    }

    if (!$alignment->isa("Bio::SimpleAlign")) {
	carp "alignment is not a Bio::SimpleAlign object\n";
	return undef;
    }

    my $features = $args{-features};
    if ($features) {
	if (ref $features ne "ARRAY"
		|| !$features->[0]
		|| !$features->[0]->isa("Bio::SeqFeatureI"))
	{
		    carp "features is not a list ref of Bio::SeqFeatureI"
			. " compliant objects\n";
		    return;
	}
    }

    # Don't create tempdir until we have to and only if it has not already
    # been created.
    if (!defined $self->{_tempdir} || !-d $self->{_tempdir}) {
	$self->{_tempdir} = tempdir("AlignConsXXXXXXXX", TMPDIR => 1,
	    				CLEANUP => 1),
    }

    #
    # Don't create tempalignfile until we have to.
    #
    $self->{_tempalignfile} = File::Temp::tempnam($self->{_tempdir}, "align"); 

    my $alignIO = Bio::AlignIO->new(
				-file	=> ">" . $self->{_tempalignfile},
				-format	=> 'fasta');
    return if !$alignIO;
    return if !$alignIO->write_aln($alignment);
    $alignIO->close();

    if ($features) {
	$self->{_tempfeatfile} = File::Temp::tempnam($self->{_tempdir}, "feat");
	return if !_write_features_as_GFF($self->{_tempfeatfile}, $features);
    }

    my $stringency = $args{-s};
    if (defined $stringency) {
        if ($stringency =~ /(.+)%/)  {
            $stringency = $1 / 100;
        }
        if ($stringency <= 0 || $stringency >= 1) {
            carp "top_pct is out of range; please specify as a number between"
                    . " 0 and 1 or as a string between 0% and 100%";
            return;
        }
    }

    my $threshold = $args{-t};
    if (defined $threshold) {
        if ($threshold =~ /(.+)%/)  {
            $threshold = $1 / 100;
        }
        if ($threshold <= 0 || $threshold > 1) {
	    carp "min_conservation is out of range; please specify as a number"
		    . " between 0 and 1 or as a string between 0% and 100%";
	    return;
        }
    }


    my %ac_params = (
		    -r => $args{-r},
		    -w => $args{-w},
		    -n => $args{-n},
		    -s => $stringency,
		    -t => $threshold,
		    -m => $args{-m},
		    -f => $args{-f},
		    -l => $args{-l},
		    -i => $self->{_tempalignfile},
		    -g => $self->{_tempfeatfile}
		);

    my @ac_args;
    while (my ($param, $value) = each %ac_params) {
	if (defined $value) {
	    push @ac_args, $param, $value;
	}
    }

    # Boolean parameters
    if ($args{-W}) {
	push @ac_args, '-W';
    }
    #
    # Now -c arg is passed in rather than computed from alignment.
    # DJA 2008/01/10
    #
    if ($args{-c}) {
	push @ac_args, '-c';
    }

    my $cmd = join " ", PROG, @ac_args;
    my @ac_output = `$cmd `;#2>&1`;
      # ^ I commented out the STDERR redirection because this module does not attempt to catch warnings.
      #   With the redirection, errors will be slurped by this module and the user will never see them.
    my $rval = $? >> 8;

    my $report;
    if ($rval == 0) {
	if (!$args{-r} || $args{-r} eq 'c') {
	    #my $conserved = ORCA::Analysis::ConservedRegionList->new();
	    my @conserved = ();
	    my $window_size;
	    my $window_inc;
	    my $top_x_pct;
	    my $min_pct_id;
	    my $comp_pct_id;
	    my $eff_pct_id;
	    foreach my $line (@ac_output) {
		chomp $line;
		if ($line
		    =~ /^\s*(\d+)\s+(\d+)\s+\d+\s+(\d+)\s+(\d+)\s+(\d+\.\d+)/)
		{
		    # old format
		    my $conserved_region
			    = ORCA::Analysis::ConservedRegion->new(
						-seq1_start	=> $1,
						-seq1_end	=> $2,
						-align_start	=> $3,
						-align_end	=> $4,
						-score		=> $5);
		    #$conserved->add_conserved_region($conserved_region);
		    push @conserved, $conserved_region;
		} elsif ($line
		    =~ /^\s*(\d+)\s+(\d+)\s+\d+\s+(\d+)\s+(\d+)\s+\d+\s+(\d+)\s+(\d+)\s+(\d+\.\d+)/)
		{
		    # new format
		    my $conserved_region
			    = ORCA::Analysis::ConservedRegion->new(
						-seq1_start	=> $1,
						-seq1_end	=> $2,
						-seq2_start	=> $3,
						-seq2_end	=> $4,
						-align_start	=> $5,
						-align_end	=> $6,
						-score		=> $7);
		    #$conserved->add_conserved_region($conserved_region);
		    push @conserved, $conserved_region;
		} elsif ($line =~ /^\s*Window Size\s*:\s*(\d+)/) {
		    $window_size = $1;
		} elsif ($line =~ /^\s*Window Increment\s*:\s*(\d+)/) {
		    $window_inc = $1;
		} elsif ($line =~ /^\s*Top X% Identities\s*:\s*(\d+\.\d+)/) {
		    $top_x_pct = $1;
		} elsif ($line =~ /^\s*Computed % Identity\s*:\s*(\d+\.\d+)/) {
		    $comp_pct_id = $1;
		} elsif ($line =~ /^\s*Minimum % Identity\s*:\s*(\d+\.\d+)/) {
		    $min_pct_id = $1;
		} elsif ($line =~ /^\s*Effective % Identity\s*:\s*(\d+\.\d+)/) {
		    $eff_pct_id = $1;
		}
	    }

	    $report = ORCA::Analysis::ConservationReport->new(
				-report_type		=> 'c',
				-alignment		=> $alignment,
				-filtered_features	=> $features,
				-conserved_regions	=> \@conserved);

	    $report->param("window_size", $window_size) if defined $window_size;
	    $report->param("window_inc", $window_inc) if defined $window_inc;
	    $report->param("top_pct", $top_x_pct) if defined $top_x_pct;
	    $report->param("comp_pct_id", $comp_pct_id) if defined $comp_pct_id;
	    $report->param("min_pct_id", $min_pct_id) if defined $min_pct_id;
	    $report->param("cutoff", $eff_pct_id) if defined $eff_pct_id;
	} elsif ($args{-r} eq 'p') {
	    my @conserved = ();
	    foreach my $line (@ac_output) {
		chomp $line;
		if ($line =~ /^\s*(\d+)\s+(\d+\.\d+)/)
		{
		    push @conserved, {
					-position	=> $1,
					-score		=> $2};
		}
	    }

	    $report = ORCA::Analysis::ConservationReport->new(
				-report_type		=> 'p',
				-alignment		=> $alignment,
				-conservation_profile	=> \@conserved);
	    $report->param("window_size", $args{-w}) if defined $args{-w};
	    $report->param("window_inc", $args{-n}) if defined $args{-n};
	}

	# only clean up if align_cons ran successfully DJA 2006/02/28
	$self->_cleanup;
    } else {
	carp "Error running $cmd\n" . join('', @ac_output) . "\n";
    }

    return $report;
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

#
# Write list of features to a GFF formatted file.
#
sub _write_features_as_GFF
{
    my ($file_name, $features) = @_;

    return if !$file_name;
    return if !$features;

    if (!open(GFF, ">$file_name")) {
	carp "Error creating temporary GFF file\n";
	return;
    }

    foreach my $feat (@$features) {
	print GFF $feat->gff_string . "\n";
    }

    close(GFF);

    return 1;
}

sub _cleanup
{
    my $self = shift;

    my $tempdir = $self->{_tempdir};
    if ($tempdir && -d $tempdir) {
    	unlink(<$tempdir/*>);
    }
}

sub DESTROY
{
    my $self = shift;

    my $tempdir = $self->{_tempdir};
    if ($tempdir && -d $tempdir) {
    	rmtree($tempdir);
    }
}

1;
