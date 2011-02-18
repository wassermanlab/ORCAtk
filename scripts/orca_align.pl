#!/usr/local/bin/perl -w

=head1 NAME

  orca_align.pl 

=head1 SYNOPSIS

  orca_align.pl -in1 seq_file1 -in2 seq_file2 [-out align_file]
  			[other options]

=head1 ARGUMENTS

  -in1		- input FastA file containing first sequence to be
    		  aligned
  -in2		- input FastA file containing second sequence to be
    		  aligned
  -out		- output file in which to write alignment
  -rm1		- sequence 1 repeat masking species; causes repeat
    		  masking to be performed against sequence 1 with this
    		  species
  -rm2		- sequence 2 repeat masking species; causes repeat
    		  masking to be performed against sequence 2 with this
    		  species
  -maxlen	- 'maximum' length of a sequence to be passed to global
    		  alignment (N-W) program - in fact as long as the
    		  product of the two sequence lengths do not exceed the
    		  square of this maximum, N-W will be called
  -format	- output format (if specified, should be either FastA or
    		  ClustalW)

  Global Alignment (Needleman-Wunsch) options.
  -gms		- nucleotide match score
  -gmms		- nucleotide mismatch score
  -ggop		- gap open penalty
  -ggxp		- gap extension penalty

  Local Alignment (Blastn) options.
  -gapped	- turn on/off gapped blast (for determining HSPs)
  -lms		- match score
  -lmms		- mismatch score
  -lmmi		- mismatch score increment (for recursive BLAST
		  anchoring)
  -lgop		- gap open penalty
  -lgxp		- gap extension penalty
  -xdrop	- X dropoff value in bits
  -expect	- expectation value
  -wordsize	- word size
  -wsd		- word size decrement (for recursive BLAST anchoring)

  Other options.
  -debug	- Turn on debugging output

=head1 DESCRIPTION

=over 2

Wrapper script for running the ORCA pair-wise DNA aligner from the command
line.

=back

=head1 REQUIREMENTS

=over 2

At least version 1.2.2 of Bioperl. A local copy of NCBI BLAST. If repeat
masking is specified, requires a local copy of RepeatMasker with the WU-Blast
option enabled.

=back

=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use Getopt::Long;
use IO::File;
use File::Spec::Functions;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Bio::AlignIO;
use Bio::Tools::Run::RepeatMasker;
use ORCA::Aligner;

#
# Some global constants
#
my $PROG = 'orca_align';	# The name of this program
my $VERSION = '1.0.6';		# Version number of this program

my $RM_PROG = 'RepeatMasker';	# Repeat masker program.
				# If the path is not given assume it's in
				# the users path (i.e. /usr/local/bin).
#
# 06/08/03
# Repeat masker options.
# The -s, -no_is and -frag 1000000 options are recommended in the maskeraid
# README file.
#
my $RM_OPTS = '-w -s -no_is -frag 1000000';

#
# ASSUMPTION:
# There appears to be no way to specify the output file name for
# RepeatMasker and no easy way to determine that output file name.
# At the time this was coded, RepeatMasker's file naming convention
# was to name the output file the same as the input file with
# the extension '.masked' appended. The assumption is made that is
# still valid.
#
my $RM_EXT = 'masked';

my $DFLT_FORMAT = 'fasta';	# Default output format. Takes any format
				# recognized by bioperl's AlignIO module.
				# If not defined, outputs 'raw' format which
				# is a simplified ClustalW-like format without
				# a header or sequence ID labels.

my @SUPPORTED_FORMATS = (	# Theoretically, any alignment format supported
				# by bioperl's AlignIO module should work, but
	'fasta',		# ORCA has only been tested with these formats
	'clustalw'		# so restrict output formats to these.
);

my $InFile1;
my $InFile2;
my $OutFile;
my $RepeatMask1;
my $RepeatMask2;
my $GAMaxLen;
my $OutputFormat;
my $GAMatchScore;
my $GAMismatchScore;
my $GAGapPenalty;
my $GAGapExtPenalty;
my $BLGapped;
my $BLExpect;
my $BLGapPenalty;
my $BLGapExtPenalty;
my $BLXDropoff;
my $BLMatchScore;
my $BLMismatchScore;
my $BLMismatchInc;
my $BLWordSize;
my $BLWordSizeDec;
my $Debug;
GetOptions(
	'debug'		=> \$Debug,
	'in1=s'		=> \$InFile1,
	'in2=s'		=> \$InFile2,
	'out=s'		=> \$OutFile,
	'rm1=s'		=> \$RepeatMask1,
	'rm2=s'		=> \$RepeatMask2,
	'format=s'	=> \$OutputFormat,
	'maxlen=i'	=> \$GAMaxLen,
	'gms=i'		=> \$GAMatchScore,
	'gmms=i'	=> \$GAMismatchScore,
	'ggop=i'	=> \$GAGapPenalty,
	'ggxp=i'	=> \$GAGapExtPenalty,
	'gapped=s'	=> \$BLGapped,
	'lms=i'		=> \$BLMatchScore,
	'lmms=i'	=> \$BLMismatchScore,
	'lmmi=i'	=> \$BLMismatchInc,
	'lgop=i'	=> \$BLGapPenalty,
	'lgxp=i'	=> \$BLGapExtPenalty,
	'expect=f'	=> \$BLExpect,
	'xdrop=f'	=> \$BLXDropoff,
	'wordsize=i'	=> \$BLWordSize,
	'wsd=i'		=> \$BLWordSizeDec);

#
# FastA input files containing sequences to be aligned
#
if (!defined $InFile1 || !defined $InFile2) {
	usage();
	exit;
}

$OutputFormat = $DFLT_FORMAT if !defined $OutputFormat;

#
# Do some initial sanity checks
#
ProgramCheck();

my $seq1 = ReadSequence($InFile1);
die "Error reading sequence 1 from FastA file $InFile1\n" if !defined $seq1;

my $seq2 = ReadSequence($InFile2);
die "Error reading sequence 2 from FastA file $InFile2\n" if !defined $seq2;

my $seq1_to_align = $seq1;
my $seq2_to_align = $seq2;

my $seq1_masked;
if (defined $RepeatMask1) {
	$seq1_masked = RepeatMask($seq1, $RepeatMask1);
	if (!defined $seq1_masked) {
		warn "Error repeat masking sequence 1:\n"
			. "file=$InFile1, species=$RepeatMask1\n"
			. " Using unmasked sequence.";
	} else {
		$seq1_to_align = $seq1_masked;
	}
}

my $seq2_masked;
if (defined $RepeatMask2) {
	my $seq2_masked = RepeatMask($seq2, $RepeatMask2);
	if (!defined $seq2_masked) {
		warn "Error repeat masking sequence 2:\n"
			. "file=$InFile2, species=$RepeatMask2\n"
			. " Using unmasked sequence.";
	} else {
		$seq2_to_align = $seq2_masked;
	}
}

my %params;
$params{'-debug'} = $Debug if defined $Debug;
$params{'-ga_max_len'} = $GAMaxLen if defined $GAMaxLen;
$params{'-ga_match'} = $GAMatchScore if defined $GAMatchScore;
$params{'-ga_mismatch'} = $GAMismatchScore if defined $GAMismatchScore;
$params{'-ga_gap_penalty'} = $GAGapPenalty if defined $GAGapPenalty;
$params{'-ga_gap_ext_penalty'} = $GAGapExtPenalty if defined $GAGapExtPenalty;
$params{'-bl_gapped'} = $BLGapped if defined $BLGapped;
$params{'-bl_match'} = $BLMatchScore if defined $BLMatchScore;
$params{'-bl_mismatch'} = $BLMismatchScore if defined $BLMismatchScore;
$params{'-bl_mismatch_inc'} = $BLMismatchInc if defined $BLMismatchInc;
$params{'-bl_gap_penalty'} = $BLGapPenalty if defined $BLGapPenalty;
$params{'-bl_gap_ext_penalty'} = $BLGapExtPenalty if defined $BLGapExtPenalty;
$params{'-bl_x_dropoff'} = $BLXDropoff if defined $BLXDropoff;
$params{'-bl_expect'} = $BLExpect if defined $BLExpect;
$params{'-bl_wordsize'} = $BLWordSize if defined $BLWordSize;
$params{'-bl_wordsize_dec'} = $BLWordSizeDec if defined $BLWordSizeDec;

my $aligner = ORCA::Aligner->new(%params);
if (!defined $aligner) {
	die "Error creating new ORCA::Aligner factory object\n";
}

my $alignment = $aligner->Align(
		'-seq1'	=> $seq1_to_align,
		'-seq2'	=> $seq2_to_align);
if (!defined $alignment) {
	die "Unable to compute alignment\n";
}

#
# Debugging only
#
$aligner->PrintAlignment() if $Debug;

OutputAlignment($alignment, $OutputFormat, $OutFile);

exit;

#
# Perform some up front checking on various program parameters, resources
# etc.
#
sub
ProgramCheck
{
	#if ($RepeatMask1 || $RepeatMask2) {
	#	unless (CheckExecutable($RM_PROG)) {
	#		die "Repeat masking program $RM_PROG not found or not"
	#			. " executable\n"
	#	}
	#}

	my $formatSupported = 0;
	foreach my $format (@SUPPORTED_FORMATS) {
		if (lc($OutputFormat) eq lc($format)) {
			$formatSupported = 1;
			last;
		}
	}
	if (!$formatSupported) {
		die "The output format specified, $OutputFormat, is not"
			. " supported\n";
	}
}

#
# Check to see if the specified program exists and is executable anywhere
# within the user's path.
#
sub
CheckExecutable
{
	my ($prog) = @_;

	my $executable = 0;

	if (-x $prog) {
		$executable = 1;
	} else {
		my $path = $ENV{PATH};
		my @pathdirs = split(':', $path);
		foreach my $dir (@pathdirs) {
			if (-x "$dir/$prog") {
				$executable = 1;
				last;
			}
		}
	}

	return $executable;
}

#
# Read sequence from a FastA file
#
sub
ReadSequence
{
	my ($file) = @_;

	my $in = Bio::SeqIO->new(
			'-file'		=> "$file",
			'-format'	=> 'fasta');

	my $seq = $in->next_seq();

	$in->close;

	return $seq;
}

#
# Read two sequences from a FastA file - no longer used
#
sub
ReadSequences
{
	my ($file) = @_;

	my $in = Bio::SeqIO->new(
			'-file'		=> "$file",
			'-format'	=> 'fasta');

	my $seq1;
	my $seq2;
	$seq1 = $in->next_seq();
	if (defined $seq1) {
		$seq2 = $in->next_seq();
	}

	$in->close;

	return ($seq1, $seq2);
}

#
# Perform repeat masking on the given file with the given species. Create a
# temporary file with the masked sequence and return the name of that file.
#
sub
RepeatMask
{
	my ($seq, $mask) = @_;

	return undef if !defined $seq;
	return undef if !defined $mask;

	my %params = (
	    		w	=> 1,	# use wu-blast
			s	=> 1,	# slow, recommended for wu-blast
			no_is	=> 1,	# recommended for wu-blast
			frag	=> 1000000	# recommended for wu-blast
		    );
	#
	# Default RepeatMasker settings are for primates so if mask
	# looks something like human or primate don't actually pass it.
	#
	if ($mask !~ /^hum|^pri|^homo|^hs/i) {  
	    	$params{$mask} = 1;
	}
	
	my $rm = Bio::Tools::Run::RepeatMasker->new(%params);
	if (!defined $rm) {
	    warn "Error setting up RepeatMasker\n";
	    return undef;
	}

	my $rm_dir = $rm->program_dir;
	my $rm_prog = $rm->program_name;
	if ($rm_dir) {
		if (!-x catfile($rm_dir, $rm_prog)) {
		    warn "Could not locate RepeatMasker or not executable\n";
		    return undef;
		}
	} else {
	    	if (!CheckExecutable($rm_prog)) {
		    warn "Could not locate RepeatMasker or not executable\n";
		    return undef;
		}
	}

	my $masked = undef;
	if ($rm->run($seq)) {
	    $masked = $rm->masked_seq();
	}

	return $masked;
}

#
# Read back the sub-alignments in the temporary alignment files, build the
# complete alignment and then write it out in the format specified.
#
sub
OutputAlignment
{
	my($alignment, $format, $file) = @_;

	return 0 if !defined $alignment;
	return 0 if !defined $format;

	my $stream;
	if (defined $file) {
		$stream = Bio::AlignIO->new('-file' => ">$file",
				       '-format' => $format);
		return 0 if !defined $stream;
	} else {
		#
		# stdout
		#
		$stream = Bio::AlignIO->new('-format' => $format);
		return 0 if !defined $stream;
	}
	
	my $rval = $stream->write_aln($alignment);

	$stream->close;

	return $rval;
}

#
# Make a directory in which to create any temporary files
#
sub
MakeTempDir
{
	my $template = sprintf("%dXXXXXX", $$);

	#
	# CLEANUP => 1 to autonatically removes dir and contents at
	# termination of program execution
	#
	my $tmpDir;
	if ($Debug) {
		$tmpDir = tempdir($template, TMPDIR => 1);
	} else {
		$tmpDir = tempdir($template, TMPDIR => 1, CLEANUP => 1);
	}

	return $tmpDir;
}

sub
usage
{
    #
    # Have to create this to get default values
    #
    my $aligner = ORCA::Aligner->new;

    print "\n$PROG $VERSION options (abbreviations allowed if unique):\n\n";
    print "-in1       input FastA file containing first sequence to"
	    . " align\n";
    print "-in2       input FastA file containing second sequence to"
	    . " align\n";
    print "-out       output sequence alignment file\n";
    print "\t\t(default = stdout)\n";
    print "-rm1       sequence 1 repeat masking species; causes repeat\n";
    print "           masking to be performed against sequence 1 with\n";
    print "           this species\n";
    print "-rm2       sequence 2 repeat masking species; causes repeat\n";
    print "           masking to be performed against sequence 2 with\n";
    print "           this species\n";
    print "-maxlen    'max' length of a sequence to pass to global (N-W)\n";
    print "           alignment - in actuality two sequences will be passed\n";
    print "           to N-W as long as the product of their lengths does\n";
    print "           not exceed maxlen squared\n";
    print "\t\t(default = " . $aligner->GAMaxLen . ")\n";
    print "-format    output format. Supported formats: ";
    printf "(%s)\n", join(', ', @SUPPORTED_FORMATS);
    print "\t\t(default = $DFLT_FORMAT)\n";
    print "-debug     turn on debugging output\n";
    print "\nGlobal Alignment (N-W) Parameters\n";
    print "-gms       nucleotide match score\n";
    print "\t\t(default = " . $aligner->GAMatchScore . ")\n";
    print "-gmms      nucleotide mismatch score\n";
    print "\t\t(default = " . $aligner->GAMismatchScore . ")\n";
    print "-ggop      penalty for opening a gap\n";
    print "\t\t(default = " . $aligner->GAGapOpenPenalty . ")\n";
    print "-ggxp      penalty for extending a gap\n";
    print "\t\t(default = " . $aligner->GAGapExtPenalty . ")\n";
    print "\nLocal Alignment (Blastn) Parameters\n";
    print "-gapped    use gapped alignment\n";
    print "\t\t(default = "
    		. ($aligner->BLGapped || "blast default value") . ")\n";
    print "-lms       nucleotide match score\n";
    print "\t\t(default = "
    		. ($aligner->BLMatchScore || "blast default value") . ")\n";
    print "-lmms      nucleotide mismatch score\n";
    print "\t\t(default = "
    		. ($aligner->BLMismatchScore || "blast default value") . ")\n";
    print "-lmmi      nucleotide mismatch incremement\n";
    print "\t\t(default = "
    		. ($aligner->BLMismatchInc || "blast default value") . ")\n";
    print "-lgop      gap open penalty\n";
    print "\t\t(default = "
    		. ($aligner->BLGapOpenPenalty || "blast default value") . ")\n";
    print "-lgxp      gap extension penalty\n";
    print "\t\t(default = "
    		. ($aligner->BLGapExtPenalty || "blast default value") . ")\n";
    print "-expect    expectation value\n";
    print "\t\t(default = "
    		. ($aligner->BLExpect || "blast default value") . ")\n";
    print "-xdrop     X dropoff value in bits\n";
    print "\t\t(default = "
    		. ($aligner->BLXDropoff || "blast default value") . ")\n";
    print "-wordsize  word size\n";
    print "\t\t(default = "
    		. ($aligner->BLWordSize || "blast default value") . ")\n";
    print "-wsd       word size decrement\n";
    print "\t\t(default = "
    		. ($aligner->BLWordSizeDec || "blast default value") . ")\n\n";
}
