=head1 NAME

ORCA::Aligner - Object for performing pair-wise DNA alignment

=head1 SYNOPSIS

    # ORCA::Aligner "factory object" creation and parameter
    # initialization (please see PARAMETERS section below for a
    # description of valid parameters):

    $orca_aligner = ORCA::Aligner->new(%params)

    # Align two sequences:

    $seqIO = Bio::SeqIO->new(
				'-file'		=> 'sequences.fa',
				'-format'	=> 'fasta');
    $seq1 = $seqIO->next_seq();
    $seq2 = $seqIO->next_seq();
    $aln = $orca_aligner->align(
				'-seq1'	=> $seq1,
				'-seq2'	=> $seq2);

=head1 DESCRIPTION

ORCA::Aligner.pm is a BioPerl compatible object for performing pair-wise DNA
sequence alignments. It calls NCBI BLAST through BioPerl's
Bio::Tools::Run::StandAloneBlast.pm module in order to find HSP's which are
then filtered for colinearity in order to build a list of anchors. It then
calls the external nwalign executable to perform Needleman-Wunsch (NW)
alignments on each of the regions bounded by these anchors. This process is
recursive, such that if any region is too long for NW, it is re-BLASTED using
looser BLAST parameters. This process stops when either the regions are short
enough for a succesful NW call or the BLAST parameters cannot be adjusted
any further at which point the two sub-sequences in the region are simply
placed side-by-side with the shorter padded with gaps at the 3' end.

=head1 PARAMETERS

The optional parameters which may be passed to the ORCA::Aligner factory
creation are as follows:

    Global Alignment (Needleman-Wunsch) parameters:
        -ga_max_len	    - (int) maximum length of a sequence to be
          			    passed to global alignment (NW)
    				    program
        -ga_match	    - (int) nucleotide match score
        -ga_mismatch	    - (int) nucleotide mismatch score
        -ga_gap_penalty	    - (int) gap open penalty
        -ga_gap_ext_penalty - (int) gap extension penalty

    Local Alignment (BLAST) parameter:
        -bl_gapped	    - (1/0) turn on/off gapped blast
        -bl_match	    - (int) match score
        -bl_mismatch	    - (int) mismatch score
        -bl_mismatch_inc    - (int) mismatch score increment
    			            (used for recursive BLAST anchoring)
        -bl_gap_penalty	    - (int) gap open penalty
        -bl_gap_ext_penalty - (int) gap extension penalty
        -bl_x_drop	    - (int) X dropoff value in bits
        -bl_expect	    - (real) expectation value
        -bl_wordsize	    - (int) word size
        -bl_wordsize_dec    - (int) word size decrement
    			            (used for recursive BLAST anchoring)
    Other parameters:
        -debug		    - (1/0) turn on/off debugging messages

=head1 DEPENDENCIES

ORCA::Aligner.pm requires at least version 1.2.2 of BioPerl, a local
installation of NCBI BLAST and an external Needleman-Wunsch global
alignment program 'nwalign' distributed with the ORCA package.

Bioperl Issues:

Unfortunately versions of Bioperl prior to v. 1.2.2 required a fix
which was made to so that bl2seq could tell which type of BLAST was run in
order for strand information to be preserved. The affected modules were
Bio::Tools::BPbl2seq.pm and Bio::AlignIO::bl2seq.pm. This fix is now fully
incorporated into version 1.2.2 of Bioperl.

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=head1 METHODS

=cut

package ORCA::Aligner;
use vars qw($VERSION);
use strict;

use Carp;
use IO::File;
use File::Spec::Functions;
use File::Basename;
use File::Copy;
use File::Path;
use File::Temp qw/ tempfile tempdir /;
use Class::Struct;
use Time::localtime;
use Bio::LocatableSeq;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Tools::BPbl2seq;
use Bio::Tools::Run::StandAloneBlast;

$VERSION = '1.1.1';

#
# This must be in sync with the definition of E_SEQLEN in nwalign.h!!!
#
use constant NW_SEQ_LEN_ERR => 130;

#
# Some package constants
#
use constant DEBUG => 0;

use constant GA_PROG => 'nwalign';	# The global alignment (N-W) program.
					# If the path is not given assume it's
					# in the users path
					# (i.e. /usr/local/bin).

use constant HSP_CHAINING_METHOD => 'LIS';	# must be either 'greedy' or
						# 'LIS' (longest increasing
						# sub-sequence)

use constant DFLT_GA_LEN => 10000;	# Default max. length of seqments to be
					# passed to the global alignment (N-W)
					# program.

use constant GA_MAX_LEN => 16383;	# Max. allowable length of seqments
					# to be passed to global alignment
					# program Note: these are only
					# 'suggestions'. The true maximum is
					# based on the product of the sequence
					# lengths not exeeding the square of
					# these values.

use constant HBS_ANCHOR_LEN => 25;	# No. of nucleotides from bounding HSPs
					# included with seqments sent to global
					# aligment program

use constant BL_PROG => 'blastn';	# Blast program

use constant BL_WORDSIZE_MIN => 7;	# Minimum Blastn word size

use constant BL_MISMATCH_MAX => -3;	# Maximum Blastn mismatch score

use constant BL_EXPECT_MAX => 10;	# Maximum Blastn expect value

#
# Default ORCA parameters for bl2seq (which differ from bl2seq defaults)
#
use constant DFLT_BL_EXPECT => 0.1;	# Expectation value
use constant DFLT_BL_EXPECT_MULT => 10;	# Expectation value multiplier

#use constant DFLT_BL_MISMATCH_SCORE => -7;# Mismatched nucleotide pair score
use constant DFLT_BL_MISMATCH_SCORE => -3;# Mismatched nucleotide pair score
#use constant DFLT_BL_MISMATCH_INC => 2;# Amount by which to increment mismatch
					# score 
use constant DFLT_BL_MISMATCH_INC => 0;	# Amount by which to increment mismatch
					# score 

use constant DFLT_BL_WORDSIZE => 15;	# Word size
use constant DFLT_BL_WORDSIZE_DEC => 4;	# Amount by which to decrememnt window
					# size
use constant DFLT_BL_FILTER => 'F';	# Filter query sequence

#
# Default parameters for global alignment (N-W) program
#
use constant DFLT_GA_MATCH_SCORE => 3;	# Matching nucleotide pair score
use constant DFLT_GA_MISMATCH_SCORE => -1;# Mismatched nucleotide pair score
use constant DFLT_GA_GAP_PENALTY => 20;	# Cost to open a gap
use constant DFLT_GA_GAP_EXT_PENALTY => 0;# Cost to extend a gap

use constant SEQ_OUT_WIDTH => 60;

#
# Structure defining an alignment consisting of a list of alignment segments
# which may be either a High Scoring Pair (HSP) produced by blast, an
# HSP bounded segment (HBS) aligned by an external global alignment program
# (Needleman-Wunsch), or a sub-alignment list. Thus this is a recursive
# structure.
#
struct Alignment => {
	seq1		=> '$',	# references to bioperl sequence objects
	seq2		=> '$',	# to be aligned
	blastParams	=> '$',	# reference to paramaters passed to blast
	blastReport	=> '$',	# blast report
	alignList	=> '$'	# reference to the list of alignment segments
};

#
# Structure defining a segment of the alignment which may be either of type
# High Scoring Pair (HSP) produced by blast, HSP bounded segment (HBS)
# aligned by an external global alignment program (Needleman-Wunsch),
# or a sub-alignment list (SAL).
#
struct AlignSegment => {
	type	=> '$',		# segment type (HSP, HBS or SAL)
	prev	=> '$',		# pointer to previous (left) segment
	next	=> '$',		# pointer to next (right) segment
	data	=> '$',		# pointer to HSP, HBR or sub-alignment
	file	=> '$'		# file containing alignment (for HSP or HBS)
};

#
# Structure defining an HSP bounded segment (HBS) which may be alignable
# by the external global alignment (N-W) program or not.
#
struct HBS => {
	alignable	=> '$',	# indicates alignability
	seg1Start	=> '$',	# start position on sequence 1
	seg1End		=> '$',	# end position on sequence 1
	seg2Start	=> '$',	# start position on sequence 2
	seg2End		=> '$',	# end position on sequence 2
	seg1LeftOL	=> '$',	# amount segment 1 overlaps previous HSP
	seg1RightOL	=> '$',	# amount segment 1 overlaps next HSP
	seg2LeftOL	=> '$',	# amount segment 2 overlaps previous HSP
	seg2RightOL	=> '$'	# amount segment 2 overlaps next HSP
};

#
# Structure defining Blast parameters
#
struct BlastParams => {
	program		=> '$',
	gapped		=> '$',
	expect		=> '$',
	gapPenalty	=> '$',
	gapExtPenalty	=> '$',
	xDropoff	=> '$',
	matchScore	=> '$',
	mismatchScore	=> '$',
	wordSize	=> '$',
	filter		=> '$'
};

my $DEBUG = DEBUG;

=head2 new

 Title	  : new
 Usage	  : $orca_aligner = ORCA::Aligner->new(%params);
 Function : Returns a new ORCA::Aligner factory object.
 Returns  : An ORCA::Aligner factory object initialized with optional
 	    parameters.
 Args	  : Optional named parameters:
	      -ga_max_len	  - (int) maximum length of a sequence to be
	      				  passed to global alignment (NW)
					  program
	      -ga_match		  - (int) nucleotide match score
	      -ga_mismatch	  - (int) nucleotide mismatch score
	      -ga_gap_penalty	  - (int) gap open penalty
	      -ga_gap_ext_penalty - (int) gap extension penalty

	      -bl_gapped	  - (1/0) turn on/off gapped blast
	      -bl_match		  - (int) match score
	      -bl_mismatch	  - (int) mismatch score
	      -bl_mismatch_inc 	  - (int) mismatch score increment
				          (used for recursive BLAST anchoring)
	      -bl_gap_penalty	  - (int) gap open penalty
	      -bl_gap_ext_penalty - (int) gap extension penalty
	      -bl_x_drop	  - (int) X dropoff value in bits
	      -bl_expect	  - (real) expectation value
	      -bl_wordsize	  - (int) word size
	      -bl_wordsize_dec 	  - (int) word size decrement
				          (used for recursive BLAST anchoring)
	      -bl_filter 	  - (T/F) filter query sequence

	      -debug		  - (1/0) turn on/off debugging messages
		
=cut

sub new
{
    my ($proto, %args) = @_;
    my $class = ref($proto) || $proto;

    my $self = {};

    #
    # Use current date information for creating temporary working
    # directory name.
    #
    my $tm = localtime;
    $self->{year} = 1900 + $tm->year;
    $self->{month} = $tm->mon + 1;
    $self->{day} = $tm->mday;

    #
    # Keep track of orientation of alignment: 1 = plus/plus, 0 = plus/minus
    #
    $self->{orient} = undef;

    bless ($self, $class);

    if (defined $args{-debug}) {
	$DEBUG = $args{-debug};
    }

    #
    # Set global alignment parameters.
    #
    if (defined $args{-ga_max_len}) {
	$self->GAMaxLen($args{-ga_max_len});
    } else {
	$self->GAMaxLen(DFLT_GA_LEN);
    }
    if (defined $args{-g_match}) {
	$self->GAMatchScore($args{-g_match});
    } else {
	$self->GAMatchScore(DFLT_GA_MATCH_SCORE);
    }
    if (defined $args{-g_mismatch}) {
	$self->GAMismatchScore($args{-g_mismatch});
    } else {
	$self->GAMismatchScore(DFLT_GA_MISMATCH_SCORE);
    }
    if (defined $args{-g_gap_penalty}) {
	$self->GAGapOpenPenalty($args{-g_gap_penalty});
    } else {
	$self->GAGapOpenPenalty(DFLT_GA_GAP_PENALTY);
    }
    if (defined $args{-g_gap_ext_penalty}) {
	$self->GAGapExtPenalty($args{-g_gap_ext_penalty});
    } else {
	$self->GAGapExtPenalty(DFLT_GA_GAP_EXT_PENALTY);
    }

    #
    # Set local alignment parameters.
    #
    if (defined $args{-bl_gapped}) {
	$self->BLGapped($args{-bl_gapped});
    }
    if (defined $args{-bl_match}) {
	$self->BLMatchScore($args{-bl_match});
    }
    if (defined $args{-bl_mismatch}) {
	$self->BLMismatchScore($args{-bl_mismatch});
    } else {
	$self->BLMismatchScore(DFLT_BL_MISMATCH_SCORE);
    }
    if (defined $args{-bl_mismatch_inc}) {
	$self->BLMismatchInc($args{-bl_mismatch_inc});
    } else {
	$self->BLMismatchInc(DFLT_BL_MISMATCH_INC);
    }
    if (defined $args{-bl_gap_penalty}) {
	$self->BLGapOpenPenalty($args{-bl_gap_penalty});
    }
    if (defined $args{-bl_gap_ext_penalty}) {
	$self->BLGapExtPenalty($args{-bl_gap_ext_penalty});
    }
    if (defined $args{-bl_x_drop}) {
	$self->BLXDropoff($args{-bl_x_drop});
    }
    if (defined $args{-bl_expect}) {
	$self->BLExpect($args{-bl_expect});
    } else {
	$self->BLExpect(DFLT_BL_EXPECT);
    }
    if (defined $args{-bl_expect_mult}) {
	$self->BLExpectMult($args{-bl_expect_mult});
    } else {
	$self->BLExpectMult(DFLT_BL_EXPECT_MULT);
    }
    if (defined $args{-bl_wordsize}) {
	$self->BLWordSize($args{-bl_wordsize});
    } else {
	$self->BLWordSize(DFLT_BL_WORDSIZE);
    }
    if (defined $args{-bl_wordsize_dec}) {
	$self->BLWordSizeDec($args{-bl_wordsize_dec});
    } else {
	$self->BLWordSizeDec(DFLT_BL_WORDSIZE_DEC);
    }
    if (defined $args{-bl_filter}) {
	$self->BLFilter($args{-bl_filter});
    } else {
	$self->BLFilter(DFLT_BL_FILTER);
    }

    #
    # Check if BLAST is available.
    #
    if (!_CheckBlast()) {
	return undef;
    }

    #
    # Set the external global alignment (N-W) program.
    #
    if (!$self->_GAProg(GA_PROG)) {
	return undef;
    }
    
    #
    # Make a temporary working directory.
    #
    if (!$self->_MakeTempDir()) {
	return undef;
    }

    return $self;
}

=head2 Align

 Title	  : Align
 Usage	  : $aln = $orca_aligner->Align('-seq1' => $seq1,
					'-seq2' => $seq2);
 Function : Perform pair-wise DNA alignment on seq1 and seq2 and return the
 	    alignment as a Bio::SimpleAlign object.
 Returns  : a Bio::SimpleAlign object
 Args	  : two Bio::Seq objects

=cut

sub align
{
    my ($self, %args) = @_;

    return $self->Align(%args);
}

sub Align
{
    my ($self, %args) = @_;

    #
    # Initialize.
    #
    $self->{orient} = undef;
    $self->{_alignment} = undef;
    $self->{_alignFiles} = undef;

    my $seq1 = $args{-seq1};
    unless (defined $seq1) {
	carp "Sequence 1 not specified\n";
	return undef;
    }

    my $seq2 = $args{-seq2};
    unless (defined $seq2) {
	carp "Sequence 2 not specified\n";
	return undef;
    }

    my $blastParams = BlastParamsCreate(
			    $self->BLGapped,
			    $self->BLExpect,
			    $self->BLGapOpenPenalty,
			    $self->BLGapExtPenalty,
			    $self->BLXDropoff,
			    $self->BLMatchScore,
			    $self->BLMismatchScore,
			    $self->BLWordSize,
			    $self->BLFilter);
    return undef if !defined $blastParams;

    print "Building top-level alignment\n" if $DEBUG;

    my $alignment = $self->_AlignmentBuild($seq1, $seq2, $blastParams);
    return undef if !defined $alignment;
    $self->{_alignment} = $alignment;

    my $alignFiles = $self->_CreateAlignmentFiles($alignment);
    return undef if !defined $alignFiles;
    $self->{_alignFiles} = $alignFiles;

    my $totalAlign = $self->_CreateTotalAlignment();

    $self->_Cleanup();

    return $totalAlign;
}

#
# Used for debugging purposes only
#
sub PrintAlignment
{
    my $self = shift;

    my $alignment = $self->{_alignment};

    if (defined $alignment) {
	printf("Orientation %s\n",
	    $self->{orient} == 1 ? 'Plus/Plus'
	    	: ($self->{orient} == -1 ? 'Plus/Minus' : '?/?'));
	AlignmentPrint($alignment);
    }
}

################################################################################
# Get/Set methods
################################################################################
sub debug
{
    my $self = shift;

    if (@_) {
	$DEBUG = shift;
    }

    return $DEBUG;
}

sub GAMaxLen
{
    my ($self, $maxLen) = @_;

    if (defined $maxLen) {
	if ($maxLen > GA_MAX_LEN || $maxLen < 0) {
	    croak "Maximum global alignment length of $maxLen is out of range.";
	}
	$self->{ga_max_len} = $maxLen;
    }

    return $self->{ga_max_len};
}

sub GAMatchScore
{
    my ($self, $matchScore) = @_;

    if (defined $matchScore) {
	if ($matchScore < 0) {
	    croak "Match score of $matchScore is out of range.";
	}
	$self->{ga_match} = $matchScore;
    }

    return $self->{ga_match};
}

sub GAMismatchScore
{
    my ($self, $mismatchScore) = @_;

    if (defined $mismatchScore) {
	if ($mismatchScore > 0) {
	    croak "Mismatch score of $mismatchScore is out of range.";
	}
	$self->{ga_mismatch} = $mismatchScore;
    }

    return $self->{ga_mismatch};
}

sub GAGapOpenPenalty
{
    my ($self, $penalty) = @_;

    if (defined $penalty) {
	if ($penalty < 0) {
	    croak "Gap open penalty of $penalty is out of range.";
	}
	$self->{ga_gap_penalty} = $penalty;
    }

    return $self->{ga_gap_penalty};
}

sub GAGapExtPenalty
{
    my ($self, $penalty) = @_;

    if (defined $penalty) {
	if ($penalty < 0) {
	    croak "Gap extension penalty of $penalty is out of range.";
	}
	$self->{ga_gap_ext_penalty} = $penalty;
    }

    return $self->{ga_gap_ext_penalty};
}

sub BLGapped
{
    my ($self, $gapped) = @_;

    if (defined $gapped) {
	$self->{bl_gapped} = $gapped;
    }

    return $self->{bl_gapped};
}

sub BLMatchScore
{
    my ($self, $matchScore) = @_;

    if (defined $matchScore) {
	if ($matchScore < 0) {
	    croak "Match score of $matchScore is out of range.";
	}
	$self->{bl_match} = $matchScore;
    }

    return $self->{bl_match};
}

sub BLMismatchScore
{
    my ($self, $mismatchScore) = @_;

    if (defined $mismatchScore) {
	if ($mismatchScore > BL_MISMATCH_MAX) {
	    croak "Mismatch score of $mismatchScore is out of range.";
	}
	$self->{bl_mismatch} = $mismatchScore;
    }

    return $self->{bl_mismatch};
}

sub BLMismatchInc
{
    my ($self, $mismatchInc) = @_;

    if (defined $mismatchInc) {
	if ($mismatchInc < 0) {
	    croak "Mismatch increment of $mismatchInc is out of range.";
	}
	$self->{bl_mismatch_inc} = $mismatchInc;
    }

    return $self->{bl_mismatch_inc};
}

sub BLGapOpenPenalty
{
    my ($self, $penalty) = @_;

    if (defined $penalty) {
	if ($penalty < 0) {
	    croak "Gap open penalty of $penalty is out of range.";
	}
	$self->{bl_gap_penalty} = $penalty;
    }

    return $self->{bl_gap_penalty};
}

sub BLGapExtPenalty
{
    my ($self, $penalty) = @_;

    if (defined $penalty) {
	if ($penalty < 0) {
	    croak "Gap extension penalty of $penalty is out of range.";
	}
	$self->{bl_gap_ext_penalty} = $penalty;
    }

    return $self->{bl_gap_ext_penalty};
}

sub BLXDropoff
{
    my ($self, $dropoff) = @_;

    if (defined $dropoff) {
	if ($dropoff < 0) {
	    croak "X dropoff value of $dropoff is out of range.";
	}
	$self->{bl_x_dropoff} = $dropoff;
    }

    return $self->{bl_x_dropoff};
}

sub BLExpect
{
    my ($self, $expect) = @_;

    if (defined $expect) {
	if ($expect <= 0) {
	    croak "Expectation value of $expect is out of range.";
	}
	$self->{bl_expect} = $expect;
    }

    return $self->{bl_expect};
}

sub BLExpectMult
{
    my ($self, $mult) = @_;

    if (defined $mult) {
	if ($mult <= 0) {
	    croak "Expectation value multiplier of $mult is out of range.";
	}
	$self->{bl_expect_mult} = $mult;
    }

    return $self->{bl_expect_mult};
}

sub BLWordSize
{
    my ($self, $wordsize) = @_;

    if (defined $wordsize) {
	if ($wordsize < BL_WORDSIZE_MIN) {
	    croak "Word size value of $wordsize is out of range.";
	}
	$self->{bl_wordsize} = $wordsize;
    }

    return $self->{bl_wordsize};
}

sub BLWordSizeDec
{
    my ($self, $wordsizeDec) = @_;

    if (defined $wordsizeDec) {
	if ($wordsizeDec < 0) {
	    croak "Word size decrement value of $wordsizeDec is out of range.";
	}
	$self->{bl_wordsize_dec} = $wordsizeDec;
    }

    return $self->{bl_wordsize_dec};
}

sub BLFilter
{
    my ($self, $filter) = @_;

    if (defined $filter) {
	if ($filter !~ /^T$/i && $filter !~ /^F$/i) {
	    croak "Filter value must be T or F";
	}
	$self->{bl_filter} = $filter;
    }

    return $self->{bl_filter};
}

################################################################################
# Private Methods
#
# The following methods are 'private' and should not be called externally.
################################################################################
sub _AlignmentBuild
{
    my ($self, $seq1, $seq2, $blastParams) = @_;

    return undef if !defined $seq1;
    return undef if !defined $seq2;
    return undef if !defined $blastParams;

    my $alignment = Alignment->new();
    return undef if !defined $alignment;

    $alignment->seq1($seq1);
    $alignment->seq2($seq2);
    $alignment->blastParams($blastParams);

    if (!$self->_AlignmentBuildSegments($alignment)) {
        carp "Could not build alignment segments\n";
	return undef;
    }

    return $alignment;
}

sub _CreateAlignmentFiles
{
    my ($self, $alignment) = @_;

    return undef if !defined $alignment;

    my @alignFiles;

    my $curSegment = $alignment->alignList;
    while ($curSegment) {
	if ($curSegment->type eq 'SAL') {
	    my $subAlignFiles = $self->_CreateAlignmentFiles($curSegment->data);
	    if (!defined $subAlignFiles) {
		return undef;
	    }
	    push @alignFiles, @$subAlignFiles;
	} else {
	    if (!$self->_CreateSegmentFile($alignment, $curSegment)) {
		return undef;
	    }
	    push @alignFiles, $curSegment->file;
	}
	$curSegment = $curSegment->next;
    }

    return \@alignFiles;
}

sub _CreateTotalAlignment
{
    my ($self) = @_;

    my $alignment = $self->{_alignment};
    return undef if !defined $alignment;

    my $alignFiles = $self->{_alignFiles};
    return undef if !defined $alignFiles;

    #
    # Build up complete sequences from temporary alignment files.
    #
    my $seq_str1 = '';
    my $seq_str2 = '';
    foreach my $file (@$alignFiles) {
	return 0 if !open(INFH, "$file");

	#
	# Append sequences from the current input file.
	#
	my $lineCount = 0;
	my $line;
	while (defined ($line = <INFH>)) {
	    chomp $line;
	    my $seqNum = $lineCount % 4;
	    if ($seqNum == 0) {
		$seq_str1 .= uc $line;
	    } elsif ($seqNum == 1) {
		$seq_str2 .= uc $line;
	    }
	    $lineCount++;
	}
	close(INFH);
    }

    my $seq_id1 = $alignment->seq1->id || 'seq1';
    my $strand1 = 0;
    my $seq_start1 = 1;
    my $seq_end1 = $alignment->seq1->length;
    #
    # If input sequence is a Bio::LocatableSeq, preserve the location info
    #
    if ($alignment->seq1->isa('Bio::LocatableSeq')) {
	$seq_start1 = $alignment->seq1->start || 1;
	$seq_end1 = $alignment->seq1->end
			|| $alignment->seq1->length - $seq_start1 + 1;
	$strand1 = $alignment->seq1->strand;
	if ($strand1) {
	    if ($strand1 eq '+' || $strand1 eq '1') {
		$strand1 = 1;
	    } elsif ($strand1 eq '-' || $strand1 eq '-1') {
		$strand1 = -1;
	    } else {
		$strand1 = 0;
	    }
	}
    }

    my $seq_id2	= $alignment->seq2->id || 'seq2';
    my $strand2 = 0;
    my $seq_start2 = 1;
    my $seq_end2 = $alignment->seq2->length;
    #
    # If input sequence is a Bio::LocatableSeq, preserve the location info
    #
    if ($alignment->seq2->isa('Bio::LocatableSeq')) {
	$seq_start2 = $alignment->seq2->start || 1;
	$seq_end2 = $alignment->seq2->end
			|| $alignment->seq2->length - $seq_start2 + 1;
	$strand2 = $alignment->seq2->strand;
	if ($strand2) {
	    if ($strand2 eq '+' || $strand2 eq '1') {
		$strand2 = 1;
	    } elsif ($strand2 eq '-' || $strand2 eq '-1') {
		$strand2 = -1;
	    } else {
		$strand2 = 0;
	    }
	}
    }

    #
    # Set strand of aligned sequences based on strand of input sequences 
    # (if set) and alignment orientation
    #
    if ($self->{orient} == 1) {
	if ($strand1) {
	    if (!$strand2) {
	    	$strand2 = $strand1;
	    }
	} else {
	    if ($strand2) {
	    	$strand1 = $strand2;
	    } else {
	    	$strand1 = 1;
	    	$strand2 = 1;
	    }
	}
    } elsif ($self->{orient} == -1) {
	if ($strand1) {
	    if (!$strand2) {
	    	$strand2 = $strand1 * -1;
	    } else {
	    	$strand2 *= -1;
	    }
	} else {
	    if ($strand2) {
	    	$strand1 = $strand2 * -1;
	    } else {
	    	$strand1 = 1;
	    	$strand2 = -1;
	    }
	}
    }

    #
    # Create locatable sequence objects.
    #
    my $aln_seq1 = Bio::LocatableSeq->new(
			-alphabet	=> "dna",
			-seq		=> $seq_str1,
			-id		=> $seq_id1,
			-start		=> $seq_start1,
			-end		=> $seq_end1,
			-strand		=> $strand1);

    #
    # This seems to confuse the Bio::SimpleAlign::displayname method!!!!!!!
    #
    #if ($self->{orient} == -1) {
    #	$seq_id2 .= " (ORCA reverse complemented)";
    #}
    my $aln_seq2 = Bio::LocatableSeq->new(
			-alphabet	=> "dna",
			-seq		=> $seq_str2,
			-id		=> $seq_id2,
			-start		=> $seq_start2,
			-end		=> $seq_end2,
			-strand		=> $strand2);

    my $aln = new Bio::SimpleAlign();
    $aln->add_seq($aln_seq1);
    $aln->add_seq($aln_seq2);

    return $aln;
}

sub _AlignmentBuildSegments
{
    my ($self, $alignment) = @_;

    return 0 if !defined $alignment;

    if (!$self->_AlignmentAddHSPs($alignment)) {
	#
	# If we couldn't get any HSPs, then unless we are at the top
	# level of the recursion, continue anyway to try to let the
	# global alignment portion compute the alignment. If at the top
	# level however, fail since there is no way to really determine
	# the optimal orientation of the alignment. Also, previously
	# when the alignment process was allowed to continue at this
	# point it (almost) always resulted in a meaningless alignment
	# in which the entire alignment was treated as a big
	# non-alignable segment. That is, a meaningless alignment was
	# created.
	#
	if (!defined $self->{orient}) {
	    return 0;
	}
    }

    if (!$self->_AlignmentAddHBSs($alignment)) {
	carp "Error adding HBSs to alignment\n";
	return 0;
    }

    return 1;
}

#
# Examine high scoring pairs (HSPs). Build a chain of HSPs which are
# colinear and retain the same orientation (+/+ or +/-) and add this
# chain to the alignment.
#
sub _AlignmentAddHSPs
{
    my ($self, $alignment) = @_;

    return 0 if !defined $alignment;

    if (HSP_CHAINING_METHOD eq 'greedy') {
	return $self->_AlignmentAddHSPsGreedy($alignment);
    } elsif (HSP_CHAINING_METHOD eq 'LIS') {
	return $self->_AlignmentAddHSPsLIS($alignment);
    }

    carp "Unknown HSP method\n";
    return 0;
}

#
# Greedy HSP chaining method. The first HSP (by definition highest scoring)
# is added to the chain. Then all other HSPs are examined in turn (highest
# to lowest scoring) and those which retain the same orientation as the initial
# HSP and which do not violate colinearity are added to the chain. This is
# simpler and more computationally efficient than the LIS method but may not
# result in the best chain of HSPs.
#
sub _AlignmentAddHSPsGreedy
{
    my ($self, $alignment) = @_;

    my $seq1 = $alignment->seq1;
    return 0 if !defined $seq1;

    my $seq2 = $alignment->seq2;
    return 0 if !defined $seq2;

    my $blastParams = $alignment->blastParams;
    return 0 if !defined $blastParams;

    my $report = undef;
    my $firstHSP = undef;
    while (defined $blastParams && !defined $firstHSP) {
	$report = $self->_RunBlast($seq1, $seq2, $blastParams);
	if (!defined $report) {
	    croak "Error running BLAST\n";
	}
	
	#
	# Get initial HSP
	#
	# At time of coding, bioperl throws
	# an exception if there are no HSPs
	# in the bl2seq report, so trap it.
	#
	eval {
	    $firstHSP = $report->next_feature;
	};

	if (defined $firstHSP) {
	    #
	    # NOTE: Blast orders HSPs from highest to lowest score.
	    # For now always add the first HSP. The result is that
	    # this first HSP determines the alignment direction
	    # (and establishes initial colinearity) for the other
	    # HSPs. The underlying assumption that this will result
	    # in the best alignment may have to be questioned down
	    # the road. If the orientation is already set
	    # (i.e. this is a subalignment), start with the first
	    # HSP which matches the existing orientation.
	    # 
	    my $orient = HSPOrient($firstHSP);
	    if (!defined $self->{orient}) {
		$self->{orient} = $orient;
	    } else {
		while ($orient != $self->{orient} && defined $firstHSP) {
		    #
		    # At time of coding, bioperl throws
		    # an exception if there are no HSPs
		    # in the bl2seq report, so trap it.
		    #
		    eval {
			    $firstHSP =
				    $report->next_feature;
		    };
		    $orient = HSPOrient($firstHSP)
			    if defined $firstHSP;
		}
	    }
	}

	if (!defined $firstHSP) {
	    $blastParams = $self->_BlastParamsAdjust($blastParams);
	}
    }

    if (!defined $firstHSP) {
	$report->close if defined $report;
	return 0;
    }

    $alignment->blastParams($blastParams);
    #$alignment->blastReport($report);

    my $hspSeg = AlignSegmentCreate('HSP', $firstHSP);
    if (!defined $hspSeg) {
	$report->close;
	return 0;
    }

    my $hspList = $hspSeg;
    while (my $hsp = $report->next_feature) {
	if (HSPOrient($hsp) == $self->{orient}) {
	    $hspList = HSPListInsert($hspList, $hsp, $self->{orient});
	}
    }
    $report->close;

    $alignment->alignList($hspList);

    return 1;
}

#
# Longest increasing sub-sequence HSP chaining method. The best HSP chain is
# found by examining all possible chains and returning the highest scoring.
# This is more computationally expensive than the greedy method but will find
# the best chain.
#
sub _AlignmentAddHSPsLIS
{
    my ($self, $alignment) = @_;

    my $seq1 = $alignment->seq1;
    return 0 if !defined $seq1;

    my $seq2 = $alignment->seq2;
    return 0 if !defined $seq2;

    my $blastParams = $alignment->blastParams;
    return 0 if !defined $blastParams;

    my $report = undef;
    my $hsp_found = 0;
    my $plus_hsp_count = 0;
    my $minus_hsp_count = 0;
    my @plus_hsps;
    my @minus_hsps;
    while (defined $blastParams && !$hsp_found) {
	$report = $self->_RunBlast($seq1, $seq2, $blastParams);
	if (!defined $report) {
	    croak "Error running BLAST\n";
	}
	
	#
	# At time of coding, bioperl throws
	# an exception if there are no HSPs
	# in the bl2seq report, so trap it.
	#
	eval {
	    while (my $hsp = $report->next_feature) {
	    	$hsp_found = 1;
		if ($hsp->hit->strand == 1) {
		     push @plus_hsps, $hsp;
		     $plus_hsp_count++;
		} elsif ($hsp->hit->strand == -1) {
		     push @minus_hsps, $hsp;
		     $minus_hsp_count++;
		}
	    }
	};

	#
	# If the overall alignment orientation has already been established
	# check to see if the HSPs return match this orientation. If not,
	# set HSP found flag to be false.
	#
	if (defined $self->{orient}) {
	    if ($self->{orient} == 1) {
	    	if (!$plus_hsp_count) {
		    $hsp_found = 0;
		}
	    } elsif ($self->{orient} == -1) {
	    	if (!$minus_hsp_count) {
		    $hsp_found = 0;
		}
	    }
	}

	if (!$hsp_found) {
	    $blastParams = $self->_BlastParamsAdjust($blastParams);
	}
    }

    $report->close if defined $report;

    return 0 if !$hsp_found;

    my $plus_chain;
    my $plus_score = 0;
    if ($plus_hsp_count && (!defined $self->{orient} || $self->{orient} == 1)) {
	($plus_chain, $plus_score) = _longest_increasing_subseq(
							    \@plus_hsps, 1);
	if ($DEBUG) {
	    print "\nPlus chain: score = $plus_score\n";
	    foreach my $hsp (@$plus_chain) {
		printf "%6d %6d %6d %6d %6d\n",
			    $hsp->query->start,
			    $hsp->query->end,
			    $hsp->hit->start,
			    $hsp->hit->end,
			    $hsp->score;
	    }
	}
    }

    my $minus_chain;
    my $minus_score = 0;
    if ($minus_hsp_count && (!defined $self->{orient} || $self->{orient} == -1))
    {
	($minus_chain, $minus_score) = _longest_increasing_subseq(
							    \@minus_hsps, -1);
	if ($DEBUG) {
	    print "\nMinus chain: score = $minus_score\n";
	    foreach my $hsp (@$minus_chain) {
		printf "%6d %6d %6d %6d %6d\n",
			    $hsp->query->start,
			    $hsp->query->end,
			    $hsp->hit->start,
			    $hsp->hit->end,
			    $hsp->score;
	    }
	}
    }

    return 0 if !$plus_score && !$minus_score;

    my $seg_list;
    if ($plus_score >= $minus_score) {
    	if (!defined $self->{orient}) {
	    $self->{orient} = 1;
	}
	
	# chains are already sorted so add HSP segments in order (do not
	# need to call HSPListInsert
	my $prev_seg;
	my $first = 1;
	foreach my $hsp (@$plus_chain) {
	    my $seg = AlignSegmentCreate('HSP', $hsp);
	    if ($first) {
	    	$seg_list = $seg;
	    	$first = 0;
	    } else {
		$prev_seg->next($seg);
	    	$seg->prev($prev_seg);
	    }
	    $prev_seg = $seg;
	}
    } else {
    	if (!defined $self->{orient}) {
	    $self->{orient} = -1;
	}
	
	# chains are already sorted so add HSP segments in order (do not
	# need to call HSPListInsert
	my $prev_seg;
	my $first = 1;
	foreach my $hsp (@$minus_chain) {
	    my $seg = AlignSegmentCreate('HSP', $hsp);
	    if ($first) {
	    	$seg_list = $seg;
	    	$first = 0;
	    } else {
		$prev_seg->next($seg);
	    	$seg->prev($prev_seg);
	    }
	    $prev_seg = $seg;
	}
    }

    $alignment->alignList($seg_list);

    return 1;
}

sub _longest_increasing_subseq
{
    my ($hsps, $strand) = @_;

    my $lisc;
    my %scpl_lab;
    
    # sort the input array if not already sorted
    @$hsps = sort {$a->query->start <=> $b->query->start} @$hsps;

    #print "\nSorted Hits:\n";
    #foreach my $hsp (@$hsps) {
    #	printf "%6d %6d %6d %6d %6d\n",
    #		    $hsp->query->start,
    #		    $hsp->query->end,
    #		    $hsp->hit->start,
    #		    $hsp->hit->end,
    #		    $hsp->score;
    #}
    
    # build the DAG
    my %SON2FATHER;	# stores the links (hash key) to father nodes
    my %FATHER2SON;	# stores the links (hash key) to children nodes
    my %SCORE;		# stores the scores of all father chains to the roots
    
    my $size = scalar @$hsps;
    for (my $i = 0; $i < $size - 1; $i++) {
	my $score = $hsps->[$i]->score;
	if ($strand == 1) {
	    my $qend1 = $hsps->[$i]->query->end;
	    my $hend1 = $hsps->[$i]->hit->end;
	    for (my $j = $i+1; $j < $size; $j++) {
		if (($hsps->[$j]->query->start > $qend1)
			&& ($hsps->[$j]->hit->start > $hend1))
		{
		    # non overlapping and non crossing hits 
		    push @{$FATHER2SON{$i}}, $j;
		    push @{$SON2FATHER{$j}}, $i;
		    push @{$SCORE{$j}}, $score;
		}
	    }
	} elsif ($strand == -1) {
	    my $qend1 = $hsps->[$i]->query->end;
	    my $hstart1 = $hsps->[$i]->hit->start;
	    for (my $j = $i+1; $j < $size; $j++) {
		if (($hsps->[$j]->query->start > $qend1)
			&& ($hsps->[$j]->hit->end < $hstart1))
		{
		    # non overlapping and non crossing hits 
		    push @{$FATHER2SON{$i}}, $j;
		    push @{$SON2FATHER{$j}}, $i;
		    push @{$SCORE{$j}}, $score;
		}
	    }
	}
    }
    
    # calculate the scores of the sub-segments
    my $glob_max = 0;
    #my $glob_end_id= "N";	# glob_end_id stores the key ($i) for the
    				# optimal chaining global score
    my $glob_end_id = undef;

    for (my $i = 0; $i < $size; $i++) {
    	my $score = $hsps->[$i]->score;
	if (exists($SON2FATHER{$i})) {
	    my $tmax;
	    my $maxi = 0;	# detect the local maximum (for 1 node)
	    for (my $j = 0; $j <= $#{$SON2FATHER{$i}}; $j++) {
	        # for all the parents
		if (($SCORE{$i}[$j] + $score) > $maxi) {
		    # on trouve le max
		    $maxi = $SCORE{$i}[$j] + $score;
		    if ($maxi > $glob_max) {
		    	$glob_max = $maxi;
			$glob_end_id = $i;
		    }
		}
	    }

	    if (exists($FATHER2SON{$i})) {
	        # if the node has children
		for (my $j = 0; $j <= $#{$FATHER2SON{$i}}; $j++) {
		    # for all the children
		    for (my $s = 0; $s <= $#{$SON2FATHER{$FATHER2SON{$i}[$j]}};					$s++)
		    {
		    	# for all the parents of the child (!)
			if ($SON2FATHER{$FATHER2SON{$i}[$j]}[$s] == $i) {
			    # le parent est le noeud en cours
			    $SCORE{$FATHER2SON{$i}[$j]}[$s] = $maxi;
			}
		    }
		}
	    }
	} else {
	    if ($score > $glob_max) {
		$glob_max = $score;
		$glob_end_id = $i;
	    }
	}
    }

    # traceback
    my $glob_start_id = $glob_end_id;
    my @path;
    push @path, $glob_end_id if defined $glob_end_id;

    while (defined $glob_start_id && exists($SON2FATHER{$glob_start_id})) {
	my $lmax = 0;
	my $vmax;
	for (my $i = 0; $i <= $#{$SON2FATHER{$glob_start_id}}; $i++) {
	    # for all the father nodes of $glob_end_id,
	    # choose those with best score
	    if ($SCORE{$glob_start_id}[$i] > $lmax) {
		$lmax = $SCORE{$glob_start_id}[$i];
		$vmax = $SON2FATHER{$glob_start_id}[$i];
		#print $vmax." ".$lmax."\n";
	    }
	}
	$glob_start_id = $vmax;
	push @path, $vmax;
    }
    
    my @chained_hsps;
    foreach my $idx (reverse @path) {
    	push @chained_hsps, $hsps->[$idx];
    }
    
    return(\@chained_hsps, $glob_max);
}

#
# Add HSP bounded segments to the alignment
#
sub _AlignmentAddHBSs
{
    my ($self, $alignment) = @_;

    return 0 if !defined $alignment;

    my $seq1 = $alignment->seq1;
    return 0 if !defined $seq1;

    my $seq2 = $alignment->seq2;
    return 0 if !defined $seq2;

    my $curElem = $alignment->alignList;
    if (!defined $curElem) {
	#
	# If blast produced no usable HSPs then it's possible that
	# the alignment list is empty. If so create a new alignment
	# without using HSPs. i.e. Use the external global alignment
	# program
	#
	my $seq1Start	= 1;
	my $seq1End	= $seq1->length;
	my $seq1Len	= $seq1->length;
	my $seq2Start	= 1;
	my $seq2End	= $seq2->length;
	my $seq2Len	= $seq2->length;

	my $segment = undef;
	if (($seq1Len * $seq2Len)
		<= ($self->{ga_max_len} * $self->{ga_max_len}))
	{
	    #
	    # Sequence lengths are within the length limits
	    # so create an alignable HBS from the complete
	    # sequences.
	    #
	    my $hbs = undef;
	    $hbs = HBSCreate(1, $seq1Start, $seq1End, $seq2Start, $seq2End,
	    			0, 0, 0, 0);
	    if (defined $hbs) {
		$segment = AlignSegmentCreate('HBS', $hbs);
	    }
	} else {
		#
		# Sequence lengths are outside the length limits
		# so create a non-alignable HBS from the complete
		# sequences.
		#
		my $hbs = undef;
		$hbs = HBSCreate(0, $seq1Start, $seq1End,
		    $seq2Start, $seq2End, 0, 0, 0, 0);
		if (defined $hbs) {
		    $segment = AlignSegmentCreate('HBS', $hbs);
		}
	}

	return 0 if !defined $segment;

	$alignment->alignList($segment);

	return 1;
    }

    return 0 if $curElem->type ne 'HSP';

    my $curHSP = $curElem->data;
    return 0 if !defined $curHSP;

    my $seq1Start	= 1;
    my $seq1End	= $seq1->length;
    my $seq2Start	= 1;
    my $seq2End	= $seq2->length;

    my $hspQStart	= $curHSP->query->start;
    my $hspQEnd	= $curHSP->query->end;
    my $hspHStart	= $curHSP->hit->start;
    my $hspHEnd	= $curHSP->hit->end;

    #
    # Get first HSP Bounded Segment (HBS).
    #
    my (
	$valid1,
	$s1Start,
	$s1OStart,
	$s1End,
	$s1OEnd,
	$s1LeftOL,
	$s1RightOL) = CalcStartSegBoundaries($seq1Start, $seq1End, $hspQStart,
								$hspQEnd, 1);
    my (
	$valid2,
	$s2Start,
	$s2OStart,
	$s2End,
	$s2OEnd,
	$s2LeftOL,
	$s2RightOL) = CalcStartSegBoundaries($seq2Start, $seq2End, $hspHStart,
						    $hspHEnd, $self->{orient});

    my $segment = $self->_SegmentBuild(
		    $alignment, $seq1, $seq2, $valid1, $valid2,
		    $s1Start, $s1End, $s1OStart, $s1OEnd, $s1LeftOL,
		    $s1RightOL, $s2Start, $s2End, $s2OStart, $s2OEnd,
		    $s2LeftOL, $s2RightOL, $self->{orient});

    if (defined $segment) {
	#
	# The new segment becomes the first of the alignment
	#
	$alignment->alignList($segment);
	$segment->next($curElem);
	$curElem->prev($segment);
    }

    #
    # Process subsequent HBSs.
    #
    my $nextElem = $curElem->next;
    while (defined $nextElem) {
	my $curHSP	= $curElem->data;
	my $curQStart	= $curHSP->query->start;
	my $curQEnd	= $curHSP->query->end;
	my $curHStart	= $curHSP->hit->start;
	my $curHEnd	= $curHSP->hit->end;

	my $nextHSP	= $nextElem->data;
	my $nextQStart	= $nextHSP->query->start;
	my $nextQEnd	= $nextHSP->query->end;
	my $nextHStart	= $nextHSP->hit->start;
	my $nextHEnd	= $nextHSP->hit->end;

	my (
	    $valid1,
	    $s1Start,
	    $s1OStart,
	    $s1End,
	    $s1OEnd,
	    $s1LeftOL,
	    $s1RightOL) = CalcMidSegBoundaries($seq1Start, $seq1End,
			    $curQStart, $curQEnd, $nextQStart, $nextQEnd, 1);
	my (
	    $valid2,
	    $s2Start,
	    $s2OStart,
	    $s2End,
	    $s2OEnd,
	    $s2LeftOL,
	    $s2RightOL) = CalcMidSegBoundaries($seq2Start, $seq2End,
			    $curHStart, $curHEnd, $nextHStart, $nextHEnd,
			    $self->{orient});

	my $segment = $self->_SegmentBuild(
			$alignment, $seq1, $seq2, $valid1,
			$valid2, $s1Start, $s1End, $s1OStart, $s1OEnd,
			$s1LeftOL, $s1RightOL, $s2Start, $s2End,
			$s2OStart, $s2OEnd, $s2LeftOL, $s2RightOL,
			$self->{orient});

	if (defined $segment) {
	    #
	    # Insert the new segment between the current and
	    # next HSPs
	    #
	    $curElem->next($segment);
	    $segment->prev($curElem);
	    $segment->next($nextElem);
	    $nextElem->prev($segment);
	}

	$curElem = $nextElem;
	$nextElem = $curElem->next;
    }

    
    #
    # Get last HSP Bounded Segment (HBS).
    #
    $curHSP		= $curElem->data;
    $hspQStart	= $curHSP->query->start;
    $hspQEnd	= $curHSP->query->end;
    $hspHStart	= $curHSP->hit->start;
    $hspHEnd	= $curHSP->hit->end;
    (
	$valid1,
	$s1Start,
	$s1OStart,
	$s1End,
	$s1OEnd,
	$s1LeftOL,
	$s1RightOL) = CalcEndSegBoundaries($seq1Start, $seq1End, $hspQStart,
						$hspQEnd, 1);
    (
	$valid2,
	$s2Start,
	$s2OStart,
	$s2End,
	$s2OEnd,
	$s2LeftOL,
	$s2RightOL) = CalcEndSegBoundaries($seq2Start, $seq2End, $hspHStart,
						$hspHEnd, $self->{orient});

    $segment = $self->_SegmentBuild(
		    $alignment, $seq1, $seq2, $valid1, $valid2,
		    $s1Start, $s1End, $s1OStart, $s1OEnd, $s1LeftOL,
		    $s1RightOL, $s2Start, $s2End, $s2OStart, $s2OEnd,
		    $s2LeftOL, $s2RightOL, $self->{orient});

    if (defined $segment) {
	#
	# The new segment becomes the last of the alignment
	#
	$curElem->next($segment);
	$segment->prev($curElem);
    }

    return 1;
}

sub _SegmentBuild
{
    my ($self, $alignment, $seq1, $seq2, $valid1, $valid2, $s1Start, $s1End,
	    $s1OStart, $s1OEnd, $s1LeftOL, $s1RightOL, $s2Start, $s2End,
	    $s2OStart, $s2OEnd, $s2LeftOL, $s2RightOL, $orient) = @_;

    my $segment = undef;

    if ($valid1) {
	my $s1OLen = $s1OEnd - $s1OStart;
	if ($valid2) {
	    my $s2OLen = $s2OEnd - $s2OStart;
	    if (($s1OLen * $s2OLen) > ($self->{ga_max_len}
			    * $self->{ga_max_len}))
	    {
		#
		# Length of at least of the segments
		# exceeds maximum for external global
		# alignment program. So attempt to build
		# a new sub-alignment for this segment
		# with adjusted Blast parameters.
		#
		my $newBlastParams = $self->_BlastParamsAdjust(
						    $alignment->blastParams);

		if (defined $newBlastParams) {
		    my $hbs1 = $seq1->trunc($s1Start, $s1End);
		    my $hbs2 = $seq2->trunc($s2Start, $s2End);
		    #
		    # Create a sub-alignment to add to the
		    # alignment.
		    #
		    if ($DEBUG) {
		    	print "Building sub-alignment $s1Start-$s1End\n";
		    }
		    my $subalign = $self->_AlignmentBuild($hbs1, $hbs2,
							    $newBlastParams);
		    if (defined $subalign) {
			    $segment = AlignSegmentCreate('SAL', $subalign);
		    }
		}

		if (!defined $segment) {
		    #
		    # We could not create a sub-aligment
		    # for this segment so just create
		    # a non-alignable segment to add
		    # to the alignment.
		    #
		    my $hbs = HBSCreate(0, $s1Start, $s1End, $s2Start, $s2End,
								0, 0, 0, 0);
		    $segment = AlignSegmentCreate('HBS', $hbs);
		}
	    } else {
		#
		# Lengths are within allowable limits for
		# external global alignment program so create
		# an HSP Bounded Segment to add to the
		# alignment.
		#
		my $hbs = HBSCreate(1, $s1OStart, $s1OEnd,
				$s2OStart, $s2OEnd,
				$s1LeftOL, $s1RightOL,
				$s2LeftOL, $s2RightOL);
		$segment = AlignSegmentCreate('HBS', $hbs);
	    }
	} else {
	    #
	    # Sequence 1 segment is valid but sequence 2 segment
	    # is not so create a non-alignable segment.
	    #
	    my $hbs = HBSCreate(0, $s1Start, $s1End, undef, undef,
			    0, 0, 0, 0);
	    $segment = AlignSegmentCreate('HBS', $hbs);
	}
    } else {
	if ($valid2) {
	    #
	    # Sequence 2 segment is valid but sequence 1 segment
	    # is not so create a non-alignable segment.
	    #
	    my $hbs = HBSCreate(0, undef, undef, $s2Start, $s2End, 0, 0, 0, 0);
	    $segment = AlignSegmentCreate('HBS', $hbs);
	}
    }

    return $segment;
}

#
# Create the alignment file for the given segment
#
sub _CreateSegmentFile
{
    my ($self, $alignment, $segment) = @_;

    return 0 if !defined $alignment;
    return 0 if !defined $segment;

    my $rval = 0;

    my $type = $segment->type;
    if ($type eq 'HSP') {
	my $hsp = $segment->data;
	if ($hsp) {
	    my $fileName = $self->_TempFileName('hsp');
	    $rval = CreateHSPAlignFile($fileName, $hsp);
	    $segment->file($fileName);
	}
    } elsif ($type eq 'HBS') {
	my $seq1 = $alignment->seq1;
	my $seq2 = $alignment->seq2;
	my $hbs = $segment->data;
	if ($hbs) {
	    my $fileName = $self->_TempFileName('hbs');
	    $rval = $self->_CreateHBSAlignFile($fileName, $hbs, $seq1, $seq2);
	    $segment->file($fileName) if $rval;
	}
    }

    return $rval;
}

#
# Create an alignment file with the given file name for the given HBS
#
sub _CreateHBSAlignFile
{
    my ($self, $fileName, $hbs, $seq1, $seq2) = @_;

    return 0 if !defined $fileName;
    return 0 if !defined $hbs;
    return 0 if !defined $seq1;
    return 0 if !defined $seq2;

    my $alignable	= $hbs->alignable;
    my $s1Start	= $hbs->seg1Start;
    my $s1End	= $hbs->seg1End;
    my $s1LeftOL	= $hbs->seg1LeftOL;
    my $s1RightOL	= $hbs->seg1RightOL;
    my $s2Start	= $hbs->seg2Start;
    my $s2End	= $hbs->seg2End;
    my $s2LeftOL	= $hbs->seg2LeftOL;
    my $s2RightOL	= $hbs->seg2RightOL;

    my $rval = 0;
    my $errcode = 0;
    my $hbsSeq1 = undef;
    my $hbsSeq2 = undef;
    if ($alignable) {
	$hbsSeq1 = $seq1->trunc($s1Start, $s1End);
	if ($self->{orient} == 1) {
	    $hbsSeq2 = $seq2->trunc($s2Start, $s2End);
	    ($rval, $errcode) = $self->_GlobalAlign(
				    $fileName, $hbsSeq1, $hbsSeq2,
				    $s1LeftOL, $s1RightOL,
				    $s2LeftOL, $s2RightOL);
	} elsif ($self->{orient} == -1) {
	    $hbsSeq2 = $seq2->trunc($s2Start, $s2End)->revcom;
	    #
	    # Note we have to reverse the left and right overlaps
	    # for the second sequence since it has been reverse
	    # complemented.
	    #
	    ($rval, $errcode) = $self->_GlobalAlign(
					    $fileName, $hbsSeq1, $hbsSeq2,
					    $s1LeftOL, $s1RightOL,
					    $s2RightOL, $s2LeftOL);
	}

	#
	# If we are unable to globally align here, make the whole
	# alignment fail unless the global alignment failed because
	# the sequences were too long, in which case create a
	# non-alignable segment and continue with the overall
	# alignment. All other errors at this point are considered
	# fatal.
	# Note that the check of whether the sequences are too long
	# for the global alignment program is performed earlier within
	# this module. However, this relies on the maximum sequence
	# length definition in this module being in sync with that of
	# the nwalign program which is a poor assumption, so we
	# consider these checks independent.
	#
	if (!$rval) {
	    if ($errcode == NW_SEQ_LEN_ERR) {
		print "Global alignment sequences too long"
			. " - creating non-alignable segment\n" if $DEBUG;
		$alignable = 0;
		$hbs->alignable(0);
		undef $hbsSeq1;
		undef $hbsSeq2;
	    } else {
		carp "Global alignment failed\n";
	    }
	}
    }
    
    if (!$alignable) {
	if (defined $s1Start && defined $s1End) {
	    $hbsSeq1 = $seq1->trunc($s1Start + $s1LeftOL,
					$s1End - $s1RightOL)->seq;
	}
	if (defined $s2Start && defined $s2End) {
	    if ($self->{orient} == 1) {
		$hbsSeq2 = $seq2->trunc($s2Start + $s2LeftOL,
					$s2End - $s2RightOL)->seq;
	    } elsif ($self->{orient} == -1) {
		$hbsSeq2 = $seq2->trunc($s2Start + $s2LeftOL,
					$s2End - $s2RightOL)->revcom->seq;
	    }
	}
	$rval = WriteNonAlignedSeqs($fileName, \$hbsSeq1, \$hbsSeq2);
    }

    return $rval;
}

#
# Call global alignment on input sequences.
#
sub _GlobalAlign
{
    my ($self, $alignFile, $seq1, $seq2, $s1LeftOL, $s1RightOL, $s2LeftOL,
	    $s2RightOL) = @_;

    return (0, 0) if !defined $alignFile;
    return (0, 0) if !defined $seq1;
    return (0, 0) if !defined $seq2;
    $s1LeftOL = 0 if !defined $s1LeftOL;
    $s1RightOL = 0 if !defined $s1RightOL;
    $s2LeftOL = 0 if !defined $s2LeftOL;
    $s2RightOL = 0 if !defined $s2RightOL;

    #
    # Create a temporary FastA file from the two input sequences
    # Use same file name as alignment file with extension changed
    # to "fa".
    #
    my $seqFile = ChangeFileExt($alignFile, 'fa');
    my $out = Bio::SeqIO->new(
		    '-file'		=> ">$seqFile",
		    '-format'	=> 'fasta');
    if (!$out) {
	carp "_GlobalAlign: error creating temp. FastA sequences file\n";
	return (0, 0);
    }

    if (!$out->write_seq($seq1)) {
	carp "_GlobalAlign: error writing sequence 1 to temp. FastA file\n";
	return (0, 0);
    }
    if (!$out->write_seq($seq2)) {
	carp "_GlobalAlign: error writing sequence 2 to temp. FastA file\n";
	return (0, 0);
    }

    $out->close;

    my $NWOutFile = ChangeFileExt($alignFile, 'nw');

    #
    # Call external Needleman-Wunsch alignment program
    #
    my ($nw_ok, $nw_rval) = $self->_NWAlign($seqFile, $self->{ga_match},
					    $self->{ga_mismatch},
					    $self->{ga_gap_penalty},
					    $self->{ga_gap_ext_penalty},
					    $NWOutFile);
    if (!$nw_ok) {
	carp "_GlobalAlign: error running N-W alignment\n";
	return (0, $nw_rval);
    }

    if (!ConvertNWFile($NWOutFile, $s1LeftOL, $s1RightOL, $s2LeftOL,
	    $s2RightOL, $alignFile))
    {
	carp "_GlobalAlign: error converting N-W file\n";
	return (0, 0);
    }

    return (1, 0);
}

#
# Call Needleman-Wunsch alignment program.
#
sub _NWAlign
{
    my ($self, $inFile, $matchScore, $mismatchScore, $gapPenalty,
	    $gapExtPenalty, $outFile) = @_;

    my $cmd = $self->{_ga_prog} . " $inFile $matchScore $mismatchScore"
		. " $gapPenalty $gapExtPenalty $outFile";
    
    my $out = `$cmd 2>&1`;
    print "_NWAlign output:\n$out\n" if $DEBUG;
 
    my $rval = ($? >> 8);
    if ($rval != 0) {
	carp "_NWAlign: error running $cmd\n$out\n";
	return (0, $rval);
    }

    return (1, 0);
}

#
# Run stand alone blast (Blastn) with the parameter given by alignment
#
sub _RunBlast
{
    my ($self, $seq1, $seq2, $blastParams) = @_;

    return undef if !defined $seq1;
    return undef if !defined $seq2;
    return undef if !defined $blastParams;

    my $paramArray = BlastParams2Array($blastParams);

    #
    # We no longer have to create an actual blast file, but since bl2seq
    # creates it's own temp file anyway, we may as well specify it and
    # keep it in our temp dir with all our other temp files.
    #
    my $blastFile = $self->_TempFileName('bls');
    push @$paramArray, ('outfile' => $blastFile);

    #
    # Necessary for later versions (>= 1.4 (1.3?)) of Bioperl to correctly
    # set the returned report type to BPbl2seq.
    #
    push @$paramArray, ('_READMETHOD' => 'BPlite');

    #
    # Setup and run bl2seq
    #
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@$paramArray);
    if (!defined $factory) {
	carp "Error setting up stand alone BLAST\n";
	return undef;
    }

    my $report = $factory->bl2seq($seq1, $seq2);
    if (!defined $report) {
	carp "Error running bl2seq. No BLAST report returned.\n";
	return undef;
    }

    return $report;
}

sub _CheckBlast
{
    my $factory = Bio::Tools::Run::StandAloneBlast->new();
    if (!$factory) {
	carp "Error checking for existence of BLAST programs\n";
	return 0;
    }

    my $exe = $factory->executable('bl2seq');
    if (!$exe) {
	carp "bl2seq program not found or not executable\n";
	return 0;
    }

    #my $exe2 = $factory->executable(BL_PROG);
    #if (!$exe2) {
    #	carp "BLAST program " . BL_PROG . " not found or not executable\n";
    #	return 0;
    #}

    return 1;
}

sub _GAProg
{
    my ($self, $gaProg) = @_;

    if (defined $gaProg) {
	if ($self->_Executable($gaProg)) {
	    $self->{_ga_prog} = $gaProg;
	} else {
	    carp "Global alignment program $gaProg not found"
		    . " or not executable\n";
	}
    }

    return $self->{_ga_prog} || undef;
}

#
# Check to see if the specified program exists and is executable anywhere
# within the user's path.
#
sub _Executable
{
    my ($self, $prog) = @_;

    my $isExecutable = 0;

    return $isExecutable if !defined $prog;

    if (-x $prog) {
	$isExecutable = 1;
    } else {
	my $path = $ENV{PATH};
	my @pathdirs = split(':', $path);
	foreach my $dir (@pathdirs) {
	    my $exe = catfile($dir, $prog);
	    if (-x $exe) {
		$isExecutable = 1;
		last;
	    }
	}
    }

    return $isExecutable;
}

#
# Make a directory in which to create all the temporary alignment and
# other files.
#
sub _MakeTempDir
{
    my $self = shift;

    my $template = sprintf("OrcaAlign%04d%02d%02dXXXXXX",
		    $self->{year}, $self->{month}, $self->{day});

    #
    # CLEANUP => 1 to autonatically removes dir and contents at
    # termination of program execution
    #
    my $tmpDir;
    if ($self->{debug}) {
	$tmpDir = tempdir($template, TMPDIR => 1);
    } else {
	$tmpDir = tempdir($template, TMPDIR => 1, CLEANUP => 1);
    }

    $self->{_tempDir} = $tmpDir;

    return $self->{_tempDir} || undef;
}

#
# Create and return a new unique temporary file name.
#
sub _TempFileName
{
    my ($self, $ext) = @_;

    my $name;

    $name = File::Temp::tempnam($self->{_tempDir}, "");
    $name .= ".$ext" if defined $ext && $ext;

    return $name;
}

sub _BlastParamsAdjust
{
    my ($self, $params) = @_;

    my $wordSize;
    my $mismatchScore;
    my $expect;

    $wordSize = $params->wordSize;
    $mismatchScore = $params->mismatchScore;
    $expect = $params->expect;

    my $newWordSize = $wordSize - $self->BLWordSizeDec;
    $newWordSize = BL_WORDSIZE_MIN if $newWordSize < BL_WORDSIZE_MIN;

    my $newMismatchScore = $mismatchScore + $self->BLMismatchInc;
    $newMismatchScore = BL_MISMATCH_MAX if $newMismatchScore > BL_MISMATCH_MAX;

    my $newExpect = $expect * $self->BLExpectMult;
    $newExpect = BL_EXPECT_MAX if $newExpect > BL_EXPECT_MAX;


    my $newParams = undef;
    if ($newWordSize < $wordSize || $newMismatchScore > $mismatchScore
    	|| $newExpect > $expect)
    {
	$newParams = BlastParamsCreate(
					$params->gapped,
					#$params->expect,
					$newExpect,
					$params->gapPenalty,
					$params->gapExtPenalty,
					$params->xDropoff,
					$params->matchScore,
					$newMismatchScore,
					$newWordSize,
					$params->filter);
    }

    return $newParams;
}

sub _Cleanup
{
    my $self = shift;

    unless ($self->{debug}) {
	#
	# Remove all files in the temporary directory.
	#
	my $tempDir = $self->{_tempDir};
	if (defined $tempDir) {
	    unlink(<$tempDir/*>);
	}
    }
}

sub DESTROY
{
    my $self = shift;

    unless ($self->{debug}) {
	#
	# Remove temporary directory/files? They are set to
	# autoclean anyway.
	#
	if (defined $self->{_tempDir}) {
		rmtree($self->{_tempDir});
	}
    }
}

################################################################################
# Helper Routines
#
# The following routines are helper routines, not methods.
################################################################################
sub BlastParamsCreate
{
    my ($gapped, $expect, $gapPenalty, $gapExtPenalty,
	$xDropoff, $matchScore, $mismatchScore,
	$wordSize, $filter) = @_;

    my $params = BlastParams->new();
    return undef if !defined $params;

    $params->program(BL_PROG);
    $params->gapped($gapped) if defined $gapped;
    $params->expect($expect) if defined $expect;
    $params->gapPenalty($gapPenalty) if defined $gapPenalty;
    $params->gapExtPenalty($gapExtPenalty) if defined $gapExtPenalty;
    $params->xDropoff($xDropoff) if defined $xDropoff;
    $params->matchScore($matchScore) if defined $matchScore;
    $params->mismatchScore($mismatchScore) if defined $mismatchScore;
    $params->wordSize($wordSize) if defined $wordSize;
    $params->filter($filter) if defined $filter;

    return $params;
}

sub BlastParams2Array
{
    my ($params) = @_;

    my @paramArray;

    push @paramArray, ('p' => $params->program)
	    if defined $params->program;
    push @paramArray, ('g' => $params->gapped)
	    if defined $params->gapped;
    push @paramArray, ('e' => $params->expect)
	    if defined $params->expect;
    push @paramArray, ('G' => $params->gapPenalty)
	    if defined $params->gapPenalty;
    push @paramArray, ('E' => $params->gapExtPenalty)
	    if defined $params->gapExtPenalty;
    push @paramArray, ('X' => $params->xDropoff)
	    if defined $params->xDropoff;
    push @paramArray, ('r' => $params->matchScore)
	    if defined $params->matchScore;
    push @paramArray, ('q' => $params->mismatchScore)
	    if defined $params->mismatchScore;
    push @paramArray, ('W' => $params->wordSize)
	    if defined $params->wordSize;
    push @paramArray, ('F' => $params->filter)
	    if defined $params->filter;

    return \@paramArray;
}

#
# Calculate the boundaries of a segment bounded by the start of a sequence
# and the first HSP.
#
sub CalcStartSegBoundaries
{
    my ($seqStart, $seqEnd, $hspStart, $hspEnd, $orient) = @_;

    my $sStart = undef;		# segment start
    my $sOStart = undef;	# segment start including overlap
    my $sEnd = undef;		# segment end
    my $sOEnd = undef;		# segment end including overlap
    my $sLeftOL = undef;	# amount segment overlaps previous HSP
    my $sRightOL = undef;	# amount segment overlaps next HSP
    my $valid = 0;
    if ($orient == 1) {
	if ($hspStart > $seqStart) {
	    $sStart = $seqStart;
	    $sOStart = $sStart;
	    $sEnd = $hspStart - 1;
	    $sOEnd = $sEnd + HBS_ANCHOR_LEN;
	    $sOEnd = $seqEnd if $sOEnd > $seqEnd;
	    $sLeftOL = 0;
	    $sRightOL = $sOEnd - $sEnd;
	    $valid = 1;
	}
    } elsif ($orient == -1) {
	if ($hspEnd < $seqEnd) {
	    $sStart = $hspEnd + 1;
	    $sOStart = $sStart - HBS_ANCHOR_LEN;
	    $sOStart = $seqStart if $sOStart < $seqStart;
	    $sEnd = $seqEnd;
	    $sOEnd = $sEnd;
	    $sLeftOL = $sStart - $sOStart;
	    $sRightOL = 0;
	    $valid = 1;
	}
    }

    return ($valid, $sStart, $sOStart, $sEnd, $sOEnd, $sLeftOL, $sRightOL);
}

#
# Calculate the boundaries of a segment bounded by the end of a sequence
# and the last HSP.
#
sub CalcEndSegBoundaries
{
    my ($seqStart, $seqEnd, $hspStart, $hspEnd, $orient) = @_;

    my $sStart = undef;		# segment start
    my $sOStart = undef;	# segment start including overlap
    my $sEnd = undef;		# segment end
    my $sOEnd = undef;		# segment end including overlap
    my $sLeftOL = undef;	# amount segment overlaps previous HSP
    my $sRightOL = undef;	# amount segment overlaps next HSP
    my $valid = 0;
    if ($orient == 1) {
	if ($hspEnd < $seqEnd) {
	    $sStart = $hspEnd + 1;
	    $sOStart = $sStart - HBS_ANCHOR_LEN;
	    $sOStart = $seqStart if $sOStart < $seqStart;
	    $sEnd = $seqEnd;
	    $sOEnd = $sEnd;
	    $sLeftOL = $sStart - $sOStart;
	    $sRightOL = 0;
	    $valid = 1;
	}
    } elsif ($orient == -1) {
	if ($hspStart > $seqStart) {
	    $sStart = $seqStart;
	    $sOStart = $sStart;
	    $sEnd = $hspStart - 1;
	    $sOEnd = $sEnd + HBS_ANCHOR_LEN;
	    $sOEnd = $seqEnd if $sOEnd > $seqEnd;
	    $sLeftOL = 0;
	    $sRightOL = $sOEnd - $sEnd;
	    $valid = 1;
	}
    }

    return ($valid, $sStart, $sOStart, $sEnd, $sOEnd, $sLeftOL, $sRightOL);
}

#
# Calculate the boundaries of a segment bounded by two HSPs in the middle
# of a sequence.
#
sub CalcMidSegBoundaries
{
    my ($seqStart, $seqEnd, $hsp1Start, $hsp1End, $hsp2Start,
	    $hsp2End, $orient) = @_;

    my $sStart = undef;		# segment start
    my $sOStart = undef;	# segment start including overlap
    my $sEnd = undef;		# segment end
    my $sOEnd = undef;		# segment end including overlap
    my $sLeftOL = undef;	# amount segment overlaps previous HSP
    my $sRightOL = undef;	# amount segment overlaps next HSP
    my $valid = 0;
    if ($orient == 1) {
	if ($hsp2Start > $hsp1End + 1) {
	    $sStart = $hsp1End + 1;
	    $sOStart = $sStart - HBS_ANCHOR_LEN;
	    $sOStart = $seqStart if $sOStart < $seqStart;
	    $sEnd = $hsp2Start - 1;
	    $sOEnd = $sEnd + HBS_ANCHOR_LEN;
	    $sOEnd = $seqEnd if $sOEnd > $seqEnd;
	    $sLeftOL = $sStart - $sOStart;
	    $sRightOL = $sOEnd - $sEnd;
	    $valid = 1;
	}
    } elsif ($orient == -1) {
	if ($hsp2End < $hsp1Start - 1) {
	    $sEnd = $hsp1Start - 1;
	    $sOEnd = $sEnd + HBS_ANCHOR_LEN;
	    $sOEnd = $seqEnd if $sOEnd > $seqEnd;
	    $sStart = $hsp2End + 1;
	    $sOStart = $sStart - HBS_ANCHOR_LEN;
	    $sOStart = $seqStart if $sOStart < $seqStart;
	    $sLeftOL = $sStart - $sOStart;
	    $sRightOL = $sOEnd - $sEnd;
	    $valid = 1;
	}
    }

    return ($valid, $sStart, $sOStart, $sEnd, $sOEnd, $sLeftOL, $sRightOL);
}

#
# Create a new HBS structure
#
sub HBSCreate
{
    my ($alignable, $seg1Start, $seg1End, $seg2Start, $seg2End, $seg1LeftOL,
	    $seg1RightOL, $seg2LeftOL, $seg2RightOL) = @_;

    my $hbs = HBS->new();
    return undef if !defined $hbs;

    $hbs->alignable($alignable);
    $hbs->seg1Start($seg1Start);
    $hbs->seg1End($seg1End);
    $hbs->seg2Start($seg2Start);
    $hbs->seg2End($seg2End);
    $hbs->seg1LeftOL($seg1LeftOL);
    $hbs->seg1RightOL($seg1RightOL);
    $hbs->seg2LeftOL($seg2LeftOL);
    $hbs->seg2RightOL($seg2RightOL);

    return $hbs;
}

#
# Create a new alignment segment of type 'type' from the given data
#
sub AlignSegmentCreate
{
    my ($type, $data) = @_;

    return undef if !defined $data || !defined $type;

    my $segment = AlignSegment->new();
    return undef if !defined $segment;

    $segment->type($type);
    $segment->data($data);

    return $segment;
}

#
# Return the orientation of an HSP
#
sub HSPOrient
{
    my ($hsp) = @_;

    return undef if !defined $hsp;

    if (!$hsp->hit->strand) {
	croak "error determining HSP orientation";
    }

    return $hsp->hit->strand;
}

#
# Compare two HSPs for colinearity.
# Assumes HSPs are already of the same orientation.
# If both the query and hit sequences of hsp1 fall completely after
# those of hsp2, then return 1
# If both the query and hit sequences hsp1 fall completely after those
# of hsp2, then return -1.
# In all other cases return 0.
#
# NOTE: In future we may have to explicitly deal with partial overlaps and so
# on.
#
sub HSPCompare
{
    my ($hsp1, $hsp2, $orient) = @_;

    return 0 if !defined $hsp1;
    return 0 if !defined $hsp2;

    my $qStart1	= $hsp1->query->start;
    my $qEnd1	= $hsp1->query->end;
    my $hStart1	= $hsp1->hit->start;
    my $hEnd1	= $hsp1->hit->end;

    my $qStart2	= $hsp2->query->start;
    my $qEnd2	= $hsp2->query->end;
    my $hStart2	= $hsp2->hit->start;
    my $hEnd2	= $hsp2->hit->end;

    if ($orient == 1) {
	#
	# Plus/Plus orientation
	#
	if ($qEnd1 < $qStart2 && $hEnd1 < $hStart2) {
	    #
	    # HSP 1 query lies completely before HSP 2 query
	    # and HSP 1 hit lies completely before HSP 2 hit
	    #
	    return -1;
	} elsif ($qStart1 > $qEnd2 && $hStart1 > $hEnd2) {
	    # 
	    # HSP 1 query lies completely after HSP 2 query
	    # and HSP 1 hit lies completely after HSP 2 hit
	    #
	    return 1;
	}
    } elsif ($orient == -1) {
	#
	# Plus/Minus (reverse complement) orientation
	#
	if ($qEnd1 < $qStart2 && $hStart1 > $hEnd2) {
	    #
	    # HSP 1 query lies completely before HSP 2 query
	    # and HSP 1 hit lies completely after HSP 2 hit
	    #
	    return -1;
	} elsif ($qStart1 > $qEnd2 && $hEnd1 < $hStart2) {
	    # 
	    # HSP 1 query lies completely after HSP 2 query
	    # and HSP 1 hit lies completely before HSP 2 hit
	    #
	    return 1;
	}
    }

    return 0;
}

#
# Insert HSP into list
#
sub HSPListInsert
{
    my ($hspList, $hsp, $orient) = @_;

    return undef if !defined $hspList;
    return $hspList if !defined $hsp;

    my $curElem = $hspList;
    my $prevElem = undef;
    my $comp = 0;
    while (defined $curElem) {
	$comp = HSPCompare($hsp, $curElem->data, $orient);
	last if $comp <= 0;
	$prevElem = $curElem;
	$curElem = $curElem->next;
    }

    if (!defined $curElem) {
	#
	# Insert new HSP after last element of list
	#
	my $segment = AlignSegmentCreate('HSP', $hsp);
	$prevElem->next($segment);
	$segment->prev($prevElem);
    } else {
	if ($comp < 0) {
	    #
	    # Insert new hsp before current element
	    #
	    my $segment = AlignSegmentCreate('HSP', $hsp);
	    $segment->next($curElem);
	    $curElem->prev($segment);
	    if (defined $prevElem) {
		$prevElem->next($segment);
		$segment->prev($prevElem);
	    } else {
		#
		# First element of list
		#
		$hspList = $segment;
	    }
	}
    }

    return $hspList;
}

#
# Used only for debugging purposes
#
sub AlignmentPrint
{
    my ($align) = @_;

    return if !defined $align;

    AlignListPrint($align->alignList);
}

#
# Used only for debugging purposes
#
sub AlignListPrint
{
    my ($alist) = @_;

    my $seg = $alist;

    while ($seg) {
	AlignSegmentPrint($seg);
	$seg = $seg->next;
    }
}

#
# Used only for debugging purposes
#
sub AlignSegmentPrint
{
    my ($seg) = @_;

    return if !defined $seg;

    my $type = $seg->type;
    my $data = $seg->data;
    my $file = $seg->file;

    if (defined $file) {
	print("${file}:\n");
    }

    if ($type eq 'HSP') {
	HSPPrint($data);
    } elsif ($type eq 'HBS') {
	HBSPrint($data);
    } elsif ($type eq 'SAL') {
	printf "BEGIN Sub-alignment\n";
	AlignmentPrint($data);
	printf "END Sub-alignment\n";
    }
}

#
# Used only for debugging purposes
#
sub HSPPrint
{
    my ($hsp) = @_;

    return if !defined $hsp;

    my $qStart	= $hsp->query->start;
    my $qEnd	= $hsp->query->end;
    my $hStart	= $hsp->hit->start;
    my $hEnd	= $hsp->hit->end;

    printf("HSP: %7d %7d %7d %7d\n", $qStart, $qEnd, $hStart, $hEnd);
}

#
# Used only for debugging purposes
#
sub HBSPrint
{
    my ($hbs) = @_;

    return if !defined $hbs;

    my $alignable	= $hbs->alignable;
    my $lol1	= $hbs->seg1LeftOL;
    my $rol1	= $hbs->seg1RightOL;
    my $lol2	= $hbs->seg2LeftOL;
    my $rol2	= $hbs->seg2RightOL;
    my $start1	= $hbs->seg1Start;
    my $end1	= $hbs->seg1End;
    my $start2	= $hbs->seg2Start;
    my $end2	= $hbs->seg2End;

    $start1	+= $lol1;
    $end1	-= $rol1;
    $start2	+= $lol2;
    $end2	-= $rol2;

    printf("HBS: %7d %7d %7d %7d %2d %2d %2d %2d %s\n",
	    $start1, $end1, $start2, $end2, $lol1, $rol1, $lol2, $rol2,
	    $alignable ? 'alignable' : 'non-alignable');
}

#
# Create an alignment file with the given file name for the given HSP
#
sub CreateHSPAlignFile
{
    my ($fileName, $hsp) = @_;

    return 0 if !defined $fileName;
    return 0 if !defined $hsp;

    my $hspSeq1 = $hsp->querySeq;
    return 0 if !defined $hspSeq1;

    my $hspSeq2 = $hsp->sbjctSeq;
    return 0 if !defined $hspSeq2;

    my $homol = $hsp->homologySeq;
    return 0 if !defined $homol;

    return WriteAlignedSeqs($fileName, \$hspSeq1, \$hspSeq2, \$homol);
}

#
# Create a new file from a N-W output file, trimming the overlapping
# parts of the sequences.
#
sub ConvertNWFile
{
    my ($NWFile, $s1LeftOL, $s1RightOL, $s2LeftOL, $s2RightOL, $outFile) = @_;

    return 0 if !defined $NWFile;
    return 0 if !defined $outFile;

    #
    # If the overlaps are all 0, simply copy the N-W file
    # Note: if we do not need to keep intermediate files (for
    # debugging or whatever) simply rename file for max. efficiencey
    #
    if (!$s1LeftOL && !$s1RightOL && !$s2LeftOL && !$s2RightOL) {
	if ($DEBUG) {
	    return copy($NWFile, $outFile);
	} else {
	    return move($NWFile, $outFile);
	}
    }

    if (!open(INFH, "$NWFile")) {
	carp "ConvertNWFile: error opening N-W file $NWFile\n";
	return 0;
    }

    my $seq1;
    my $seq2;
    my $homol;
    my $lineCount = 0;
    my $line;
    while (defined ($line = <INFH>)) {
	chomp $line;
	my $seqNum = $lineCount % 4;
	if ($seqNum == 0) {
		$seq1 .= $line;
	} elsif ($seqNum == 1) {
		$seq2 .= $line;
	} elsif ($seqNum == 2) {
		$homol .= $line;
	}
	$lineCount++;
    }
    close(INFH);

    ReconcileHBSBoundaries(\$seq1, \$seq2, \$homol, $s1LeftOL, $s1RightOL,
	    $s2LeftOL, $s2RightOL);

    if (!open(OUTFH, ">$outFile")) {
	carp "ConvertNWFile: error creating output file $outFile\n";
	return 0;
    }

    my $curPos = 0;
    my $size = length $seq1;
    while ($curPos < $size) {
	my $subSeq1 = substr($seq1, $curPos, SEQ_OUT_WIDTH);
	my $subSeq2 = substr($seq2, $curPos, SEQ_OUT_WIDTH);
	my $subSeq3 = substr($homol, $curPos, SEQ_OUT_WIDTH);
	my $length1 = length $subSeq1;
	my $length2 = length $subSeq2;
	if ($length1 != $length2) {
	    carp "mismatched alignment\n"
		    . "file = $NWFile\n"
		    . "length1 = $length1\n"
		    . "length2 = $length2\n"
		    . "subSeq1:\n$subSeq1\n\n"
		    . "subSeq2:\n$subSeq2\n\n"
		    . "seq1:\n$seq1\n\n"
		    . "seq2:\n$seq2\n\n";
	    close(OUTFH);
	    return 0;
	}
	print OUTFH "$subSeq1\n";
	print OUTFH "$subSeq2\n";
	print OUTFH "$subSeq3\n";
	print OUTFH "\n";
	$curPos += SEQ_OUT_WIDTH;
    }
    close(OUTFH);

    return 1;
}

#
# Adjust the boundaries of the HBS by removing the regions of overlap
# with the bounding HSPs and adding gap characters as necessary.
#
sub ReconcileHBSBoundaries
{
    my ($seq1, $seq2, $homol, $s1LeftOL, $s1RightOL, $s2LeftOL, $s2RightOL)
	    = @_;

    return if !defined $seq1 || !defined $seq2 || !defined $homol;
    $s1LeftOL	= 0 if !defined $s1LeftOL;
    $s1RightOL	= 0 if !defined $s1RightOL;
    $s2LeftOL	= 0 if !defined $s2LeftOL;
    $s2RightOL	= 0 if !defined $s2RightOL;


    return if $s1LeftOL == 0 && $s1RightOL == 0 && $s2LeftOL == 0
	    && $s2RightOL == 0;

    #
    # XXX
    # Is this more efficient than manipulating the sequence strings
    # directly with subseq (which requires string copying), chop etc?
    #
    my @a1 = unpack('c*', $$seq1);
    my @a2 = unpack('c*', $$seq2);

    my ($lc1, $rc1) = RemoveOverlaps(\@a1, $s1LeftOL, $s1RightOL);
    my ($lc2, $rc2) = RemoveOverlaps(\@a2, $s2LeftOL, $s2RightOL);

    my $lcDiff = $lc1 - $lc2;
    my $rcDiff = $rc1 - $rc2;

    #
    # Add gaps back onto the left hand side according to the difference
    # in characters removed between the two sequences
    #
    my $dash = ord('-');
    if ($lcDiff < 0) {
	my $count = 0;
	while ($count > $lcDiff) {
	    unshift @a2, $dash;
	    $count--;
	}
    } elsif ($lcDiff > 0) {
	my $count = 0;
	while ($count < $lcDiff) {
	    unshift @a1, $dash;
	    $count++;
	}
    }

    #
    # Add gaps back onto the right hand side according to the difference
    # in characters removed between the two sequences
    #
    if ($rcDiff < 0) {
	my $count = 0;
	while ($count > $rcDiff) {
	    push @a2, $dash;
	    $count--;
	}
    } elsif ($rcDiff > 0) {
	my $count = 0;
	while ($count < $rcDiff) {
	    push @a1, $dash;
	    $count++;
	}
    }

    $$seq1 = pack('c*', @a1);
    $$seq2 = pack('c*', @a2);

    #
    # Amounts to chop off homology sequence are the greater of the
    # corresponding amounts chopped off the two sequences.
    #
    my $hlc = $lcDiff > 0 ? $lc1 : $lc2;
    my $hrc = $rcDiff > 0 ? $rc1 : $rc2;

    #
    # The difference in the number chopped off each sequence have to 
    # be added back to the homology sequence as non-homologous chars.
    #
    my $hLen = (length $$homol) - $hlc - $hrc;
    $$homol = ' ' x abs($lcDiff) . substr($$homol, $hlc, $hLen)
	    . ' ' x abs($rcDiff);
}

#
# Remove the number of nucleotides from the left and right sides of sequence
# as specified by lol and rol. The numbers specified do not include any gap
# characters which fall in the area being pruned. Add the number of gap
# characters to the counts of characters actually removed and return these
# counts.
#
sub RemoveOverlaps
{
    my ($seq, $lol, $rol) = @_;

    return (0, 0) if !defined $seq;
    $lol = 0 if !defined $lol;
    $rol = 0 if !defined $rol;
    return (0, 0) if $lol == 0 && $rol == 0;

    my $dash = ord('-');

    my $leftCut = 0;
    my $rightCut = 0;

    #
    # Remove nucleotide characters from the left of sequence
    #
    my $curChar;
    my $toCut = $lol;
    while ($leftCut < $toCut) {
	$curChar = shift @$seq;
	#
	# If a gap is encountered increment the number to cut
	#
	$toCut++ if $curChar eq $dash;
	$leftCut++;
    }
    #
    # Remove any remaining gaps on the left
    #
    do {
	$curChar = shift @$seq;
	$leftCut++;
    } until ($curChar ne $dash);
    #
    # We've taken off one too many characters to put the last one
    # back on
    #
    unshift @$seq, $curChar;
    $leftCut--;

    $toCut = $rol;
    while ($rightCut < $toCut) {
	$curChar = pop @$seq;
	#
	# If a gap is encountered increment the number to cut
	#
	$toCut++ if $curChar eq $dash;
	$rightCut++;
    }
    #
    # Remove any remaining gaps on the right
    #
    do {
	$curChar = pop @$seq;
	$rightCut++;
    } until ($curChar ne $dash);
    #
    # We've taken off one too many characters to put the last one
    # back on
    #
    push @$seq, $curChar;
    $rightCut--;

    return ($leftCut, $rightCut);
}

#
# Write two aligned sequences to a file.
# Takes a file name, two sequence and one homology string.
#
sub WriteAlignedSeqs
{
    my ($file, $seq1, $seq2, $homol) = @_;

    return 0 if !defined $file || !defined $seq1 || !defined $seq2
	    || !defined $homol;

    return 0 if !open(FH, ">$file");

    my $len1 = length $$seq1;
    my $len2 = length $$seq2;
    if ($len1 != $len2) {
	close(FH);
	return 0;
    }

    my $curPos = 0;
    my $inc = SEQ_OUT_WIDTH;
    while ($curPos < $len1) {
	print FH substr($$seq1, $curPos, $inc) . "\n";
	print FH substr($$seq2, $curPos, $inc) . "\n";
	print FH substr($$homol, $curPos, $inc) . "\n";
	print FH "\n";
	$curPos += $inc;
    }

    close(FH);

    return 1;
}

#
# Output two non-aligned sequencess side by side.
# Pad shorter sequence with gaps as appropriate to ensure sequences are
# of equal length
#
sub WriteNonAlignedSeqs
{
    my($file, $seq1, $seq2) = @_;

    return 0 if !defined $file;
    return 0 if !defined $seq1;
    return 0 if !defined $seq2;

    return 0 if !open(OUTFH, ">$file");

    #
    # If one of the sequences is not defined, create the sequence
    # as a string of all gaps
    #
    my $seq1Len;
    my $seq2Len;
    if (!defined $$seq1) {
	$seq1Len = $seq2Len = length $$seq2;
	$$seq1 = '-' x $seq1Len;
    } elsif (!defined $$seq2) {
	$seq2Len = $seq1Len = length $$seq1;
	$$seq2 = '-' x $seq2Len;
    } else {
	$seq1Len = length $$seq1;
	$seq2Len = length $$seq2;
    }

    #
    # If the lengths of the sequences are different, pad the shorter
    # sequences with the gap character to the length of the longer
    # sequence.
    #
    my $lenDiff = $seq1Len - $seq2Len;
    my $outLen;
    if ($lenDiff >= 0) {
	$$seq2 .= '-' x $lenDiff;
	$outLen = $seq1Len;
    } elsif ($lenDiff < 0) {
	$$seq1 .= '-' x -$lenDiff;
	$outLen = $seq2Len;
    }

    #
    # Since we have no significant alignment between the two sequences
    # create a blank homology string (for compatibility with the alignment
    # files).
    #
    my $homol = ' ' x SEQ_OUT_WIDTH;

    #
    # Write sequences to output file SEQ_OUT_WIDTH chars at a time.
    #
    my $curPos = 0;
    while ($curPos < $outLen) {
	my $subSeq1 = substr($$seq1, $curPos, SEQ_OUT_WIDTH);
	my $subSeq2 = substr($$seq2, $curPos, SEQ_OUT_WIDTH);
	my $subSeq3 = substr($homol, 0, length $subSeq1);

	print OUTFH "$subSeq1\n";
	print OUTFH "$subSeq2\n";
	print OUTFH "$subSeq3\n";
	print OUTFH "\n";

	$curPos += SEQ_OUT_WIDTH;
    }
    close(OUTFH);

    return 1;
}

sub ChangeFileExt
{
    my ($fileName, $newExt) = @_;

    my $newName;

    my ($base, $dir, $ext) = fileparse($fileName, qr{\..*});

    $newName = "${dir}${base}.$newExt";

    return $newName;
}

1;
