#!/usr/local/bin/perl -w

=head1 NAME

enriched_tfbss.pl

=head1 SYNOPSIS

  enriched_tfbss.pl -i infile

=head1 ARGUMENTS

Arguments switches may be abbreviated where unique.

   -i infile	= File containing input sequences (from PAZAR)
   -o outfile	= File to which results are written
   -l logfile	= File to which logging info and debug messages are written
   -D		= Turn on debugging output

=head1 DESCRIPTION

Read input sequences and separate into long (regulatory regions) and short
(binding sites). Map short to long sequences.  If more than 1 long seq per
gene, combine if close enough.

for each (combined) long sequence:
	extend out to ~1kb (avoiding exons if possible) 
	extract coords of orthologous seq
	fetch orthologous seq
	align with orthologous seq

for each window length, L <- 10..20, L += 2
    for each ID score threshold T <- L/2..L, T++
	n_good_nt <- 0
	n_bad_nt <- 0
	n_good_tf <- 0
	n_bad_tf <- 0
    	for each alignment, A
	    seq <- A.seq1
	    sequence nucleotide, nt(0..seq.length - 1).good <- false
	    for each nucleotide position, i <- 0..seq.length - L (i++)
	        compute score, s of window, w of length L => s(w(i))
		if s(w(i)) >= T
		    w(i).good <- true
		    nt(i..i+L-1).good <- true
		else
		    w(i).good <- false

	    for each i <- 0..seq.length - 1 (i++)
	        if (nt(i).good)
		    n_good_nt++
		else
		    n_bad_nt++

	    for each tfbs, t
		t.good <- true
	        for i <- t.start..t.end, i++
		    if nt(i).good = false
		        t.good <- false
			break

		if t.good
		    n_good_tf++
		else
		    n_bad_tf++

	print L, T, (n_good_tf / n_good_nt) / (n_bad_tf / n_bad_tf)


=head1 AUTHOR

  David Arenillas
  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  E-mail: dave@cmmt.ubc.ca

=cut

use strict;

use lib '/usr/local/src/ensembl-48/ensembl/modules';

use POSIX qw(ceil floor);
use Getopt::Long;
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);
use Log::Dispatch::File;
use File::Temp qw/ :POSIX /;
use Bio::SeqIO;
use Bio::LocatableSeq;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Orca;
use ConservationAnalysis::ConservationAnalysis;

use constant DEBUG		=> 0;
use constant OPT_SEQ_LEN	=> 1000;
use constant MIN_WIN_LEN	=> 10;
use constant MAX_WIN_LEN	=> 20;
use constant WIN_LEN_STEP	=> 2;
use constant LOG_FILE 		=> 'enriched_tfbss.log';

use constant ENSEMBL_DB_HOST		=> 'napa.cmmt.ubc.ca';
use constant ENSEMBL_HS_DB_NAME		=> 'homo_sapiens_core_48_36j';
use constant ENSEMBL_MM_DB_NAME		=> 'mus_musculus_core_48_37a';
use constant ENSEMBL_DB_USER		=> 'ensembl_r';
use constant ENSEMBL_DB_PASS		=> '';

use constant UCSC_DB_HOST		=> 'genome-mysql.cse.ucsc.edu';
use constant UCSC_HG_DB_NAME		=> 'hg18';
use constant UCSC_MM_DB_NAME		=> 'mm9';
use constant UCSC_DB_USER		=> 'genome';
use constant UCSC_DB_PASS		=> '';

use constant LIFT_OVER_EXE		=> '/usr/local/bin/liftOver.linux.i386';
use constant LIFT_OVER_HG_MM		=> '/space/data/UCSCLiftOver/hg18ToMm9.over.chain';
use constant LIFT_OVER_MM_HG		=> '/space/data/UCSCLiftOver/mm9ToHg18.over.chain';
use constant LIFT_OVER_MIN_MATCH	=> 0.25;

my $in_file;
my $out_file;
my $log_file = LOG_FILE;
my $debug = DEBUG;
GetOptions(
    'i=s'	=> \$in_file,
    'o=s'	=> \$out_file,
    'l=s'	=> \$log_file,
    'D'		=> \$debug
);

if (!$in_file) {
    pod2usage(
	    -msg	=> "No input positive set CRM region file specified",
	    -verbose	=> 1);
}

my $logger = get_logger();
if ($debug) {
    $logger->level($DEBUG);
} else {
    $logger->level($INFO);
}
my $appender = Log::Log4perl::Appender->new("Log::Dispatch::File",
					    filename    => $log_file,
					    mode        => "write");
#my $layout = Log::Log4perl::Layout::PatternLayout->new("%d %M:%L %p: %m%n");
my $layout = Log::Log4perl::Layout::PatternLayout->new("%p: %m%n");
$appender->layout($layout);
$logger->add_appender($appender);

my $start_time = time;
my $localtime = localtime($start_time);

$logger->info("Date: $localtime");
$logger->info("");

my $ehs_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    			-host		=> ENSEMBL_DB_HOST,
			-user		=> ENSEMBL_DB_USER,
			-dbname		=> ENSEMBL_HS_DB_NAME,
			-pass		=> ENSEMBL_DB_PASS,
			-species	=> 'human');
if (!$ehs_db) {
    $logger->logdie("connecting to human Ensembl DB");
}

my $emm_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    			-host		=> ENSEMBL_DB_HOST,
			-user		=> ENSEMBL_DB_USER,
			-dbname		=> ENSEMBL_MM_DB_NAME,
			-pass		=> ENSEMBL_DB_PASS,
			-species	=> 'mouse');
if (!$emm_db) {
    $logger->logdie("connecting to mouse Ensembl DB");
}

my $hs_slice_adaptor = $ehs_db->get_SliceAdaptor;
if (!$hs_slice_adaptor) {
    $logger->logdie("getting human SliceAdaptor");
}

my $mm_slice_adaptor = $emm_db->get_SliceAdaptor;
if (!$mm_slice_adaptor) {
    $logger->logdie("getting mouse SliceAdaptor");
}

my $uhg_dsn = sprintf("dbi:mysql:database=%s;host=%s",
			    UCSC_HG_DB_NAME, UCSC_DB_HOST);
my $uhg_db = DBI->connect($uhg_dsn, UCSC_DB_USER, UCSC_DB_PASS);
if (!$uhg_db) {
    $logger->logdie("connecting to UCSC human DB " . UCSC_HG_DB_NAME);
}

my $umm_dsn = sprintf("dbi:mysql:database=%s;host=%s",
			    UCSC_MM_DB_NAME, UCSC_DB_HOST);
my $umm_db = DBI->connect($umm_dsn, UCSC_DB_USER, UCSC_DB_PASS);
if (!$umm_db) {
    $logger->logdie("connecting to UCSC human DB " . UCSC_MM_DB_NAME);
}

my ($prom_seqs, $tfbs_seqs) = read_seqs($in_file);
$prom_seqs = expand_seqs($prom_seqs);
# XXX combine sequences if overlapping?
my $prom_tfbs_map = map_tfbss_to_promoters($prom_seqs, $tfbs_seqs);
#$prom_seqs = fetch_alignments($prom_seqs);
$prom_seqs = lift_over($prom_seqs);

$logger->info("Orthologs:");
foreach my $pkey (keys %$prom_seqs) {
    my $prom = $prom_seqs->{$pkey};

    $logger->info(sprintf("%s", $prom->{-pazar_id}));
    $logger->info(sprintf("%s: chr%s:%d-%d",
    				$prom->{-species},
    				$prom->{-chr},
    				$prom->{-start},
    				$prom->{-end}));
    $logger->info(sprintf("%s: chr%s:%d-%d",
    				$prom->{-species2},
    				$prom->{-chr2},
    				$prom->{-start2},
    				$prom->{-end2}));
}

$prom_seqs = create_alignments($prom_seqs);

if ($out_file) {
    if (!open(OFH, ">$out_file")) { 
    	$logger->error("Opening output file $out_file\n");
	open(OFH, ">-");
    }
} else {
    open(OFH, ">-");
}

printf OFH "wlen\tscore\t%%id\tgood_nt\tbad_nt\tgood_tf\tbad_tf\tratio\tpseu+1\tpseu+.1\n";
my $wlen = MIN_WIN_LEN;
while ($wlen <= MAX_WIN_LEN) {
    $logger->info("Window length = $wlen");
    foreach my $wscore (int($wlen/2)..$wlen) {
        my $thresh = $wscore / $wlen;
	$logger->info("ID threshold = $thresh");

	my $n_good_nt = 0;
	my $n_bad_nt = 0;
	my $n_good_tf = 0;
	my $n_bad_tf = 0;
    	foreach my $pkey (keys %$prom_seqs) {
	    my $prom = $prom_seqs->{$pkey};

	    $logger->info("Promoter region = " . $prom->{-pazar_id});

	    if (!$prom->{-alignment}) {
	    	$logger->warn("No alignment - skipping");
		next;
	    }

	    my $ca = ConservationAnalysis::ConservationAnalysis->new(
			base_seq                => $prom->{-seq_obj},
			comparison_seq          => $prom->{-seq_obj2},
			alignment               => $prom->{-alignment}
		    );
	    if (!$ca) {
	    	$logger->error("Could not create ConservationAnalysis object");
		next;
	    }

	    $ca->param('window_size', $wlen);
	    $ca->param('min_filtered_cr_length', $wlen);

	    if (!$ca->compute_conservation_profile(position_type => 's')) {
	    	$logger->logdie("Could not compute conservation profile");
	    }

	    my $seq1 = $prom->{-seq};
	    my $seq1_len = length $seq1;

	    my @nts = map(0, 0 .. ($seq1_len - 1));
	    			
	    my $cp = $ca->conservation_profile;
	    foreach my $pos (0 .. $seq1_len - $wlen) {
		if ($cp->[$pos]->{score} >= $thresh) {
		    # can't use slice as lvalue?
		    #$nts[$pos .. $pos + $wlen - 1] = 1;
		    foreach my $i ($pos .. $pos + $wlen - 1) {
		    	$nts[$i] = 1;
		    }
		}
	    }

	    foreach my $i (0 .. $seq1_len - 1) {
	        if ($nts[$i]) {
		    $n_good_nt++
		} else {
		    $n_bad_nt++
		}
	    }

	    foreach my $tkey (@{$prom_tfbs_map->{$pkey}}) {
		my $tfbs = $tfbs_seqs->{$tkey};
		my $tfbs_good = 1;
	        foreach my $i ($tfbs->{-start} .. $tfbs->{-end}) {
		    if (!$nts[$i - $prom->{-start}]) {
		        $tfbs_good = 0;
			last;
		    }
		}

		if ($tfbs_good) {
		    $n_good_tf++;
		} else {
		    $n_bad_tf++;
		}
	    }
	}

	my $ratio;
	my $pseudo_ratio1;
	my $pseudo_ratio2;
	if ($n_good_nt == 0 || $n_bad_nt == 0) {
	    $ratio = "NaN";
	    $pseudo_ratio1 = "NaN";
	    $pseudo_ratio2 = "NaN";
	} else {
	    if ($n_bad_tf == 0) {
		$ratio = "NaN";
	    } else {
		$ratio = sprintf "%.3f",
		    ($n_good_tf / $n_good_nt) / ($n_bad_tf / $n_bad_nt);
	    }
	    $pseudo_ratio1 = sprintf "%.3f",
		(($n_good_tf + 1) / $n_good_nt)
			/ (($n_bad_tf + 1) / $n_bad_nt);
	    $pseudo_ratio2 = sprintf "%.3f",
		(($n_good_tf + 0.1) / $n_good_nt)
			/ (($n_bad_tf + 0.1) / $n_bad_nt);
	}

	printf OFH "%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n",
		    $wlen, $wscore, $thresh,
		    $n_good_nt, $n_bad_nt,
		    $n_good_tf, $n_bad_tf,
		    $ratio, $pseudo_ratio1, $pseudo_ratio2;
    }

    $wlen += WIN_LEN_STEP;
}
close(OFH);

sleep 1;

exit;

sub read_seqs
{
    my ($file) = @_;

    open(INFH, "$in_file")
	|| die "Error opening input PAZAR sequence file $in_file\n";

    my %lseqs;	# long seqs
    my %sseqs;	# short seqs (binding sites)

    my $line = <INFH>;	# read past header
    while ($line = <INFH>) {
    	chomp $line;
	my ($pazar_id, $seq, $chr, $start, $end, $strand, $species,
		$pazar_gene_id, $ensembl_id, $gene_name) = split "\t", $line;
		

    	if ($end - $start > 30) {
	    $lseqs{$pazar_id} = {
				    -pazar_id => $pazar_id,
				    -seq => $seq,
				    -chr => $chr,
				    -start => $start,
				    -end => $end,
				    -strand => $strand,
				    -species => $species,
				    -pazar_gene_id => $pazar_gene_id,
				    -ensembl_id => $ensembl_id,
				    -gene_name => $gene_name};
	} else {
	    $sseqs{$pazar_id} = {
				    -pazar_id => $pazar_id,
				    -seq => $seq,
				    -chr => $chr,
				    -start => $start,
				    -end => $end,
				    -strand => $strand,
				    -species => $species,
				    -pazar_gene_id => $pazar_gene_id,
				    -ensembl_id => $ensembl_id,
				    -gene_name => $gene_name};
	}
    }

    return (\%lseqs, \%sseqs);
}

sub map_tfbss_to_promoters
{
    my ($proms, $tfbss) = @_;

    my %map;
    foreach my $pkey (keys %$proms) {
    	my $prom = $proms->{$pkey};
	foreach my $tkey (keys %$tfbss) {
	    my $tfbs = $tfbss->{$tkey};
	    if ($prom->{-pazar_gene_id} eq $tfbs->{-pazar_gene_id}
	    	&& $tfbs->{-start} >= $prom->{-start}
	    	&& $tfbs->{-end} <= $prom->{-end})
	    {
	    	push @{$map{$prom->{-pazar_id}}}, $tfbs->{-pazar_id};
	    }
	}
    }

    return \%map;
}

sub expand_seqs
{
    my ($seqs) = @_;

    foreach my $key (keys %$seqs) {
	my $seq = $seqs->{$key};

	my $pazar_id = $seq->{-pazar_id};
	my $species = $seq->{-species};
    	my $chr = $seq->{-chr};
    	my $start = $seq->{-start};
    	my $end = $seq->{-end};
	my $len = $end - $start + 1;
	my $xlen = OPT_SEQ_LEN - $len;

	my $nstart = $start - $xlen;
	my $nend = $end + $xlen;
	my $nlen = $nend - $nstart + 1;

	$logger->info("Original sequence $species chr$chr:$start-$end\n");

	my $slice;
	if ($species eq 'HOMO SAPIENS') {
	    $slice = $hs_slice_adaptor->fetch_by_region('chromosome', $chr,
	    						    $nstart, $nend);
	} elsif ($species eq 'MUS MUSCULUS') {
	    $slice = $mm_slice_adaptor->fetch_by_region('chromosome', $chr,
	    						    $nstart, $nend);
	} else {
	    $logger->logdie("Illegal species, $species");
	}

	my $exons = $slice->get_all_Exons;
	my $mid = $nlen / 2;
	foreach my $exon (@$exons) {
	    my $ex_rstart = $exon->start;	# relative start
	    my $ex_rend = $exon->end;		# relative end
	    my $ex_astart = $slice->start + $ex_rstart - 1;	# abs. start
	    my $ex_aend = $slice->start + $ex_rend - 1;		# abs. end

	    if ($ex_aend >= $nstart && $ex_rstart < $mid && $ex_aend < $start)
	    {
	    	$nstart = $ex_aend + 1;
	    } elsif ($ex_astart < $nend && $ex_rend > $mid
			&& $ex_astart > $end)
	    {
	    	$nend = $ex_astart - 1;
	    }
	}
	$nlen = $nend - $nstart + 1;

	$logger->info(
		"Expanded sequence $species chr$chr:$nstart-$nend ($nlen)");

	if ($nlen < $len) {
	    $logger->error("New $pazar_id sequence length, $nlen, is below"
	    		    . " original length, $len");
	} elsif ($nlen < 1000) {
	    $logger->info("New $pazar_id sequence length, $nlen, is below"
	    		    . " optimal length");
	} else {
	    # XXX trim to optimal length
	    $logger->info("New $pazar_id sequence length, $nlen, is above"
			    . " optimal length");
	}

	if ($species eq 'HOMO SAPIENS') {
	    $slice = $hs_slice_adaptor->fetch_by_region('chromosome', $chr,
	    						    $nstart, $nend);
	} elsif ($species eq 'MUS MUSCULUS') {
	    $slice = $mm_slice_adaptor->fetch_by_region('chromosome', $chr,
	    						    $nstart, $nend);
	}

	$seq->{-start} = $nstart;
	$seq->{-end} = $nend;
	$seq->{-seq} = $slice->seq;
    }

    return $seqs;
}

sub fetch_alignments
{
    my ($proms) = @_;

    foreach my $pkey (keys %$proms) {
    	my $prom = $proms->{$pkey};

    	my $species1 = $prom->{-species};
    	my $chr1 = $prom->{-chr};
    	my $start1 = $prom->{-start};
    	my $end1 = $prom->{-end};
    	my $strand1 = $prom->{-strand};

	my ($chr2,
	    $start2,
	    $end2) = compute_aligned_seq_bounds(
					    $species1, $chr1, $start1, $end1);

	if (!$chr2) {
	    $logger->error("Could not compute orthologous species seq bounds");
	    next;
	}

	my $species2;
	my $slice;
	if ($species1 eq 'HOMO SAPIENS') {
	    $species2 = 'MUS MUSCULUS';
	    $slice = $mm_slice_adaptor->fetch_by_region('chromosome', $chr2,
	    						    $start2, $end2);
	} elsif ($species1 eq 'MUS MUSCULUS') {
	    $species2 = 'HOMO SAPIENS';
	    $slice = $hs_slice_adaptor->fetch_by_region('chromosome', $chr2,
	    						    $start2, $end2);
	} else {
	    $logger->logdie("Illegal species $species1");
	}

	if (!$slice) {
	    $logger->error("Could not fetch orthologous species slice");
	    next;
	}

	my $seq_obj1 = Bio::LocatableSeq->new(
				    -display_id	=> $prom->{-pazar_id},
				    -start	=> $start1,
				    -end	=> $end1,
				    -strand	=> $strand1 eq '-' ? -1 : 1,
				    -alphabet	=> 'dna',
				    -seq	=> $prom->{-seq});

	my $seq_obj2 = Bio::LocatableSeq->new(
				    -display_id	=> $prom->{-pazar_id} . '_orth',
				    -start	=> $start2,
				    -end	=> $end2,
				    -strand	=> 1,
				    -alphabet	=> 'dna',
				    -seq	=> $slice->seq);

	my $orca = Orca->new();
	my $aln = $orca->align(-seq1 => $seq_obj1, -seq2 => $seq_obj2);

	if (!$aln) {
	    $logger->error("Could not create alignment");
	    next;
	}

	$prom->{-chr2} = $chr2;
	$prom->{-start2} = $start2;
	$prom->{-end2} = $end2;
	$prom->{-strand2} = '+';
	$prom->{-seq2} = $slice->seq;
	$prom->{-seq_obj} = $seq_obj1;
	$prom->{-seq_obj2} = $seq_obj2;
	$prom->{-alignment} = $aln;
    }

    return $proms;
}

sub lift_over
{
    my ($proms) = @_;

    my %spec_proms;
    foreach my $pkey (keys %$proms) {
    	my $prom = $proms->{$pkey};

	push @{$spec_proms{$prom->{-species}}}, $prom;
    }

    foreach my $spec (keys %spec_proms) {
    	my ($fh, $bed_file) = tmpnam();

	foreach my $prom (@{$spec_proms{$spec}}) {
	    printf $fh "chr%s\t%d\t%d\n",
	    		$prom->{-chr}, $prom->{-start}, $prom->{-end};
	}
	close $fh;

	my $out_file = tmpnam();

	my $cmd = sprintf
		    "%s $bed_file %s $out_file /tmp/unMapped -minMatch=%s",
				LIFT_OVER_EXE,
				$spec eq 'HOMO SAPIENS'
					? LIFT_OVER_HG_MM : LIFT_OVER_MM_HG,
				LIFT_OVER_MIN_MATCH;

    	system $cmd;

	open(IFH, "<$out_file");
	my $idx = 0;
	while (my $line = <IFH>) {
	    chomp $line;
	    my ($chr, $start, $end) = split /\s+/, $line;

	    my $prom = $spec_proms{$spec}[$idx]; 
	    $prom->{-species2} = $spec eq 'HOMO SAPIENS'
					? 'MUS MUSCULUS' : 'HOMO SAPIENS';
	    if ($chr =~ /chr(\S+)/) {
	    	$chr = $1;
	    }
	    $prom->{-chr2} = $chr;
	    $prom->{-start2} = $start;
	    $prom->{-end2} = $end;

	    $idx++;
	}
	close IFH;
    }

    return $proms;
}

sub create_alignments
{
    my ($proms) = @_;

    foreach my $pkey (keys %$proms) {
    	my $prom = $proms->{$pkey};

    	my $species1 = $prom->{-species};
    	my $chr1 = $prom->{-chr};
    	my $start1 = $prom->{-start};
    	my $end1 = $prom->{-end};
    	my $strand1 = $prom->{-strand};

    	my $species2 = $prom->{-species2};
    	my $chr2 = $prom->{-chr2};
    	my $start2 = $prom->{-start2};
    	my $end2 = $prom->{-end2};

	my $slice;
	if ($species1 eq 'HOMO SAPIENS') {
	    $slice = $mm_slice_adaptor->fetch_by_region('chromosome', $chr2,
	    						    $start2, $end2);
	} elsif ($species1 eq 'MUS MUSCULUS') {
	    $slice = $hs_slice_adaptor->fetch_by_region('chromosome', $chr2,
	    						    $start2, $end2);
	} else {
	    $logger->logdie("Illegal species $species1");
	}

	if (!$slice) {
	    $logger->error("Could not fetch orthologous species slice");
	    next;
	}

	my $seq_obj1 = Bio::LocatableSeq->new(
			    -display_id	=> $prom->{-pazar_id},
			    -start	=> $start1,
			    -end	=> $end1,
			    -strand	=> $strand1 eq '-' ? -1 : 1,
			    -alphabet	=> 'dna',
			    -seq	=> $prom->{-seq});

	my $seq_obj2 = Bio::LocatableSeq->new(
			    -display_id	=> $prom->{-pazar_id} . '_orth',
			    -start	=> $start2,
			    -end	=> $end2,
			    -strand	=> 1,
			    -alphabet	=> 'dna',
			    -seq	=> $slice->seq);

	my $orca = Orca->new();
	my $aln = $orca->align(-seq1 => $seq_obj1, -seq2 => $seq_obj2);

	if (!$aln) {
	    $logger->error("Could not create alignment");
	    next;
	}

	$prom->{-seq2} = $slice->seq;
	$prom->{-seq_obj} = $seq_obj1;
	$prom->{-seq_obj2} = $seq_obj2;
	$prom->{-alignment} = $aln;
    }


    return $proms;
}

#
# Compute aligned sequence bounds using UCSC 'net<species>' table
#
sub compute_aligned_seq_bounds
{
    my ($species, $chr, $start, $end) = @_;

    $logger->info(
	"Computing ortholog alignment bounds for $species chr$chr:$start-$end");

    # first look for level 1 alignment that spans entire target region
    my $species2;
    my $sth;
    if ($species eq 'HOMO SAPIENS') {
    	$species2 = 'MUS MUSCULUS';
	my $sql = qq{select tname, tstart, tend, qname, qstart, qend, strand
		    from netMm9
		    where tname = 'chr$chr'
		    and tstart <= $start and tend >= $end
		    and level = 1};
	$sth = $uhg_db->prepare($sql);
    } elsif ($species eq 'MUS MUSCULUS') {
    	$species2 = 'HOMO SAPIENS';
	my $sql = qq{select tname, tstart, tend, qname, qstart, qend, strand
		    from netHg18
		    where tname = 'chr$chr'
		    and tstart <= $start and tend >= $end
		    and level = 1};
	$sth = $umm_db->prepare($sql);
    } else {
	$logger->logdie("Unknown species $species");
    }
    
    if (!$sth) {
    	$logger->logdie("Preparing UCSC alignment query - " . $DBI::errstr);
    }

    if (!$sth->execute()) {
    	$logger->logdie("Executing UCSC alignment query - " . $DBI::errstr);
    }

    my $qname;
    my $qstart;
    my $qend;
    my $strand;
    my $count = 0;
    while (my @row = $sth->fetchrow_array()) {
    	$count++;
	if ($count > 1) {
	    $logger->error("More than one level 1 alignment spanning region");
	    return;
	}
	$qname = $row[3];
	$qstart = $row[4];
	$qend = $row[5];
	$strand = $row[6];
    }

    if ($species eq 'HOMO SAPIENS') {
	my $sql = qq{select tname, tstart, tend, qname, qstart, qend, strand
		    from netMm9
		    where tname = 'chr$chr'
		    and qname = '$qname'
		    and tstart <= $end and tend >= $start
		    and level = 2};
	$sth = $uhg_db->prepare($sql);
    } elsif ($species eq 'MUS MUSCULUS') {
	my $sql = qq{select tname, tstart, tend, qname, qstart, qend, strand
		    from netHg18
		    where tname = 'chr$chr'
		    and qname = '$qname'
		    and tstart <= $end and tend >= $start
		    and level = 2};
	$sth = $umm_db->prepare($sql);
    }
    
    if (!$sth) {
    	$logger->logdie("Preparing UCSC alignment query - " . $DBI::errstr);
    }

    if (!$sth->execute()) {
    	$logger->logdie("Executing UCSC alignment query - " . $DBI::errstr);
    }

    my $min_tstart = 999999999;
    my $min_tend = 0; 
    my $min_qstart = 999999999;
    my $min_qend = 0; 
    my $start_qstrand;
    my $end_qstrand;
    $count = 0;
    while (my @row = $sth->fetchrow_array()) {
    	$count++;

	my $tname	= $row[0];
	my $tstart	= $row[1];
	my $tend	= $row[2];
	my $qstart	= $row[4];
	my $qend	= $row[5];
	my $strand	= $row[6];

	if ($strand eq '+') {
	    if ($tstart < $min_tstart) {
	    	$min_tstart = $tstart;
	    	$min_qstart = $qstart;
		$start_qstrand = $strand;
	    }
	    if ($tend > $min_tend) {
	    	$min_tend = $tend;
	    	$min_qend = $qend;
		$end_qstrand = $strand;
	    }
	} elsif ($strand eq '-') {
	    if ($tstart < $min_tstart) {
	    	$min_tstart = $tstart;
	    	$min_qend = $qend;
		$end_qstrand = $strand;
	    }
	    if ($tend > $min_tend) {
	    	$min_tend = $tend;
	    	$min_qstart = $qstart;
		$start_qstrand = $strand;
	    }
	} else {
	    $logger->logdie("Unknown strand $strand");
	}
    }

    if ($count == 0) {
    	$logger->error("No level 2 alignments found");
	return;
    }

    if ($min_qstart == 999999999) {
    	$logger->error("Could not determine min. query start");
	return;
    }

    if ($min_qend == 0) {
    	$logger->error("Could not determine min. query end");
	return;
    }

    if ($min_qstart > $min_qend) {
    	$logger->error("Inconsistent min. query start and end");
	return;
    }

    # XXX is this really an error?
    if ($start_qstrand ne $end_qstrand) {
    	$logger->error("Start and end strands do not match");
	return;
    }

    # extend ends
    if ($min_tstart > $start) {
    	my $diff = $min_tstart - $start;
	if ($start_qstrand eq '+') {
	    $min_qstart -= $diff * 2;
	} else {
	    $min_qend += $diff * 2;
	}
    }

    # extend ends
    if ($min_tend < $end) {
    	my $diff = $end - $min_tend;
	if ($end_qstrand eq '+') {
	    $min_qend += $diff * 2;
	} else {
	    $min_qstart -= $diff * 2;
	}
    }

    $logger->info(
	"Ortholog alignment bounds: $species2 $qname:$min_qstart-$min_qend");

    if ($qname =~ /^chr(\S+)/) {
    	$qname = $1;
    }

    return ($qname, $min_qstart, $min_qend);
}
