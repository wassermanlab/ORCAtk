=head1 NAME

ORCA::Web::Application - The main ORCA web application module

=head1 SYNOPSIS

  use ORCA::Web::Application;

  my $app = ORCA::Web::Application->new();
  $app->run;

=head1 AUTHOR

  David Arenillas (dave@cmmt.ubc.ca)

=head1 COPYRIGHT

  Wasserman Lab
  Centre for Molecular Medicine and Therapeutics
  University of British Columbia

  Distributed under the terms of the GNU General Public License (GPL)

=cut

package ORCA::Web::Application;

use strict;

use base 'CGI::Application';

use ORCA::Web::Options;

#use lib ORCA_LIB;
use lib ENSEMBL_LIB;
use lib ENSEMBL_COMPARA_LIB;

use Data::Dumper;    # for debugging only
use Template;
use CGI;

use DBI;
use GD;
use File::Temp qw/ tempfile /;
use Bio::LocatableSeq;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;
use Bio::SeqIO;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Tools::Est2Genome;
use Bio::Graphics::Panel;
use TFBS::DB::JASPAR5;
use TFBS::MatrixSet;
use TFBS::Matrix::PFM;
use ORCA::Web::State;
use ORCA::Analysis::PhastCons;
use ORCA::Analysis::Pairwise;
use ORCA::Graphics::PhastCons;
use ORCA::Graphics::Pairwise;

use CGI::Carp qw(carpout);    # fatalsToBrowser;

$CGI::POST_MAX = MAX_FILE_SIZE;

my $log_file;
if ($ENV{'USER'}) {
    $log_file = sprintf "%s/ORCAtk_%s.log", ABS_LOG_PATH, $ENV{'USER'};
} else {
    $log_file = sprintf "%s/ORCAtk.log", ABS_LOG_PATH;
}

open(LOG, ">>$log_file") || die "Could not open log file $log_file.\n";

carpout(\*LOG);

use constant DEBUG => 1;

sub setup
{
    my $self = shift;

    #printf STDERR "[%s] setup\n", scalar localtime(time);

    #$Carp::Verbose = 1;

    $self->{'errors'} = [];

    $self->mode_param('rm');
    $self->start_mode('home');
    $self->run_modes(
        'home'                     => 'analysis_start',
        'analysis_start'           => 'analysis_start',
        'new_analysis'             => 'analysis_start',
        'select_seq_input_method'  => 'select_seq_input_method',
        'select_seqs_paste_upload' => 'select_seqs_paste_upload',
        'select_seq1_coords'       => 'select_seq1_coords',
        'select_seq2_coords'       => 'select_seq2_coords',
        'select_seq1_name'         => 'select_seq1_name',
        'select_gene'              => 'select_gene',
        'select_gene1'             => 'select_gene1',
        'select_ortholog'          => 'select_ortholog',
        'select_transcripts'       => 'select_transcripts',
        'select_transcript1'       => 'select_transcript1',
        'select_ca_params'         => 'select_ca_params',
        'select_tfs'               => 'select_tfs',
        'results'                  => 'results',
        'error'                    => 'error',
        'warning'                  => 'warning',
        'contact'                  => 'contact'
    );

    my $q         = $self->query;
    my $rm        = $q->param('rm');
    my $from      = $q->param("from");
    my $sid       = $q->param('sid');

    $rm = "analysis_start" if !$rm;

    #printf STDERR "\nrm = $rm\n";
    #printf STDERR "from = $from\n" if $from;

    my $new_analysis = 0;
    if (  !$sid
        || $rm eq 'new_analysis'
        || $rm eq 'analysis_start'
        || $rm eq 'home'
        || $rm eq ''
        || !defined $rm)
    {
        $new_analysis = 1;
    }
    #printf STDERR "new_analysis = $new_analysis\n";

    if ($new_analysis) {
        #
        # If we have an existing session, expire it (delete the temp file).
        #
        if ($sid) {
            #my $state = ORCA::Web::State->new(__Fn => _session_tmp_file($sid));
            #$state->expire();
            $self->{-state} = undef;
        }

        #
        # Create a new session ID and associated state object.
        #
        $sid = $$ . time;
        my $state = ORCA::Web::State->new(
            -sid    => $sid,
            __Fn    => _session_tmp_file($sid)
        );

        $state = $self->initialize_state($state);

        $self->{-state}     = $state;
        $self->{-errors}    = undef;
        $self->{-warnings}  = undef;
    } else {
        #
        # Continuing with existing analysis
        #

        #
        # Check for outstanding errors
        #
        if ($q->cgi_error()) {
            $self->add_error("CGI error: " . $q->cgi_error());
            $q->cgi_error(undef);
        }

        if ($self->errors) {
            return $self->error();
        }

        #
        # Retrieve the current state.
        #
        my $state = ORCA::Web::State->new(__Fn => _session_tmp_file($sid));
        if (!$state) {
            return $self->error("Error retrieving state for session $sid");
        }
        $self->{-state} = $state;

        #
        # Check "from" param and execute "event driven" routines here
        #
        if ($from) {
            if ($from eq "select_seqs_paste_upload") {
                if (!$self->seqs_selected_paste_upload) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_seq1_coords") {
                if (!$self->seq1_selected_coords) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_seq2_coords") {
                if (!$self->seq2_selected_coords) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_seq1_name") {
                if (!$self->seq1_selected_name) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_gene") {
                if (!$self->gene_selected) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_gene1") {
                if (!$self->gene1_selected) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_ortholog") {
                if (!$self->ortholog_selected) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_transcripts") {
                if (!$self->transcripts_selected) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_transcript1") {
                if (!$self->transcript1_selected) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_ca_params") {
                #printf STDERR "calling ca_params_selected\n";
                if (!$self->ca_params_selected) {
                    $q->param('rm', 'error');
                }
            } elsif ($from eq "select_tfs") {
                if (!$self->tfs_selected) {
                    $q->param('rm', 'error');
                }
            }
        }
    }
}

sub teardown
{
    my $self = shift;

    #
    # Save current state
    #
    my $state = $self->state();
    if ($state) {
        $state->dumper->Purity(1);
        $state->dumper->Deepcopy(1);
        $state->commit();
    }

    $self->clean_tempfiles;
}

sub analysis_start
{
    my $self = shift;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Start Analysis',
        sid              => $self->state->sid,
        var_template     => "analysis_start.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_seq_input_method
{
    my $self = shift;

    my $state = $self->state;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Select Sequence Input Method',
        sid              => $state->sid,
        var_template     => "select_seq_input_method.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_seqs_paste_upload
{
    my $self = shift;

    #printf STDERR "[%s] select_seqs_paste_upload\n", scalar localtime(time);

    my $q = $self->query;

    my $seq1 = $q->param('seq1');

    my $state = $self->state;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Paste/Upload Sequences)',
        sid              => $state->sid,
        in_seq1          => $seq1,
        #in_seq_file1        => $state->seq_file1,
        #in_seq2             => $state->seq2,
        #in_seq_file2        => $state->seq_file2,
        #in_cdna             => $state->cdna,
        #in_cdna_file        => $state->cdna_file,
        var_template => "select_seqs_paste_upload.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_seq1_coords
{
    my $self = shift;

    #printf STDERR "[%s] select_seq1_coords\n", scalar localtime(time);

    my $state = $self->state;
    my $q     = $self->query;

    my $latin_name    = $q->param('species');
    my $in_seq_chr1   = $q->param('chr');
    my $in_seq_start1 = $q->param('start');
    my $in_seq_end1   = $q->param('end');
    my $in_species1;

    if ($latin_name) {
        $in_species1 =
            $state->species_latin_to_common_names()->{lc $latin_name};
        $state->species1($in_species1);
    } else {
        $in_species1 = $state->species1;
    }
    #printf STDERR "select_seq1_coords: in_species1 = $in_species1\n";
    #printf STDERR "select_seq1_coords: in_seq_chr1 = $in_seq_chr1\n";
    #printf STDERR "select_seq1_coords: in_seq_start1 = $in_seq_start1\n";
    #printf STDERR "select_seq1_coords: in_seq_end1 = $in_seq_end1\n";

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        species          => $state->species_common_names,
        ensembl_dbs      => $state->species_ensembl_dbs,
        title            => 'ORCAtk: Select Sequence by Genomic Coordinates)',
        sid              => $state->sid,

        in_species1   => $in_species1,
        in_seq_chr1   => $in_seq_chr1,
        in_seq_start1 => $in_seq_start1,
        in_seq_end1   => $in_seq_end1,

        var_template => "select_seq1_coords.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_seq2_coords
{
    my $self = shift;

    #printf STDERR "[%s] select_seq2_coords\n", scalar localtime(time);

    my $state = $self->state;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        species          => $state->species_common_names,
        ensembl_dbs      => $state->species_ensembl_dbs,
        title            => 'ORCAtk: Select Orthologoud Sequence by Genomic'
            . ' Coordinates)',
        sid => $state->sid,

        in_species2   => $state->species2,
        in_seq_chr2   => $state->seq_chr2,
        in_seq_start2 => $state->seq_start2,
        in_seq_end2   => $state->seq_end2,

        var_template => "select_seq2_coords.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_seq1_name
{
    my $self = shift;

    #printf STDERR "[%s] select_seq1_name\n", scalar localtime(time);

    my $state = $self->state;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        species          => $state->species_common_names,
        ensembl_dbs      => $state->species_ensembl_dbs,
        title            => 'ORCAtk: Select Gene Name/Symbol)',
        sid              => $state->sid,
        in_gene_name     => $state->gene_name1,
        in_species1      => $state->species1,
        var_template     => "select_seq1_name.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_gene
{
    my $self = shift;

    #printf STDERR "[%s] select_gene\n", scalar localtime(time);

    my $state = $self->state;

    my $species1   = $state->species1;
    my $species2   = $state->species2;
    my $gene_name1 = $state->gene_name1;
    my $wc         = $state->use_wildcards;

    my $genes =
        $self->fetch_ensembl_genes_by_name($species1, $gene_name1, $wc);
    if (!$genes || !$genes->[0]) {
        return $self->warning("There were no $species1 genes found with"
                . " name/symbol '$gene_name1'");
    }

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        species          => $state->species_common_names,
        ensembl_dbs      => $state->species_ensembl_dbs,
        title            => 'ORCAtk: Select Gene',
        sid              => $state->sid,
        in_gene_list     => $genes,
        species1         => $species1,
        in_species2      => $species2,
        var_template     => "select_gene.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_gene1
{
    my $self = shift;

    #printf STDERR "[%s] select_gene1\n", scalar localtime(time);

    my $q = $self->query;

    my $ensembl_id1 = $q->param('ensembl_id');

    my $state = $self->state;

    my $species1 = $state->species1;
    #printf STDERR "species1 = $species1\n";
    my $species2   = $state->species2;
    my $gene_name1 = $state->gene_name1;
    my $wc         = $state->use_wildcards;

    # Entry from PAZAR - latin species name passed
    my $latin_name = $q->param('species');
    #printf STDERR "species1 latin name = $latin_name\n";
    if ($latin_name) {
        #printf STDERR "Latin name = $latin_name\n";
        $species1 = $state->species_latin_to_common_names()->{lc $latin_name};
        if (!$species1) {
            return $self->error(
                "Could not determine species common name from latin name $latin_name"
            );
        }
        #printf STDERR "species1 = $species1\n";
        $state->species1($species1);
    }

    my $genes;
    if ($ensembl_id1) {
        my $gene =
            $self->fetch_ensembl_gene_by_stable_id($species1, $ensembl_id1);
        if (!$gene) {
            return $self->warning(
                "Could not retrieve Ensembl $species1 gene '$ensembl_id1'");
        }
        push @$genes, $gene;
    } else {
        #printf STDERR "species1 = $species1; gene_name1 = $gene_name1\n";
        $genes =
            $self->fetch_ensembl_genes_by_name($species1, $gene_name1, $wc);
        if (!$genes || !$genes->[0]) {
            return $self->warning(
                "There were no $species1 genes found with name/symbol '$gene_name1'"
            );
        }
    }

    my $specref = $state->species_common_names;
    my @species = @$specref;
    unshift @species, 'Select orthologous species';

    my $phastcons_species = SPECIES_TRACK_NAMES;
    my $supports_phastcons = 0;
    if ($phastcons_species->{$species1}) {
        $supports_phastcons = 1;
    }

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Select Gene',
        sid              => $state->sid,
        species          => \@species,
        ensembl_dbs      => $state->species_ensembl_dbs,
        in_gene_list     => $genes,
        in_species1      => $species1,
        supports_phastcons => $supports_phastcons,
        var_template     => "select_gene1.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_ortholog
{
    my $self = shift;

    #printf STDERR "[%s] select_ortholog\n", scalar localtime(time);

    my $state = $self->state;

    my $species1 = $state->species1;
    my $species2 = $state->species2;
#printf STDERR "select_ortholog: species1 = $species1; species2 = $species2\n";
    my $gene_name1 = $state->gene_name1;

    # get orthologs
    my $genes = $self->fetch_ensembl_orthologs_compara();
    if (!$genes || !$genes->[0]) {
        return $self->warning(
            "There were no $species2 genes found orthologous to $species1 gene '$gene_name1'"
        );
    }

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Select Ortholog',
        sid              => $state->sid,
        species2         => $species2,
        in_gene_list     => $genes,
        var_template     => "select_ortholog.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_transcripts
{
    my $self = shift;

    #printf STDERR "[%s] select_transcripts\n", scalar localtime(time);

    my $state = $self->state;
    my $sid   = $state->sid;

    my $species1    = $state->species1;
    my $species2    = $state->species2;
    my $ensembl_id1 = $state->ensembl_gene_id1;
    my $ensembl_id2 = $state->ensembl_gene_id2;

    my $transcripts1 =
        $self->fetch_ensembl_transcripts($species1, $ensembl_id1);
    my $transcripts2 =
        $self->fetch_ensembl_transcripts($species2, $ensembl_id2);

    #
    # Flip orientation of graph it seq1 transcripts are on - strand.
    #
    #if ($transcripts1->[0]->strand == -1) {
    #    $state->flip_graph(1);
    #}

    my $trans_panel1 =
        $self->create_transcript_panel($transcripts1, $species1);
    my $trans_panel2 =
        $self->create_transcript_panel($transcripts2, $species2);

    my $panel_file1     = "${sid}_trans_panel1.png";
    my $panel_file2     = "${sid}_trans_panel2.png";
    my $panel_rel_path1 = REL_TMP_PATH . "/$panel_file1";
    my $panel_rel_path2 = REL_TMP_PATH . "/$panel_file2";
    my $panel_abs_path1 = ABS_TMP_PATH . "/$panel_file1";
    my $panel_abs_path2 = ABS_TMP_PATH . "/$panel_file2";

    if (!open(FH, ">$panel_abs_path1")) {
        return $self->error("Could not create transcript 1 panel file");
    }
    print FH $trans_panel1->gd->png;
    close FH;

    if (!open(FH, ">$panel_abs_path2")) {
        return $self->error("Could not create transcript 2 panel file");
    }
    print FH $trans_panel2->gd->png;
    close FH;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Select Transcripts',
        sid              => $state->sid,
        species1         => $species1,
        species2         => $species2,
        transcripts1     => $transcripts1,
        transcripts2     => $transcripts2,
        panel_file1      => $panel_rel_path1,
        panel_file2      => $panel_rel_path2,
        in_upstream_bp   => $state->seq_upstream_bp,
        in_downstream_bp => $state->seq_downstream_bp,
        in_down_rel_to   => $state->seq_down_rel_to,
        var_template     => "select_transcripts.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_transcript1
{
    my $self = shift;

    #printf STDERR "[%s] select_transcript1\n", scalar localtime(time);

    my $state = $self->state;
    my $sid   = $state->sid;

    my $q = $self->query;

    my $species1    = $q->param('species')    || $state->species1;
    my $ensembl_id1 = $q->param('ensembl_id') || $state->ensembl_gene_id1;

    my $transcripts1 =
        $self->fetch_ensembl_transcripts($species1, $ensembl_id1);

    #
    # Flip orientation of graph it seq1 transcripts are on - strand.
    #
    #if ($transcripts1->[0]->strand == -1) {
    #    $state->flip_graph(1);
    #}

    my $trans_panel1 =
        $self->create_transcript_panel($transcripts1, $species1);

    my $panel_file1     = "${sid}_trans_panel1.png";
    my $panel_rel_path1 = REL_TMP_PATH . "/$panel_file1";
    my $panel_abs_path1 = ABS_TMP_PATH . "/$panel_file1";

    if (!open(FH, ">$panel_abs_path1")) {
        return $self->error("Could not create transcript panel file");
    }
    print FH $trans_panel1->gd->png;
    close FH;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Select Transcript',
        sid              => $state->sid,
        species1         => $species1,
        transcripts1     => $transcripts1,
        panel_file1      => $panel_rel_path1,
        in_upstream_bp   => $state->seq_upstream_bp,
        in_downstream_bp => $state->seq_downstream_bp,
        in_down_rel_to   => $state->seq_down_rel_to,
        var_template     => "select_transcript1.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_ca_params
{
    my $self = shift;

    #printf STDERR "[%s] select_ca_params\n", scalar localtime(time);

    my $state    = $self->state;
    my $species1 = $state->species1;

    if ($state->analysis_type eq 'pairwise') {
        # top percentile take precedence
        $state->ca_min_conservation(CA_MIN_CONS_PAIRWISE);
        $state->ca_min_cr_length(CA_MIN_CR_LEN_PAIRWISE);
    } elsif ($state->analysis_type eq 'phastcons') {
        $state->ca_min_conservation(CA_MIN_CONS_PHASTCONS);
        $state->ca_min_cr_length(CA_MIN_CR_LEN_PHASTCONS);
    } else {
        return $self->error("Analysis type is undefined!");
    }

    my $vars = {
        abs_html_path       => ABS_HTML_PATH,
        rel_html_path       => REL_HTML_PATH,
        abs_cgi_bin_path    => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path    => REL_CGI_BIN_PATH,
        version             => VERSION,
        title               => 'ORCAtk: Select Conservation Parameters',
        sid                 => $state->sid,
        analysis_type       => $state->analysis_type,
        in_top_percentile   => $state->ca_top_percentile,
        in_min_conservation => $state->ca_min_conservation,
        in_filter_exons     => $state->ca_filter_exons,
        in_window_size      => $state->ca_window_size,
        in_min_cr_length    => $state->ca_min_cr_length,
        var_template        => "select_ca_params.html"
    };

    return $self->process_template('master.html', $vars);
}

sub select_tfs
{
    my $self = shift;

    #printf STDERR "[%s] select_tfs\n", scalar localtime(time);

    my $state = $self->state;

    my %in_core_tf_ids_selected;
    if ($state->tf_core_ids) {
        foreach my $tf_id (@{$state->tf_core_ids}) {
            $in_core_tf_ids_selected{$tf_id} = 1;
        }
    }

    #my %in_phylo_tf_ids;
    #if ($state->tf_phylo_ids) {
    #    foreach my $tf_id (@{$state->tf_phylo_ids}) {
    #        $in_phylo_tf_ids{$tf_id} = 1;
    #    }
    #}

    #my %in_fam_tf_ids;
    #if ($state->tf_fam_ids) {
    #    foreach my $tf_id (@{$state->tf_fam_ids}) {
    #        $in_fam_tf_ids{$tf_id} = 1;
    #    }
    #}

    my %in_tax_groups_selected;
    if ($state->tf_tax_groups) {
        foreach my $tax_group (@{$state->tf_tax_groups}) {
            $in_tax_groups_selected{$tax_group} = 1;
        }
    }

    my $in_tf_matrix_text = $state->tf_matrix_paste_text;
    my $in_tf_matrix_file = $state->tf_matrix_file;

    #
    # Now the only TF set you can select from individually will be the
    # JASPAR CORE vertebrate set. DJA 10/01/21
    #
    my $core_vert_tf_set  = $self->fetch_matrix_set(
        -collection => 'CORE',
        -tax_group  => 'vertebrates'
    );

    my @core_vert_tf_list;
    my $iter = $core_vert_tf_set->Iterator();
    while (my $matrix = $iter->next()) {
        push @core_vert_tf_list, $matrix;
    }

    # Sort by name in a case insensitive fashion
    @core_vert_tf_list = sort {
        uc($a->name) cmp uc($b->name)
    } @core_vert_tf_list;

    my $vars = {
        abs_html_path           => ABS_HTML_PATH,
        rel_html_path           => REL_HTML_PATH,
        abs_cgi_bin_path        => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path        => REL_CGI_BIN_PATH,
        tf_collections          => JASPAR_COLLECTIONS,
        tax_groups              => JASPAR_TAX_GROUPS,
        version                 => VERSION,
        title                   => 'ORCAtk: Select TFBS Profiles',
        sid                     => $state->sid,
        core_vert_tf_list       => \@core_vert_tf_list,
        in_tf_select_method     => $state->tf_select_method,
        in_core_tf_ids_selected => \%in_core_tf_ids_selected,
        in_tax_groups_selected  => \%in_tax_groups_selected,
        in_tf_collection        => $state->tf_collection,
        in_tf_matrix_file       => $in_tf_matrix_file,
        in_tf_matrix_text       => $in_tf_matrix_text,
        in_tf_min_ic            => $state->tf_min_ic || TF_MIN_IC,
        in_tf_threshold         => $state->tf_threshold,
        in_tf_search_start      => $state->tf_search_start,
        in_tf_search_end        => $state->tf_search_end,
        in_tf_filter_sites      => $state->tf_filter_sites,
        var_template            => "select_tfs.html"
    };

    return $self->process_template('master.html', $vars);
}

sub results
{
    my $self = shift;

    my $state = $self->state;

    my $tfs_selected = $state->tfs_selected();

    my $matrix_set;
    if ($tfs_selected) {
        my $tf_select_method = $state->tf_select_method();

        printf STDERR "\n[%s] results: tf_select_method = $tf_select_method\n",
            scalar localtime(time);

        if ($tf_select_method eq 'paste') {
            $matrix_set = $self->parse_tf_matrix_text(
                $state->tf_matrix_paste_text()
            );

            if (!$matrix_set) {
                return $self->error(
                    "Parsing TFBS profile matrices from pasted text");
            }
        } elsif ($tf_select_method eq 'upload') {

            #printf STDERR "matrix file text:\n%s\n",
            #    $state->tf_matrix_file_text();

            $matrix_set = $self->parse_tf_matrix_text(
                $state->tf_matrix_file_text()
            );

            if (!$matrix_set) {
                return $self->error(
                    "Parsing TFBS profile matrices from uploaded file"
                );
            }
        } else {
            my %matrix_args = (
                -collection     => $state->tf_collection(),
                -tax_group      => $state->tf_tax_groups(),
                -min_ic         => $state->tf_min_ic(),
                -ID             => $state->tf_core_ids(),
                -matrixtype     => 'PWM'
            );

            #printf STDERR "[%s] results matrix_args = \n%s\n\n",
            #    scalar localtime(time),
            #    Data::Dumper::Dumper(%matrix_args);

            $matrix_set = $self->fetch_matrix_set(%matrix_args);

            if (!$matrix_set) {
                return $self->error(
                    sprintf(
                        "Reading TFBS profile matrices from %s with:\n"
                        . "-collection  => %s\n"
                        . "-tax_group   => %s\n"
                        . "-min_ic      => %s\n"
                        . "-ID          => %s\n\n",
                        JASPAR_DB_NAME,
                         $state->tf_collection() || '',
                         $state->tf_tax_groups()
                            ? join(', ', @{$state->tf_tax_groups()})
                            : '',
                         $state->tf_min_ic() || '',
                         $state->tf_core_ids() 
                            ? join(', ', @{$state->tf_core_ids()})
                            : ''
                     )
                );
            }
        }
    }

    $state->tf_matrix_set($matrix_set);

    #if ($state->tf_matrix_set()) {
    #    printf STDERR "\n[%s] results: matrix_set:\n%s\n\n",
    #        scalar localtime(time),
    #        #Data::Dumper::Dumper($state->tf_matrix_set());
    #        Data::Dumper::Dumper($self->state->tf_matrix_set());
    #}

    my $analysis_type = $state->analysis_type;

    if ($analysis_type eq 'pairwise') {
        $self->results_pairwise();
    } elsif ($analysis_type eq 'phastcons') {
        $self->results_phastcons();
    }
}

sub results_pairwise
{
    my $self = shift;

    #printf STDERR "[%s] results_pairwise\n", scalar localtime(time);

    my $q = $self->query;

    my $state = $self->state;

    #printf STDERR "[%s] Results state:\n%s\n",
    #        scalar localtime(time), Data::Dumper::Dumper($state);

    my $pid             = $state->pid;
    my $sid             = $state->sid;
    my $seq_obj1        = $state->seq_obj1;
    my $masked_seq_obj1 = $state->masked_seq_obj1;
    my $seq_obj2        = $state->seq_obj2;
    my $masked_seq_obj2 = $state->masked_seq_obj2;
    my $seq_exons1      = $state->seq_exons1;
    my $seq_cpgs1       = $state->seq_cpgs1;
    my $seq_chr1        = $state->seq_chr1;
    my $seq_start1      = $state->seq_start1;
    my $seq_end1        = $state->seq_end1;
    my $ucsc_db         = $state->ucsc_db;
    my $track_name      = $state->track_name;

    if (!$seq_obj1) {
        return $self->warning(
            "Cannot perform analysis - no sequence entered");
    }

#
# Create a unique file ID for output files. Cannot just use sid as it's possible
# to do multiple analyses with the same sid. DJA 08/10/08
#
    my @ltime   = localtime(time);
    my $datestr = sprintf "%02d%02d%02d%02d%02d%02d",
        $ltime[5] % 100, $ltime[4] + 1, $ltime[3], $ltime[2],
        $ltime[1], $ltime[0];
    my $ufid  = $pid . $datestr;
    my $fbase = "ORCA${ufid}";

    my $alignment;
    my $ca;
    my $phca;
    my $seq1_rel_path;
    my $seq2_rel_path;
    my $masked_seq1_rel_path;
    my $masked_seq2_rel_path;
    my $aln_rel_path;
    my $css_rel_path;
    my $cr_rel_path;
    my $tfbss;
    my $tfbs_rel_path;
    my $plot_rel_path;
    my $arch_rel_path;
    my $ucsc_rel_path;
    my $cons_graph;
    my @tf_sites;

    my $graph_start = $seq_obj1->start;
    my $graph_end   = $seq_obj1->end;

    #printf STDERR "seq1 strand = %d\n", $seq_obj1->strand;
    #printf STDERR "seq2 strand = %d\n", $seq_obj2->strand;

    $ca = ORCA::Analysis::Pairwise->new(
        -chr                   => $seq_chr1,
        -start                 => $graph_start,
        -end                   => $graph_end,
        -base_seq              => $seq_obj1,
        -comparison_seq        => $seq_obj2,
        -masked_base_seq       => $masked_seq_obj1,
        -masked_comparison_seq => $masked_seq_obj2,
        -base_seq_exons        => $seq_exons1,
    );

    if (!$ca) {
        return $self->error("Could not set up conservation analysis");
    }

    $alignment = $ca->compute_alignment();

    if (!$alignment) {
        return $self->error("Could not align sequences");
    }

#printf STDERR "aligned seq1 strand = %d\n", $alignment->get_seq_by_pos(1)->strand;
#printf STDERR "aligned seq2 strand = %d\n", $alignment->get_seq_by_pos(2)->strand;

    my %ccp_params;
    #
    # Explicitly set position_type => 'c' to centre conservation profile
    # (to make consistent with orcatk stand-alone script). DJA 080409
    #
    $ccp_params{'-position_type'} = 'c';
    $ccp_params{'-window_size'}   = $state->ca_window_size
        if $state->ca_window_size;
    $ccp_params{'-window_inc'}    = 1;

    if (!$ca->compute_conservation_profile(%ccp_params)) {
        return $self->error("Could not compute conservation profile");
    }

    my %ccr_params;
    $ccr_params{'-top_pct'} = $state->ca_top_percentile / 100
        if $state->ca_top_percentile;
    $ccr_params{'-min_conservation'} = $state->ca_min_conservation / 100
        if $state->ca_min_conservation;
    $ccr_params{'-filter_exons'}           = $state->ca_filter_exons;
    $ccr_params{'-min_filtered_cr_length'} = $state->ca_min_cr_length
        if $state->ca_min_cr_length;

    if (!$ca->compute_conserved_regions(%ccr_params)) {
        # not an error, there could be no siginificant conservation
        $self->add_warning("No conserved regions found");
    }

    my $matrix_set = $state->tf_matrix_set();
    if ($matrix_set) {
        my $tf_search_start = $state->tf_search_start || 1;
        my $tf_search_end = $state->tf_search_end || $seq_obj1->length();

        $tfbss = $ca->compute_conserved_tfbss(
            -matrix_set     => $matrix_set,
            -tfbs_threshold => $state->tf_threshold
                               ? $state->tf_threshold . '%'
                               : undef,
            -min_tfbs_cr_overlap      => MIN_TFBS_CONSERVATION_OVERLAP,
            -filter_overlapping_sites => $state->tf_filter_sites,
            -start                    => $tf_search_start,
            -end                      => $tf_search_end
        );
    }

    $cons_graph = ORCA::Graphics::Pairwise->new(
        -analysis       => $ca,
        -cpg_islands    => $seq_cpgs1,
        -flip           => $self->state->flip_graph
    );

    if (!$cons_graph || !$cons_graph->{-gd_image}) {
        return $self->error("Unable to set up ORCAtk analysis graph");
    }

 #my ($plot_fh, $plot_file) = tempfile(DIR => ABS_TMP_PATH, SUFFIX => '.png');
 #if (!$plot_fh) {
 #    return $self->error("Unable to create output plot file");
 #}
    my $plot_file = "$fbase.png";
    $plot_rel_path = REL_TMP_PATH . "/$plot_file";
    my $plot_abs_path = ABS_TMP_PATH . "/$plot_file";
    if (!open(PLOT, ">$plot_abs_path")) {
        return $self->error("Could not create PNG plot file $plot_abs_path");
    }
    print PLOT $cons_graph->{-gd_image}->png;
    close PLOT;

    #
    # Write sequences/alignment/conserved regions etc.  out to files
    # (to be archived in results)
    #
    my $seq1_fname = "${fbase}_seq1.fa";
    $seq1_rel_path = REL_TMP_PATH . "/$seq1_fname";
    my $seq1_abs_path = ABS_TMP_PATH . "/$seq1_fname";
    my $seqIO         = Bio::SeqIO->new(
        -file   => ">$seq1_abs_path",
        -format => "fasta"
    );
    $seqIO->write_seq($seq_obj1);
    $seqIO->close;

    if ($masked_seq_obj1) {
        my $masked_seq1_fname = "${fbase}_seq1_masked.fa";
        $masked_seq1_rel_path = REL_TMP_PATH . "/$masked_seq1_fname";
        my $masked_seq1_abs_path = ABS_TMP_PATH . "/$masked_seq1_fname";
        $seqIO = Bio::SeqIO->new(
            -file   => ">$masked_seq1_abs_path",
            -format => "fasta"
        );
        $seqIO->write_seq($masked_seq_obj1);
        $seqIO->close;
    }

    # write seq2
    my $seq2_fname    = "${fbase}_seq2.fa";
    $seq2_rel_path = REL_TMP_PATH . "/$seq2_fname";
    my $seq2_abs_path = ABS_TMP_PATH . "/$seq2_fname";
    $seqIO = Bio::SeqIO->new(
        -file   => ">$seq2_abs_path",
        -format => "fasta"
    );
    $seqIO->write_seq($seq_obj2);
    $seqIO->close;

    if ($masked_seq_obj2) {
        my $masked_seq2_fname    = "${fbase}_seq2_masked.fa";
        my $masked_seq2_rel_path = REL_TMP_PATH . "/$masked_seq2_fname";
        my $masked_seq2_abs_path = ABS_TMP_PATH . "/$masked_seq2_fname";
        $seqIO = Bio::SeqIO->new(
            -file   => ">$masked_seq2_abs_path",
            -format => "fasta"
        );
        $seqIO->write_seq($masked_seq_obj2);
        $seqIO->close;
    }

    # write alignment
    my $aln_fname = "${fbase}_alignment.txt";
    $aln_rel_path = REL_TMP_PATH . "/$aln_fname";
    my $aln_abs_path = ABS_TMP_PATH . "/$aln_fname";
    my $alnIO        = Bio::AlignIO->new(
        -file   => ">$aln_abs_path",
        -format => ALIGNMENT_FORMAT
    );
    $alnIO->write_aln($alignment);
    $alnIO->close;

    # write conserved regions report out to file
    my $cr_file = "${fbase}_conserved_regions.txt";
    $self->write_conserved_regions_report($ca, $cr_file);
    $cr_rel_path = REL_TMP_PATH . "/$cr_file";

    my $css_file = "${fbase}_conserved_subsequences.txt";
    $self->write_conserved_subsequences($ca, $css_file);
    $css_rel_path = REL_TMP_PATH . "/$css_file";

    # write TFBS pairs
    my $tfbs_file = "${fbase}_TFBSs.txt";
    $self->write_tf_site_pairs($tfbss, $tfbs_file);
    $tfbs_rel_path = REL_TMP_PATH . "/$tfbs_file";

    #
    # Archive all result files
    #
    my $arch_file = "${fbase}.tar.gz";
    $arch_rel_path = REL_TMP_PATH . "/$arch_file";
    my $arch_abs_path = ABS_TMP_PATH . "/$arch_file";
    my $ABS_TMP_PATH  = ABS_TMP_PATH;
    `cd $ABS_TMP_PATH; gtar czf $arch_file ${fbase}* --exclude $arch_file`;

    #
    # Create UCSC track file (if chromosomal location is known, i.e. seqs
    # were retrieved from Ensembl rather than uploaded/pasted
    #
    my $ucsc_url;
    if ($seq_chr1) {
        my $track = $ca->ucsc_track;

        if ($track) {
            my $ucsc_file = "${fbase}_UCSC_track.txt";

            $ucsc_rel_path = REL_TMP_PATH . "/$ucsc_file";
            my $ucsc_abs_path = ABS_TMP_PATH . "/$ucsc_file";

            unless (open(UCSC, ">$ucsc_abs_path")) {
                $self->add_warning(
                    "Could not create UCSC track file $ucsc_file"
                );
            } else {
                print UCSC "$track\n";

                close(UCSC);

                $ucsc_url =
                    sprintf "http://genome.ucsc.edu/cgi-bin/hgTracks?org=%s&db=%s&position=%s:%d-%d&hgt.customText=%s/tmp/%s",
                    $state->species1,
                    $state->species_ucsc_dbs->{$state->species1},
                    $seq_chr1,
                    $graph_start,
                    $graph_end,
                    ORCA_URL,
                    $ucsc_file;
            }
        } else {
            $self->add_warning("Could not create UCSC track");
        }
    }

    my $vars = {
        abs_html_path          => ABS_HTML_PATH,
        rel_html_path          => REL_HTML_PATH,
        abs_cgi_bin_path       => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path       => REL_CGI_BIN_PATH,
        version                => VERSION,
        title                  => 'ORCAtk: Results',
        sid                    => $sid,
        in_ca_top_percentile   => $state->ca_top_percentile,
        in_ca_min_conservation => $state->ca_min_conservation,
        in_ca_window_size      => $state->ca_window_size,
        in_ca_min_cr_length    => $state->ca_min_cr_length,
        in_ca_filter_exons     => $state->ca_filter_exons,
        in_tf_threshold        => $state->tf_threshold,
        in_tf_search_start     => $state->tf_search_start,
        in_tf_search_end       => $state->tf_search_end,
        in_tf_filter_sites     => $state->tf_filter_sites,
        in_flip_graph          => $state->flip_graph,
        plot_file              => $plot_rel_path,
        aln_file               => $aln_rel_path,
        cr_file                => $cr_rel_path,
        css_file               => $css_rel_path,
        tfbs_file              => $tfbs_rel_path,
        arch_file              => $arch_rel_path,
        ucsc_file              => $ucsc_rel_path,
        ucsc_url               => $ucsc_url,
        var_template           => "results.html"
    };

    return $self->process_template('master.html', $vars);
}

sub results_phastcons
{
    my $self = shift;

    #printf STDERR "[%s] results_phastcons\n", scalar localtime(time);

    my $q = $self->query;

    my $state = $self->state;

    #printf STDERR "[%s] Results state:\n%s\n",
    #        scalar localtime(time), Data::Dumper::Dumper($state);

    my $pid             = $state->pid;
    my $sid             = $state->sid;
    my $seq_obj1        = $state->seq_obj1;
    my $masked_seq_obj1 = $state->masked_seq_obj1;
    my $seq_obj2        = $state->seq_obj2;
    my $masked_seq_obj2 = $state->masked_seq_obj2;
    my $seq_exons1      = $state->seq_exons1;
    my $seq_cpgs1       = $state->seq_cpgs1;
    my $seq_chr1        = $state->seq_chr1;
    my $seq_start1      = $state->seq_start1;
    my $seq_end1        = $state->seq_end1;
    my $ucsc_db         = $state->ucsc_db;
    my $track_name      = $state->track_name;

    if (!$seq_obj1) {
        return $self->warning(
            "Cannot perform analysis - no sequence entered");
    }

#
# Create a unique file ID for output files. Cannot just use sid as it's possible
# to do multiple analyses with the same sid. DJA 08/10/08
#
    my @ltime   = localtime(time);
    my $datestr = sprintf "%02d%02d%02d%02d%02d%02d",
        $ltime[5] % 100, $ltime[4] + 1, $ltime[3], $ltime[2],
        $ltime[1], $ltime[0];
    my $ufid  = $pid . $datestr;
    my $fbase = "ORCA${ufid}";

    my $alignment;
    my $ca;
    my $phca;
    my $seq1_rel_path;
    my $seq2_rel_path;
    my $masked_seq1_rel_path;
    my $masked_seq2_rel_path;
    my $aln_rel_path;
    my $css_rel_path;
    my $cr_rel_path;
    my $tfbs_rel_path;
    my $plot_rel_path;
    my $arch_rel_path;
    my $ucsc_rel_path;
    my $cons_graph;

    my $graph_start = $seq_obj1->start;
    my $graph_end   = $seq_obj1->end;

    # phastCons analysis
    $phca = ORCA::Analysis::PhastCons->new(
        -seq      => $seq_obj1,
        -exons    => $seq_exons1,
        #-exons    => \@exons,
        -db       => $ucsc_db,
        -track    => $track_name,
        -chr      => $seq_chr1,
        -start    => $seq_start1,
        -end      => $seq_end1
    );

    if (!$phca) {
        return $self->error("Could not set up phastCons analysis");
    }

    $phca->compute_conservation_profile();

    #printf STDERR "phastCons conservation_profile:\n"
    #    . Data::Dumper::Dumper($phca->conservation_profile) . "\n";

    my %ccr_params;
    $ccr_params{'-min_conservation'} = $state->ca_min_conservation / 100
        if $state->ca_min_conservation;
    $ccr_params{'-filter_exons'} = $state->ca_filter_exons;
    $ccr_params{'-min_length'}   = $state->ca_min_cr_length
        if $state->ca_min_cr_length;
    if (!$phca->compute_conserved_regions(%ccr_params)) {
        # not an error, there could be no significant conservation
        $self->add_warning("No conserved regions found");
    }

    #printf STDERR "results_phastcons conserved regions:\n%s\n",
    #    Data::Dumper::Dumper($phca->conserved_regions);

    my $matrix_set = $state->tf_matrix_set();
    my $tfbss;
    if ($matrix_set) {
        #printf STDERR "\n[%s] results_phastcons: matrix_set:\n%s\n\n",
        #    scalar localtime(time),
        #    Data::Dumper::Dumper($state->tf_matrix_set());

        my $tf_search_start = $state->tf_search_start || 1;
        my $tf_search_end = $state->tf_search_end || $seq_obj1->length();
        # XXX TF search start/end specified in relative coords - convert
        # to chromosomal
        #$tf_search_start += $seq_start1 - 1;
        #$tf_search_end += $seq_start1 - 1;

        $tfbss = $phca->compute_conserved_tfbss(
            -matrix_set                 => $matrix_set,
            -min_tfbs_score             => $state->tf_threshold
                                            ? $state->tf_threshold . '%'
                                            : undef,
            -min_tfbs_cr_overlap        => MIN_TFBS_CONSERVATION_OVERLAP,
            -filter_overlapping_tfbss   => $state->tf_filter_sites,
            -start                      => $tf_search_start,
            -end                        => $tf_search_end
        );
    }

    #printf STDERR "\nresults_phastcons: phca:\n%s\n\n",
    #    Data::Dumper::Dumper($phca),

    #printf STDERR "results_phastcons: conserved TFBSs:\n%s\n",
    #    Data::Dumper::Dumper($tfbss);


    #printf STDERR "results_phastcons: CpGs:\n%s\n",
    #    Data::Dumper::Dumper($seq_cpgs1);

    #
    # XXX convert CpGs to chromosomal coords
    #
    #my @graph_cpgs;
    #foreach my $cpg (@$seq_cpgs1) {
    #    my $gcpg = Bio::SeqFeature::Generic->new(
    #        -primary_tag => $cpg->primary_tag,
    #        -source_tag  => $cpg->source_tag,
    #        -strand      => $cpg->strand
    #    );
    #    $gcpg->start($cpg->start + $graph_start - 1);
    #    $gcpg->end($cpg->end + $graph_start - 1);
    #    push @graph_cpgs, $gcpg;
    #}
    #print STDERR "graph_cpgs:\n" . Data::Dumper::Dumper(@graph_cpgs) . "\n";

    $cons_graph = ORCA::Graphics::PhastCons->new(
        -analysis       => $phca,
        -cpg_islands    => $seq_cpgs1,
        #-cpg_islands    => @graph_cpgs  ? \@graph_cpgs  : undef,
        -flip           => $self->state->flip_graph
    );

    if (!$cons_graph || !$cons_graph->{-gd_image}) {
        return $self->error("Unable to set up ORCAtk analysis graph");
    }

    #my ($plot_fh, $plot_file) = tempfile(DIR => ABS_TMP_PATH, SUFFIX => '.png');
    #if (!$plot_fh) {
    #    return $self->error("Unable to create output plot file");
    #}
    my $plot_file = "$fbase.png";
    $plot_rel_path = REL_TMP_PATH . "/$plot_file";
    my $plot_abs_path = ABS_TMP_PATH . "/$plot_file";
    if (!open(PLOT, ">$plot_abs_path")) {
        return $self->error("Could not create PNG plot file $plot_abs_path");
    }
    print PLOT $cons_graph->{-gd_image}->png;
    close PLOT;

    #
    # Write sequences/alignment/conserved regions etc.  out to files
    # (to be archived in results)
    #
    my $seq1_fname = "${fbase}_seq1.fa";
    $seq1_rel_path = REL_TMP_PATH . "/$seq1_fname";
    my $seq1_abs_path = ABS_TMP_PATH . "/$seq1_fname";
    my $seqIO         = Bio::SeqIO->new(
        -file   => ">$seq1_abs_path",
        -format => "fasta"
    );
    $seqIO->write_seq($seq_obj1);
    $seqIO->close;

    if ($masked_seq_obj1) {
        my $masked_seq1_fname = "${fbase}_seq1_masked.fa";
        $masked_seq1_rel_path = REL_TMP_PATH . "/$masked_seq1_fname";
        my $masked_seq1_abs_path = ABS_TMP_PATH . "/$masked_seq1_fname";
        $seqIO = Bio::SeqIO->new(
            -file   => ">$masked_seq1_abs_path",
            -format => "fasta"
        );
        $seqIO->write_seq($masked_seq_obj1);
        $seqIO->close;
    }

    # write conserved regions report out to file
    my $cr_file = "${fbase}_conserved_regions.txt";
    $self->write_conserved_regions_report($phca, $cr_file);
    $cr_rel_path = REL_TMP_PATH . "/$cr_file";

    my $css_file = "${fbase}_conserved_subsequences.txt";
    $self->write_conserved_subsequences($phca, $css_file);
    $css_rel_path = REL_TMP_PATH . "/$css_file";

    # write TFBSs
    my $tfbs_file = "${fbase}_TFBSs.txt";
    $self->write_conserved_tfbss($phca, $tfbs_file);
    $tfbs_rel_path = REL_TMP_PATH . "/$tfbs_file";

    #
    # Archive all result files
    #
    my $arch_file = "${fbase}.tar.gz";
    $arch_rel_path = REL_TMP_PATH . "/$arch_file";
    my $arch_abs_path = ABS_TMP_PATH . "/$arch_file";
    my $ABS_TMP_PATH  = ABS_TMP_PATH;
    `cd $ABS_TMP_PATH; gtar czf $arch_file ${fbase}* --exclude $arch_file`;

    #
    # Create UCSC track file (if chromosomal location is known, i.e. seqs
    # were retrieved from Ensembl rather than uploaded/pasted
    #
    my $ucsc_url;

    my $track = $phca->ucsc_track();

    my $ucsc_file = "${fbase}_UCSC_track.txt";
    $ucsc_rel_path = REL_TMP_PATH . "/$ucsc_file";
    my $ucsc_abs_path = ABS_TMP_PATH . "/$ucsc_file";
    if (!open(UCSC, ">$ucsc_abs_path")) {
        return $self->error(
            "Could not create UCSC track file $ucsc_file");
    }

    print UCSC "$track\n";

    close(UCSC);

    if (   $seq_chr1 !~ /scaffold/i
        && $seq_chr1 !~ /contig/i
        && $seq_chr1 !~ /ultra/i
        && $seq_chr1 !~ /super/i)
    {
        $seq_chr1 = "chr$seq_chr1";
    }

    $ucsc_url =
        sprintf "http://genome.ucsc.edu/cgi-bin/hgTracks?org=%s&db=%s&position=%s:%d-%d&hgt.customText=%s/tmp/%s",
        $state->species1,
        $state->species_ucsc_dbs->{$state->species1},
        $seq_chr1,
        $graph_start,
        $graph_end,
        ORCA_URL,
        $ucsc_file;

    my $vars = {
        abs_html_path          => ABS_HTML_PATH,
        rel_html_path          => REL_HTML_PATH,
        abs_cgi_bin_path       => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path       => REL_CGI_BIN_PATH,
        version                => VERSION,
        title                  => 'ORCAtk: Results',
        sid                    => $sid,
        in_ca_top_percentile   => $state->ca_top_percentile,
        in_ca_min_conservation => $state->ca_min_conservation,
        in_ca_window_size      => $state->ca_window_size,
        in_ca_min_cr_length    => $state->ca_min_cr_length,
        in_ca_filter_exons     => $state->ca_filter_exons,
        in_tf_threshold        => $state->tf_threshold,
        in_tf_search_start     => $state->tf_search_start,
        in_tf_search_end       => $state->tf_search_end,
        in_tf_filter_sites     => $state->tf_filter_sites,
        in_flip_graph          => $state->flip_graph,
        plot_file              => $plot_rel_path,
        aln_file               => $aln_rel_path,
        cr_file                => $cr_rel_path,
        css_file               => $css_rel_path,
        tfbs_file              => $tfbs_rel_path,
        arch_file              => $arch_rel_path,
        ucsc_file              => $ucsc_rel_path,
        ucsc_url               => $ucsc_url,
        var_template           => "results.html"
    };

    return $self->process_template('master.html', $vars);
}

sub seqs_selected_paste_upload
{
    my $self = shift;

   #printf STDERR "[%s] seqs_selected_paste_upload\n", scalar localtime(time);

    my $q = $self->query;

    my $seq_file1 = $q->param("seq_file1");
    my $seq_file2 = $q->param("seq_file2");
    my $gff_file  = $q->param("gff_file");
    my $cdna_file = $q->param("cdna_file");

    my $state = $self->state;

    my $sid      = $state->sid;
    my $species1 = $state->species1;

    my $seq1 = "";
    foreach my $seq_line ($q->param("seq1")) {
        $seq1 .= $seq_line;
    }

    my $seq2 = "";
    foreach my $seq_line ($q->param("seq2")) {
        $seq2 .= $seq_line;
    }

    my $gff_text = "";
    foreach my $gff_line ($q->param("gff_text")) {
        $gff_text .= $gff_line;
    }

    #printf STDERR "\nGFF:\n$gff_text\n\n" if $gff_text;

    my $cdna_text = "";
    foreach my $seq_line ($q->param("cdna")) {
        $cdna_text .= $seq_line;
    }

    my ($seq_select_method1, $seq_obj1, $masked_seq_obj1);
    if ($seq1) {
        $seq_select_method1 = 'paste';

        my @seq_lines = split /\n/, $seq1;
        my $seq_id;
        my $seq_str = "";
        foreach my $seq_line (@seq_lines) {
            chomp($seq_line);
            # if any fasta header, use as sequence ID
            if ($seq_line =~ /^>(.*)/) {
                $seq_id = $1;
                $seq_id =~ s///;
            } else {
                # strip all non-nucleotide, non-masking characters
                $seq_line =~ s/[^ACTGactgNXnx]//g;
                $seq_str .= $seq_line;
            }
        }

        if (!$seq_str) {
            $self->add_error("Sequence 1 not pasted correctly.");
            return;
        }

        my $seq_end1 = length $seq_str;
        $seq_obj1 = Bio::LocatableSeq->new(
            -display_id => $seq_id || "seq1",
            -start      => 1,
            -end        => $seq_end1,
            -strand     => 1,
            -alphabet   => 'dna',
            -seq        => $seq_str
        );
    } elsif ($seq_file1) {
        $seq_select_method1 = 'upload';

        my $seq_fh = $q->upload('seq_file1');
        #$seq_obj1 = $self->upload_seq($seq_fh, $seq_rc1, 'seq1');
        $seq_obj1 = $self->upload_seq($seq_fh, 0, 'seq1');
        if (!$seq_obj1) {
            $self->add_error("Could not upload sequence 1 file");
            return;
        }
    }

    if (!$seq_select_method1 || !$seq_obj1) {
        $self->add_error("Sequence 1 not entered correctly. Please paste the"
                . " sequence directly into the area provided"
                . " OR upload it from a file.");
        return;
    }

    my ($seq_select_method2, $seq_obj2, $masked_seq_obj2);
    if ($seq2) {
        $seq_select_method2 = 'paste';

        my @seq_lines = split /\n/, $seq2;
        my $seq_id;
        my $seq_str = "";
        foreach my $seq_line (@seq_lines) {
            chomp($seq_line);
            # if any fasta header, use as sequence ID
            if ($seq_line =~ /^>(.*)/) {
                $seq_id = $1;
                $seq_id =~ s///;
            } else {
                # strip all non-nucleotide, non-masking characters
                $seq_line =~ s/[^ACTGactgNXnx]//g;
                $seq_str .= $seq_line;
            }
        }

        if (!$seq_str) {
            $self->add_error("Sequence 2 not pasted correctly.");
            return;
        }

        my $seq_end2 = length $seq_str;
        $seq_obj2 = Bio::LocatableSeq->new(
            -display_id => $seq_id || "seq2",
            -start      => 1,
            -end        => $seq_end2,
            -strand     => 1,
            -alphabet   => 'dna',
            -seq        => $seq_str
        );
    } elsif ($seq_file2) {
        $seq_select_method2 = 'upload';

        my $seq_fh = $q->upload('seq_file2');
        $seq_obj2 = $self->upload_seq($seq_fh, 0, 'seq2');
        if (!$seq_obj2) {
            $self->add_error("Could not upload sequence 2 file");
            return;
        }
    }

    if (!$seq_select_method2 || !$seq_obj2) {
        $self->add_error("Sequence 2 not entered correctly. Please paste the"
                . " sequence directly into the area provided"
                . " OR upload from a file.");
        return;
    }

    my @seq_exons1;
    my ($exon_select_method, $cdna_obj, $gff_obj, $gff_exons);
    if ($gff_text) {
        $exon_select_method = "gff_paste";
    } elsif ($gff_file) {
        $exon_select_method = "gff_file";

        my $gff_fh = $q->upload('gff_file');
        $gff_exons = $self->upload_gff($gff_fh);
        if (!$gff_exons) {
            $self->add_error("Could not upload exon GFF file");
            return;
        }
        @seq_exons1 = @$gff_exons;
    } elsif ($cdna_text) {
        $exon_select_method = "cdna_paste";

        my @seq_lines = split /\n/, $cdna_text;
        my $seq_id;
        my $seq_str = "";
        foreach my $seq_line (@seq_lines) {
            chomp($seq_line);
            # if any fasta header, use as sequence ID
            if ($seq_line =~ /^>(.*)/) {
                $seq_id = $1;
                $seq_id =~ s///;
            } else {
                # strip all non-nucleotide, non-masking characters
                $seq_line =~ s/[^ACTGactgNXnx]//g;
                $seq_str .= $seq_line;
            }
        }

        if (!$seq_str) {
            $self->add_error("cDNA sequence not pasted correctly.");
            return;
        }

        my $seq_end1 = length $seq_str;
        $cdna_obj = Bio::LocatableSeq->new(
            -display_id => $seq_id || "seq1",
            -start      => 1,
            -end        => $seq_end1,
            -strand     => 1,
            -alphabet   => 'dna',
            -seq        => $seq_str
        );
    } elsif ($cdna_file) {
        $exon_select_method = 'cdna_upload';

        my $cdna_fh = $q->upload('cdna_file');
        $cdna_obj = $self->upload_seq($cdna_fh, 0, 'cDNA');
        if (!$cdna_obj) {
            $self->add_error("Could not upload cDNA file");
            return;
        }
    }

    if ($cdna_obj) {
        # Run est2genome and create exon objects
        # XXX this should be separated out and put in it's own module under
        # ORCAtk/lib/ORCA/Run
        my $seq1_tmpfile = ABS_TMP_PATH . "/${sid}_seq1.fa";
        my $cdna_tmpfile = ABS_TMP_PATH . "/${sid}_cdna.fa";
        my $est_outfile  = ABS_TMP_PATH . "/${sid}_est2genome.out";
        my $seqIO        = Bio::SeqIO->new(
            -file   => ">$seq1_tmpfile",
            -format => 'fasta'
        );
        if (!$seqIO->write_seq($seq_obj1)) {
            $self->add_error("Writing temporary sequence file $seq1_tmpfile");
            return;
        }
        $seqIO->close;
        $seqIO = Bio::SeqIO->new(
            -file   => ">$cdna_tmpfile",
            -format => 'fasta'
        );
        if (!$seqIO->write_seq($cdna_obj)) {
            $self->add_error("Writing temporary sequence file $cdna_tmpfile");
            return;
        }
        $seqIO->close;

        my $out =
            `/usr/local/bin/est2genome -est $cdna_tmpfile -genome $seq1_tmpfile -outfile $est_outfile -space 500 2>&1`;

        my $rval = $? >> 8;
        if ($rval != 0) {
            $self->add_error("Running est2genome:\n$out\n");
            return;
        }

        my $featiter = Bio::Tools::Est2Genome->new(-file => $est_outfile);
        my $feats = $featiter->parse_next_gene();
        foreach my $f (@$feats) {
            if (lc $f->primary_tag eq 'exon') {
                push @seq_exons1, $f;
            }
        }
    }

    #    if ($state->analysis_type eq 'phastcons') {
    #    my $track_name = SPECIES_TRACK_NAMES->{$species1};
    #    if ($track_name) {
    #        $state->track_name($track_name);
    #        my $ucsc_db = $state->species_ucsc_dbs->{$species1};
    #        $state->ucsc_db($ucsc_db);
    #    } else {
    #        $self->add_error("PhastCons analysis is not available for this"
    #                . " species. Please enter coordinates for an"
    #                . " orthologous species sequence to perform"
    #                . " pairwise analysis");
    #        return;
    #    }
    #    }

    $state->analysis_type("pairwise");
    $state->seq_select_method1($seq_select_method1);
    $state->seq1($seq1);
    $state->seq_file1($seq_file1);
    $state->seq_obj1($seq_obj1);
    $state->masked_seq_obj1($masked_seq_obj1);

    $state->seq_select_method2($seq_select_method2);
    $state->seq2($seq2);
    $state->seq_file2($seq_file2);
    $state->seq_obj2($seq_obj2);
    $state->masked_seq_obj2($masked_seq_obj2);

    $state->exon_select_method($exon_select_method);
    $state->seq_exons1(@seq_exons1 ? \@seq_exons1 : undef);
    $state->cdna_obj($cdna_obj) if $cdna_obj;

    return 1;
}

sub seq1_selected_coords
{
    my $self = shift;

    #printf STDERR "[%s] seq1_selected_coords\n", scalar localtime(time);

    my $q = $self->query;

    my $species1      = $q->param("species1");
    my $seq_chr1      = $q->param("seq_chr1");
    my $seq_start1    = $q->param("seq_start1");
    my $seq_end1      = $q->param("seq_end1");
    my $analysis_type = $q->param("analysis_type");

    my $state = $self->state;

    if (!$species1 || !$seq_chr1 || !$seq_start1 || !$seq_end1) {
        $self->add_error("Sequence 1 coordinates entered incorrectly");
        return;
    }

    if ($analysis_type eq 'phastcons') {
        my $track_name = SPECIES_TRACK_NAMES->{$species1};
        if ($track_name) {
            $state->track_name($track_name);
            my $ucsc_db = $state->species_ucsc_dbs->{$species1};
            $state->ucsc_db($ucsc_db);
        } else {
            $self->add_error("PhastCons analysis is not available for"
                    . " $species1. Please enter coordinates for an"
                    . " orthologous species sequence to perform"
                    . " pairwise analysis");
            return;
        }
    }

    $seq_chr1   =~ s/ //g;
    $seq_start1 =~ s/[, ]//g;
    $seq_end1   =~ s/[, ]//g;

    if ($seq_end1 <= $seq_start1) {
        $self->add_error(
            "Sequence 1 start position must be less than end position");
        return;
    }

    if ($seq_end1 - $seq_start1 + 1 > MAX_SEQ_LEN) {
        $self->add_error(
            "Sequence 1 region is too large, please select a smaller region");
        return;
    }

    my ($seq_select_method1, $seq_obj1, $masked_seq_obj1, $seq_exons1,
        $seq_cpgs1);

    $seq_select_method1 = 'coords';

    ($seq_obj1, $masked_seq_obj1, $seq_exons1, $seq_cpgs1) =
        $self->fetch_ensembl_seqs(
            $species1, $seq_chr1, $seq_start1, $seq_end1,
            #$seq_rc1);
            0
        );

    if (!$seq_select_method1 || !$seq_obj1) {
        $self->add_error("Could not retrieve sequence 1 from Ensembl."
                . " Please make sure the Ensembl DB and the"
                . " chromosomal coordinates are entered correctly");
        return;
    }

    $self->state->seq_select_method1($seq_select_method1);
    $self->state->species1($species1);
    $self->state->seq_chr1($seq_chr1);
    $self->state->seq_start1($seq_start1);
    $self->state->seq_end1($seq_end1);
    $self->state->seq_obj1($seq_obj1);
    $self->state->masked_seq_obj1($masked_seq_obj1);
    $self->state->seq_exons1($seq_exons1);
    $self->state->seq_cpgs1($seq_cpgs1);
    $self->state->analysis_type($analysis_type);

    return 1;
}

sub seq2_selected_coords
{
    my $self = shift;

    #printf STDERR "[%s] seq2_selected_coords\n", scalar localtime(time);

    my $q = $self->query;

    my $species2   = $q->param("species2")   || undef;
    my $seq_chr2   = $q->param("seq_chr2")   || undef;
    my $seq_start2 = $q->param("seq_start2") || undef;
    my $seq_end2   = $q->param("seq_end2")   || undef;

    if (!$species2 || !$seq_chr2 || !$seq_start2 || !$seq_end2) {
        $self->add_error("Sequence 2 coordinates entered incorrectly");
        return;
    }

    $seq_chr2   =~ s/ //g    if $seq_chr2;
    $seq_start2 =~ s/[, ]//g if $seq_start2;
    $seq_end2   =~ s/[, ]//g if $seq_end2;

    if ($seq_end2 <= $seq_start2) {
        $self->add_error(
            "Sequence 2 start position must be less than end position");
        return;
    }

    if ($seq_end2 - $seq_start2 + 1 > MAX_SEQ_LEN) {
        $self->add_error(
            "Sequence 2 region is too large, please select a smaller region");
        return;
    }

    my ($seq_select_method2, $seq_obj2, $masked_seq_obj2);

    $seq_select_method2 = 'coords';

    ($seq_obj2, $masked_seq_obj2) =
        $self->fetch_ensembl_seqs(
            $species2, $seq_chr2, $seq_start2, $seq_end2);

    if (!$seq_obj2) {
        $self->add_error("Could not retrieve sequence 2 from Ensembl."
                . " Please make sure the Ensembl DB and the"
                . " chromosomal coordinates are entered correctly");
        return;
    }

    $self->state->seq_select_method2($seq_select_method2);
    $self->state->species2($species2);
    $self->state->seq_chr2($seq_chr2);
    $self->state->seq_start2($seq_start2);
    $self->state->seq_end2($seq_end2);
    $self->state->seq_obj2($seq_obj2);
    $self->state->masked_seq_obj2($masked_seq_obj2);

    return 1;
}

sub seq1_selected_name
{
    my $self = shift;

    printf STDERR "[%s] seq1_selected_name\n", scalar localtime(time);

    my $q = $self->query;

    my $species1  = $q->param("species1");
    my $gene_name = $q->param("gene_name");
    my $wc        = $q->param("use_wildcards");

    printf STDERR "[%s] species1 = $species1\n", scalar localtime(time);

    $self->state->seq_select_method1('name');
    # XXX ???
    #if ($species1 eq $self->state->species2) {
    #    $self->state->species2($self->state->species1);
    #}
    $self->state->species1($species1);
    $self->state->gene_name1($gene_name);
    $self->state->use_wildcards(defined $wc && $wc eq "on" ? 1 : 0);

    return 1;
}

sub gene_selected
{
    my $self = shift;

    #printf STDERR "[%s] gene_selected\n", scalar localtime(time);

    my $q = $self->query;

    # get species 1 gene selection

    my $ensembl_id1 = $q->param("ensembl_id1");
    my $species2    = $q->param("species2");

    if ($species2 eq $self->state->species1) {
        $self->add_warning(
            "Orthologous sequence species cannot be the same as species 1");
        return;
    }

    $self->state->ensembl_gene_id1($ensembl_id1);
    $self->state->species2($species2);

    return 1;
}

sub gene1_selected
{
    my $self = shift;

    #printf STDERR "[%s] gene1_selected\n", scalar localtime(time);

    my $q = $self->query;

    my $rm            = $q->param('rm');
    my $analysis_type = $q->param("analysis_type");
    my $ensembl_id1   = $q->param("ensembl_id1");

    my $state = $self->state;

    my $species1 = $state->species1;

    if ($analysis_type eq 'pairwise') {
        # pairwise analysis - check orthologous species input
        my $species2 = $q->param("species2");

        if ($species2 eq 'Select orthologous species') {
            $self->add_error("You must select an orthologous species for"
                    . " pairwise analysis");
            return;
        }

        if ($species2 eq $species1) {
            $self->add_error("Orthologous species cannot be the same as"
                    . " analysis species");
            return;
        }

        $self->state->species2($species2);
    } elsif ($analysis_type eq 'phastcons') {
        my $track_name = SPECIES_TRACK_NAMES->{$species1};
        if ($track_name) {
            $state->track_name($track_name);
            my $ucsc_db = $state->species_ucsc_dbs->{$species1};
            $state->ucsc_db($ucsc_db);
        } else {
            $self->add_error("PhastCons analysis is not available for"
                    . " $species1. Please enter coordinates for an"
                    . " orthologous species sequence to perform"
                    . " pairwise analysis");
            return;
        }
    } else {
        $self->add_error("Analysis type is undefined!");
        return;
    }

    $state->ensembl_gene_id1($ensembl_id1);
    $state->analysis_type($analysis_type);

    return 1;
}

sub ortholog_selected
{
    my $self = shift;

    #printf STDERR "[%s] ortholog_selected\n", scalar localtime(time);

    my $q = $self->query;

    my $ensembl_id2 = $q->param("ensembl_id2");

    $self->state->ensembl_gene_id2($ensembl_id2);

    return 1;
}

sub transcripts_selected
{
    my $self = shift;

    #printf STDERR "[%s] transcripts_selected\n", scalar localtime(time);

    my $q = $self->query;

    my $trans_id1     = $q->param("trans_id1");
    my $trans_id2     = $q->param("trans_id2");
    my $upstream_bp   = $q->param("upstream_bp");
    my $downstream_bp = $q->param("downstream_bp");
    my $down_rel_to   = $q->param("down_rel_to");

    my $state = $self->state;

    my $species1 = $state->species1;
    my $species2 = $state->species2;
    my $gene_id1 = $state->ensembl_gene_id1;
    my $gene_id2 = $state->ensembl_gene_id2;

    # get sequence objects based on transcripts
    my ($seq_obj1, $masked_seq_obj1, $seq_exons1, $seq_cpgs1, $seq_chr1) =
        $self->fetch_ensembl_transcript_sequences($species1, $trans_id1,
        $upstream_bp, $downstream_bp, $down_rel_to);
    my ($seq_obj2, $masked_seq_obj2, $seq_exons2, $seq_cpgs2, $seq_chr2) =
        $self->fetch_ensembl_transcript_sequences($species2, $trans_id2,
        $upstream_bp, $downstream_bp, $down_rel_to);

    if ($state->analysis_type eq 'phastcons') {
        my $track_name = SPECIES_TRACK_NAMES->{$species1};
        if ($track_name) {
            $state->track_name($track_name);
            my $ucsc_db = $state->species_ucsc_dbs->{$species1};
            $state->ucsc_db($ucsc_db);
        } else {
            $self->add_error("PhastCons analysis is not available for this"
                    . " species. Please enter coordinates for an"
                    . " orthologous species sequence to perform"
                    . " pairwise analysis");
            return;
        }
    }

    $state->seq_obj1($seq_obj1);
    $state->masked_seq_obj1($masked_seq_obj1);
    $state->seq_chr1($seq_chr1);
    $state->seq_start1($seq_obj1->start);
    $state->seq_end1($seq_obj1->end);
    $state->seq_exons1($seq_exons1);
    $state->seq_cpgs1($seq_cpgs1);
    $state->seq_down_rel_to($down_rel_to);

    $state->seq_obj2($seq_obj2);
    $state->masked_seq_obj2($masked_seq_obj2);
    $state->seq_chr2($seq_chr2);
    $state->seq_start2($seq_obj2->start);
    $state->seq_end2($seq_obj2->end);

    return 1;
}

sub transcript1_selected
{
    my $self = shift;

    #printf STDERR "[%s] transcript1_selected\n", scalar localtime(time);

    my $state = $self->state;

    my $q = $self->query;

    my $trans_id1     = $q->param("trans_id1");
    my $upstream_bp   = $q->param("upstream_bp");
    my $downstream_bp = $q->param("downstream_bp");
    my $down_rel_to   = $q->param("down_rel_to");

    my $species1 = $state->species1;
    my $gene_id1 = $state->ensembl_gene_id1;

    # get sequence objects based on transcripts
    my ($seq_obj1, $masked_seq_obj1, $seq_exons1, $seq_cpgs1, $seq_chr1) =
        $self->fetch_ensembl_transcript_sequences($species1, $trans_id1,
        $upstream_bp, $downstream_bp, $down_rel_to);

    #printf STDERR "[%s] sequence from transcript:\n%s\n",
    #    scalar localtime(time), Data::Dumper::Dumper($seq_obj1);

    #printf STDERR "[%s] exons from transcript:\n%s\n",
    #    scalar localtime(time), Data::Dumper::Dumper($seq_exons1);

    $self->state->seq_obj1($seq_obj1);
    $self->state->masked_seq_obj1($masked_seq_obj1);
    $self->state->seq_exons1($seq_exons1);
    $self->state->seq_cpgs1($seq_cpgs1);
    $self->state->seq_chr1($seq_chr1);
    $self->state->seq_start1($seq_obj1->start);
    $state->seq_down_rel_to($down_rel_to);

    return 1;
}

sub ca_params_selected
{
    my $self = shift;

    #printf STDERR "[%s] ca_params_selected\n", scalar localtime(time);

    my $q = $self->query;

    $self->state->ca_top_percentile($q->param("top_percentile"));
    $self->state->ca_min_conservation($q->param("min_conservation"));
    $self->state->ca_window_size($q->param("window_size"));
    $self->state->ca_filter_exons(
        ($q->param("filter_exons") && $q->param("filter_exons") eq 'on')
        ? 1
        : 0);
    $self->state->flip_graph(
        ($q->param("flip_graph") && $q->param("flip_graph") eq 'on') ? 1 : 0);
    $self->state->ca_min_cr_length($q->param("min_cr_length"));

    #printf STDERR "ca_params_selected: flip_graph = %s\n",
    #                        $q->param('flip_graph');

    #printf STDERR "[%s] state:\n%s\n",
    #                scalar localtime(time),
    #                Data::Dumper::Dumper($self->state);

    return 1;
}

sub tfs_selected
{
    my $self = shift;

    #printf STDERR "[%s] tfs_selected\n", scalar localtime(time);

    my $q = $self->query;

    my $select_method = $self->state->tf_select_method(
        $q->param('tf_select_method')
    );

    my $matrix_set;
    #if ($select_method eq 'min_ic') {
    #    # selected JASPAR CORE vertebrate TFs by min. IC
    #    $self->state->tf_collection('CORE');
    #    $self->state->tf_tax_groups(['vertebrates']);
    #    $self->state->tf_min_ic($q->param('min_ic'));
    #    $self->state->tf_core_ids(undef);
    #    $self->state->tf_matrix_file(undef);
    #    $self->state->tf_matrix_paste_text(undef);
    #    $self->state->tf_matrix_set(undef);
    if ($select_method eq 'tax_group') {
        # selected JASPAR CORE TFs by taxonomic supergroup
        $self->state->tf_collection('CORE');

        my @tax_groups;
        foreach my $tax_group ($q->param('core_tax_groups')) {
            push @tax_groups, $tax_group;
        }

        $self->state->tf_tax_groups(undef);
        if (!@tax_groups) {
            $self->add_error("No JASPAR CORE taxonomic supergroups defined");
            return;
        }
        $self->state->tf_tax_groups(\@tax_groups);

        $self->state->tf_min_ic($q->param('min_ic'));
        $self->state->tf_core_ids(undef);
        $self->state->tf_matrix_file(undef);
        $self->state->tf_matrix_paste_text(undef);
        $self->state->tf_matrix_set(undef);
    } elsif ($select_method eq 'specific') {
        # selected specific JASPAR CORE (vertebrate) TFs
        $self->state->tf_collection('CORE');

        my @tf_ids;
        foreach my $tf_id ($q->param('core_tf_ids')) {
            push @tf_ids, $tf_id;
        }

        $self->state->tf_core_ids(undef);
        if (!@tf_ids) {
            $self->add_error("No JASPAR CORE vertebrate TFs selected");
            return;
        }
        $self->state->tf_core_ids(\@tf_ids);

        $self->state->tf_min_ic(undef);
        $self->state->tf_tax_groups(undef);
        $self->state->tf_matrix_file(undef);
        $self->state->tf_matrix_paste_text(undef);
        $self->state->tf_matrix_set(undef);
    } elsif ($select_method eq 'upload') {
        # uploaded user-defined TFBS profile(s)
        my $tf_file = $q->param('tf_matrix_file');
        my $tf_fh   = $q->upload('tf_matrix_file');

        my $matrix_text = $self->upload_tf_matrix_file($tf_fh);
        if (!$matrix_text) {
            $self->add_error("Uploading TFBS profile matrices");
            return;
        }

        printf STDERR "matrix file text:\n$matrix_text\n";

        $self->state->tf_matrix_file($tf_file);
        $self->state->tf_matrix_file_text($matrix_text);

        $self->state->tf_min_ic(undef);
        $self->state->tf_tax_groups(undef);
        $self->state->tf_core_ids(undef);
        $self->state->tf_matrix_paste_text(undef);
    } elsif ($select_method eq 'paste') {
        my $matrix_text = '';
        foreach my $matrix_line ($q->param('tf_matrix_text')) {
            $matrix_text .= $matrix_line;
        }

        $self->state->tf_matrix_paste_text($matrix_text);

        $self->state->tf_min_ic(undef);
        $self->state->tf_tax_groups(undef);
        $self->state->tf_core_ids(undef);
        $self->state->tf_matrix_file(undef);
        $self->state->tf_matrix_set(undef);
    } else {
        my $jaspar_collections = JASPAR_COLLECTIONS;
        my $collection_matched = 0;
        foreach my $collection (@$jaspar_collections) {
            if ($select_method eq $collection) {
                $self->state->tf_collection($collection);

                $collection_matched = 1;

                last;
            }
        }

        if (!$collection_matched) {
            $self->add_error("Could not determine TF selection method");
            return;
        }

        $self->state->tf_min_ic(undef);
        $self->state->tf_tax_groups(undef);
        $self->state->tf_core_ids(undef);
        $self->state->tf_matrix_file(undef);
        $self->state->tf_matrix_paste_text(undef);
        $self->state->tf_matrix_set(undef);
    }

    $self->state->tf_threshold($q->param('tf_threshold'));

    my $tf_search_start = $q->param('tf_search_start');
    if (defined $tf_search_start && $tf_search_start ne '') {
        if ($tf_search_start < 1) {
            $self->add_warning(
                "TF search start ($tf_search_start) < 1;" . " setting to 1");
            $tf_search_start = 1;
        } elsif (defined $self->state->seq_obj1) {
            my $seqlen = $self->state->seq_obj1->length;
            if ($tf_search_start > $seqlen) {
                $self->add_warning("TF search start ($tf_search_start)"
                        . " > sequence 1 length; setting to 1");
                $tf_search_start = 1;
            }
        }
    }

    my $tf_search_end = $q->param('tf_search_end');
    if (defined $tf_search_end && $tf_search_end ne '') {
        my $seqlen;
        if (defined $self->state->seq_obj1) {
            $seqlen = $self->state->seq_obj1->length;
            if ($tf_search_end > $seqlen) {
                $self->add_warning("TF search end ($tf_search_end) > sequence 1"
                        . " length; setting to $seqlen");
                $tf_search_end = $seqlen;
            }
        }

        if ($tf_search_end < 1) {
            if ($seqlen) {
                $self->add_warning("TF search end ($tf_search_end) < 1;"
                        . " setting to $seqlen");
                $tf_search_end = $seqlen;
            } else {
                $self->add_warning("TF search end < 1");
                $tf_search_end = undef;
            }
        }
    }

    if (   $tf_search_start
        && $tf_search_end
        && $tf_search_start > $tf_search_end)
    {
        $self->add_warning("TF search start ($tf_search_start)"
                . " > TF search end ($tf_search_end); switching");
        ($tf_search_start, $tf_search_end) =
            ($tf_search_end, $tf_search_start);
    }

    $self->state->tf_search_start($tf_search_start);
    $self->state->tf_search_end($tf_search_end);

    my $fos = 0;
    if ($q->param('tf_filter_sites') && $q->param('tf_filter_sites') eq 'on')
    {
        $fos = 1;
    }
    $self->state->tf_filter_sites($fos);

    $self->state->tfs_selected(1);

    return 1;
}

sub write_conserved_regions_report
{
    my ($self, $ao, $cr_file) = @_;

    return if !$ao;

    my $sid = $self->state->sid;
    my $analysis_type = $self->state->analysis_type;

    my $cr_abs_path = ABS_TMP_PATH . "/$cr_file";
    if (!open(FH, ">$cr_abs_path")) {
        $self->add_error(
            "Could not create conserved regions report file $cr_abs_path");
        return;
    }

    if ($analysis_type eq 'pairwise') {
        my $report = $ao->conserved_regions_report;
        if (!$report) {
            $self->add_warning("No conserved regions report");
            close(FH);
            return;
        }

        printf FH "Window size         : %d bp\n",
            $report->param('window_size');
        printf FH "Window increment    : %d bp\n",
            $report->param('window_inc');
        printf FH "Top X%% identities   : %.1f\n", $report->param('top_pct')
            if defined $report->param('top_pct');
        printf FH "Minimum %% identity  : %.1f\n",
            $report->param('min_pct_id')
            if defined $report->param('min_pct_id');
        printf FH "Computed %% identity : %.1f\n",
            $report->param('comp_pct_id')
            if defined $report->param('comp_pct_id');
        printf FH "Effective %% identity: %.1f\n", $report->param('cutoff')
            if defined $report->param('cutoff');
        print FH "\n";
        print FH
            " start1\t   end1\tlength1\t start2\t   end2\tlength2\t aln_st\taln_end\t    \%id\n";

        my $crs = $report->conserved_regions;
        if (!$crs || !@$crs) {
            $self->add_warning(
                "Conserved regions report - no conserved regions");
        }

        foreach my $cr (@$crs) {
            printf FH "%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t%7.1f\n",
                $cr->seq1_start,
                $cr->seq1_end,
                $cr->seq1_end - $cr->seq1_start + 1,
                $cr->seq2_start,
                $cr->seq2_end,
                $cr->seq2_end - $cr->seq2_start + 1,
                $cr->align_start,
                $cr->align_end,
                $cr->score * 100;
        }
    } elsif ($analysis_type eq 'phastcons') {
        # XXX convert to chromosomal coords?
        my $crs = $ao->conserved_regions;

        printf FH "Minimum %% identity  : %.1f\n", $ao->min_conservation * 100
                if defined $ao->min_conservation;
        print FH "\n";
        print FH
            "    start\t      end\t length\t    \%id\n";

        if (!$crs || !@$crs) {
            $self->add_warning(
                "Conserved regions report - no conserved regions");
        }

        foreach my $cr (@$crs) {
            printf FH "%9d\t%9d\t%7d\t%7.1f\n",
                $cr->start,
                $cr->end,
                $cr->end - $cr->start + 1,
                $cr->score * 100;
        }
    }
    close FH;
}

sub write_conserved_subsequences
{
    my ($self, $ao, $css_file) = @_;

    return if !$ao;

    my $sid = $self->state->sid;
    my $analysis_type = $self->state->analysis_type;

    my $css = $ao->compute_conserved_subsequences;

    if (!$css) {
        $self->add_warning("No conserved sub-sequences");
        return;
    }

    my $css_abs_path = ABS_TMP_PATH . "/$css_file";
    if (!open(FH, ">$css_abs_path")) {
        $self->add_error("Could not create conserved sequences report file"
                . " $css_abs_path");
        return;
    }

    foreach my $css (@$css) {
        printf FH ">%s\n%s\n", $css->display_id, $css->seq;
    }

    close FH;
}

sub write_conserved_tfbss
{
    my ($self, $phca, $tfbs_file) = @_;

    return if !$phca;

    my $start = $phca->start;

    my $tfbss = $phca->conserved_tfbss;
    return if !$tfbss;

    my $sid           = $self->state->sid;
    my $tfbs_abs_path = ABS_TMP_PATH . "/$tfbs_file";
    if (!open(TFBS, ">$tfbs_abs_path")) {
        $self->add_error("Could not create TFBS file $tfbs_abs_path");
        return;
    }

    printf TFBS
        "TF\t  Cons.\t  Start\t    End\t Strand\t  Score\t%%Score\tSeq.\n";

    foreach my $site (@$tfbss) {
        my $strand;

        if ($site->strand eq '+' || $site->strand == 1) {
            $strand = '+';
        } elsif ($site->strand eq '-' || $site->strand == -1) {
            $strand = '-';
        }

        printf TFBS "%s\t%7.1f\t%7d\t%7d\t%7s\t%7.3f\t%7.1f\t%s\n",
            $site->pattern->name,
            ($site->get_tag_values('conservation'))[0] * 100,
            $site->start,
            $site->end,
            #$site->start + $start - 1,
            #$site->end + $start - 1,
            $strand,
            $site->score,
            $site->rel_score * 100,
            $site->seq->seq;
    }

    close TFBS;
}

sub write_tf_site_pairs
{
    my ($self, $tfbss, $tfbs_file) = @_;

    return if !$tfbss;

    my $sid           = $self->state->sid;
    my $tfbs_abs_path = ABS_TMP_PATH . "/$tfbs_file";
    if (!open(TFBS, ">$tfbs_abs_path")) {
        $self->add_error("Could not create TFBS file $tfbs_abs_path");
        return;
    }

    printf TFBS
        "TF\t  Cons.\t Start1\t   End1\tStrand1\t Score1\t%%Score1\tSeq1"
        . "\t Start2\t   End2\tStrand2\t Score2\t%%Score2\tSeq2\n";

    foreach my $site_pair (@$tfbss) {
        my $site1 = $site_pair->site1;
        my $site2 = $site_pair->site2;

        my $strand1;
        my $strand2;
        if ($site1->strand eq '+' || $site1->strand == 1) {
            $strand1 = '+';
        } elsif ($site1->strand eq '-' || $site1->strand == -1) {
            $strand1 = '-';
        }
        if ($site2->strand eq '+' || $site2->strand == 1) {
            $strand2 = '+';
        } elsif ($site2->strand eq '-' || $site2->strand == -1) {
            $strand2 = '-';
        }

        printf TFBS "%s\t%7.1f\t%7d\t%7d\t%7s\t%7.3f\t%7.1f\t%s",
            $site1->pattern->name,
            ($site1->get_tag_values('conservation'))[0] * 100,
            $site1->start,
            $site1->end,
            $strand1,
            $site1->score,
            $site1->rel_score * 100,
            $site1->seq->seq;
        printf TFBS "\t%7d\t%7d\t%7s\t%7.3f\t%7.1f\t%s\n",
            $site2->start,
            $site2->end,
            $strand2,
            $site2->score,
            $site2->rel_score * 100,
            $site2->seq->seq;
    }

    close TFBS;
}

sub fetch_matrix_set
{
    my ($self, %args) = @_;

    if (!$args{-matrixtype}) {
        $args{-matrixtype} = 'PWM';
    }

    #my $db =
    #    TFBS::DB::JASPAR4->connect("dbi:mysql:$db_name:" . JASPAR_DB_HOST,
    #    JASPAR_DB_USER, JASPAR_DB_PASS);
    my $db_name = JASPAR_DB_NAME;
    my $db = TFBS::DB::JASPAR5->connect(
            "dbi:mysql:$db_name:" . JASPAR_DB_HOST,
            JASPAR_DB_USER,
            JASPAR_DB_PASS
    );

    if (!$db) {
        $self->add_error("Could not connect to JASPAR DB $db_name");
        return;
    }

    my $matrix_set;
    #if ($min_ic) {
    #    $matrix_set = $db->get_MatrixSet(
    #        -matrixtype => 'PWM',
    #        -collection => $collection,
    #        -min_ic     => $min_ic
    #    );
    #} elsif ($ids) {
    #    $matrix_set = $db->get_MatrixSet(
    #        -matrixtype => 'PWM',
    #        -collection => $collection,
    #        -ID         => $ids
    #    );
    #} else {
    #    $matrix_set = $db->get_MatrixSet(
    #        -matrixtype => 'PWM',
    #        -collection => $collection
    #    );
    #}

    $matrix_set = $db->get_MatrixSet(%args);

    return $matrix_set;
}

sub fetch_ensembl_seqs
{
    my ($self, $species, $chr, $start, $end, $rc) = @_;

    $species = lc $species;

    my $coord_sys_name = 'chromosome';
    if ($chr =~ /genescaffold/i || $chr =~ /gene_scaffold/i) {
        $coord_sys_name = 'genescaffold';
    } elsif ($chr =~ /scaffold/i) {
        $coord_sys_name = 'scaffold';
    } elsif ($chr =~ /ultra/i) {
        $coord_sys_name = 'ultracontig';
    } elsif ($chr =~ /super/i) {
        $coord_sys_name = 'supercontig';
    } elsif ($chr =~ /contig/i) {
        $coord_sys_name = 'contig';
    }

    #printf STDERR
    #    "species = $species; coord_sys = $coord_sys_name; position = $chr:$start-$end\n";

    my $db_name = $self->state->species_ensembl_dbs->{$species};
    if (!$db_name) {
        $self->add_error(
            "Could not determine Ensembl DB name for species $species");
        return;
    }

    #my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    #    -host    => ENSEMBL_DB_HOST,
    #    -port    => ENSEMBL_DB_PORT,
    #    -user    => ENSEMBL_DB_USER,
    #    -pass    => ENSEMBL_DB_PASS,
    #    -dbname  => $db_name,
    #    -species => $species,
    #    -driver  => 'mysql'
    #);
    my $db = $self->get_ensembl_db_adaptor($species, 'core');

    if (!$db) {
        $self->add_error("Could not connect to Ensembl DB $db_name");
        return;
    }

    my $sa = $db->get_SliceAdaptor;
    if (!$sa) {
        $self->add_error("Could not get Ensembl slice adaptor");
        return;
    }

    my $slice = $sa->fetch_by_region($coord_sys_name, $chr, $start, $end);
    if (!$slice) {
        return $self->error("Could not fetch Ensembl $coord_sys_name"
                . " slice $chr:$start-$end");
    }

    my $seq = Bio::LocatableSeq->new(
        -primary_id => $slice->primary_id,
        -display_id => $slice->display_id,
        -alphabet   => 'dna',
        -seq        => $slice->seq,
        -start      => $start,
        -end        => $end,
        -strand     => 1
    );

    my $mseq = Bio::LocatableSeq->new(
        -primary_id => $slice->primary_id,
        -display_id => $slice->display_id,
        -alphabet   => 'dna',
        -seq        => $slice->get_repeatmasked_seq->seq,
        -start      => $start,
        -end        => $end,
        -strand     => 1
    );

    #printf STDERR "seq1 strand = %d\n", $seq->strand;
    if ($rc) {
        #print STDERR "reverse complementing seq1\n";
        $seq = $seq->revcom;
        #$seq->strand(-1);
        if (defined $mseq) {
            $mseq = $mseq->revcom;
            #$mseq->strand(-1);
        }
    }
    #printf STDERR "seq1 strand = %d\n", $seq->strand;

    # get exons
    my $exons = _get_exons_from_slice($slice, $self->state->analysis_type);
    $exons = _combine_exons($exons);
    if ($rc) {
        _revcom_feature_coords($exons, $seq);
    }

    #printf STDERR "[%s] fetch_ensembl_seqs exons\n:%s\n\n",
    #    scalar localtime(time), Data::Dumper::Dumper($exons);

    # get CpG islands
    my $cpgs = _get_cpg_islands_from_slice($slice);
    if ($rc) {
        _revcom_feature_coords($cpgs, $seq);
    }

    return ($seq, $mseq, $exons, $cpgs, $chr);
}

sub fetch_ensembl_gene_by_stable_id
{
    my ($self, $species, $stable_id) = @_;

    $species = lc $species;

    my $db_name = $self->state->species_ensembl_dbs->{$species};
    #printf STDERR "species, db: $species, $db_name\n";
    if (!$db_name) {
        $self->add_error(
            "Could not determine Ensembl DB name for species $species");
        return;
    }

    #my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    #    -host    => ENSEMBL_DB_HOST,
    #    -port    => ENSEMBL_DB_PORT,
    #    -user    => ENSEMBL_DB_USER,
    #    -pass    => ENSEMBL_DB_PASS,
    #    -dbname  => $db_name,
    #    -driver  => 'mysql',
    #    -species => $species
    #);
    my $db = $self->get_ensembl_db_adaptor($species, 'core');

    if (!$db) {
        $self->add_error("Could not connect to Ensembl DB $db_name");
        return;
    }

    my $ga = $db->get_GeneAdaptor;
    if (!$ga) {
        $self->add_error("Could not get Ensembl gene adaptor");
        return;
    }

    my $gene = $ga->fetch_by_stable_id($stable_id);

    return $gene;
}

sub fetch_ensembl_genes_by_name
{
    my ($self, $species, $name, $wc) = @_;

    #printf STDERR "fetch_ensembl_genes_by_name\n";

    $species = lc $species;

    #printf STDERR "species = $species; gene name = $name\n";

    #printf STDERR "state:\n%s\n", Data::Dumper::Dumper($self->state);

    my $db_name = $self->state->species_ensembl_dbs->{$species};
    #printf STDERR "species db = $db_name\n";

    if (!$db_name) {
        $self->add_error(
            "Could not determine Ensembl DB name for species $species");
        return;
    }

    #my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    #    -host    => ENSEMBL_DB_HOST,
    #    -port    => ENSEMBL_DB_PORT,
    #    -user    => ENSEMBL_DB_USER,
    #    -pass    => ENSEMBL_DB_PASS,
    #    -dbname  => $db_name,
    #    -driver  => 'mysql',
    #    -species => $species
    #);
    my $db = $self->get_ensembl_db_adaptor($species, 'core');

    if (!$db) {
        $self->add_error("Could not connect to Ensembl DB $db_name");
        return;
    }

    my $ga = $db->get_GeneAdaptor;
    if (!$ga) {
        $self->add_error("Could not get Ensembl gene adaptor");
        return;
    }

    # XXX use wildcards. This doesn't work as the Ensembl API code does not
    # allow for wildcards. Could replicate Ensembl code but this is kind of
    # ugly.
    #if ($wc) {
    #    $name = "\%$name\%";
    #}
    #print STDERR "Fetching genes with name $name\n";
    my $genes = $ga->fetch_all_by_external_name($name);

    return $genes;
}

#sub fetch_ensembl_orthologs_mart
#{
#    my ($self) = @_;
#
#    my $state    = $self->state;
#    my $species1 = $state->species1;
#    my $species2 = $state->species2;
##printf STDERR "fetch_ensembl_orthologs_mart: species1 = $species1; species2 = $species2\n";
#    my $ensembl_id1 = $state->ensembl_gene_id1;
#
#    $species1 = lc $species1;
#    $species2 = lc $species2;
#    my $table = '';
#
#    my $abbrev1 = $state->species_latin_abbrevs->{$species1};
#    if (!$abbrev1) {
#        $self->add_error("Species 1 unknown fetching orthologs");
#        return;
#    }
#
#    my $abbrev2 = $state->species_latin_abbrevs->{$species2};
#    if (!$abbrev2) {
#        $self->add_error("Species 2 unknown fetching orthologs");
#        return;
#    }
#
#    # mart 49 version
#    #$table = sprintf(ENSEMBL_MART_HOMOL_TBL_TEMPLATE, $abbrev1, $abbrev2);
#
#    # mart 54 version
#    $table = sprintf(ENSEMBL_MART_HOMOL_TBL_TEMPLATE,
#        $abbrev1, ucfirst(substr($abbrev2, 0, 4)));
#
#    my $mart_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#        -host   => ENSEMBL_MART_DB_HOST,
#        -port   => ENSEMBL_MART_DB_PORT,
#        -user   => ENSEMBL_MART_DB_USER,
#        -pass   => ENSEMBL_MART_DB_PASS,
#        -dbname => ENSEMBL_MART_DB,
#        -driver => 'mysql'
#    );
#
#    if (!$mart_db) {
#        $self->add_error("Could not connect to Ensembl Mart DB");
#        return;
#    }
#
#    my $db_name = $self->state->species_ensembl_dbs->{$species2};
#    if (!$db_name) {
#        $self->add_error(
#            "Could not determine Ensembl DB name for species $species2");
#        return;
#    }
#
#    my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#        -host    => ENSEMBL_DB_HOST,
#        -port    => ENSEMBL_DB_PORT,
#        -user    => ENSEMBL_DB_USER,
#        -pass    => ENSEMBL_DB_PASS,
#        -dbname  => $db_name,
#        -species => $species2,
#        -driver  => 'mysql'
#    );
#
#    if (!$db) {
#        $self->add_error("Could not connect to Ensembl DB $db_name");
#        return;
#    }
#
#    my $ga = $db->get_GeneAdaptor;
#    if (!$ga) {
#        $self->add_error("Could not get Ensembl gene adaptor");
#        return;
#    }
#
#    # mart 49 version
#    #my $sql = qq{select distinct homol_stable_id
#    #            from $table
#    #            where homol_stable_id is not NULL
#    #            and gene_stable_id = '$ensembl_id1'};
#
#    # mart 54 version
#    my $sql = qq{select distinct stable_id_4016_r2
#                from $table
#                where stable_id_4016_r2 is not NULL
#                and stable_id_4016 = '$ensembl_id1'};
#
#    my $sth = $mart_db->dbc->prepare($sql);
#    if (!$sth) {
#        $self->add_error("Preparing ortholog select statement");
#        return;
#    }
#
#    if (!$sth->execute) {
#        $self->add_error("Executing ortholog select statement");
#        return;
#    }
#
#    my @genes;
#    while (my @row = $sth->fetchrow_array) {
#        push @genes, $ga->fetch_by_stable_id($row[0]);
#    }
#
#    return @genes ? \@genes : undef;
#}

sub fetch_ensembl_orthologs_compara
{
    my ($self) = @_;

    my $state    = $self->state;
    my $species1 = $state->species1;
    my $species2 = $state->species2;
#printf STDERR "fetch_ensembl_orthologs_compara: species1 = $species1; species2 = $species2\n";
    my $ensembl_id1 = $state->ensembl_gene_id1;

    $species1 = lc $species1;
    $species2 = lc $species2;
    my $table = '';

    my $latin_name2 = $state->species_latin_names->{$species2};
    if (!$latin_name2) {
        $self->add_error("Species 2 unknown latin name fetching orthologs");
        return;
    }

    my $latin_name2a = $latin_name2;

    # e.g. convert 'mus musculus' to 'mus_musculus'
    $latin_name2a =~ s/ /_/;

    #Bio::EnsEMBL::Registry->load_registry_from_db(
    #    -host => ENSEMBL_DB_HOST,
    #    -user => ENSEMBL_DB_USER,
    #    -port => ENSEMBL_DB_PORT
    #);
    $self->load_ensembl_registry();


    #my $member_adaptor =
    #    Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Member');

    #my $db_name = $self->state->species_ensembl_dbs->{$species2};
    #if (!$db_name) {
    #    $self->add_error(
    #        "Could not determine Ensembl DB name for species $species2");
    #    return;
    #}

    #my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    #    -host    => ENSEMBL_DB_HOST,
    #    -port    => ENSEMBL_DB_PORT,
    #    -user    => ENSEMBL_DB_USER,
    #    -pass    => ENSEMBL_DB_PASS,
    #    -dbname  => $db_name,
    #    -species => $species2,
    #    -driver  => 'mysql'
    #);
    #my $db = $self->get_ensembl_db_adaptor($species2, 'core');

    #if (!$db) {
    #    $self->add_error("Could not connect to Ensembl DB $db_name");
    #    return;
    #}

    #my $ga = $db->get_GeneAdaptor;
    #if (!$ga) {
    #    $self->add_error("Could not get Ensembl gene adaptor");
    #    return;
    #}

    #my $member =
    #    $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',
    #    $ensembl_id1);

    #my $homology_adaptor =
    #    Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');

    #printf STDERR "latin_name2 = $latin_name2\n";
    #my $homologies =
    #    $homology_adaptor->fetch_all_by_Member_paired_species($member,
    #    ucfirst $latin_name2, ['ENSEMBL_ORTHOLOGUES']);

    #my @orthologs;
    #foreach my $homology (@{$homologies}) {
    #    foreach my $member_attribute (@{$homology->get_all_Member_Attribute})
    #    {
    #        my ($member, $attribute) = @{$member_attribute};
    #        #printf STDERR "genome_db->name = %s; latin_name2 = %s\n",
    #        #    $member->genome_db->name, $latin_name2;
    #        if (lc $member->genome_db->name eq lc $latin_name2) {
    #            push @orthologs, $ga->fetch_by_stable_id($member->stable_id);
    #        }
    #    }
    #}


    my $db = $self->get_ensembl_db_adaptor($species1, 'core');

    my $ga = $db->get_GeneAdaptor;

    my $ens_gene = $ga->fetch_by_stable_id($ensembl_id1);

    my $gene_homol_species = $ens_gene->get_all_homologous_Genes();

    my @orthologs;
    foreach my $ghs (@$gene_homol_species) {
        my $ens_gene2   = $ghs->[0];
        my $ghs_homol   = $ghs->[1];
        my $ghs_spec2   = lc $ghs->[2];

        #printf STDERR "ghs_spec2: $ghs_spec2\n";
        if (   $ghs_spec2 eq lc $species2
            || $ghs_spec2 eq lc $latin_name2
            || $ghs_spec2 eq lc $latin_name2a
        ) {
            push @orthologs, $ens_gene2;
        }
    }


    return @orthologs ? \@orthologs : undef;
}

sub fetch_ensembl_transcripts
{
    my ($self, $species, $ens_gene_id) = @_;

    $species = lc $species;

    my $db_name = $self->state->species_ensembl_dbs->{$species};
    if (!$db_name) {
        $self->add_error(
            "Could not determine Ensembl DB name for species $species");
        return;
    }

    #my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    #    -host    => ENSEMBL_DB_HOST,
    #    -port    => ENSEMBL_DB_PORT,
    #    -user    => ENSEMBL_DB_USER,
    #    -pass    => ENSEMBL_DB_PASS,
    #    -dbname  => $db_name,
    #    -species => $species,
    #    -driver  => 'mysql'
    #);
    my $db = $self->get_ensembl_db_adaptor($species, 'core');

    if (!$db) {
        $self->add_error("Could not connect to Ensembl DB $db_name");
        return;
    }

    my $ga = $db->get_GeneAdaptor;
    if (!$ga) {
        $self->add_error("Could not get Ensembl gene adaptor");
        return;
    }

    my $gene = $ga->fetch_by_stable_id($ens_gene_id);
    if (!$gene) {
        $self->add_error(
            "Could not get fetch Ensembl $species gene $ens_gene_id");
        return;
    }

    my $transcripts = $gene->get_all_Transcripts;
    if (!$transcripts) {
        $self->add_error("Could not get fetch Ensembl transcripts for $species"
                . " gene $ens_gene_id");
        return;
    }

    return $transcripts;
}

sub fetch_ensembl_exons_from_transcript
{
    my ($self, $trans, $species) = @_;

    $species = lc $species;

    my $db_name = $self->state->species_ensembl_dbs->{$species};
    if (!$db_name) {
        $self->add_error(
            "Could not determine Ensembl DB name for species $species");
        return;
    }

    #my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    #    -host    => ENSEMBL_DB_HOST,
    #    -port    => ENSEMBL_DB_PORT,
    #    -user    => ENSEMBL_DB_USER,
    #    -pass    => ENSEMBL_DB_PASS,
    #    -dbname  => $db_name,
    #    -species => $species,
    #    -driver  => 'mysql'
    #);
    my $db = $self->get_ensembl_db_adaptor($species, 'core');

    if (!$db) {
        $self->add_error("Could not connect to Ensembl DB $db_name");
        return;
    }

    my $ea = $db->get_ExonAdaptor;
    if (!$ea) {
        $self->add_error("Could not get Ensembl exon adaptor");
        return;
    }

    my $exons = $ea->fetch_all_by_Transcript($trans);

    return $exons;
}

sub fetch_ensembl_transcript_sequences
{
    #my ($self, $species, $gene_id, $trans_id, $upstream_bp,
    my ($self, $species, $trans_id, $upstream_bp, $downstream_bp,
        $down_rel_to)
        = @_;

    $species = lc $species;

    my $db_name = $self->state->species_ensembl_dbs->{$species};
    if (!$db_name) {
        $self->add_error(
            "Could not determine Ensembl DB name for species $species");
        return;
    }

    #my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    #    -host    => ENSEMBL_DB_HOST,
    #    -port    => ENSEMBL_DB_PORT,
    #    -user    => ENSEMBL_DB_USER,
    #    -pass    => ENSEMBL_DB_PASS,
    #    -dbname  => $db_name,
    #    -species => $species,
    #    -driver  => 'mysql'
    #);
    my $db = $self->get_ensembl_db_adaptor($species, 'core');

    if (!$db) {
        $self->add_error("Could not connect to Ensembl DB $db_name");
        return;
    }

    #my $ga = $db->get_GeneAdaptor;
    #if (!$ga) {
    #    $self->add_error("Could not get Ensembl gene adaptor");
    #    return;
    #}

    #my $gene = $ga->fetch_by_stable_id($gene_id);

    my $ta = $db->get_TranscriptAdaptor;
    if (!$ta) {
        $self->add_error("Could not get Ensembl transcript adaptor");
        return;
    }

    #printf STDERR "fetching transcript $trans_id\n";
    my $trans = $ta->fetch_by_stable_id($trans_id);
    if (!$trans) {
        $self->add_error("Fetching Ensembl $species transcript $trans_id");
        return;
    }

    my $coord_sys_name = $trans->slice->coord_system->name;
    #printf STDERR "coord_sys_name = $coord_sys_name\n";
    my $chr    = $trans->slice->seq_region_name;
    my $strand = $trans->strand;
    my $start  = $trans->start;
    my $end    = $trans->end;

    if ($strand == 1) {
        my $tss = $start;

        if (defined $upstream_bp) {
            $start = $tss - $upstream_bp;
        }

        if (defined $downstream_bp) {
            if ($down_rel_to eq 'end') {
                $end += $downstream_bp;
            } elsif ($down_rel_to eq 'tss') {
                $end = $tss + $downstream_bp - 1;
            }
        }
    } elsif ($strand == -1) {
        my $tss = $end;

        if (defined $upstream_bp) {
            $end = $tss + $upstream_bp;
        }

        if (defined $downstream_bp) {
            if ($down_rel_to eq 'end') {
                $start -= $downstream_bp;
            } elsif ($down_rel_to eq 'tss') {
                $start = $tss - $downstream_bp + 1;
            }
        }
    }

    return $self->fetch_ensembl_seqs($species, $chr, $start, $end, 0);
}

sub create_transcript_panel
{
    my ($self, $transcripts, $species) = @_;

    #printf STDERR "transcripts:\n%s\n", Data::Dumper::Dumper($transcripts);

    my ($start, $end) = (1e21, 0);
    foreach my $trans (@$transcripts) {
        $start = $trans->start if $trans->start < $start;
        $end   = $trans->end   if $trans->end > $end,;
    }
    my $size = $end - $start;

    ($start, $end) = ($start - int($size / 8), $end + int($size / 8));

    my $panel = Bio::Graphics::Panel->new(
        -width    => TRANSCRIPT_PANEL_WIDTH,
        -start    => $start,
        -stop     => $end,
        -flip     => ($transcripts->[0]->strand == -1 ? 1 : 0),
        -pad_left => 50
    );

    # ruler
    $panel->add_track(
        "arrow" => [
            Bio::SeqFeature::Generic->new(
                -start => $start,
                -end   => $end
            )
        ],
        -fgcolor => "#7f7f7f",
        -tick    => 2
    );

    # draw transcript tracks
    foreach my $ens_trans (@$transcripts) {
        #if (!defined $trans->{'_trans_exon_array'}) {
        #    carp "transcript exon array not defined\n";
        #}

        #if (!defined $trans->adaptor()) {
        #    carp "transcript adapaptor not defined\n";
        #}

        # Convert to Bioperl exons (SeqFeatures)
        my $transcript = _ensembl_to_bio_transcript($ens_trans);
        #printf STDERR "transcript:\n%s\n", Data::Dumper::Dumper($transcript);

        $panel->add_track(
            transcript2  => $transcript,
            -bgcolor     => TRANSCRIPT_BGCOLOUR,
            -fgcolor     => TRANSCRIPT_FGCOLOUR,
            -connector   => "hat",
            -height      => 15,
            -label       => 0,
            -description => 0
        );
    }

    # write transcript IDs manually
    my %written;    #each one only once
    my @boxes = $panel->boxes;
    my $black = $panel->gd->colorAllocate(0, 0, 0);
    foreach my $box (@boxes) {
        my $acc = $box->[0]->display_name;
        if (defined($acc) && !$written{$acc}) {
            $panel->gd->string(GD::gdSmallFont(), 10, $box->[2], $acc,
                $black);
            $written{$acc} = 1;    # write each only once
        }
    }

    return $panel;
}

sub fetch_species_names_and_dbs
{
    my $self = shift;

    my $dsn = sprintf "dbi:mysql:host=%s;database=%s;",
        ENSEMBL_DB_HOST, ENSEMBL_DATABASES_DB_NAME;

    my $dbh = DBI->connect($dsn, ENSEMBL_DB_USER, ENSEMBL_DB_PASS);
    if (!$dbh) {
        $self->add_error("Could not connect to species databases DB");
        return;
    }

    my $ensembl_lib = ENSEMBL_LIB;
    my $ensembl_version;
    if ($ensembl_lib =~ /ensembl-(\w+)/) {
        $ensembl_version = $1;
    }
    #printf STDERR "ensembl_version = $ensembl_version\n";

    my $sql;
    if ($ensembl_version eq 'current') {
        $sql = qq{select sn.common_name, ed.organism, ed.db_name
                from ens_dbnames ed, species_names sn
                where ed.organism = sn.latin_name
                and ed.current = "yes"};
    } else {
        $sql = qq{select sn.common_name, ed.organism, ed.db_name
                from ens_dbnames ed, species_names sn
                where ed.organism = sn.latin_name
                and ed.db_name like '%_core_$ensembl_version%'};
    }

    my $sth = $dbh->prepare($sql);
    if (!$sth) {
        $self->add_error(
            "Could not prepare SQL for species databases DB\n" . "$sql\n");
        return;
    }

    if (!$sth->execute()) {
        $self->add_error(
            "Could not execute SQL for species databases DB\n" . "$sql\n");
        return;
    }

    my @common_names;
    my %latin_names;
    my %latin_to_common_names;
    my %latin_abbrevs;
    my %ensembl_dbs;
    my %ucsc_dbs;
    my $species_ucsc_db_names = SPECIES_UCSC_DB_NAMES;
    while (my @row = $sth->fetchrow_array()) {
        my $common_name = lc $row[0];
        my $latin_name  = lc $row[1];
        my $ens_db_name = $row[2];

        push @common_names, $common_name;
        $latin_names{$common_name}          = $latin_name;
        $latin_to_common_names{$latin_name} = $common_name;
        if ($latin_name =~ /^(\w)\w+\s+(\w+)/) {
            $latin_abbrevs{$common_name} = "$1$2";
        }
        $ensembl_dbs{$common_name} = $ens_db_name;

        my $ucsc_db_name = $species_ucsc_db_names->{$common_name};
        $ucsc_dbs{$common_name} = $ucsc_db_name if $ucsc_db_name;
    }

    @common_names = sort @common_names;

    return (
        @common_names          ? \@common_names          : undef,
        %latin_names           ? \%latin_names           : undef,
        %latin_to_common_names ? \%latin_to_common_names : undef,
        %latin_abbrevs         ? \%latin_abbrevs         : undef,
        %ensembl_dbs           ? \%ensembl_dbs           : undef,
        %ucsc_dbs              ? \%ucsc_dbs              : undef
    );
}

sub _ensembl_to_bio_transcript
{
    my ($ens_trans, $tag) = @_;

    #printf STDERR "_ensembl_to_bio_transcript ensembl exons:\n%s\n",
    #            Data::Dumper::Dumper($ens_trans);

    return if !$ens_trans;

    my $trans = Bio::SeqFeature::Gene::Transcript->new(
        -source_tag => 'ensembl',
        #-source_tag    => $ens_trans->display_id,
        -display_name => $ens_trans->external_name,
        -start        => $ens_trans->start,
        -end          => $ens_trans->end,
        -strand       => $ens_trans->strand
    );

    my $ens_exons = $ens_trans->get_all_Exons;
    foreach my $ens_exon (@$ens_exons) {
        $trans->add_exon(
            Bio::SeqFeature::Gene::Exon->new(
                -primary_tag  => 'exon',
                -source_tag   => 'ensembl',
                -display_name => $ens_exon->display_id,
                -start        => $ens_exon->start,
                -end          => $ens_exon->end,
                -strand       => $ens_exon->strand
            )
        );
    }

    return $trans;
}

sub upload_seq
{
    my ($self, $fh, $rc, $seq_name) = @_;

    my $seq_id;
    my $seq = "";
    while (my $line = <$fh>) {
        chomp($line);
        # if any fasta header, use as sequence ID
        if ($line =~ /^>(.*)/) {
            $seq_id = $1;
            $seq_id =~ s///;
        } else {
            # strip all non-nucleotide, non-masking characters
            $seq =~ s/[^ACTGactgNXnx]//g;
            $seq .= $line;
        }
    }

    return if !$seq;

    if ($rc) {
        $seq =~ tr/actgACTG/tgacTGAC/;
        $seq = reverse $seq;
    }

    my $seq_obj = Bio::LocatableSeq->new(
        -display_id => $seq_id || $seq_name,
        -start      => 1,
        -end        => length $seq,
        -strand     => 1,
        -alphabet   => 'dna',
        -seq        => $seq
    );

    return $seq_obj;
}

sub upload_gff
{
    my ($self, $fh) = @_;

    my $gffio = Bio::Tools::GFF->new(-fh => $fh, -gff_version => 2);

    my @exons;
    if ($gffio) {
        while (my $feature = $gffio->next_feature) {
            push @exons, $feature;
        }
    } else {
        $self->add_error("Could not open exon GFF file\n");
    }

    return @exons ? \@exons : undef;
}

sub upload_tf_matrix_file
{
    my ($self, $fh) = @_;

    return if !$fh;

    my $text;
    while (my $line = <$fh>) {
        $text .= $line;
    }
    close($fh);

    return $text;
}

sub parse_tf_matrix_text
{
    my ($self, $text) = @_;

    return if !$text;

    my $matrix_set    = TFBS::MatrixSet->new();
    my $name          = '';
    my $matrix_string = '';
    my $line_count    = 0;
    my $matrix_count  = 0;
    my @lines         = split /\n/, $text;
    foreach my $line (@lines) {
        chomp $line;
        next if !$line;

        if ($line =~ /^>\s*(\S+)/) {
            $name = $1;
        } else {
            if ($line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/) {
                # line of the form: A [ # # # ... # ]
                $matrix_string .= "$1\n";
            } elsif ($line =~ /^\s*\d+/) {
                # line of the form: # # # ... #
                $matrix_string .= "$line\n";
            } else {
                next;
            }
            $line_count++;

            if ($line_count == 4) {
                my $pfm = TFBS::Matrix::PFM->new(
                    -matrixstring => $matrix_string,
                    -name         => $name,
                    -ID           => $name
                );

                $matrix_set->add_Matrix($pfm->to_PWM);

                $matrix_count++;

                $line_count    = 0;
                $name          = '';
                $matrix_string = '';
            }
        }
    }

    return $matrix_count ? $matrix_set : undef;
}

sub ensembl_registry_loaded
{
    my $self = shift;

    if (@_) {
        $self->{-ensembl_registry_loaded} = shift;;
    }

    return $self->{-ensembl_registry_loaded};
}

sub load_ensembl_registry
{
    my $self = shift;

    unless ($self->ensembl_registry_loaded()) {
        Bio::EnsEMBL::Registry->load_registry_from_db(
            -host    => ENSEMBL_DB_HOST,
            -user    => ENSEMBL_DB_USER,
            -pass    => ENSEMBL_DB_PASS
        );

        $self->ensembl_registry_loaded(1);
    }

    return $self->ensembl_registry_loaded();
}

sub get_ensembl_db_adaptor
{
    my ($self, $species, $type) = @_;

    $self->load_ensembl_registry();

    if (lc $type eq 'core') {
        my $latin_name = $self->state->species_latin_names->{$species};

        $species = $latin_name if $latin_name;
    }

    #printf STDERR "get_ensembl_db_adaptor: species = $species\n";

    my $db = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $type);

    return $db;
}

sub _get_exons_from_slice
{
    my ($slice, $analysis_type) = @_;

    my $genes_in_slice  = $slice->get_all_Genes;
    my $slice_chr_start = $slice->start;
    my $slice_chr_end   = $slice->end;
    my $slice_start     = 1;
    my $slice_end       = $slice_chr_end - $slice_chr_start + 1;

    my @exon_list;
    foreach my $gene (@$genes_in_slice) {
        my $transcripts = $gene->get_all_Transcripts;
        foreach my $transcript (@$transcripts) {
            my $exons = $transcript->get_all_Exons;
            foreach my $exon (@$exons) {
                #
                # Check transcript exons fall within slice and convert from
                # Ensembl Exons to Bioperl Exons.
                #
                if ($exon->end >= $slice_start && $exon->start <= $slice_end) {
                    #
                    # The analysis objects expecte exons to be passed with
                    # chromosomal start/ends not coords relative to transcript
                    # start.
                    #
                    my $new_exon = Bio::SeqFeature::Gene::Exon->new(
                        -primary_tag => 'exon',
                        -source_tag  => 'ensembl',
                        -start       => $slice->start + $exon->start - 1,
                        -end         => $slice->start + $exon->end - 1,
                        -strand      => $exon->strand
                    );

                    push(@exon_list, $new_exon);
                }
            }
        }
    }

    return \@exon_list;
}

sub _get_cpg_islands_from_slice
{
    my ($slice) = @_;

    my $cpg_islands = $slice->get_all_SimpleFeatures('CpG');

    my $slice_start = $slice->start;
    my $slice_end   = $slice->end;
    my @cpgs;
    foreach my $cpg (@$cpg_islands) {
        my $new_cpg = Bio::SeqFeature::Generic->new(
            -primary_tag => 'CpG',
            -source_tag  => 'ensembl',
            -start       => $cpg->start,
            -end         => $cpg->end,
            -strand      => $cpg->strand
        );

        #$new_cpg->display_name($gene_name);
        push(@cpgs, $new_cpg);
    }

    return @cpgs ? \@cpgs : undef;
}

sub _combine_exons
{
    my ($exons) = @_;

    return undef if !$exons || !@$exons;

    my $combined_exons = [];
    foreach my $exon (@$exons) {
        push @$combined_exons, $exon;
    }

    $combined_exons = _sort_and_combine_exons($combined_exons);

    return $combined_exons;
}

sub _sort_and_combine_exons
{
    my ($exon_list) = @_;

    @$exon_list = sort {$a->start <=> $b->end} @$exon_list;
    my $num_exons = scalar @$exon_list;
    for (my $i = 0; $i < $num_exons; $i++) {
        my $exon1 = $exon_list->[$i] if exists($exon_list->[$i]);
        if ($exon1) {
            for (my $j = $i + 1; $j < $num_exons; $j++) {
                my $exon2 = $exon_list->[$j] if exists($exon_list->[$j]);
                if ($exon2) {
                    if (_combine_exon($exon1, $exon2)) {
                       #printf STDERR "Combining exons " . $exon1->start . "-"
                       #    . $exon1->end . " and " . $exon2->start
                       #    . "-" . $exon2->end . "\n";
                        if ($exon2->start < $exon1->start) {
                            $exon1->start($exon2->start);
                        }
                        if ($exon2->end > $exon1->end) {
                            $exon1->end($exon2->end);
                        }
                        #
                        # If the second exon is bigger, the new exon gets
                        # it's name (might makes right :).
                        #
                        if (   $exon2->display_name
                            && $exon2->length > $exon1->length)
                        {
                            $exon1->display_name($exon2->display_name);
                        }
                        delete $exon_list->[$j];
                    } else {
                        last;
                    }
                }
            }
        }
    }

    my @unique_exons;
    foreach my $exon (@$exon_list) {
        if (defined $exon) {
            push @unique_exons, $exon;
        }
    }

    undef @$exon_list;

    return \@unique_exons;
}

sub _combine_exon
{
    my ($exon1, $exon2) = @_;

    my $combine = 1;
    $combine = 0
        if $exon1->start > $exon2->end + 1
            || $exon1->end < $exon2->start - 1;

    return $combine;
}

sub _revcom_feature_coords
{
    my ($features, $seq) = @_;

    return if !$features || !@$features;

    my $seq_len = $seq->length;
    foreach my $feat (@$features) {
        my $start = $seq_len - $feat->end + 1;
        my $end   = $seq_len - $feat->start + 1;
        $feat->start($start);
        $feat->end($end);
    }

    #
    # Re-sort them (purely for aesthetic reasons)
    #
    @$features = sort {$a->start <=> $b->end} @$features;
}

sub state
{
    my $self = shift;

    if (@_) {
        $self->{-state} = shift;
    }

    return $self->{-state};
}

#
# High-level error routine. Call low level add_error routine with current
# error and output all current errors to HTML error template.
#
sub error
{
    my ($self, $error) = @_;

    $self->add_error($error) if $error;

    my $errors = $self->errors;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Error',
        # state may not be defined if we are reporting an error from
        # loading state object from file DJA 10/01/22
        sid              => $self->state ? $self->state->sid : undef,
        errors           => $errors,
        var_template     => "error.html"
    };

    my $output = $self->process_template('master.html', $vars);
    return $output;
}

#
# High-level warning routine. Call low level add_warning routine with current
# warning and output all current warnings to HTML warning template.
#
sub warning
{
    my ($self, $warning) = @_;

    $self->add_warning($warning) if $warning;

    my $warnings = $self->warnings;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Warning',
        sid              => $self->state->sid,
        warnings         => $warnings,
        var_template     => "warning.html"
    };

    my $output = $self->process_template('master.html', $vars);

    return $output;
}

sub contact
{
    my ($self) = @_;

    my $vars = {
        abs_html_path    => ABS_HTML_PATH,
        rel_html_path    => REL_HTML_PATH,
        abs_cgi_bin_path => ABS_CGI_BIN_PATH,
        rel_cgi_bin_path => REL_CGI_BIN_PATH,
        version          => VERSION,
        title            => 'ORCAtk: Contact',
        sid              => $self->state->sid,
        contact_name     => CONTACT_NAME,
        contact_email    => CONTACT_EMAIL,
        var_template     => "contact.html"
    };

    my $output = $self->process_template('master.html', $vars);

    return $output;
}

sub process_template
{
    my ($self, $template_name, $vars) = @_;

    my $config = {
        ABSOLUTE     => 1,
        INCLUDE_PATH => ABS_TEMPLATE_PATH . "/", # or list ref
        INTERPOLATE  => 1,                       # expand "$var" in plain text
        POST_CHOMP   => 1,                       # cleanup whitespace
             #PRE_PROCESS    => 'header',         # prefix each template
        EVAL_PERL => 1    # evaluate Perl code blocks
    };

    my $string   = '';
    my $template = Template->new($config);
    my $input    = ABS_TEMPLATE_PATH . "/$template_name";
    $template->process($input, $vars, \$string) || die $template->error();

    return $string;
}

sub errors
{
    my ($self, $errors) = @_;

    if ($errors) {
        $self->{-errors} = $errors;
    }

    return $self->{-errors};
}

sub warnings
{
    my ($self, $warnings) = @_;

    if ($warnings) {
        $self->{-warnings} = $warnings;
    }

    return $self->{-warnings};
}

sub clean_tempfiles
{
    my $self = shift;

    my @tempfiles = glob(ABS_TMP_PATH . "/*");
    foreach my $file (@tempfiles) {
        unlink $file if -M $file > CLEAN_TEMPFILES_OLDER_THAN;
    }
}

#
# Add latest error to internal error list and output to stderr (log file).
#
sub add_error
{
    my ($self, $error) = @_;

    return if !$error;

    printf STDERR "[%s] ERROR: $error\n", scalar localtime(time);

    my $errors = $self->errors;
    push @$errors, $error;

    $self->errors($errors);
}

#
# Add latest warning to internal warning list and output to stderr (log file).
#
sub add_warning
{
    my ($self, $warning) = @_;

    return if !$warning;

    printf STDERR "[%s] Warning: $warning\n", scalar localtime(time);

    my $warnings = $self->warnings;
    push @$warnings, "$warning";

    $self->warnings($warnings);
}

sub initialize_state
{
    my ($self, $state) = @_;

    return if !$state;

    #
    # Set default/initial state values
    #
    $state->analysis_type(undef);
    $state->species1(SPECIES1);
    $state->species2(SPECIES2);
    $state->gene_name1(undef);
    $state->gene_name2(undef);
    $state->ensembl_gene_id1(undef);
    $state->ensembl_gene_id2(undef);
    $state->ensembl_transcript_id1(undef);
    $state->ensembl_transcript_id2(undef);

    # Sequence information
    $state->seq_select_method1(undef);
    $state->seq1(undef);
    $state->seq_db1(undef);    # not actually used
    $state->seq_chr1(undef);
    $state->seq_start1(undef);
    $state->seq_end1(undef);
    $state->seq_exons1(undef);
    $state->seq_cpgs1(undef);
    $state->seq_rc1(undef);

    $state->seq_select_method2(undef);
    $state->seq2(undef);
    $state->seq_db2(undef);    # not actually used
    $state->seq_chr2(undef);
    $state->seq_start2(undef);
    $state->seq_end2(undef);
    $state->seq_exons2(undef);

    $state->cdna_select_method(undef);
    $state->cdna(undef);
    $state->cdna_file(undef);
    $state->seq_upstream_bp(UPSTREAM_BP);
    $state->seq_downstream_bp(DOWNSTREAM_BP);
    $state->seq_down_rel_to(DOWN_REL_TO);

    # Conservation parameters
    $state->ca_top_percentile(CA_TOP_PERCENTILE);
    $state->ca_min_conservation(undef);
    $state->ca_filter_exons(CA_FILTER_EXONS);
    $state->ca_window_size(CA_WINDOW_SIZE);
    $state->ca_min_cr_length(undef);

    # TFBS parameters
    $state->tfs_selected(0);
    $state->tf_select_method(TF_SELECT_METHOD);
    $state->tf_collection(undef);
    $state->tf_min_ic(TF_MIN_IC);
    $state->tf_threshold(TF_THRESHOLD);
    $state->tf_filter_sites(TF_FILTER_SITES);
    $state->tf_core_ids(undef);
    $state->tf_tax_groups(undef);
    $state->tf_phylo_ids(undef);
    $state->tf_fam_ids(undef);
    $state->tf_matrix_file(undef);
    $state->tf_matrix_file_text(undef);
    $state->tf_matrix_paste_text(undef);
    $state->tf_search_start(undef);
    $state->tf_search_end(undef);

    # Computed/compound items
    $state->seq_obj1(undef);
    $state->masked_seq_obj1(undef);
    $state->seq_obj2(undef);
    $state->masked_seq_obj2(undef);
    $state->cdna_obj(undef);

    $state->flip_graph(0);

    my ($common_names,
        $latin_names, 
        $latin_to_common_names,
        $latin_abbrevs,
        $ensembl_dbs,
        $ucsc_dbs) = $self->fetch_species_names_and_dbs();

    #printf STDERR "latin_names:\n" . Data::Dumper::Dumper($latin_names) . "\n";
    #printf STDERR "latin_to_common_names:\n"
    #            . Data::Dumper::Dumper($latin_to_common_names) . "\n";
    #printf STDERR "latin_abbrevs:\n"
    #            . Data::Dumper::Dumper($latin_abbrevs) . "\n";
    #printf STDERR "ensembl_dbs:\n" . Data::Dumper::Dumper($ensembl_dbs) . "\n";
    #printf STDERR "ucsc_dbs:\n" . Data::Dumper::Dumper($ucsc_dbs) . "\n";

    $state->species_common_names($common_names);
    $state->species_latin_names($latin_names);
    $state->species_latin_to_common_names($latin_to_common_names);
    $state->species_latin_abbrevs($latin_abbrevs);
    $state->species_ensembl_dbs($ensembl_dbs);
    $state->species_ucsc_dbs($ucsc_dbs);


    return $state;
}

sub _session_tmp_file
{
    my $sid = shift;

    return sprintf("%s/$sid", ABS_TMP_PATH);
}

1;
