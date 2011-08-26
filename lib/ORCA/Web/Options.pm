# Version
use constant VERSION => '1.0.0';

# Contact info
use constant CONTACT_NAME   => 'David Arenillas';
use constant CONTACT_EMAIL  => 'dave@cmmt.ubc.ca';

use constant ORCA_HOME  => '/disk1/ORCAtk';

# Libs
use constant ORCA_LIB   => ORCA_HOME . '/lib';

#
# We don't want to automatically update to current Ensembl libs/DBs in case
# the genome assembly becomes out of sync with that of the current UCSC
# phastCons DB
#
#use constant ENSEMBL_LIB            => '/usr/local/src/ensembl-current/ensembl/modules';
#use constant ENSEMBL_COMPARA_LIB    => '/usr/local/src/ensembl-current/ensembl-compara/modules';
use constant ENSEMBL_LIB            => '/usr/local/src/ensembl-62/ensembl/modules';
use constant ENSEMBL_COMPARA_LIB    => '/usr/local/src/ensembl-62/ensembl-compara/modules';

# Environment variables
use constant LD_LIBRARY_PATH    => '/usr/local/mysql/lib/mysql';

# Paths and URLs
use constant ABS_HTML_PATH      => '/usr/local/apache/htdocs/ORCAtk';
use constant ABS_CGI_BIN_PATH   => '/usr/local/apache/cgi-bin/ORCAtk';
use constant ABS_TEMPLATE_PATH  => '/usr/local/apache/htdocs/ORCAtk/templates';
use constant ABS_TMP_PATH       => '/usr/local/apache/htdocs/ORCAtk/tmp';
use constant REL_TMP_PATH       => '/ORCAtk/tmp';
use constant REL_CGI_BIN_PATH   => '/cgi-bin/ORCAtk';
use constant REL_HTML_PATH      => '/ORCAtk';
use constant ABS_LOG_PATH       => ORCA_HOME . '/log';

# Constants
use constant MAX_FILE_SIZE  => 2000000;    # limit upload file sizes
use constant MAX_SEQ_LEN    => 1000001;    # limit sequence lengths

use constant ENSEMBL_DATABASES_DB_NAME => "ensembl_databases";

# XXX No longer using ensembl_databases DB for this!
use constant SPECIES_UCSC_DB_NAMES => {
    'human' => 'hg19',
    'mouse' => 'mm9',
    'rat'   => 'rn4',
    'fly'   => 'dm3',
    #'worm'  => 'ce6',
    'yeast' => 'sacCer2'
};

use constant SPECIES_TRACK_NAMES => {
    'human' => 'phastCons46wayPlacental',
    'mouse' => 'phastCons30wayPlacental',
    'rat'   => 'phastCons9way',
    'fly'   => 'phastCons15way',
    #'worm'  => 'phastCons6way',
    'yeast' => 'phastCons7way'
};

#use constant ENSEMBL_MART_HOMOL_TBL_TEMPLATE
#                                    => '%s_gene_ensembl__homolog_%s__dm';

#use constant ENSEMBL_DB_HOST	    => 'ensembldb.ensembl.org';
#use constant ENSEMBL_DB_PORT	    => '5306';
#use constant ENSEMBL_DB_USER	    => 'anonymous';
use constant ENSEMBL_DB_HOST => 'vm2.cmmt.ubc.ca';
use constant ENSEMBL_DB_PORT => '3306';
use constant ENSEMBL_DB_USER => 'ensembl_r';
use constant ENSEMBL_DB_PASS => '';

#use constant ENSEMBL_MART_DB_HOST	=> 'martdb.ensembl.org';
#use constant ENSEMBL_MART_DB_PORT	=> '5316';
#use constant ENSEMBL_MART_DB_USER	=> 'anonymous';
#use constant ENSEMBL_MART_DB		=> 'ensembl_mart_49';
#use constant ENSEMBL_MART_DB_HOST	=> 'vm2.cmmt.ubc.ca';
#use constant ENSEMBL_MART_DB_PORT	=> '3306';
#use constant ENSEMBL_MART_DB_USER	=> 'ensembl_r';
#use constant ENSEMBL_MART_DB		=> 'ensembl_mart_56';
#use constant ENSEMBL_MART_DB_PASS	=> '';

use constant JASPAR_DB_HOST => 'vm5.cmmt.ubc.ca';
use constant JASPAR_DB_NAME => 'JASPAR_2010';
use constant JASPAR_DB_USER => 'jaspar_r';
use constant JASPAR_DB_PASS => '';

# Collections OTHER than CORE
use constant JASPAR_COLLECTIONS =>
    ['PHYLOFACTS', 'FAM', 'PBM', 'PBM_HOMEO', 'PBM_HLH'];

use constant JASPAR_TAX_GROUPS =>
    ['vertebrates', 'insects', 'nematodes', 'plants', 'fungi'];

# Default values
use constant SPECIES1          => 'human';
use constant SPECIES2          => 'mouse';
use constant CA_TOP_PERCENTILE => 10.0;
# Min. conservation default is only defined for multi-species (phastCons)
# analysis. For pairwise analysis CA_TOP_PERCENTILE takes precedence (although
# a min. conservation can still be specified by user.
use constant CA_MIN_CONS_PAIRWISE  => undef;
use constant CA_MIN_CONS_PHASTCONS => 70.0;
use constant CA_FILTER_EXONS       => 1;
use constant CA_WINDOW_SIZE        => 100;
# Multi-species (phastCons) conserved regions tend to me much smaller
use constant CA_MIN_CR_LEN_PAIRWISE  => 50;
use constant CA_MIN_CR_LEN_PHASTCONS => 20;

use constant TF_THRESHOLD     => '80.0';
use constant TF_MIN_IC        => 10;
use constant TF_FILTER_SITES  => 1;
use constant TF_SELECT_METHOD => 'min_ic';

use constant ALIGNMENT_FORMAT => 'clustalw';

use constant UPSTREAM_BP                   => 5000;
use constant DOWNSTREAM_BP                 => 0;
use constant DOWN_REL_TO                   => 'end';
use constant MIN_TFBS_CONSERVATION_OVERLAP => 1;

# Transcript graphical display params
use constant TRANSCRIPT_PANEL_WIDTH => 400;
use constant TRANSCRIPT_BGCOLOUR    => "cyan";
use constant TRANSCRIPT_FGCOLOUR    => "cyan";

# Delete temporary files older than this number of days
use constant CLEAN_TEMPFILES_OLDER_THAN => 3;

1;
