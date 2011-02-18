# Version
use constant VERSION		=> '1.0.0';

# Contact info
use constant CONTACT_NAME	=> 'David Arenillas';
use constant CONTACT_EMAIL	=> 'dave@cmmt.ubc.ca';

# Libs
use constant ORCA_LIB		=> '/space/devel/ORCAtk/lib';
# Use current Ensembl API
use constant ENSEMBL_LIB	=> '/usr/local/src/ensembl-current/ensembl/modules';

# Environment variables
use constant LD_LIBRARY_PATH	=> '/usr/local/mysql/lib/mysql';

# Paths and URLs
use constant ABS_HTML_PATH	=> '/usr/local/apache/htdocs/ORCAtk';
use constant ABS_CGI_BIN_PATH	=> '/usr/local/apache/cgi-bin/ORCAtk';
use constant ABS_TEMPLATE_PATH	=> '/usr/local/apache/htdocs/ORCAtk/templates';
use constant ABS_TMP_PATH	=> '/usr/local/apache/htdocs/ORCAtk/tmp';
use constant REL_TMP_PATH	=> '/ORCAtk/tmp';
use constant REL_CGI_BIN_PATH	=> '/cgi-bin/ORCAtk';
use constant REL_HTML_PATH	=> '/ORCAtk';
use constant ABS_LOG_PATH	=> '/tmp';

# Constants
use constant MAX_FILE_SIZE	=> 2000000;	# limit upload file sizes
use constant MAX_SEQ_LEN	=> 1000001;	# limit sequence lengths

use constant ENSEMBL_DATABASES_DB_NAME	=> "ensembl_databases";

# No longer used - fetch from ensembl_databases db
use constant SPECIES		=> ['human','mouse','armadillo','cat',
				    'chicken','chimp','cow','dog','elephant',
				    'fly','frog','fugu','guinea pig','horse',
				    'mosquito','opossum','platypus','rabbit',
				    'rat','shrew','worm','yeast','zebrafish'];

use constant SPECIES_ABBREVS	=> {'human'	=> 'hsapiens',
				    'mouse'	=> 'mmusculus',
				    'armadillo'	=> 'dnovemcinctus',
				    'chicken'	=> 'ggallus',
				    'dog'	=> 'cfamiliaris',
				    'horse'	=> 'ecaballus',
				    'rabbit'	=> 'ocuniculus',
				    'rat'	=> 'rnorvegicus'
				   };

use constant SPECIES_LATIN_NAMES => {'human'	=> 'Homo sapiens',
				    'mouse'	=> 'Mus musculus',
				    'armadillo'	=> 'Dasypus novemcinctus',
				    'chicken'	=> 'Gallus gallus',
				    'dog'	=> 'Canis familiaris',
				    'horse'	=> 'Equus caballus',
				    'rabbit'	=> 'Oryctolagus cuniculus',
				    'rat'	=> 'Rattus norvegicus'
				   };

# No longer used - fetch from ensembl_databases db
use constant SPECIES_ENSEMBL_DBS => {
			    'human'	=> 'homo_sapiens_core_51_36m',
			    'mouse'	=> 'mus_musculus_core_51_37d',
			    'armadillo'	=> 'dasypus_novemcinctus_core_51_1g',
			    'chicken'	=> 'gallus_gallus_core_51_2i',
			    'dog'	=> 'canis_familiaris_core_51_2i',
			    'horse'	=> 'equus_caballus_core_51_2a',
			    'rabbit'	=> 'oryctolagus_cuniculus_core_51_1g',
			    'rat'	=> 'rattus_norvegicus_core_51_34t'
			};

# No longer used - fetch from ensembl_databases db
use constant SPECIES_UCSC_DBS => {
				    'human'	=> 'hg18',
				    'mouse'	=> 'mm9'
				};

use constant SPECIES_TRACK_NAMES => {
				    'human'	=> 'phastCons28way',
				    'mouse'	=> 'phastCons30way'
				};

#use constant ENSEMBL_DB_HOST	=> 'ensembldb.ensembl.org';
#use constant ENSEMBL_DB_PORT	=> '5306';
#use constant ENSEMBL_DB_USER	=> 'anonymous';
use constant ENSEMBL_DB_HOST	=> 'vm5.cmmt.ubc.ca';
use constant ENSEMBL_DB_PORT	=> '3306';
use constant ENSEMBL_DB_USER	=> 'ensembl_r';
use constant ENSEMBL_DB_PASS	=> '';

use constant ENSEMBL_COMPARA_DB_HOST	=> 'vm5.cmmt.ubc.ca';
use constant ENSEMBL_COMPARA_DB_PORT	=> '3306';
use constant ENSEMBL_COMPARA_DB_USER	=> 'ensembl_r';
use constant ENSEMBL_COMPARA_DB_NAME	=> 'ensembl_compara_51';
use constant ENSEMBL_COMPARA_DB_PASS	=> '';

use constant JASPAR_DB_HOST	=> 'vm5.cmmt.ubc.ca';
use constant JASPAR_DB_USER	=> 'jaspar_r';
use constant JASPAR_DB_PASS	=> '';
use constant JASPAR_CORE_DB	=> 'JASPAR_CORE_2008';
use constant JASPAR_PHYLO_DB	=> 'JASPAR_PHYLOFACTS_2008';
use constant JASPAR_FAM_DB	=> 'JASPAR_FAM_2008';
use constant JASPAR_DBS		=> [JASPAR_CORE_DB,
				    JASPAR_PHYLO_DB,
				    JASPAR_FAM_DB];

# Default values
use constant SPECIES1			=> 'human';
use constant SPECIES2			=> 'mouse';
use constant CA_TOP_PERCENTILE		=> 10.0;
# Min. conservation default is only defined for multi-species (phastCons)
# analysis. For pairwise analysis CA_TOP_PERCENTILE takes precedence (although
# a min. conservation can still be specified by user.
use constant CA_MIN_CONS_PAIRWISE	=> undef;
use constant CA_MIN_CONS_PHASTCONS	=> 70.0;
use constant CA_FILTER_EXONS		=> 1;
use constant CA_WINDOW_SIZE		=> 100;
# Multi-species (phastCons) conserved regions tend to me much smaller
use constant CA_MIN_CR_LEN_PAIRWISE	=> 50;
use constant CA_MIN_CR_LEN_PHASTCONS	=> 20;

use constant TF_THRESHOLD		=> '80.0';
use constant TF_MIN_IC			=> 10;
use constant TF_FILTER_SITES		=> 1;
use constant TF_SELECT_METHOD		=> 'min_ic';

use constant ALIGNMENT_FORMAT		=> 'clustalw';

use constant UPSTREAM_BP		=> 5000;
use constant DOWNSTREAM_BP		=> 0;
use constant DOWN_REL_TO		=> 'end';
use constant MIN_TFBS_CONSERVATION_OVERLAP	=> 1;

# Transcript graphical display params
use constant TRANSCRIPT_PANEL_WIDTH	=> 400;
use constant TRANSCRIPT_BGCOLOUR	=> "cyan";
use constant TRANSCRIPT_FGCOLOUR	=> "cyan";

# Delete temporary files older than this number of days
use constant CLEAN_TEMPFILES_OLDER_THAN	=> 3;

1;
