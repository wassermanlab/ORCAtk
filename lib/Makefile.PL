use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
my $DISTNAME = 'ORCA';
my $VERSION = '0.0.1';
WriteMakefile(
	'NAME'		=> 'ORCA',
	'DISTNAME'	=> $DISTNAME,
	'VERSION'	=> $VERSION,
	'dist'      =>	{
		COMPRESS     => 'gzip -9f',
		SUFFIX       => '.gz', 
		DIST_DEFAULT => 'all tardist',
	},
);
