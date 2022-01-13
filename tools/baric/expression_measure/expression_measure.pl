#!/usr/bin/perl

use FindBin;
use lib "$FindBin::RealBin";

use strict;
use warnings;
use ParamParser;
use GeneralBioinfo;
use General;
use File::Basename;
use Data::Dumper;
use Appli;
use File::Copy;
use Cwd;

my $DEBUG = 0;

MAIN:
{
	my $o_appli = &SetAppli();
	my $o_param = New ParamParser( 'APPLI', $o_appli );

	if ( $o_param->IsDefined('galaxy') )
	{
		print $o_appli->GetGalaxyXml();
		exit 0;
	}
	
	my $count_result = $o_param->Get('count');
	my $rpkm_result  = $o_param->Get('rpkm');

	Initialization($o_param);
	PrepareReference($o_param) if ( $o_param->IsDefined('annotation') );

	# build command line
	my $cmd = 'lipm_rnaseq.pl';
	
	# bidouille mate_dist (comp Appli) => mate-dist (lipm_rnaseq)
	if ( $o_param->IsDefined('mate_dist') )
	{
		my $mate_dist = $o_param->Get('mate_dist');
		$o_param->Delete('mate_dist');
		$o_param->Set('mate-dist', $mate_dist);
	}
	if ( $o_param->IsDefined('skip_rdna') )
	{
		my $skip_rdna = $o_param->Get('skip_rdna');
		$o_param->Delete('skip_rdna');
		$o_param->Set('skip-rdna', $skip_rdna);
	}

	foreach my $param ( $o_param->GetKeys() )
	{
		$cmd .= " --$param " . $o_param->Get($param) if ( $o_param->IsDefined($param) );
	}
	print STDERR "$cmd --ncpus 16 \n" if($DEBUG);
	system($cmd);

	# move results to Galaxy files
	my $outdir  = $o_param->Get('outdir');
	my @a_count = <$outdir/lipm_rnaseq.summary.count.*>;
	my @a_rpkm  = <$outdir/lipm_rnaseq.summary.rpkm.*>;

	if ( @a_count == 0 || @a_rpkm == 0 )
	{
		print STDERR "Error occured: $!\n";
		exit 2;
	}
	
	my $rpkm_ori  =  shift(@a_rpkm);
	move( $rpkm_ori,  $rpkm_result )  or die "Copy failed: $!";
	my $count_ori = shift(@a_count);
	move( $count_ori, $count_result ) or die "Copy failed: $!";
	
	exit 0;
}

=head2 procedure PrepareReference

 Title        : PrepareReference
 Usage        : PrepareReference($o_param)
 Prerequisite : ParamParser
 Procedure    : Build Glint index, search rRNA and create Genome info file
 Args         : $o_param - object - parameters
 Access       : public
 TODO         : none
 Globals      : none

=cut

sub PrepareReference
{
	my $o_param = shift;

	my $workdir = $o_param->Get('outdir') . "/$$/";
	mkdir($workdir);
	$ENV{REFERENCE_REPOSITORY} = $o_param->Get('outdir');

	my %h_files = (
					annotation => $workdir . "$$.gff3",
					reference  => $workdir . "$$",
					rdna       => $workdir . "$$.rdna.gff3",
					rnammer    => $workdir . "$$.rdna.gff",
					glintinfo  => $workdir . "$$.he",
					geninfo    => $workdir . "$$.genome",
	);

	$o_param->AssertFileExists('annotation');
	$o_param->AssertFileExists('reference');

	# link file to workdir
	symlink( $o_param->Get('annotation'), $h_files{annotation} );
	symlink( $o_param->Get('reference'),  $h_files{reference} );

	# link or search rRNA
	system("touch $h_files{rdna}");
#	if ( $o_param->IsDefined('rdna') )
#	{
#		$o_param->AssertFileExists('rdna');
#		symlink( $o_param->Get('rdna'), $h_files{rdna} );
#	}
#	else
#	{
#		my $cmd_rnammer = "rnammer -S " . $o_param->Get('kingdom') . " -gff $h_files{rnammer} $h_files{reference}";
#		print "$cmd_rnammer\n" if($DEBUG);
#		system($cmd_rnammer);
#		ConvertGFF( $h_files{rnammer}, $h_files{rdna} );
#		unlink( $h_files{rnammer} );
#	}

	# glint index
	my $cmd_glint = "glint index $h_files{reference}";
	system($cmd_glint);

	# file genome
	BuildGenomeInfo( $h_files{reference}, $h_files{geninfo} );

	# clean $o_param
	$o_param->Delete('annotation');
	$o_param->Delete('rdna');
	$o_param->Delete('kingdom');
	$o_param->Set( 'reference', $$ );

	return;
}

=head2 procedure BuildGenomeInfo

 Title        : BuildGenomeInfo
 Usage        : BuildGenomeInfo($genome_file, $genome_info_filename)
 Prerequisite : GeneralBioinfo, General
 Procedure    : Build file with molecule accession and lenght for glint
 				Chr1 30427671
 				Chr2 19698289
 				
 Args         : $genome_file - string - filepath of fasta genome file
 				$genome_info_filename - string - filepath of genome info file
 Access       : public
 TODO         : none
 Globals      : none

=cut

sub BuildGenomeInfo
{
	my ( $genome_file, $genome_info ) = @_;

	my %h_seqinfo = ();
	GeneralBioinfo::FastaToHash( $genome_file, \%h_seqinfo, 1 );

	my $f_out = General::GetStreamOut($genome_info);
	foreach my $acc ( keys(%h_seqinfo) )
	{
		print $f_out "$acc $h_seqinfo{$acc}->{len}\n";
	}
	$f_out->close;

	return;
}

=head2 procedure ConvertGFF

 Title        : ConvertGFF (from EugenePP (Erika))
 Usage        : ConvertGFF($infile, $outfile)
 Prerequisite : General
 Procedure    : Convert GFF2 result of RNAmmer to GFF3
 Args         : $infile - string - file path of RNAmmer GFF result
 				$outfile - string - file path of RNAmmer GFF3
 Access       : public
 TODO         : 
 Globals      : none

=cut

sub ConvertGFF
{
	my ( $infile, $outfile ) = @_;

	my $f_in  = General::GetStreamIn($infile);
	my $f_out = General::GetStreamOut($outfile);

	my %h_uniq_copy = ();

	while ( my $line = <$f_in> )
	{
		chomp($line);
		next if ( $line =~ /^#/ || $line eq '' );

		#Mt35Chr5        RNAmmer-1.2     rRNA    25254354        25254468        50.3    +       .       8s_rRNA
		my @a_ = split( /\t/, $line );
		$a_[0] =~ s/ /_/;

		# keep one copy by chromosome
		next if ( exists $h_uniq_copy{"$a_[0]$a_[$#a_]"} );

		my ( $beg, $end, $strand ) = ( $a_[3], $a_[4], $a_[6] );
		my $id = "$a_[0]-rRNA-$a_[$#a_]-$beg-$end";

		print $f_out "$a_[0]\t$a_[1]\t$a_[2]\t$beg\t$end\t$a_[5]\t$strand\t.\tID=$id;Name=$id;Ontology_term=SO:0000252;\n";
		$h_uniq_copy{"$a_[0]$a_[$#a_]"} = 1;
	}
	$f_in->close;
	$f_out->close;

	return;
}

=head2 procedure Initialization

 Title        : Initialization
 Usage        : Initialization($o_param)
 Prerequisite : File::Basename
 Procedure    : Prepare parameters for lipm_rnaseq.pl
 Args         : $o_param - object - parameters
 Access       : private
 TODO         : none
 Globals      : none

=cut

sub Initialization
{
	my $o_param = shift;

	# create workdir
	my $workdir = getcwd() . "/rnaseq";
	mkdir($workdir);
	$o_param->Set( 'outdir', $workdir );

	# add extension
	# make symbolic link to working directory (managed by galaxy)
	my @a_readfiles = split( /,/, $o_param->Get('readfiles') );
	my @a_filenames = split( /,/, $o_param->Get('filenames') );
	my @a_fastqfiles = ();
	my %h_filename_uniq = ();
	
	for(my $index = 0; $index <= scalar(@a_readfiles); $index++)
	{
		my $input = $a_readfiles[$index];
		my $filename = $a_filenames[$index];
		# skip empty
		next if ( !defined $input || $input =~ m/^\s*$/ );
		# trim leading space
		$input =~ s/^\s*//;
		$filename =~ s/^\s*//;
		# remove extension if exists
		$filename =~ s/\.fastq(\.gz)?$//;
		# replace space by '_'
		$filename =~ s/[^\w\d\.-:]+/_/g;
		
		if($o_param->Get('orientation') =~ m/pe/) 
		{
			$filename .= ($index%2)?'_2':'_1';
		}
		
		my $extension = ( -B $input ) ? '.fastq.gz' : '.fastq';
		my $fastq_filename = "$workdir/${filename}$extension";
		symlink( $input, $fastq_filename );
		push( @a_fastqfiles, $fastq_filename );

		# duplicate filename
		if( exists $h_filename_uniq{$fastq_filename} )
		{
			print STDERR "Duplicate file: ${filename}$extension\n";
			exit 2;
		}
		$h_filename_uniq{$fastq_filename} = 1;
	}
	$o_param->Set( 'readfiles', join( ',', @a_fastqfiles ) );
	$o_param->Delete('filenames');

	# remove wrapper parameter
	$o_param->Delete('rpkm');
	$o_param->Delete('count');

	# reference builtin
	if ( !$o_param->IsDefined('annotation') )
	{
		$o_param->Delete('rdna');
		$o_param->Delete('kingdom');
	}

	if ( $o_param->IsDefined('repository') )
	{
		$ENV{REFERENCE_REPOSITORY} = $o_param->Get('repository');
		$o_param->Delete('repository');
	}

	return;
}

sub SetAppli
{

	my %h_program_description = (
								  'name'     => 'LIPM RNAseq',
								  'descr'    => 'LIPM RNAseq',
								  'authors'  => 'jerome.gouzy@toulouse.inra.fr',
								  'category' => 'RNASeq'
	);

	my %h_program_inputs = (
							 'readfiles' => {
											  'label'       => 'Input reads file(s)',
											  'descr'       => 'Reads file(s)',
											  'type'        => 'String',
											  'type_galaxy' => 'fastq',
											  'card'        => '1,n'
							 },
							 'reference' => {
											  'label'       => 'Reference genome file',
											  'descr'       => 'Reference genome',
											  'type'        => 'String',
											  'type_galaxy' => 'fasta',
							 },
							 'filenames' => {
											  'label'       => 'filenames',
											  'descr'       => 'filenames',
											  'type'        => 'String',
											  'type_galaxy' => 'text',
							 }
							 
	);

	my %h_program_outputs = (
							  'rpkm' => {
										  'descr'       => 'Output RPKM result file',
										  'type'        => 'String',
										  'type_galaxy' => 'tabular'
							  },
							  'count' => {
										   'descr'       => 'Output Count result file',
										   'type'        => 'String',
										   'type_galaxy' => 'tabular'
							  }
	);

	my %h_program_parameters = (
								 'annotation' => {
												   'label'       => 'Genome annotation',
												   'descr'       => 'Genome annotation file (GFF3)',
												   'type'        => 'String',
												   'type_galaxy' => 'gff3'
								 },
								 'rdna' => {
											 'label'       => 'Ribosomal annotation',
											 'descr'       => 'Ribosomal RNA annotation file (GFF3)',
											 'type'        => 'String',
											 'type_galaxy' => 'gff3'
								 },
								 'skip_rdna' => {
											 'label'       => 'Skip rRNA step',
											 'descr'       => '',
											 'type'        => 'boolean',
											 'type_galaxy' => 'boolean'
								 },
								 'kingdom' => {
												'label' => 'Kingdom of genome',
												'descr' => 'Kingdom of genome',
												'type'  => 'Choice',
												'enum'  => {
															'bac' => 'Bacteria',
															'euk' => 'Eukaryote',
															'arc' => 'Archaea',
												},
												'default'     => 'bac',
												'type_galaxy' => 'select'
								 },
								 'orientation' => {
													'label' => 'Library type',
													'descr' => 'Library type (Single end, Oriented Single-end, Pair-end, Oriented Pair-end)',
													'type'  => 'Choice',
													'enum'  => {
																'se'  => 'Single end',
																'ose' => 'Oriented Single-end',
																'pe'  => 'Pair-end',
																'ope' => 'Oriented Pair-end'
													},
													'default'     => 'se',
													'type_galaxy' => 'select'
								 },
								 'small' => {
											  'label'       => 'Small inserts analysis',
											  'descr'       => 'Small inserts analysis. Change mapping algorithm.',
											  'type'        => 'boolean',
											  'type_galaxy' => 'boolean'
								 },
								 'mmis' => {
											 'label'       => 'Maximum number of mismatches',
											 'descr'       => 'Maximum number of mismatches',
											 'type'        => 'integer',
											 'type_galaxy' => 'integer',
											 'min'         => 0
								 },
								 'lmin' => {
											 'label'       => 'Minimal hit length',
											 'descr'       => 'Minimal hit length',
											 'type'        => 'integer',
											 'type_galaxy' => 'integer',
											 'default'     => 18,
											 'min'         => 0
								 },
								 'mate_dist' => {
												  'label'       => 'Maximal distance between paired reads',
												  'descr'       => 'Maximal distance between paired reads',
												  'type'        => 'integer',
												  'type_galaxy' => 'integer',
												  'min'         => 0
								 },
								 'oriented_end' => {
													 'label' => 'Oriented end',
													 'descr' => 'Specification of the insert end which gives the orientation of the transcript (depends on protocols). WARNING: the oriented end can be 2 even on oriented single end (the one that is not in the data)',
													 'type'  => 'Choice',
													 'enum'  => {
																 '1' => 'Reads file 1',
																 '2' => 'Reads file 2'
													 },
													 'default'     => '1',
													 'type_galaxy' => 'select'
								 },
								 'nospecific' => {
												   'label'       => 'Use no specific mapped reads',
												   'descr'       => 'Considering all mapped read/pairs (default: considering unambigously mapped reads/pairs only)',
												   'type'        => 'boolean',
												   'type_galaxy' => 'boolean'
								 },
								 'level' => {
											  'label' => 'Expression reported for',
											  'descr' => 'Select biological objects (gff3 type) for which the expression is reported',
											  'type'  => 'Choice',
											  'enum'  => {
														  'gene' => 'Gene',
														  'RNA'  => 'RNA'
											  },
											  'default'     => 'gene',
											  'type_galaxy' => 'select'
								 },
								 'add_annotations' => {
														'label'       => 'Add annotation',
														'descr'       => 'Add annotations found in the repository if available',
														'type'        => 'boolean',
														'type_galaxy' => 'boolean'
								 },
								 'repository' => {
													'label'	=> 'Repository path',
													'label'	=> 'Repository path',
													'type'	=> 'String'
								 },
								 'minimum_overlap_as_a_fraction_of_read_hit' => {
																				  'label'       => 'Minimum overlap',
																				  'descr'       => 'Minimum overlap required as a fraction of the mapped read (intersecBed parameter -f)',
																				  'type'        => 'float',
																				  'type_galaxy' => 'float',
																				  'default'     => 1,
																				  'min'         => 0,
																				  'max'         => 1
								 }
	);

	my $o_appli = New Appli(
							 -general => \%h_program_description,
							 -inputs  => \%h_program_inputs,
							 -outputs => \%h_program_outputs,
							 -params  => \%h_program_parameters
	);

	return $o_appli;
}
