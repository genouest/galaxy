#!/usr/bin/env perl -w
use strict;
use Bio::SearchIO;
use FileHandle;
use Parallel::ForkManager;

#@uthor: Bernhard GSCHLOESSL
#date: 06/10/2012
#last update: 27/04/2021
##getOrthologHitRatio(..) calculates a percentage of query alignment similarity (referred to total hit length)to a specific hit in order to determine if query is ortholog of hit.
#and draws a histogram of OHRs (R script written by F Dorkeld).

&main(@ARGV);

sub main{
  my $xmlBlastOutFile=shift(@_);
  my $threads=4;#number of threads to be used
  my $blastType;
  my @xmlFileList;
  my $pm;
  #fuse output files
  my $fuseOHRcmd="cat ";
  my $rmCmd="rm -f ";
  my $globalOHRoutf="OHRandBlastFeatures.csv";#name for Galaxy output file
  my $globalOHRhistogramf="OHR_distribution.pdf";#histogram of OHR distribution

  $blastType=&getBlastType($xmlBlastOutFile);
  $blastType=lc($blastType);

  @xmlFileList=($xmlBlastOutFile);
  system("echo \"QueryID;HitAccNr;HitFunction;HitBitScore;HitEvalue;QueryALNstart;QueryALNend;nonGapALNlength;Identity;OHR\" > $globalOHRoutf");

  #if several threads are requested -> split XML file in $thread subfiles -> parallelise on $threads threads
  my $queryCount=`grep -e \"<Iteration>\" $xmlBlastOutFile -c`;
  chomp($queryCount);
  $queryCount+=0;


  if($threads>1){
    #retrieve number of sequences per file respecting the requested number of threads
    my $nrOfSeqsPerXMLfile=($queryCount%$threads == 0) ? $queryCount/$threads : int($queryCount/$threads)+1;#add 1 sequence if there is a mod
    #split files into $threads XML files
    @xmlFileList=&splitXmlBLASTfile($nrOfSeqsPerXMLfile,$xmlBlastOutFile,$blastType);
  }

  $pm = Parallel::ForkManager->new($threads);

  #parallelise processing of XML file
  DATA_LOOP:

  foreach my $partialXMLf (@xmlFileList){
    chomp($partialXMLf);

    # Forks and returns the pid for the child:
    my $pid = $pm->start and next DATA_LOOP;

    &getOrthologHitRatio($partialXMLf,$blastType);
    $pm->finish; # Terminates the child process
  }#\treat each partialXMLf file OR the entire file if $threads==1
  $pm->wait_all_children;

  foreach my $partialXMLf (@xmlFileList){
    chomp($partialXMLf);
    $fuseOHRcmd.=$partialXMLf."*.ohr ";
    $rmCmd.=$partialXMLf."* ";
  }
  system($fuseOHRcmd." >> ".$globalOHRoutf);
  system("$rmCmd");

  &drawHistogramWithR($globalOHRoutf,$globalOHRhistogramf);
}#\main(..)


#NEW: -> OHR is calculated on first HSP of first HIT (is ranged by best score)
#Implementation of method suggested by O. Neil et al. 2010
#Ortholog Hit Ratio (OHR) is the ""percentage of an ortholog"" found in a unigene by blast search.
#This percentage is determined by dividing the number of non-gap characters of the query alignment
#by the total length of the subject (hit) sequence.
#the transcript/gene couple is retained which shows the best bitscore and evalue (in case that there is the same bitscore) ->
#if there are multiple best reciprocal hits (score and evalue are identic) all this hits will
#be kept (later bi-directional blast analysis will filter)
sub getOrthologHitRatio{
  #print "Do getOrthologHitRatio(..)\n";
  my $hitTXT="hit";
  my $queryTXT="query";
  my $infile=$_[0];
  my $blastType=$_[1];
  my $alnReferred2=(!defined($_[1])) ? $hitTXT : $queryTXT;#take either hit or query ungapped aln length into account
  my $blastOutputFormat="blastxml";
  my $fhOHR;
  my $allBestOHRf=$infile."_referredTO".$alnReferred2."seq.ohr";

  if($blastType=~m/tblast/){
	die "Check/Update if subject lengths agree OR if the length ratio has to be adapted similar to e.g. blastx!!\n";
  }

  # one can also request that the parser does NOT keep the XML data in memory
  # by using the tempfile initialization flag (-tempfile => 1).
  #Bio::SearchIO::blast
  my $in = new Bio::SearchIO(-tempfile => 1,
			    -format => $blastOutputFormat,
			    -file   => $infile );#$ARGV[0] blastxml,blasttable

  if(!defined($fhOHR = FileHandle->new(">$allBestOHRf"))){ die "Cannot write file $allBestOHRf!\n"; }

  #treat query one by one
  while(my $result=$in->next_result){#Bio::Search::Result::ResultI
    my $nrHits=$result->num_hits;
    my $qname=$result->query_name();
    my $qdesc=$result->query_description();
    my $qLen=$result->query_length();#length of total query sequence

    $qdesc=~s/(^gi\|){0,1}([a-zA-Z0-9]+)(.*$)/$2/g;#keep only the GI geneID
    $qdesc=~s/;/,/g;

    #in case that there is no querys description (e.g. for selfcreated (non-GenBank) FASTA sequences) use the query name as qdesc
    if($qdesc eq ""){
      $qdesc=$qname;
    }
    #if there is at least one hit
    if($nrHits>0){
      my $hit = $result->next_hit;#Bio::Search::Hit::HitI
      my $bestHitAcc=$hit->accession();
      my $bestHitName=$hit->name();#calculate OHR on best hit only
      my $hsp = $hit->next_hsp;#Bio::Search::HSP::HSPI
      my $hdesc=$hit->description;
      $hdesc=~s/;/,/g;
      #retrieve HSP with best OHR
      my $hitLen=$hit->length();#length of total hit sequence
      #NEW refer always to aln based on query sequence
      my $actHspTypeAlnLen=$hsp->length($queryTXT);#NEW refer always to aln based on query sequence ##OLD $hsp->length($alnReferred2);#$hsp->length('query') #alignment length without gaps referred to query OR hit
      #in case of blastn and blastp the query and subject lengths refer to the same sequence type
      #in case of blastx/blastp/tblastx the query length is in nt whereas the subject(->hit) is in aa => if the alignment type of interest is the hitlen the query has to be adapted to aa format (divided by 3)
      if(($blastType!~m/blastn/)&&($alnReferred2 eq $hitTXT)){
	      $actHspTypeAlnLen/=3;
      }

      #By default $actHspTypeAlnLen is now the ungapped query aln lentgth in order to calculate the OHR
      #(->OHR > 1 indicates that query is better reconstituted than hit sequence) -> for this use also "hit" option (->divide by hitLen)
      my $actOHR=($alnReferred2 eq $queryTXT) ? $actHspTypeAlnLen/$qLen : $actHspTypeAlnLen/$hitLen ;#taken from O'Neil et al., 2010: ATTENTION: OHR refers here only to the longest HSP of a hit!
      my $actIdentity=$hsp->percent_identity();#returns calculated percent identity for the given HSP
      my $actScore = $hsp->bits();#BIT score for this HSP
      my $actEvalue = $hsp->evalue();
      my $hspQueryStart=$hsp->start($queryTXT);
      my $hspQueryEnd=$hsp->end($queryTXT);
      $fhOHR->print($qdesc.";".$bestHitAcc.";".$hdesc.";".$actScore.";".$actEvalue.";".$hspQueryStart.";".$hspQueryEnd.";".$actHspTypeAlnLen.";".$actIdentity.";".$actOHR."\n");
    }#\if there are hits
  }#\treat each query
  $fhOHR->close;
}#\getOrthologHitRatio(..)



#this function is baseq on awk code originally written by Laurent Manchon (lmanchon@univ-montp2.fr )
sub splitXmlBLASTfile{
  my $nrOfSeqsPerXMLfile=$_[0];# max is number of sequences per output file
  my $xmlBlastOutFile=$_[1];
  my $blastType=(defined($_[2])) ? $_[2] : "blastx";
  my $queryCount=0;
  my @xmlFiles=();
  my $xmlPrefix=$xmlBlastOutFile;
  $xmlPrefix=~s/(.+)(\.xml)/$1/g;
  my $actXMLoutf;
  my $outfCount=0;#count of XML outfiles
  my $write=0;#boolean indicating if a line should be written or not
  my $blastOutputVersion=`grep  -e  \"<BlastOutput_version>\" $xmlBlastOutFile -m1 `;
  chomp($blastOutputVersion);

  my $outSuffix=".xml";
  my $blastXMLbegin="<?xml version=\"1.0\"?>\n<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n";
  $blastXMLbegin.="<BlastOutput>\n<BlastOutput_program>$blastType</BlastOutput_program>\n$blastOutputVersion\n";
  $blastXMLbegin.="<BlastOutput_reference></BlastOutput_reference>\n<BlastOutput_db>/home/data/blastdb/nr</BlastOutput_db>\n<BlastOutput_query-ID>lcl|1_0</BlastOutput_query-ID>\n";
  $blastXMLbegin.="<BlastOutput_query-def></BlastOutput_query-def>\n<BlastOutput_query-len></BlastOutput_query-len>\n<BlastOutput_param>\n<Parameters>\n<Parameters_matrix>BLOSUM62</Parameters_matrix>\n";
  $blastXMLbegin.="<Parameters_expect>0.1</Parameters_expect>\n<Parameters_gap-open>11</Parameters_gap-open>\n<Parameters_gap-extend>1</Parameters_gap-extend>\n<Parameters_filter>F</Parameters_filter>\n";
  $blastXMLbegin.="</Parameters>\n</BlastOutput_param>\n<BlastOutput_iterations>\n";
  my $blastXMLend="</BlastOutput_iterations>\n</BlastOutput>";


  open(RH,"< $xmlBlastOutFile") || die "Cannot read file $xmlBlastOutFile!\n";
  while(my $line=<RH>){
    #start of a new query
    if($line=~m/<Iteration>/){
      $queryCount++;
      $write=1;

      #set new output XML file
      if(($queryCount==1)||($queryCount==$nrOfSeqsPerXMLfile+1)){
	$outfCount++;

	#write XML end tags to last file
	if($queryCount==$nrOfSeqsPerXMLfile+1){
	  open(WH,">>$actXMLoutf") || die "Cannot write $actXMLoutf!\n";
	  print WH $blastXMLend;
	  close WH;
	  $queryCount=1;#reset query count
	}
	#set new XML output file
	$actXMLoutf=$xmlPrefix."_".$outfCount.$outSuffix;
	push(@xmlFiles,$actXMLoutf);

	open(WH,">$actXMLoutf") || die "Cannot write $actXMLoutf!\n";
	print WH $blastXMLbegin;
	close WH;
      }
    }

    if($write){
      open(WH,">>$actXMLoutf") || die "Cannot write $actXMLoutf!\n";
      print WH $line;
      close WH;
    }

    if($line=~m/<\/Iteration>/){
      $write=0;
    }
  }#\while reading input XML file
  close RH;

  #write last XML end tags
  open(WH,">>$actXMLoutf") || die "Cannot write $actXMLoutf!\n";
  print WH $blastXMLbegin;
  close WH;

  return @xmlFiles;
}#\splitXmlBLASTfile(..)


#check blast type in XML BLAST output file
sub getBlastType{
    my $xmlBlastOutFile=$_[0];
    my $blastType=`grep -e \"<BlastOutput_program>\" $xmlBlastOutFile -m1 `;
    chomp($blastType);

    $blastType=~s/(<BlastOutput_program>)([a-zA-Z]+)(<\/BlastOutput_program>)/$2/g;
    #print "The blast type is $blastType!\n";

    return $blastType;
}#\getBlastType(..)


#R script contributed by Franck Dorkeld
sub drawHistogramWithR{
  my $inf=shift(@_);
  my $outPDF=shift(@_);
  my $rScript="all\<-read.csv(file=\"$inf\",header=TRUE,sep=\";\"); maxOHR\<-round(max(all\$OHR),digits=1)+0.2; bins\<-seq(0.0,maxOHR,by=0.01);pdf(\"$outPDF\");hist(all\$OHR,breaks=bins,xlim=c(0,maxOHR),xlab=\"ortholog hit ratio\",ylab=\"transcript count\",main=\"OHR distribution\");dev.off();";

  system("R -e \'$rScript\'");
}#\drawHistogramWithR(..)

