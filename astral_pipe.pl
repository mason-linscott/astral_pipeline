#! /usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;
use Cwd;
use Cwd 'abs_path';
use List::Util 'shuffle';
use File::Path;
use File::Spec;
#use File::Copy qw( move );
#use File::Copy::Recursive qw( dirmove );

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( 'ab:c:g:hi:l:m:M:o:p:P:rs:S:u', \%opts );

if( $opts{u} ){
  &updates;
  die "Exiting program because the list of updates was requested.\n\n";
}

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# declare variables
my $dir = getcwd(); # get current working directory
chomp( $dir );
my $raxml = "raxml"; # name of directory to be created for holding raxml output
my $raxmldir = $dir . "/" . $raxml; # full path to phylip file directory
my $phylip = "phylip"; # directory to be created for holding phylip files
my $phylipdir = $dir . "/" . $phylip; # full path to phylip file directory
my $grouped_loci = "grouped_loci"; # directory to be created for holding grouped phylip files
my $grouped_locidir = $dir . "/" . $grouped_loci; # full path to folder holding grouped phylip files
my $astral = "astral"; # directory to be created for holding astral files
my $astraldir = $dir . "/" . $astral; # full path to folder holding astral files
my $treeout = "best_trees.tre"; # file to hold concatenated trees
my $treeoutdir = join( '/', $astraldir, $treeout );
my $bs = "bs_path.txt"; # file to hold paths to bootstrap trees
my $bspath = join( '/', $astraldir, $bs ); # path to file holding paths to bootstrap trees

# parse the command line
my( $loci, $map, $proc, $groupsize, $astralonly, $boot, $constrain, $version, $remove, $memory, $pis, $snps, $model ) = &parsecom( \%opts );

my $mem = &setmem($memory);

# check if the specified astral version is installed
my $astralpath = system( "which $version | cat" ) ==0 or die "\n$version is not installed: $!\nMake sure you specify the full name of the executable jar file for astral.\nFor example, astral.4.10.12.jar.\n\n";
chomp( $astralpath );

# insert a loop to run only astral on previously generated trees and kill program
if( $astralonly == 1 ){
	my $astraloutdir = &runastral( $treeoutdir, $astraldir, $astralpath, $boot, $bspath, $version, $mem );
	&astralparse( $astraldir, $astraloutdir );
	die "\n\nExiting program because pipeline was called to run Astral only.\n\n";
}

# create input / output directories
mkpath( [$raxml], 1, 0755 ) or die "Couldn't make directory to hold raxml files in $dir: $!";
mkpath( [$phylip], 1, 0755 ) or die "Couldn't make directory to hold phylip files in $dir: $!";
mkpath( [$grouped_loci], 1, 0755 ) or die "Couldn't make directory to hold grouped locus files in $dir: $!";
mkpath( [$astral], 1, 0755 ) or die "Couldn't make directory to hold astral input and output in $dir: $!";

# parse the species map
my $indhash = &parsemap( $map, $astraldir );

# capture number of individuals
my $count = scalar keys %$indhash;

# take list of loci from pyrad and convert it to series of phylip files
my $locuscount = &pyloci( $loci, $indhash, $count, $phylipdir, $snps, $pis );

# randomly group loci together
my $full = &grouploci( $locuscount, $groupsize, $phylipdir, $grouped_locidir, $count );

# open directory containing phylip files
opendir( WD, $grouped_locidir ) or die "Can't open $grouped_locidir: $!\n";

# read contents of directory
my @contents = readdir( WD );

# remove . and .. from list
shift @contents;
shift @contents;

# constrain raxml to only run on a subset of files
if( $constrain > 0 ){
	print "My constraint = ", $constrain, "\n";
	@contents = shuffle @contents;
	# select out first N files for running trees
	@contents = splice( @contents, 0, $constrain );
	
	print "RAxML will run on:\n";
	foreach my $item( @contents ){
		print $item, "\n";
	}
}

# double the number of bootstraps for RAxML so the -g option can be used in astral
$boot = ($boot*2);

# run raxml on each locus and store tree from each
foreach my $file( @contents ){
	if( $file =~ /\.phy$/ ){
		if($remove == 1){
			$file = &removeMissing("$grouped_locidir/$file");		
		}
		&runraxml( $proc, $file, $raxmldir, $grouped_locidir, $boot, $model );
		next;
	}
}

# close the directory
closedir( WD );

# concatenate trees from raxml runs into single file as input to astral
&treecat( $raxmldir, $full, $astraldir, $treeoutdir );

# generate file of locations of bootstrap trees
if( $boot > 0 ){
	&bspaths( $raxmldir, $astraldir, $bspath );
}

# return the number of bootstraps to normal so the -g option can be used in astral
$boot = ($boot/2);

# run astral
my $astraloutdir = &runastral( $treeoutdir, $astraldir, $astralpath, $boot, $bspath, $version, $mem );

# parse astral trees
&astralparse( $astraldir, $astraloutdir );


exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
  print "\nastral_pipe.pl is a perl script developed by Steven Michael Mussmann\n";
  print "This script executes the astral pipeline using the .loci file from pyRAD and a species map text file as input.\n\n";
  print "To report bugs send an email to mussmann\@email.uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  print "Program Options:\n";
  print "\t\t[ -a | -b | -c | -g | -l | -m | -M | -p | -P | -r | -s | -S | -v ]\n\n";
  
  print "\t-a:\tUse this flag if you want to re-run astral on a pre-existing dataset generated by this pipeline.\n";
  print "\t\tBy default this is turned off.\n\n";
  
  print "\t-b:\tUse this flag to specify the number of bootstraps to be conducted in both RAxML and astral.\n";
  print "\t\tWARNING: This number will be multiplied by the number of processors.\n";
  print "\t\tFor example, if you are running on a 32 processor machine, specify a value of either 3 or 4 to do approximately 100 bootstraps.\n";
  print "\t\tDefault value: 0.\n\n";
  
  print "\t-c:\tUse this flag to constrain the number of trees to be built by RAxML.\n";
  print "\t\tBy default this pipeline will run on all loci provided to it. Provide this flag the maximum number of trees you wish to generate in RAxML.\n\n";
  
  print "\t-g:\tUse this flag to randomly concatenate groups of loci.\n";
  print "\t\tThe default value is 50 loci.\n\n";
  
  print "\t-l:\tUse this flag to specify the input .loci file from pyRAD.\n";
  print "\t\tThe program will die if no input file is provided.\n\n";
  
  print "\t-m:\tUse this flag to specify your species map file.\n";
  print "\t\tThis is a two-column tab-delimited file with your individual identifier in the left column and species name in the right column.\n";
  print "\t\tThe program will die if no species map file is provided.\n\n";
  
  print "\t-M:\tUse this flag to specify the amount of memory available to ASTRAL.\n";
  print "\t\tThe default value is 4GB.  You may need to make adjustments to make ASTRAL run properly on your system\n\n";

  print "\t-p:\tUse this flag to specify the number of processors you wish to use for running RAxML.\n";
  print "\t\tThe program will die if no number is provided.\n\n";
  
  print "\t-P:\tUse this flag to specify a minimum number of parsimony informative sites that must be present in a locus.\n";
  print "\t\tThe default value is 0 sites (no loci will be thrown out).\n\n";
  
  print "\t-r:\tUse this flag to remove drop individuals from grouped loci files when their sequences consist only of Ns\n";
  print "\t\tBy default this option is turned off.\n\n";

  print "\t-s:\tSNP filter.  Use this filter to specify a minimum number of SNPs that must be contained in a locus for it to be retained.\n";
  print "\t\tThe default value is 2 SNPs.\n\n";

  print "\t-S:\tSequence evolution model for RAxML.  Use this to set the model of sequence evolution.\n";
  print "\t\tThe default value is GTRGAMMA.\n\n";
  
  print "\t-v:\tUse this flag to specify the version of astral you wish to run.\n";
  print "\t\tastral.4.10.12.jar is used by default.  Feature is not yet fully implemented.\n\n";
  
}

#####################################################################################################
# subroutine to print a list of updates
sub updates{
  
  print "\nastral_pipe.pl is a perl script developed by Steven Michael Mussmann\n";
  print "This script executes the astral pipeline using the .loci file from pyRAD and a species map text file as input.\n\n";
  print "To report bugs send an email to mussmann\@email.uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  print "Recent Updates:\n";
  
  print "\t2016-April-01: \tPipeline was updated to include the -v flag.\n\t\t\t\tThis allows you to use a version of astral other than the one that is hard-coded into the pipeline (version 4.10.12)\n\t\t\t\tWARNING: Not fully implemented yet.\n\n";
  
}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
  my( $params ) =  @_;
  my %opts = %$params;
  
  # set default values for command line arguments
  # if this value = 1, only astral will be run on pre-existing raxml output
  my $astralonly = 0;
  if( $opts{a} ){
	$astralonly = 1;
  }

  my $remove = 0;
  if( $opts{r} ){
	$remove = 1;
  }
 
  my $pis = $opts{P} || "0";
   
  my $boot = $opts{b} || "0"; # sets the number of bootstraps
  
  my $constrain = $opts{c} || "0"; # sets constraint on number of trees to generate
  #$constrain += 1;
  my $groupsize = $opts{g} || "50"; #sets the number of loci to be grouped together randomly
  
  my $loci = $opts{l} or die "\nNo input file specified\n\n"; # used to specify input file name.  Program will die if no file name is provided.

  my $map = $opts{m} or die "\nNo species map file specified\n\n"; # used to specify input species map file.  Program will die if no file name is provided.
  
  my $proc = $opts{p} or die "\nNumber of processors not specified.\n\n"; # used to specify number of processors for parallel steps
  
  my $version = $opts{v} || "astral.4.10.12.jar"; # sets the version of astral to be used.  Default is astral.4.10.12.
  
  my $memory = $opts{M} || "4"; # sets the amount of memory available to ASTRAL

  my $snp = $opts{s} || "2"; # sets the SNP filter

  my $model = $opts{S} || "GTRGAMMA"; #sets the model used by RAxML
  # calculate number of bootstraps to be done
  $boot = $boot * $proc;
  
  return( $loci, $map, $proc, $groupsize, $astralonly, $boot, $constrain, $version, $remove, $memory, $pis, $snp, $model );

}

#####################################################################################################
# subroutine to convert loci file from pyrad to individual aligned phylip files

sub pyloci{
	# requires list of individuals - this will be used to make hash keys
	my( $loci, $indhash, $count, $phylipdir, $snps, $pis ) = @_;
	
	# print status update
	print "\nNow writing phylip files for individual loci\n";
	
	# declare variables
	my @temp; # temporary array to hold all loci
	my @genes; # array to hold all loci for all individuals
	
	# read the input loci file
	open( LOCI, $loci ) or die "Can't open $loci: $!\n";
	while( my $line = <LOCI> ){
		# read each line
		if( $line =~ /^>/ ){
		chomp( $line );
			# push individual onto temporary array
			push( @temp, $line );
		}elsif( $line =~ /^\/\// ){
			my @tokens = split( /\s+/, $line );
			my $locus_pis = 0;
			my $locus_snps = 0;
			foreach my $token( @tokens ){
				if( $token =~ /\-/ ){
					$locus_snps++;
				}elsif( $token =~ /\*/){
					$locus_pis++;
				}
			}
			# print $locus_pis, "\n";
			# print $locus_snps, "\n";
			if( ($locus_snps >= $snps) && ($locus_pis >= $pis) ){
				# print "TRUE", "\n";
				# combine all individuals into comma separated string
				my $string = join( ',', @temp );
				# push string onto @genes
				push( @genes, $string );
			}
			print "\n";
			# clear temporary array
			@temp = ();
		}
	}
	close LOCI;
	
	my $locuscount = scalar( @genes ); # get the number of loci recovered
	
	# for each locus write a phylip file
	my $counter = 1;
	foreach my $locus( @genes ){
		# strip > from beginning of individual names
		$locus =~ s/>//g;
		
		# declare variables
		my @present; # array to hold names of individuals with data present at locus
		my %hash; # temporary hash to hold names of individuals as key values for each locus
		
		# create new output file name
		my $out = join( '.', "locus", $counter, "phy" );
		
		my $outdir = join( '/', $phylipdir, $out );
		
		# open the output file for writing
		open( OUT, '>', $outdir ) or die "Can't write $outdir: $!\n";
		#print "Writing $out\n";
		
		# write header for phylip file
		print OUT $count;
		print OUT "\t";
		
		
		# split each individual
		my @inds = split( /,/, $locus );
		
		# get length of alignment
		my @getlength = split( /\s+/, $inds[0] );
		my $length = length($getlength[1]);
		
		# write the remainder of the header to the phylip file
		print OUT $length;
		print OUT "\n";
		
		# loop through all individuals
		foreach my $ind (@inds ){
			my @indtemp = split( /\s+/, $ind );
			
			# push name onto array
			push( @present, $indtemp[0] );
			
			#write output to file
			print OUT $indtemp[0];
			print OUT "\t";
			print OUT $indtemp[1];
			print OUT "\n";
		}
		
		#convert names to hash
		@hash{@present} = 1;
		
		# compare the two hashes
		foreach( keys %$indhash ){
			if( exists $hash{$_} ){
				next;
			}else{
				print OUT $_;
				print OUT "\t";
				print OUT "N" x $length;
				print OUT "\n";
			}
		}
		
		# close the output file
		close OUT;
		$counter++;
	}
	
	return( $locuscount );
}

#####################################################################################################
# subroutine to remove sequences with no data
sub removeMissing{
	my( $file ) = @_;
	
	#subroutine variables
	my @lines; #holds lines from the input file
	my @kept; #holds lines to keey from the input file	
	my $keptfile = "$file.kept";
	
	open( PHY, $file ) or die "Can't open $file: $!\n";
	
	#read the lines from the file
	while( my $line = <PHY> ){
		chomp( $line );
		push( @lines, $line );
	}
	close PHY;

	#get the length of the alignment
	my @temp = split( /\t/, $lines[0] );
	my $length = $temp[1];
	
	#make string of Ns of same length as alignment
	my $Nstring = ("N" x $length);
	
	#print $Nstring, "\n";
	
	# get rid of first line
	shift @lines;
	
	# check which lines to keep
	foreach my $sequence( @lines ){
		my @tempseq = split( /\t/, $sequence );
		if($tempseq[1] ne $Nstring ){
			push(@kept, $sequence);
		}
	}

	#write the kept sequences to the output file
	open( KEPT, '>', $keptfile ) or die "Can't open $file: $!\n";
	
	print KEPT scalar(@kept), "\t", $length, "\n";

	foreach my $seq( @kept ){
		print KEPT $seq, "\n";
	}	


	close KEPT;

	my @path = split(/\//, $keptfile);

	$keptfile = pop(@path);

	return $keptfile;
	
}

#####################################################################################################
# subroutine to parse species map
sub parsemap{
	my( $map, $astraldir ) = @_;
	
	# declare variables
	my @inds; # array to hold names of individuals
	my %indhash; #hash to hold names of individuals as key values
	my %spechash; # hash to hold groups of species
	my $astral_specmap = "astral_species_map.txt"; # text file to hold species map input to astral
	my $astral_specmapdir = join( '/', $astraldir, $astral_specmap );
	
	# read the input map file
	open( MAP, $map ) or die "Can't open $map: $!\n";
	
	# go through each line and get the names of individuals
	while( my $line = <MAP> ){
		chomp( $line );
		my @temp = split( /\t/, $line );
		push( @inds, $temp[0] );
		my $ind = $temp[0];
		my $species = $temp[1];
		push( @{$spechash{$species}}, $ind );
	}
	
	# close the map file
	close MAP;

	# make individual hash
	@indhash{@inds} = 1;

	# uncomment next line to view %spechash
	#print Dumper( \%spechash );
	
	# open the new species map file for writing
	open( NEWMAP, '>', $astral_specmapdir ) or die "Can't open $astral_specmapdir: $!\n";
	
	# print the contents of the %spechash to new file
	foreach my $key( sort keys %spechash ){
		# print the species name and colon to the new file
		print NEWMAP $key, ":";
		my @temp; # array to temporarily hold all names of individuals belonging to the species
		# push all the individual names to an array
		foreach my $ind( @{$spechash{$key}} ){
			push( @temp, $ind );
		}
		# join the names of all individuals (taxa) together by commas
		my $taxa = join( ",", @temp );
		# print the string of taxa and a newline character to the new species map file
		print NEWMAP $taxa, "\n";
	}
	
	# close new species map file
	close NEWMAP;
	
	# return the indhash and $astral_specmap
	return \%indhash;
	
}

#####################################################################################################
# subroutine to randomly group sets of loci
sub grouploci{
	my( $max, $groupsize, $phylipdir, $grouped_locidir, $count ) = @_;

	# provide status update
	print "\nNow writing phylip files for grouped loci\n";
	
	#generate array of numbers 1 to N
	my @numarray = (1..$max);
	
	#shuffle array of numbers
	@numarray = shuffle( @numarray );
	
	# evaluate the number of times the groupsize can go into the number of loci in the file, also find the remainder
	my $full = int($max/$groupsize);
	my $remainder = $max%$groupsize;

	# 
	for( my $i=0; $i<$full; $i++ ){
		# declare variables available to this loop
		my %temphash; # declare a temporary hash of arrays to hold loci that will be in each new file
		my $groupname = join( '.', "group", $i, "phy"); # declare a numbered groupname
		my $groupnamedir = join( '/', $grouped_locidir, $groupname );
		my $alignlength = 0;
		for( my $j=0; $j<$groupsize; $j++ ){
			my $temp = shift( @numarray );
			# construct the name of the file to be opened
			my $input = join( '.', "locus", $temp, "phy" );
			# append the path to the name of the file to be opened
			my $inputdir = join( '/', $phylipdir, $input );
			# open the single locus file for reading
			open( INFILE, $inputdir ) or die "Can't open $inputdir: $!\n";
			# read the single locus contents into hash
			while( my $line = <INFILE> ){
				if( $line =~ /^\d+\s(\d+)$/ ){
					$alignlength+=$1;
					next;
				}else{
					chomp( $line );
					my @temp = split( /\t/, $line );
					push( @{$temphash{$temp[0]}}, $temp[1] );
				}
			}
			close INFILE;
		}
		# print contents of hash to new phylip file
		open( OUTPUT, '>', $groupnamedir ) or die "can't open $groupnamedir: $!\n";
		# print header
		print OUTPUT $count, "\t", $alignlength, "\n";
		
		# print each individual to output file
		foreach my $ind( sort keys %temphash ){
			print OUTPUT $ind, "\t";
			foreach my $locus( @{$temphash{$ind}} ){
				print OUTPUT $locus;
			}
			print OUTPUT "\n";
		}
		close OUTPUT;
		
		# uncomment next line to view hash contents
		#print Dumper \%temphash;
	}
	return $full;
}

#####################################################################################################
# subroutine to run raxml
sub runraxml{
	my( $proc, $file, $raxmldir, $grouped_locidir, $boot, $model ) = @_;

	# declare variables
	my $name; # hold locus name to be fed to raxml command
	
	# get locus name
	if( $file =~ /\w+\.(\d+)\.phy/ ){
		$name = $1;
	}
	
	# join file name and path to file
	my $filepath = join( '/', $grouped_locidir, $file );

	# system call to run raxml
	system( "mpirun -np $proc -machinefile $ENV{'PBS_NODEFILE'} raxmlSSE3 -p 1234 -s $filepath -n $name -m $model -w $raxmldir -# $proc" ); # == 0 or die "ERROR: RAxML exited with non-zero status: $?";
	
	# run bootstrapping if requested
	if( $boot > 0 ){
		$name = join( '_', $name, "boot" );
		system( "mpirun -np $proc -machinefile $ENV{'PBS_NODEFILE'} raxmlSSE3 -x 1234 -p 1234 -# $boot -s $filepath -n $name -m $model -w $raxmldir" )
	}
	
}

#####################################################################################################
# subroutine to concatenate trees produced by RAxML
sub treecat{
	my( $raxmldir, $full, $astraldir, $treeoutdir ) = @_;

	# declare variables
	my @trees; # array to hold trees from RAxML output
	
	# open directory containing raxml output files
	opendir( RAX, $raxmldir ) or die "Can't open $raxmldir: $!\n";

	# read contents of directory
	my @unsorted = readdir( RAX );
	
	my @contents = sort( @unsorted );
	
	print Dumper( \@contents );

	# run raxml on each locus and store tree from each
	foreach my $file( @contents ){
		my $filedir = join( "/", $raxmldir, $file );
		if( $file =~ /^RAxML_bestTree\.\d+$/ ){
			open( TREE, $filedir ) or die "Can't open $filedir: $!\n";
			while( my $line = <TREE> ){
				chomp( $line );
				push( @trees, $line );
			}
			close( TREE );
		}
	}
	
	# write trees to new concatenated file
	open( OUT, '>', $treeoutdir ) or die "Can't open $treeout: $!\n";
	foreach my $tree( @trees ){
		print OUT $tree, "\n";
	}
	close( OUT );

	# close the directory
	closedir( RAX );	
}

#####################################################################################################
# subroutine to concatenate trees produced by RAxML
sub bspaths{
	my( $raxmldir, $astraldir, $bspath ) = @_;
	
	# open file for writing
	open( BS, '>', $bspath ) or die "Can't open $bspath: $!\n";
	
	# open RAxML directory to find bootstrap trees
	opendir( RAX, $raxmldir ) or die "Can't open $raxmldir: $!\n";
	
	# read the contents of the directory
	my @unsorted = readdir( RAX );
	
	my @contents = sort( @unsorted );
	
	print Dumper( \@contents );
	
	# find the file names that contain bootstrap trees
	foreach my $file( @contents ){
		chomp( $file );
		if( $file =~ /^RAxML_bootstrap\.\d+_boot$/ ){
			my $filedir = join( "/", $raxmldir, $file );
			print BS $filedir, "\n";
		}
	}
	# close RAxML directory
	closedir( RAX );
	
	# close file for writing
	close BS;
	
}

#####################################################################################################
# subroutine to run astral
sub runastral{
	my( $treeoutdir, $astraldir, $astralpath, $boot, $bspath, $version, $mem ) = @_;

	# declare variables
	my $astralout = "astral_output.tre"; # name of file to hold output tree from astral
	my $astraloutdir = join( '/', $astraldir, $astralout ); # concatenate output directory to astral output file name
	
	print $astralpath, "\n";
	
	# system call to run astral
if( $boot > 0 ){
		system( "java $mem -jar /home/mussmann/local/src/ASTRAL/astral.4.10.12.jar -i $treeoutdir -o $astraloutdir -b $bspath -r $boot -g" ); # == 0 or die "ERROR: astral exited with non-zero status: $?";
	}else{
		system( "java $mem -jar /home/mussmann/local/src/ASTRAL/astral.4.10.12.jar -i $treeoutdir -o $astraloutdir" ); # == 0 or die "ERROR: astral exited with non-zero status: $?";
	}
	
	return $astraloutdir;
}

#####################################################################################################
# subroutine to take astral output and split it into more appropriately named files
sub astralparse{
	my( $astraldir, $astraloutdir ) = @_;

	# declare variables
	my @trees; # array to hold trees from astral
	my $besttree = join( '/', $astraldir, "astral_best_tree.tre" );
	my $bootcon = join( '/', $astraldir, "astral_bootstrap_consensus.tre" );
	my $boottrees = join( '/', $astraldir, "astral_bootstrap_trees.tre");
	
	# read astral output into array
	open( TREES, $astraloutdir ) or die "Can't open $astraloutdir: $!\n";
	while( my $line = <TREES> ){
		chomp( $line );
		push( @trees, $line );
	}
	close TREES;
	
	# remove last item of array (astral best tree) and print to new file
	my $best = pop( @trees );
	open( BEST, '>', $besttree ) or die "Can't open $besttree: $!\n";
	print BEST $best, "\n";
	close BEST;
	
	# remove next to last item of array (astral bootstrap consensus tree) and print to new file
	my $aboot = pop( @trees );
	open( BOOT, '>', $bootcon ) or die "Can't open $bootcon: $!\n";
	print BOOT $aboot, "\n";
	close BOOT;
	
	# print remaining trees (bootstrap trees) to new file
	open( BOOTS, '>', $boottrees ) or die "Can't open $boottrees: $!\n";
	foreach my $tree( @trees ){
		print BOOTS $tree, "\n";
	}
	close BOOTS;
	
	
}
####################################################################################################u#
# subroutine to set the memory parameter to be passed to java
sub setmem{
	my( $memory ) = @_;

	my $newmem = $memory * 1000;
	my $mem = join( '', "-Xmx", $newmem, "M" );

	return $mem;

}

#####################################################################################################
