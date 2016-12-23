#!/usr/bin/perl
use strict;
use warnings;
#use IO::Handle;
use IO::File;
use IO::Zlib;
use PerlIO::gzip;
use Data::Dumper;
use Getopt::Long;
use Carp;

##
## This program splits a FASTQ/FASTA file into several smaller files,
## Based on barcode matching.
##
## run with "--help" for usage information
##
## Assaf Gordon <gordon@cshl.edu> , 11sep2008

# Forward declarations
sub load_barcode_seq ($$);
sub parse_command_line ;
sub match_sequences ;
sub mismatch_count($$) ;
sub print_results;
sub open_and_detect_input_format ($);
sub read_record;
sub write_record($$);
sub usage();

# Global flags and arguments, 
# Set by command line argumens
my $barcode_seq;
my $barcode_index;
my $reads_file;
my %hash_bc;
my $barcodes_at_eol = 0 ;
my $barcodes_at_bol = 0 ;
my $exact_match = 0 ;
my $allow_partial_overlap = 0;
my $allowed_mismatches = 1;
my $newfile_suffix = '';
my $newfile_prefix  ;
my $quiet = 0 ;
my $debug = 0 ;
my $fastq_format = 1;
my $RealBio_barcode="ATCTCGTA";
my $outer_barcode_mismatch=2;

# Global variables 
# Populated by 'create_output_files'
my %filenames;
my %files_handle;
my %reads;
my %counts = ( 'unmatched-unmatched' => 0 );
my $barcodes_length;
my @barcodes;
my $fh_f;
my $fh_r;


# The Four lines per record in FASTQ format.
# (when using FASTA format, only the first two are used)
my $seq_namef;
my $seq_namer;
my $seq_basesf;
my $seq_basesr;
my $seq_namef2;
my $seq_namer2;
my $seq_qualitiesf;
my $seq_qualitiesr;


#
# Start of Program
#
parse_command_line ;
print STDERR "loading barcodes..." if($debug);
load_barcode_seq ( $barcode_seq, $barcode_index) ;
print STDERR "done\nread raw_data list:$reads_file\n" if($debug);
open_and_detect_input_format( $reads_file );

#match_sequences ;

print_results unless $quiet;

foreach my $k(keys %files_handle)
{
	my ($fh, $rh) = @{$files_handle{$k}};
	$fh->close if $fh->opened;
	$rh->close if $rh->opened;
}
#
# End of program
#

sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions (
		"bcseq=s" => \$barcode_seq,
		"bcindex=s" => \$barcode_index,
		"obc=s" => \$RealBio_barcode,
		"obcm=i"=> \$outer_barcode_mismatch,
		"reads=s" => \$reads_file, 
		"eol"  => \$barcodes_at_eol,
		"bol"  => \$barcodes_at_bol,
		"exact" => \$exact_match,
		"prefix=s" => \$newfile_prefix,
		"suffix=s" => \$newfile_suffix,
		"quiet" => \$quiet, 
		"partial=i" => \$allow_partial_overlap,
		"debug" => \$debug,
		"mismatches=i" => \$allowed_mismatches,
		"help" => \$help
	) ;

	usage() if ($help);

	die "Error: barcode file not specified (use '--bcindex [FILENAME]')\n" unless defined $barcode_seq;
	die "Error: reads file list not specified (use '--reads [FILENAME]')\n" unless defined $reads_file;
	die "Error: prefix path/filename not specified (use '--prefix [PATH]')\n" unless defined $newfile_prefix;

	if ($barcodes_at_bol == $barcodes_at_eol) {
		die "Error: can't specify both --eol & --bol\n" if $barcodes_at_eol;
		die "Error: must specify either --eol or --bol\n" ;
	}

	die "Error: invalid for value partial matches (valid values are 0 or greater)\n" if $allow_partial_overlap<0;

	$allowed_mismatches = 0 if $exact_match;

	die "Error: invalid value for mismatches (valid values are 0 or more)\n" if ($allowed_mismatches<0);

	die "Error: partial overlap value ($allow_partial_overlap) bigger than " . 
	"max. allowed mismatches ($allowed_mismatches)\n" if ($allow_partial_overlap > $allowed_mismatches);


	exit unless $result;
}



#
# Read the barcode file
#
sub load_barcode_seq ($$) {
	my $bcseq = shift or croak "Missing barcode seq file name";
	my $bcindex;#bcindex
	if($barcode_index ne "")
	{
		$bcindex=shift or croak "Missing barcode-sample information file";
	}

	open BCFILE,"<$bcseq" or die "Error: failed to open barcode file ($bcseq)\n";
	while (<BCFILE>) {
		next if m/^#/;
		chomp;
		my ($ident, $barcode) = split ;

		$barcode = uc($barcode);

		# Sanity checks on the barcodes
		die "Error: bad data at barcode file ($bcseq) line $.\n" unless defined $barcode;
		die "Error: bad barcode value ($barcode) at barcode file ($bcseq) line $.\n"
		unless $barcode =~ m/^[AGCT]+$/;

		die "Error: bad identifier value ($ident) at barcode file ($bcseq) line $. (must be alphanumeric)\n" 
		unless $ident =~ m/^\w+$/;

		die "Error: badcode($ident, $barcode) is shorter or equal to maximum number of " .
		"mismatches ($allowed_mismatches). This makes no sense. Specify fewer  mismatches.\n" 
		if length($barcode)<=$allowed_mismatches;

		$barcodes_length = length($barcode) unless defined $barcodes_length;
		die "Error: found barcodes in different lengths. this feature is not supported yet.\n" 
		unless $barcodes_length == length($barcode);
		my $ident_rev = $ident;
		push @barcodes, [$ident, $barcode];
		$hash_bc{$barcode} = $ident;

		if ($allow_partial_overlap>0) {
			foreach my $i (1 .. $allow_partial_overlap) {
				substr $barcode, ($barcodes_at_bol)?0:-1, 1, '';
				push @barcodes, [$ident, $barcode];
				$hash_bc{$barcode} = $ident;
			}
		}
	}
	close BCFILE;
	if($bcindex ne ""){
		&create_output_files_new($bcindex);
	}else{
		&create_output_files;
	}

	if ($debug) {
		print STDERR "barcode\tsequence\n";
		foreach my $barcoderef (@barcodes) {
			my ($ident, $seq) = @{$barcoderef};
			print STDERR $ident,"\t", $seq ,"\n";
		}
	}
}

# Create one output file for each barcode combination.
sub create_output_files_new{
	my $bcidx = shift;
	my ($new_filenamef,$new_filenamer, $single);
	open BCI,$bcidx or die $!;
	while(<BCI>) 
	{
		chomp;
		#forward	reverse	sample
		#F24     R24     A5
		my @a = split /\t/, $_;
		$new_filenamef = $newfile_prefix . "/". $a[0]."-" . $a[1] . "_1." . $newfile_suffix . ".gz"; 
		$new_filenamer = $newfile_prefix . "/". $a[0]."-" . $a[1] . "_2." . $newfile_suffix . ".gz"; 
		print STDERR "\t......create output file: $new_filenamef and $new_filenamer\n" if $debug;
		push(@{$filenames{"$a[0]-$a[1]"}}, $new_filenamef);
		push(@{$filenames{"$a[0]-$a[1]"}}, $new_filenamer);
		open my $filef, ">:gzip", "$new_filenamef" or die "Error: failed to create output file ($new_filenamef)\n"; 
		open my $filer, ">:gzip", "$new_filenamer" or die "Error: failed to create output file ($new_filenamer)\n"; 
		print STDERR "handle:$a[0]-$a[1]\n";
		push(@{$files_handle{"$a[0]-$a[1]"}}, $filef);
		push(@{$files_handle{"$a[0]-$a[1]"}}, $filer);
	}
	close BCI;
}

# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub create_output_files {
	my %barcodes = map { $_->[0] => 1 } @barcodes; #generate a uniq list of barcode identifiers;
	#$barcodes{'unmatched'} = 1 ;
	foreach my $ident (keys %barcodes) 
	{
		my ($new_filenamef,$new_filenamer, $single);
		foreach my $ident2 (keys %barcodes) 
		{
			$ident2=~s/F/R/;
			$new_filenamef = $newfile_prefix . "/". $ident."-" . $ident2 . "_1." . $newfile_suffix . ".gz"; 
			$new_filenamer = $newfile_prefix . "/". $ident."-" . $ident2 . "_2." . $newfile_suffix . ".gz"; 
			print STDERR "\t......create output file: $new_filenamef and $new_filenamer\n" if $debug;
			push(@{$filenames{"$ident-$ident2"}}, $new_filenamef);
			push(@{$filenames{"$ident-$ident2"}}, $new_filenamer);
			open my $filef, ">:gzip", "$new_filenamef" or die "Error: failed to create output file ($new_filenamef)\n"; 
			open my $filer, ">:gzip", "$new_filenamer" or die "Error: failed to create output file ($new_filenamer)\n"; 
			#my $filef = IO::Zlib->new("$new_filenamef", "wb9") or die "Error: failed to create output file ($new_filenamef)\n";
			#my $filer = IO::Zlib->new("$new_filenamer", "wb9") or die "Error: failed to create output file ($new_filenamer)\n";
			print STDERR "handle:$ident-$ident2\n";
			push(@{$files_handle{"$ident-$ident2"}}, $filef);
			push(@{$files_handle{"$ident-$ident2"}}, $filer);
		}

	}
}

sub match_sequences {

	#reset counters
	foreach my $ident ( keys %filenames ) {
		$counts{$ident} = 0;
	}

	# Read file FASTQ file
	# split accotding to barcodes
	my $n = 0;
	my $total = 0;
	while ( read_record ) {
		$n++;
		if($RealBio_barcode ne ""){
			my $outer_barcode=(split /:/,$seq_namef)[-1];
			chomp $outer_barcode;
			my $c = length($RealBio_barcode)-(($RealBio_barcode^$outer_barcode)=~tr/\0/\0/);
			if($c>$outer_barcode_mismatch)
			{
				next;
			}
		}
		chomp $seq_basesf;chomp $seq_basesr;

		#print STDERR "sequence $seq_basesf: \n" if $debug;
		#print STDERR "sequence $seq_basesr: \n" if $debug;

		my ($barcode_readf, $barcode_readr);
		my ($barcode_read_mismatchf,$barcode_read_mismatchr);
		my $best_barcode_mismatches_count_f = $barcodes_length;
		my $best_barcode_ident_f = undef;
		my $best_barcode_mismatches_count_r = $barcodes_length;
		my $best_barcode_ident_r = undef;

		if ($barcodes_at_bol) {
			$barcode_readf = substr $seq_basesf, 0, $barcodes_length;
			$barcode_readr = substr $seq_basesr, 0, $barcodes_length;
		} else {
			$barcode_readf = substr $seq_basesf, - $barcodes_length;
			$barcode_readr = substr $seq_basesr, - $barcodes_length;
		}

		#quick mode or no mismatch allowed
		if(exists $hash_bc{$barcode_readf})
		{
			$best_barcode_ident_f = $hash_bc{$barcode_readf};
			$best_barcode_mismatches_count_f = 0;
			#print STDERR "E\t$n\t$barcode_readf\t$best_barcode_ident_f\n";
		}elsif($allowed_mismatches>0){
			#for mismatch detecting function
			#$mm_f,$mm_r,$min_mm_f,$min_mm_r,$bc_f,$bc_r
			&Try_all_barcodes(\$barcode_readf,\$barcode_read_mismatchf,\$best_barcode_mismatches_count_f,\$best_barcode_ident_f);
			#print STDERR "M\t$n\t$barcode_readf\t$best_barcode_ident_f\n";
		}
		if(exists $hash_bc{$barcode_readr})
		{
			$best_barcode_ident_r = $hash_bc{$barcode_readr};
			$best_barcode_mismatches_count_r = 0;
			$best_barcode_ident_r =~ s/F/R/;
			#print STDERR "E\t$n\t$barcode_readr\t$best_barcode_ident_r\n";
		}elsif($allowed_mismatches>0){
			#for mismatch detecting function
			&Try_all_barcodes(\$barcode_readr,\$barcode_read_mismatchr,\$best_barcode_mismatches_count_r,\$best_barcode_ident_r);
			$best_barcode_ident_r =~ s/F/R/;
			#print STDERR "E\t$n\t$barcode_readr\t$best_barcode_ident_r\n";
		}
		#$best_barcode_ident_f = 'unmatched' if ( (!defined $best_barcode_ident_f) || $best_barcode_mismatches_count_f > $allowed_mismatches) ;
		#$best_barcode_ident_r = 'unmatched' if ( (!defined $best_barcode_ident_r) || $best_barcode_mismatches_count_r > $allowed_mismatches) ;
		next if ( (!defined $best_barcode_ident_f) || $best_barcode_mismatches_count_f > $allowed_mismatches) ;
		next if ( (!defined $best_barcode_ident_r) || $best_barcode_mismatches_count_r > $allowed_mismatches) ;

		#print STDERR "sequence $seq_basesf matched barcode: $best_barcode_ident_f\n" if $debug;
		#print STDERR "sequence $seq_basesr matched barcode: $best_barcode_ident_r\n" if $debug;
		my $key = "$best_barcode_ident_f-$best_barcode_ident_r";
		if(exists $files_handle{$key})
		{
			$counts{$key}++;
			#get the file associated with the matched barcode.
			#(note: there's also a file associated with 'unmatched' barcode)
			push(@{$reads{$key}{'r1'}}, $seq_namef);
			push(@{$reads{$key}{'r1'}}, "$seq_basesf\n");
			push(@{$reads{$key}{'r1'}}, $seq_namef2);
			push(@{$reads{$key}{'r1'}}, $seq_qualitiesf);
			push(@{$reads{$key}{'r2'}}, $seq_namer);
			push(@{$reads{$key}{'r2'}}, "$seq_basesr\n");
			push(@{$reads{$key}{'r2'}}, $seq_namer2);
			push(@{$reads{$key}{'r2'}}, $seq_qualitiesr);
		}
		if($n>=100000)
		{
			foreach my $k(keys %reads)
			{
				write_record2($k);
			}
			%reads = ();
			$total += $n;
			$n = 0;
			print STDERR  "\t......split $total reads done......\n" if $debug;
		}

	}#while end

	foreach my $k(keys %reads)
	{
		write_record2($k);
	}
	%reads = ();
	$total += $n;
	$n = 0;
	print STDERR  "\t......split $total reads done......\n" if $debug;
}

#Try all barcodes, find the one with the lowest mismatch count
sub Try_all_barcodes{
	my ($seq,$mm,$min_mm,$bc)=@_;

	foreach my $barcoderef (@barcodes) {
		my ($ident, $barcode) = @{$barcoderef};

		# Get DNA fragment (in the length of the barcodes)
		# The barcode will be tested only against this fragment
		# (no point in testing the barcode against the whole sequence)

		$$mm = mismatch_count($$seq, $barcode); 

		# if this is a partial match, add the non-overlap as a mismatch
		# (partial barcodes are shorter than the length of the original barcodes)
		$$mm += ($barcodes_length - length($barcode)); 

		if ( $$mm < $$min_mm ) {
			$$min_mm = $$mm ;
			$$bc = $ident ;
		}
	}
}
#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }



sub print_results
{
	print "Barcode\tCount\tLocation\n";
	my $total = 0 ;
	foreach my $ident (sort keys %counts) {
		print $ident, "\t", $counts{$ident},"\t",join("\t", @{$filenames{$ident}}),"\n";
		$total += $counts{$ident};
	}
	print "total\t",$total,"\n";
}


sub read_record
{
	$seq_namef = $fh_f->getline();
	$seq_namer = $fh_r->getline();
	#$seq_namef = <$fh_f>;
	#$seq_namer = <$fh_r>;
	return undef unless defined $seq_namef; # End of file?
	return undef unless defined $seq_namer; # End of file?

	$seq_basesf = $fh_f->getline();
	$seq_basesr = $fh_r->getline();
	#$seq_basesf = <$fh_f>;
	#$seq_basesr = <$fh_r>;
	die "Error: bad input file, expecting line with sequences\n" unless defined $seq_basesf;
	die "Error: bad input file, expecting line with sequences\n" unless defined $seq_basesr;

	# If using FASTQ format, read two more lines
	if ($fastq_format) {
		$seq_namef2 = $fh_f->getline();
		$seq_namer2 = $fh_r->getline();
		#$seq_namef2 = <$fh_f>;
		#$seq_namer2 = <$fh_r>;
		die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_namef2;
		die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_namer2;

		$seq_qualitiesf = $fh_f->getline();
		$seq_qualitiesr = $fh_r->getline();
		#$seq_qualitiesf = <$fh_f>;
		#$seq_qualitiesr = <$fh_r>;
		die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualitiesf;
		die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualitiesr;
	}
	return 1;
}

sub write_record($$)
{
	my ($file_f, $file_r) = @_;

	croak "Bad file handle" unless defined $file_f;
	croak "Bad file handle" unless defined $file_r;

	print $file_f $seq_namef;
	print $file_f $seq_basesf,"\n";
	print $file_r $seq_namer;
	print $file_r $seq_basesr,"\n";

	#if using FASTQ format, write two more lines
	if ($fastq_format) {
		print $file_f $seq_namef2;
		print $file_f $seq_qualitiesf;
		print $file_r $seq_namer2;
		print $file_r $seq_qualitiesr;
	}
}

sub write_record2($)
{
	my $k = shift;
	my ($file_f, $file_r) = @{$files_handle{$k}};
	#print STDERR "write to $k\t$file_f\t$file_r\n";
	croak "Bad file handle" unless defined $file_f;
	croak "Bad file handle" unless defined $file_r;
	my $read1 = join("", @{$reads{$k}{'r1'}});
	my $read2 = join("", @{$reads{$k}{'r2'}});
	print $file_f $read1;
	print $file_r $read2;
}

sub open_and_detect_input_format ($)
{
	my $list = shift;

#	$input_file_io  = new IO::Handle;
#	die "Failed to open STDIN " unless $input_file_io->fdopen(fileno(STDIN),"r");
	#
#	# Get the first characeter, and push it back
#	my $first_char = $input_file_io->getc();
#	$input_file_io->ungetc(ord $first_char);
	#
#	if ($first_char eq '>') {
#		# FASTA format
#		$fastq_format = 0 ;
#		print STDERR "Detected FASTA format\n" if $debug;
#	} elsif ($first_char eq '@') {
#		# FASTQ format
#		$fastq_format = 1;
#		print STDERR "Detected FASTQ format\n" if $debug;
#	} else {
#		die "Error: unknown file format. First character = '$first_char' (expecting > or \@)\n";
#	}
	open LL, $list or die $!;
	while(<LL>)
	{
		chomp;
		my @a =split /[:;,\s+\t]/, $_;
		if($a[0] =~ /.gz$/)
		{
			$fh_f = IO::Zlib->new($a[0], "rb");
			$fh_r = IO::Zlib->new($a[1], "rb");
		}else{
			$fh_f = IO::File->new($a[0], "r");
			$fh_r = IO::File->new($a[1], "r");
		}
		print STDERR "\topen: $a[0]\n\topen: $a[1]\n";
		$fastq_format = 1;
		match_sequences;
		$fh_f->close;
		$fh_r->close;
	}
}

sub usage()
{
	print<<EOF;

This program reads FASTA/FASTQ file and splits it into several smaller files,
Based on barcode matching.
FASTA/FASTQ data is read from STDIN (format is auto-detected.)
Output files will be writen to disk.
Summary will be printed to STDOUT.

usage: $0 --bcseq FILE --bcindex FILE --reads FILE --prefix PREFIX [--suffix SUFFIX] [--bol|--eol] 
		 [--mismatches N] [--exact] [--partial N] [--help] [--quiet] [--debug]

Arguments:

--bcseq	FILE	- File of barcodes name. (see explanation below.)
--bcindex FILE - File of forward-reverse barcodes combination information
--obc	STR	- outer barcode sequence.(defaut: ATCTCGTA)
--obcm	INT	- max number of mismatches allowed in outer barcode.(default:2)
--reads FILE	- reads list name. (see explantation below.)
--prefix PREFIX	- File prefix. will be added to the output files. Can be used
		  to specify output directories.
--suffix SUFFIX	- File suffix (optional). Can be used to specify file
		  extensions.
--bol		- Try to match barcodes at the BEGINNING of sequences.
		  (What biologists would call the 5' end, and programmers
		  would call index 0.)
--eol		- Try to match barcodes at the END of sequences.
		  (What biologists would call the 3' end, and programmers
		  would call the end of the string.)
		  NOTE: one of --bol, --eol must be specified, but not both.
--mismatches N	- Max. number of mismatches allowed. default is 1.
--exact		- Same as '--mismatches 0'. If both --exact and --mismatches 
		  are specified, '--exact' takes precedence.
--partial N	- Allow partial overlap of barcodes. (see explanation below.)
		  (Default is not partial matching)
--quiet		- Don't print counts and summary at the end of the run.
		  (Default is to print.)
--debug		- Print lots of useless debug information to STDERR.
--help		- This helpful help screen.

Example (Assuming 'reads.list' is a FASTQ file, 'barcodes.txt' is 
the barcodes index file):

   \$perl $0 --bcseq barcodes.txt --bcindex sample.barcode.txt --reads reads.list --bol --mismatches 2 \\
	--prefix /tmp/bla_ --suffix ".fq"

Readslist file format
---------------------
Reads list file is a simple text file. Each line should contain Read1 and Read2 as follow format:
fq1	fq2
or, 
fq1 fq2
or, 
fq1,fq2
or
fq1:fq2
or,
fq1;fq2

Barcode index file format
-------------------
Barcode files are simple text files. Each line should contain an identifier 
(descriptive name for the barcode), and the barcode itself (A/C/G/T), 
separated by a TAB character. Example:

	#This line is a comment (starts with a 'number' sign)
	BC1 GATCT
	BC2 ATCGT
	BC3 GTGAT
	BC4 TGTCT

For each barcode, a new FASTQ file will be created (with the barcode's 
identifier as part of the file name). Sequences matching the barcode 
will be stored in the appropriate file.

Running the above example (assuming "barcodes.txt" contains the above 
barcodes), will create the following files:
	/tmp/bla_BC1.txt
	/tmp/bla_BC2.txt
	/tmp/bla_BC3.txt
	/tmp/bla_BC4.txt
	/tmp/bla_unmatched.txt
The 'unmatched' file will contain all sequences that didn't match any barcode.

bcfile
----------------
**Barcode information for sample marking 
--bcfile sample.bc.txt
>cat sample.bc.txt
>F23 R28 熊大
>F23 R26 熊二
>F23 R24 光头强
>F22 R22 BigDog

Barcode matching
----------------

** Without partial matching:

Count mismatches between the FASTA/Q sequences and the barcodes.
The barcode which matched with the lowest mismatches count (providing the
count is small or equal to '--mismatches N') 'gets' the sequences.

Example (using the above barcodes):
Input Sequence:
	GATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG

Matching with '--bol --mismatches 1':
   GATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
   GATCT (1 mismatch, BC1)
   ATCGT (4 mismatches, BC2)
   GTGAT (3 mismatches, BC3)
   TGTCT (3 mismatches, BC4)

This sequence will be classified as 'BC1' (it has the lowest mismatch count).
If '--exact' or '--mismatches 0' were specified, this sequence would be 
classified as 'unmatched' (because, although BC1 had the lowest mismatch count,
it is above the maximum allowed mismatches).

Matching with '--eol' (end of line) does the same, but from the other side
of the sequence.

** With partial matching (very similar to indels):

Same as above, with the following addition: barcodes are also checked for
partial overlap (number of allowed non-overlapping bases is '--partial N').

Example:
Input sequence is ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
(Same as above, but note the missing 'G' at the beginning.)

Matching (without partial overlapping) against BC1 yields 4 mismatches:
   ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
   GATCT (4 mismatches)

Partial overlapping would also try the following match:
   -ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
   GATCT (1 mismatch)

Note: scoring counts a missing base as a mismatch, so the final
mismatch count is 2 (1 'real' mismatch, 1 'missing base' mismatch).
If running with '--mismatches 2' (meaning allowing upto 2 mismatches) - this 
seqeunce will be classified as BC1.

EOF

	exit 1;
}
