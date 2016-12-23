use strict;
use warnings;

my $r1=shift;
my $r2=shift;
my $id=shift;
my $pre=shift;

my %hash;

open ID, $id or die $!;
while(<ID>)
{
	chomp;
	my @a = split;
	$hash{$a[0]} = "";
}
close ID;

if($r1=~/.gz$/){
	open RA, "gzip -dc $r1|" or die $!;
	open RB, "gzip -dc $r2|" or die $!;
}else{
	open RA, $r1 or die $!;
	open RB, $r2 or die $!;
}
open OA, ">$pre\_r1.fastq" or die $!;
open OB, ">$pre\_r2.fastq" or die $!;
while(<RA>)
{
	my @a = split /\s+/, $_;
	my ($id, $seq, $sign, $qual);
	if(exists $hash{$a[0]})
	{
		$id = $_;$seq = <RA>;$sign=<RA>;$qual=<RA>;
		print OA $id, $seq, $sign, $qual;
		$id = <RB>;$seq = <RB>;$sign=<RB>;$qual=<RB>;
		print OB $id, $seq, $sign, $qual;
	}else{
		<RA>,  <RA>,  <RA>;
		<RB>,  <RB>,  <RB>,  <RB>;
	}
}
close RA;close RB;close OA;close OB;
