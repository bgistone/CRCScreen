use strict;
use warnings;

if(@ARGV<2)
{
	print "perl $0 <genotype.list> <output>	<sampleinfo> <date>\n";
	exit;
}

open SS,$ARGV[2] or die $!;
my %sinfo=();
while(<SS>)
{
	chomp;
	my @a = split /\t/,$_;
	#$sinfo{$a[0]} = $a[1];
	$sinfo{"$a[0]-$a[1]"} = $a[2];
}
close SS;

open OO, ">$ARGV[1]" or die $!;

my @dict=('total', 'normal', 'pathogenic', '12S', '12R', '12C', '12D', '12A', '12V', '13R', '13C', '13D', '13V', '14I', 'uncountable');
my @sample=();
my %h=();
my %info=();

open LL, $ARGV[0] or die $!;
while(<LL>)
{
	chomp;
	my @a = split;
	next unless(exists $sinfo{$a[0]});
	push(@sample, $a[0]);
	open FF, $a[1] or die $!;
	while(<FF>)
	{
		chomp;
		my @b = split /\t/, $_;
		#print STDERR "$b[0]\t$a[0]\t$b[1]\t$b[2]\n";
		$h{$b[0]}{$a[0]}{'C'} = $b[1];
		$h{$b[0]}{$a[0]}{'R'} = $b[2];
		$info{$b[0]} = $b[3];
	}
	close FF;
}
close LL;

for(my $i=0;$i<@dict;$i++)
{
	print OO "#$dict[$i]:$info{$dict[$i]}\n";
}
print OO "\nID\tSample\tDate";
for(my $i=0;$i<@dict;$i++)
{
	print OO "\t$dict[$i]_count";
}
for(my $i=0;$i<@dict;$i++)
{
	print OO "\t$dict[$i]_ratio";
}
print OO "\n";
for(my $x=0;$x<@sample;$x++)
{
	print OO "$sample[$x]\t$sinfo{$sample[$x]}\t$ARGV[3]";;
	for (my $i=0;$i<@dict;$i++ )
	{
		if(exists $h{$dict[$i]}{$sample[$x]}){
			print OO "\t$h{$dict[$i]}{$sample[$x]}{'C'}";
		}else{
			print OO "\t0";
		}
	}
	for(my $j=0;$j<@dict;$j++)
	{
		if(exists $h{$dict[$j]}{$sample[$x]}){
			print OO "\t$h{$dict[$j]}{$sample[$x]}{'R'}%";
		}else{
			print OO "\t0.00%";
		}
	}
	print OO "\n";
}
