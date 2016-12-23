use strict;
use warnings;
my $inputs=shift;

my %hash = ();
#my @samples = ();
#my @types = ();
my %types = ();
my %samples = ();
open IN,$inputs or die $!;
while(<IN>)
{
	chomp;
	my @a = split /\s+/,$_;
	#sample	type	path
	&load($a[0],$a[1],$a[2]);
}
close IN;
print "Chr\tLocation";
#for(my $i=0;$i<@types;$i++)
foreach my $t(sort{$a<=>$b} keys %types)
{
	#for(my $j=0;$j<@samples;$j++)
	foreach my $s(sort{$a<=>$b} keys %samples)
	{
		#print "\t$types[$i]\_$samples[$j]\_M\t$types[$i]\_$samples[$j]\_UM\t$types[$i]\_$samples[$j]\_Context";
		print "\t$t\_$s\_M\t$t\_$s\_UM\t$t\_$s\_Context";
	}
}
print "\n";

foreach my $loc(keys %hash)
{
	print "$loc";
	#for my $t(@types)
	foreach my $t(sort{$a<=>$b} keys %samples)
	{
		my %hash_tmp_a = ();
		if(exists $hash{$loc}{$t}){
		
			%hash_tmp_a = %{$hash{$loc}{$t}};
		}else{
			foreach my $s(@samples)
			{
				print "\t0\t0\t0.00";
			}
			next;
		}
		#for my $s(@samples)
		foreach my $s(sort{$a<=>$b} keys %samples)
		{
			my %hash_tmp_b = ();
			if(exists $hash_tmp_a{$s})
			{
				%hash_tmp_b = %{$hash_tmp_a{$s}};
				my $methylated = 0;
				my $unmethylated = 0;
				if (exists $hash_tmp_b{"M"})
				{
					$methylated = $hash_tmp_b{"M"};
				}
				if(exists $hash_tmp_b{"UM"}){
					$unmethylated = $hash_tmp_b{"UM"};
				}

				my $context = 0;
				if($methylated <= 100 && $unmethylated <= 100)
				{
					$context = 0.00;
				}else{
					$context = $methylated / ($methylated + $unmethylated) * 100;
				}
				print "\t$methylated\t$unmethylated\t$context";
			}else{
				print "\t0\t0\t0.00";
			}
		}
	}
	print "\n";
}

sub load{
	my ($sample,$type,$path) = @_;
	#push(@types,$type);
	#push(@samples,$sample);
	$types{$type} = 1;
	$samples{$sample} = 1;
	print STDERR "loading $path\n";
	open PH,$path or die $!;
	<PH>;
	while(<PH>)
	{
		my @t = split;
		my $sign = "M";
		if($t[1] eq "-")
		{
			$sign = "UM";
		}
		$hash{"$t[2]\t$t[3]"}{$type}{$sample}{$sign}++;
	}
	close PH;
}
