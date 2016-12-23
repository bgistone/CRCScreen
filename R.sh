less /data_center_07/User/lixf/RealMed/2016-8-18-p039/genotype/F2-F2/F2-F2_filter_1.fastq|grep "^@"|cut -d ':' -f 10|perl -lane '$h{$F[0]}+=1;END{foreach my $k(keys %h){print "$k\t$h{$k}";}}'>bc.stat
