dir=$1 #input dir
bci=$2 #bcindex
ob=$3  #outer barcode sequence
pre=`basename $dir`
if  test -z "$dir" 
then
	echo "USAGE: sh pip.sh <direction_or_project>	[outer barcode seq]"
	echo example:
	echo "       sh pip.sh /data_center_07/User/lixf/RealMed/2016-4-28 TTAGGC"

	echo exit......
	exit
fi

if [ ! -d $dir/split_data ];then
	mkdir $dir/split_data
fi
if [ ! -d $dir/genotype ];then
	mkdir $dir/genotype
fi

ls $dir/raw*/*gz > $dir/reads.list
ls $dir/raw*/*fastq >>$dir/reads.list
sed -i 'N;s/gz\n/gz /' $dir/reads.list
sed -i 'N;s/fastq$\n/fastq /' $dir/reads.list
echo preparation finished
echo 
echo spliting reads.....
if [ -z "$ob" ]; then
	echo $ob is null
	#echo /usr/bin/perl /home/lixf/RNA_module_V1.0/Filter/fastx_bin/fastx_barcode_splitter_PE_v2.pl --bcfile /data_center_07/User/lixf/RealMed/bin/barcode.txt --reads $dir/reads.list --bol --mismatches 1 --partial 1 --prefix $dir/split_data --suffix fq >$dir/step1.split.sh
	echo /usr/bin/perl /home/lixf/RNA_module_V1.0/Filter/fastx_bin/fastx_barcode_splitter_PE_v3.pl --bcseq /data_center_07/User/lixf/RealMed/bin/barcode.txt --bcindex $bci --reads $dir/reads.list --bol --mismatches 1 --partial 1 --prefix $dir/split_data --suffix fq >$dir/step1.split.sh
else
	echo $ob is not null
	#echo /usr/bin/perl /home/lixf/RNA_module_V1.0/Filter/fastx_bin/fastx_barcode_splitter_PE_v2.pl --bcfile /data_center_07/User/lixf/RealMed/bin/barcode.txt --reads $dir/reads.list --obc $ob --bol --mismatches 1 --partial 1 --prefix $dir/split_data --suffix fq >$dir/step1.split.sh
	echo /usr/bin/perl /home/lixf/RNA_module_V1.0/Filter/fastx_bin/fastx_barcode_splitter_PE_v3.pl --bcseq /data_center_07/User/lixf/RealMed/bin/barcode.txt --bcindex $bci --reads $dir/reads.list --bol --mismatches 1 --partial 1 --prefix $dir/split_data --suffix fq --obc $ob >$dir/step1.split.sh
fi

#/usr/bin/perl /home/lixf/bin/tools/qsub-sge.pl --resource vf=100m --lines 1 $dir/step1.split.sh
echo spliting done

echo =========================
echo extracting KRAS genes data
ls $dir/split_data/F*-R*1.fq.gz|while read r1; 
do 

	r2=`echo $r1|sed 's/_1.fq/_2.fq/'`; 
	s=${r1##*\/}; 
	s=${s%%_*}; 
	outdir=$dir/genotype/$s
	if [ ! -d $outdir ] ;
	then
		mkdir $outdir
	fi

	if [ -s "$r1" ] && [ -s "$r2" ] ;
	then
		echo /home/lixf/Python-2.7.11/python /data_center_07/User/lixf/RealMed/bin/extract_Reads.py -f $r1 -r $r2 -o $outdir/$s -c 20; 
	fi
done >$dir/step2.extract.sh
/usr/bin/perl /home/lixf/bin/tools/qsub-sge.pl --resource vf=100m --lines 2 $dir/step2.extract.sh
echo =========================
echo call snv
find $dir/genotype -name "*_1.fastq"|while read r1;rm bad
do 
	size=`stat --printf "%s" $r1`
	#if [ $size -lt 1000000 ]; then
	#	continue
	#fi
	r2=`echo $r1|sed 's/_1.fastq/_2.fastq/'`; 
	s=`basename $r1`; 
	s1=${s%_*}
	s2=${s%%_*};

	outdir=$dir/genotype/$s2
	#if [ ! -d $outdir ] ;
	#then
	#	mkdir $outdir
	#fi
	echo /home/lixf/Python-2.7.11/python /data_center_07/User/lixf/RealMed/bin/AfterQC/after.py -1 $r1 -2 $r2 -f 0 -t 0 -q 10 -u 20 -p 20 -n 5 --barcode=false 
	read1 = $outdir/genotype/$s2/good/$s1\_1.good.fq
	read2 = $outdir/genotype/$s2/good/$s1\_2.good.fq
	#echo /home/lixf/Python-2.7.11/python /data_center_07/User/lixf/RealMed/bin/genotype_simple_v2.1beta.py -f $read1 -r $read2 -o $outdir/$s -c 20
	echo /home/lixf/Python-2.7.11/python /data_center_07/User/lixf/RealMed/bin/genotype_simple_v2.2.py -f $read1 -r $read2 -o $outdir/$s -s 1000 -c 20
done >$dir/step3.genotype.sh

/usr/bin/perl /home/lixf/bin/tools/qsub-sge.pl --resource vf=100m --lines 2 --maxjob 50  $dir/step3.genotype.sh
echo extracting and genotyping done
echo 
echo counting...
for i in $dir/genotype/*/*xls; do d=${i%\/*}; s=${d##*\/}; echo $s $i; done > $dir/$pre.genotype.list
ymd=`date "+%Y-%m-%d"`
#YaMaDei? 
/usr/bin/perl /data_center_07/User/lixf/RealMed/bin/collect_all_genotype.pl  $dir/$pre.genotype.list $dir/$pre.genotype.all.xls	$bci $ymd
echo counting done 
echo ============================All Job Finished====================================
