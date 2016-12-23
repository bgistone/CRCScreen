dir=$1 #input dir
bci=$2 #bcindex
#ob=$3  #outer barcode sequence
pre=`basename $dir`
if  test -z "$dir" 
then
	echo "USAGE: sh pip.sh <direction_or_project>	<barcode>"
	echo example:
	echo "       sh pip.sh /data_center_07/User/lixf/RealMed/2016-4-28 NEB3.bc.txt"

	echo exit......
	exit
fi

if [ ! -d $dir/split_data ];then
	mkdir $dir/split_data
fi
if [ ! -d $dir/QC ];then
	mkdir $dir/QC
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
echo QC
ls $dir/split_data/F*-R*1.fq.gz|while read r1; 
do 

	r2=`echo $r1|sed 's/_1.fq/_2.fq/'`; 
	s=${r1##*\/}; 
	s=${s%%_*}; 
	outdir=$dir/QC/$s
	if [ ! -d $outdir ] ;
	then
		mkdir $outdir
	fi

	if [ -s "$r1" ] && [ -s "$r2" ] ;
	then
		echo /home/lixf/Python-2.7.11/python /data_center_07/User/lixf/RealMed/bin/AfterQC/after.py -1 $r1 -2 $r2 -f 10 -t 10 -q 10 -u 20 -p 20 -n 5 --barcode=false -g $outdir; 
	fi
done >$dir/step2.QC.sh
/usr/bin/perl /home/lixf/bin/tools/qsub-sge.pl --resource vf=100m --lines 2 $dir/step2.QC.sh
echo =========================
echo call methylation
find $dir/QC -name "*_1.fq"|while read r1;
do 
	size=`stat --printf "%s" $r1`
	#if [ $size -lt 1000000 ]; then
	#	continue
	#fi
	r2=`echo $r1|sed 's/_1.fq/_2.fq/'`; 
	s=`basename $r1`; 
	s=${s%%_*};

	outdir=$dir/methylation/$s
	if [ ! -d $outdir ] ;
	then
		mkdir $outdir
	fi
	echo /data_center_07/User/lixf/RealMed/bin/RRBS/bsmap-2.90/bsmap -a $r1 -b $r2 -d /data_center_07/User/lixf/RealMed/bin/genome/ucsc.hg19.fasta -v 10 -s 12 -q 25 -o $outdir/$s.sam
	echo python /data_center_07/User/lixf/RealMed/bin/RRBS/bsmap-2.90/methratio.py -o $outdir/$s.methylation.xls -d /data_center_07/User/lixf/RealMed/bin/genome/ucsc.hg19.fasta $outdir/$s.sam
done >$dir/step3.methylation.sh

/usr/bin/perl /home/lixf/bin/tools/qsub-sge.pl --resource vf=1G --lines 2 --maxjob 50  $dir/step3.methylation.sh
echo methylaton analysis done
echo 
echo counting...
for i in $dir/methylation/*/*xls; do d=${i%\/*}; s=${d##*\/}; echo $s $i; done > $dir/$pre.methylation.list
ymd=`date "+%Y-%m-%d"`
perl /home/miyq/code/methyratio_summary.pl $dir/$pre.methylation.list >$dir/$pre.methylation.xls
echo counting done 
echo ============================All Job Finished====================================
