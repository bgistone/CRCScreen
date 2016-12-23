$s=shift;
$s= reverse $s;
$s=~tr/ATCG/TAGC/;
print "$s\n";
