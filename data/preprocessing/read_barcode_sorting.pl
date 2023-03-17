#!/usr/bin/env perl -w
use strict;

my @data;
open (FILE, "< data/preprocessing/primers.txt") or die;
while (defined (my $l = <FILE>)) {
	chomp $l;
	push @data, $l;
}
close FILE;

my %barcodes;  my $j=1;
for (my $i=0; $i < @data; $i+=2) {
	my $exp = $j;
	my @d = split/\t/,$data[$i];
	my @r = split/\D+/,$d[0];
	my $replicate = join("","R",$r[1]);
	my @d2 = split/\t/,$data[$i+1];
	my $info = join("_",$replicate,uc($d[1]),$d[3]);
	my $bc_r = reverse($d2[4]);
	$bc_r =~ tr/AGCT/TCGA/;
	my $bc = join("\t",join("",$d[4],"AGTGAT"),$bc_r); #print "$bc\n";
	$barcodes{$bc} = $info;
	$j++;
}

my @k = keys(%barcodes);

foreach my $k  (keys %barcodes) {
	#foreach my $k2 (keys % {$barcodes{$k}} ) {
		#print "$k $barcodes{$k}\n";
	#}
}

my %sgRNA; my @sgRNAs;
open (FILE, "< data/preprocessing/guides.txt") or die;
while (defined (my $j= <FILE>)) { 
	chomp $j;
	my @d= split/\t/,$j;
	push @sgRNAs, $d[1];
	$sgRNA{$d[1]} = $d[1];
}

	
	
my @reads;
open (DATA, "< data/preprocessing/merged_read_file.fastq") or die;
while (defined ( my $l = <DATA> )) { 
	chomp $l;
	push @reads,$l
}
close DATA;

my %counts; 

my $left = "TGATAGAGATACTGAGCACG";
my $right = "GTTTTAGAGCTAGA";


for (my $i=0; $i < @reads; $i +=4) {
	my $guide20 = substr($reads[$i+1],38,20);
	my $guide21 = substr($reads[$i+1],38,21);
	if (exists($sgRNA{$guide20})) {
		my $lbc = substr($reads[$i+1],4,12);
		my $rbc = substr($reads[$i+1],93,12);
		my $bcs = join("\t",join("",$lbc,"AGTGAT"),$rbc);
		if (exists($barcodes{$bcs})) {
			my $exp = $barcodes{$bcs};
			$counts{$exp}{$guide20}++;
		}
	}
	if (exists($sgRNA{$guide21})) {
		my $lbc = substr($reads[$i+1],4,12);
		my $rbc = substr($reads[$i+1],94,12);
		my $bcs = join("\t",join("",$lbc,"AGTGAT"),$rbc);
		if (exists($barcodes{$bcs})) {
			my $exp = $barcodes{$bcs};
			$counts{$exp}{$guide21}++;
		}
	}
}
	


my @bc = sort {$a cmp $b} keys (%counts);
print "exp\t";
foreach (@sgRNAs) {
	print "$_\t";
}
print "\n";
foreach (@bc) {
	print "$_\t";
	for (my $i=0; $i < @sgRNAs; $i++) {
		my $num=0;
		$num = $counts{$_}{$sgRNAs[$i]} if $counts{$_}{$sgRNAs[$i]};
		print "$num\t";
		}
	print "\n";
}
	
my @k = keys(%barcodes);

foreach my $k  (keys %barcodes) {
	print "$k $barcodes{$k}\n";
}	

=cut

exit;	