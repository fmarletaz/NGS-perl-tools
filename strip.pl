#! /usr/bin/perl
use warnings;
use strict;

die "Usage: strip.pl <fq_1> <fq_2>\n" unless @ARGV==2;

my $file1=shift @ARGV;
my $file2=shift @ARGV;

#cutoff is number of reads to be checked to find the missing pair 
my $cutoff=50;
#open(my $fq1,"<:gzip",$file1) or die $!;
#open(my $fq2,"<:gzip",$file2) or die $!;

open(my $fq1,"<",$file1) or die $!;
open(my $fq2,"<",$file2) or die $!;

(my $name1)=($file1=~/(.+)\.\w+$/);
(my $name2)=($file2=~/(.+)\.\w+$/);

open(my $pair1,">","$name1\.pe.fq");
open(my $pair2,">","$name2\.pe.fq");
open(my $single1,">","$name1\.se.fq");
open(my $single2,">","$name2\.se.fq");

my $old_eol=$/;
local $/="@";
my $ct=0; my $line=0; 
while(!eof($fq1) and !eof($fq2)) {
	my $read1=<$fq1>;
	my $read2=<$fq2>;
	
	chomp $read1;
	chomp $read2;
	
	if($read1=~/^\s*$/) { next; }
	if($read2=~/^\s*$/) { next; }
	
	#print "£ $read1\n";
	#print "£ $read2\n";
	
	$line++;
	
	
	(my $id1)=($read1=~/^([^\#]+)\#/);
	(my $id2)=($read2=~/^([^\#]+)\#/);
	
	die "problem with $id1!!\n\n\>$read2\<\n\n" unless(defined($id2));
	

	if($line%500000==0) {
		
		print ($line/1000000)"M\:\t$id1\t$id2 \n"; 
	}
	#if((tell $fq1)%1==0) { print "."; }

	if($id1 eq $id2) { 
		print $pair1 "@".$read1;
		print $pair2 "@".$read2;
	}
	
	
	else {
		my $pos1=tell $fq1;
		my $pos2=tell $fq2;
		my $nbm=0; my $pf='0';
		my @reads1=$read1; my @reads2=$read2;
		while($pf eq '0') {
			$nbm++;
			
			my $nread1=<$fq1>;
			my $nread2=<$fq2>;

			chomp $nread1;
			chomp $nread2;
			
			if($nread1=~/^\s*$/) { next; }
			if($nread2=~/^\s*$/) { next; }

			
			(my $nid1)=($nread1=~/^([^\#]+)\#/);
			(my $nid2)=($nread2=~/^([^\#]+)\#/);
			
			#print "N1/ $nid1\nN2/ $nid2\n\n";
			
			if($nid1 eq $id2) {
				print $single1 join("@",@reads1);
				print $pair2 "@".$read2;
				print $pair1 "@".$nread1;
				seek $fq2,$pos2,0;
				$pf='1';
				$ct+=scalar @reads1
			}
			elsif($nid2 eq $id1) {
				print $single2 join("@",@reads2);
				print $pair1 "@".$read1;
				print $pair2 "@".$nread2;
				seek $fq1,$pos1,0;
				$pf='1';
				$ct+=scalar @reads2
			}
			elsif($nbm==$cutoff) {
				print $single1 "@".$read1;
				print $single2 "@".$read2;
				seek $fq1,($pos1+1),0;
				seek $fq2,($pos2+1),0;
				$pf='1';
				$ct+=2;
				print "Case of cross-deletions, moving forward...\n";
				print "Possibly involved reads: $id1 , $id2\n";
			}
			push @reads1,$nread1;
			push @reads2,$nread2;
		}
	}
}

print "$ct singletons!\n";



