#!/usr/bin/perl -w

use strict;
use warnings;

#NOTES
#PURPOSE: Assist Dr. Renthal's lab in analyizing data output from AlphaFold and identify whuch amino acids are fully exposed to solvent.
#Program was written by Alvaro Ortola Tortosa to aid the research of Dr. Renthal. Full citation below
#ORCID: https://orcid.org/0009-0007-4907-5512

#Initial arguments
my $in = $ARGV[0];
my $areafile = $ARGV[1];


#Help sections
if ($#ARGV !=1){print "ERROR: Two files are necessary, use perl AminoArranger.pl -h for help.\n" and exit;}

if (($in or $areafile) eq -h){
print "\nUsage: perl AminoArranger.pl <CERS.txt> <CERS.area.txt>
Ensure that the CERS file is input first and that the area file is second PROGRAM WILL FAIL OTHERWISE
Current program only works with one input pair at a time, rename any files generated after generation
so they aren't overwritten.\n" and exit;}

#Opening necessary files
open(my $IF, '<',$in) or die "$in cannot be opened.\n";
open(my $BM1,'>', <BenchMark1.txt>); #Extracts all relevant data from file 1
open(my $BM2,'>', <BenchMark2.txt>); #Extracts area data and adds it up
open(my $BM3,'>', <BenchMark3.txt>); #Does the division step

#This obtains any lines that start with ATOM or TER and prints into a new file

print "Extracting information from $in.\n";
#Step 1: Extract any lines that begin with ATOM or TERM and place into new file
my $line = <$IF>; 
foreach $line (<$IF>){
	my @line = split /\s+/,$line;
	if ($line[0] eq "ATOM"){
		print $BM1 $line;
}elsif ($line[0] eq "TER"){
	print $BM1 $line;
}}
close $BM1; #Benchmark files are generated throughout for proofchecking/debugging and data transparency.
print "Done!\n";
print "Placing values into a new file for processing.\n";



#Step 2: Place values into a new file for processing
#Reopen benchmark file, open area file and create a comparison for the analsysis
open($BM1,'<', <BenchMark1.txt>) or die "$in cannot be opened.\n";
open(my $comp,'<', <BenchMark1.txt>) or die "$in cannot be opened.\n";
open (my $AF, '<',$areafile) or die "$areafile cannot be opened.\n";

#Need a dummy line to remove header and shift second array down 1
my $dummy = <$AF>;
my $dummy2 = $comp;

while(defined((my $line1 = <$BM1>)) and (defined((my $line2 = <$AF>))) and (defined((my $line3 = <$comp>)))){
my @LINE_BM1 = split(/\s+/, $line1);
my @LINE_AF = split(/\s+/, $line2);
#Extract what I need and write into a new file 
my $out_val = $LINE_AF[2];
my $in_val = $LINE_BM1[5];
my $AminoAcid = $LINE_BM1[3];

print $BM2 "$in_val\t$out_val\t$AminoAcid\n";
}
#File needs a sacrificial line for some reason, including it here
print $BM2 "9999999999999\t0\tEND\n"; #Will be automatically omitted from any calculations
close $BM2;
print "Done!\n";



print "Adding up solvent exposure area for each atom.\n";
#Step 3: Add up all the values in the second column that share the same value for the first column
open($BM2,'<', <BenchMark2.txt>);

#Set up total for later
my $total = 0;
my $count2 = 1;

foreach my $line (<$BM2>){
my @LINE = split(/\s+/, $line);
if ($count2 eq $LINE[0]){
	$total += $LINE [1];
}
elsif($count2 ne $LINE[0]){
	$count2++;
	my $line_num = $count2-1;
	print $BM3 "$line_num\t$total\n";
#Where it gets added to a new file
	$total = 0;
	$total += $LINE [1];
}else{}
}
close $BM2;
close $BM3;
print "Done!\n";



print "Associating amino acids with correct values.\n";
#Step 4: Get the AA in line with the correct values
open($BM2,'<', <BenchMark2.txt>); #Open BM2 file for the AA column
open($BM3,'<', <BenchMark3.txt>); #Open BM3 file for the data
open(my $BM4,'>', <BenchMark4.txt>); #Open our output file for writing

#Yet another counter
my $count3 = 1;

#This is able to get each of the amino acids in the correct spot
my $line3 = <$BM2>;
foreach $line3 (<$BM2>){
my @LINE3 = split(/\s/, $line3);
if ($LINE3[0] == $count3){
my $data = &GetData();
print $BM4 "$data\t$LINE3[2]\n";
$count3++}
else{}
}

#Subroutine to pull the data
sub GetData{
	while (my $data = <$BM3>){
		chomp $data;
return $data}}
close $BM2;
close $BM3;
close $BM4;




print "Done!\n";
#STEP 5: Comparison of data with fully solvent exposed amino acids\n";

print "Performing C\comparison with fully solvent exposed amino acids.\n";
open($BM4,'<', <BenchMark4.txt>); #Open BM3 file for the data
open(my $BM5,'>', <BenchMark5.txt>); #Open our output file for writing
print $BM5 "#Atom_Num\tValue\tAmino_ID\tExternal_AA?\n";
open($BM5,'>>', <BenchMark5.txt>); #Open our output file for writing
#Reference table for next analysis
my $ALA_Ao = 67.04;
my $ARG_Ao = 144.24;
my $ASN_Ao = 89.16;
my $ASP_Ao = 86.91;
my $CYS_Ao = 84.18;
my $GLN_Ao = 108.05;
my $GLU_Ao = 105.72;
my $GLY_Ao = 50.78;
my $HIS_Ao = 107.31;
my $ILE_Ao = 107.75;
my $LEU_Ao = 108;
my $LYS_Ao = 131.2;
my $MET_Ao = 118.62;
my $PHE_Ao = 114.92;
my $PRO_Ao = 82.91;
my $SER_Ao = 73.05;
my $THR_Ao = 84.92;
my $TRP_Ao = 119.89;
my $TYR_Ao = 123.25;
my $VAL_Ao = 91.95;

#Subroutine to match amino acid with appropiate Ao value is here
my $AA_ID = "WORD";
sub Get_Ao{
if($AA_ID eq "ALA"){
	return $ALA_Ao;
}elsif($AA_ID eq "ARG"){
	return $ARG_Ao;
}elsif($AA_ID eq "ASN"){
	return $ASN_Ao;
}elsif($AA_ID eq "ASP"){
	return $ASP_Ao;
}elsif($AA_ID eq "CYS"){
	return $CYS_Ao;
}elsif($AA_ID eq "GLN"){
	return $GLN_Ao;
}elsif($AA_ID eq "GLU"){
	return $GLU_Ao;
}elsif($AA_ID eq "GLY"){
	return $GLY_Ao;
}elsif($AA_ID eq "HIS"){
	return $HIS_Ao;
}elsif($AA_ID eq "ILE"){
	return $ILE_Ao;
}elsif($AA_ID eq "LEU"){
	return $LEU_Ao;
}elsif($AA_ID eq "LYS"){
	return $LYS_Ao;
}elsif($AA_ID eq "MET"){
	return $MET_Ao;
}elsif($AA_ID eq "PHE"){
	return $PHE_Ao;
}elsif($AA_ID eq "PRO"){
	return $PRO_Ao;
}elsif($AA_ID eq "SER"){
	return $SER_Ao;
}elsif($AA_ID eq "THR"){
	return $THR_Ao;
}elsif($AA_ID eq "TRP"){
	return $TRP_Ao;
}elsif($AA_ID eq "TYR"){
	return $TYR_Ao;
}elsif($AA_ID eq "VAL"){
	return $VAL_Ao;
}else{
return "EMPTY";
}}

#open file again
#split things by column
foreach my $line_Ao (<$BM4>){
my @LINE_Ao = split(/\s+/, $line_Ao);
$AA_ID = $LINE_Ao[2];
my $Ao = &Get_Ao($LINE_Ao[2]);
my $Value = $LINE_Ao[1];
my $num = $LINE_Ao[0];

#calculation step
my $outval = $Value/$Ao;
if ($outval > 0.5){
print $BM5 "$num\t$outval\t$AA_ID\tYES\n"}
else{
print $BM5 "$num\t$outval\t$AA_ID\n"
}}



print "Done!\n\nCheck benchmark files to review data or for data transparency.\nNOTE: Benchmark 5 is final output.\n";



exit;