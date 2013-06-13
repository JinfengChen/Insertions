#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);

GetOptions (\%opt,"mPing:s","ref:s","help");


my $help=<<USAGE;
perl $0 --mPing --ref
--mPing: sam alignment of read mapped to mPing
--ref:   reference genome
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

##read mPing sam and store alignment in hash as read name with key and alignment 1 and 2 as value, read->[read1,read2]
print "Step1: Read mPing Sam\n";
my $refmping=readsam($opt{mPing});
##get all mPing associated reads and map them to genome using stampy
##cut hte part of mPing in read and then mapp to 
print "Step2: Map read to genome\n";
#### status: store information for pair, readname->read1or2->0,1,2: 0 mean read in mping, 1 mean read cover partial mping, 2 mean read not cover mping
my $status=map2genome($opt{ref},$refmping);
##read reference sam and store alignment in hash same as mPing

print "Step3: Read genome sam\n";
my $refgenome1=readsam2("HEG4_2.3.mPing.stampy.sam");
my $refgenome2=readsam2("HEG4_2.3.mPing.bwa.sam");
my $refgenome3=readsam2("HEG4_2.3.mPing.stampy.clean_mPing.sam");
my $refgenome4=readsam2("HEG4_2.3.mPing.bwa.clean_mPing.sam");

print "Step4: Find insertions\n";
my %insert;
my %summary;
foreach(keys %$refmping){
   $summary{"01.Total"}++;
   my $read=$_;
   print ">$read\n$refmping->{$read}->[0]\n$refmping->{$read}->[1]\n";
   if ($refgenome1->{$read}){
       print "Stampy:\n$refgenome1->{$read}->[0]\n$refgenome1->{$read}->[1]\n";
   }
   if ($refgenome2->{$read}){
       print "BWA:\n$refgenome2->{$read}->[0]\n$refgenome2->{$read}->[1]\n";
   } 
   if ($refgenome3->{$read}){
       print "Stampy clean:\n$refgenome3->{$read}->[0]\n$refgenome3->{$read}->[1]\n";
   }
   if ($refgenome4->{$read}){
       print "BWA clean:\n$refgenome4->{$read}->[0]\n$refgenome4->{$read}->[1]\n";
   }
   print "\n\n";
   if ($refgenome3->{$read}){
      #my ($chr1,$breakpoint1,$startonreads1,$read1,$strand1); ##chromosome, breakpoint on chromosome, start position on reads, read sequence
      #my ($chr2,$breakpoint2,$startonreads2,$read2,$strand2);
      #status: store information for pair, readname->read1or2->(status,mping in read) status 0,1,2: 0 mean read in mping, 1 mean read cover partial mping, 2 mean read not cover mping. mping in read 5 or 3: 5 mean mping in 5 primer of reads and 3 mean in 3 primer of reads
      if ($status->{$read}->{1}->[0]==0 and $status->{$read}->{2}->[0]==0){     ###both reads in mPing, no effect
         $summary{"02.Pair_in_mPing"}++;
      }elsif($status->{$read}->{1}->[0]==0 and $status->{$read}->{2}->[0]==2){  ###one in mPing, one on flanking, support reads
         $summary{"03.Support_mPing"}++;
      }elsif($status->{$read}->{2}->[0]==0 and $status->{$read}->{1}->[0]==2){  ###one in mPing, one on flanking, support reads
         $summary{"03.Support_mPing"}++;
      }elsif($status->{$read}->{1}->[0]==1 and $status->{$read}->{2}->[0]==1){  ###both reads cover boundary
         $summary{"04.Pair_on_Boundary"}++;
         my ($chr1,$breakpoint1,$startonreads1,$read1,$strand1)=analyze($refgenome3->{$read}->[0],$status->{$read}->{1}->[1]);
         push @{$insert{$chr1}{$breakpoint1}},[$startonreads1,$read1,$strand1,$read,1,$status->{$read}->{1}->[1]] unless ($chr1 eq "NA");
         my ($chr2,$breakpoint2,$startonreads2,$read2,$strand2)=analyze($refgenome3->{$read}->[1],$status->{$read}->{2}->[1]);
         push @{$insert{$chr2}{$breakpoint2}},[$startonreads2,$read2,$strand2,$read,2,$status->{$read}->{1}->[1]] unless ($chr2 eq "NA");
      }elsif($status->{$read}->{1}->[0]==1 and $status->{$read}->{2}->[0]==0){  ###one in mPing, one on boundary
         $summary{"05.One_mPing_one_Boundary"}++;
      }elsif($status->{$read}->{1}->[0]==0 and $status->{$read}->{2}->[0]==1){  ###one in mPing, one on boundary
         $summary{"05.One_mPing_one_Boundary"}++;
      }elsif($status->{$read}->{1}->[0]==1 and $status->{$read}->{2}->[0]==2){  ###one on boundary, one on flanking
         $summary{"06.One_Flanking_one_Boundary"}++;
      }elsif($status->{$read}->{1}->[0]==2 and $status->{$read}->{2}->[0]==1){  ###one on boundary, one on flanking
         $summary{"06.One_Flanking_one_Boundary"}++;
      }
   }
}

print "Step5: Summary:\n";
foreach(sort keys %summary){
   print "$_\t$summary{$_}\n";
}

my $refseq=getfastaseq($opt{ref});
foreach my $chr (sort keys %insert){
   foreach my $bp (sort {$a <=> $b} keys %{$insert{$chr}}){
     my $ref=substr($refseq->{$chr},$bp-201,500);
     $ref=~tr/actg/ACTG/;
     #my $ref="N" x 400;
     print ">$chr $bp\n$ref\n";
     my @align=@{$insert{$chr}{$bp}};
     for(my $i=0;$i<@align;$i++){
        if (($align[$i]->[5] == 5 and $align[$i]->[2] == 1) or ($align[$i]->[5] == 3 and $align[$i]->[2] == 0)){ ##mping in 5 and read formard mapped and 3,0
           my $read=" " x 200;
           $sequence=sprintf("%-90s",$align[$i]->[1]);
           $read.="$sequence  $align[$i]->[2]  $align[$i]->[3]  $align[$i]->[4]  $align[$i]->[5]\n";
           print "$read\n";
        }else{
           my $read=sprintf("%200s",$align[$i]->[1]);
           my $temp=" " x 90;
           $read.="$temp  $align[$i]->[2]  $align[$i]->[3]  $align[$i]->[4]  $align[$i]->[5]\n";
           print "$read\n";
        }
     }
   }
}



############################
sub map2genome
{
my ($ref,$mping)=@_;
$read1="HEG4_2.3.mPing.p1.fq";
$read2="HEG4_2.3.mPing.p2.fq";
`/rhome/cjinfeng/software/tools/seqtk-master/seqtk subseq HEG4_2.3_p1.fq HEG4_2.3.mPing.reads.list > $read1` unless (-e $read1);
`/rhome/cjinfeng/software/tools/seqtk-master/seqtk subseq HEG4_2.3_p2.fq HEG4_2.3.mPing.reads.list > $read2` unless (-e $read2);

##remove part that match to mping
my $status=remove($read1,$read2,$mping) unless (-e "HEG4_2.3.mPing.p1.clean_mPing.fq" and -e "HEG4_2.3.mPing.p2.clean_mPing.fq");

##map to genome
$read3="HEG4_2.3.mPing.p1.clean_mPing.fq";
$read4="HEG4_2.3.mPing.p2.clean_mPing.fq";

$prefix="HEG4_2.3.mPing";
$index="rice_MSU7";


print "Mapping with bwa\n";
`/opt/tyler/bin/bwa aln -l 20 $ref $read1 > $read1.sai` unless (-e "$read1.sai");
`/opt/tyler/bin/bwa aln -l 20 $ref $read2 > $read2.sai` unless (-e "$read2.sai");
`/opt/tyler/bin/bwa aln -l 20 $ref $read3 > $read3.sai` unless (-e "$read3.sai");
`/opt/tyler/bin/bwa aln -l 20 $ref $read4 > $read4.sai` unless (-e "$read4.sai");
`/opt/tyler/bin/bwa sampe $ref $read1.sai $read2.sai $read1 $read2 > $prefix.bwa.sam` unless (-e "$prefix.bwa.sam");
`/usr/local/bin/samtools view -Sb -o $prefix.bwa.bam $prefix.bwa.sam` unless (-e "$prefix.bwa.bam") ;
`/opt/tyler/bin/bwa sampe $ref $read3.sai $read4.sai $read3 $read4 > $prefix.bwa.clean_mPing.sam` unless (-e "$prefix.bwa.clean_mPing.sam");
`/usr/local/bin/samtools view -Sb -o $prefix.bwa.clean_mPing.bam $prefix.bwa.clean_mPing.sam` unless (-e "$prefix.bwa.clean_mPing.bam");

print "Mapping with stampy\n";
`/opt/stampy/1.0.21-py2.7/stampy.py --species=rice --assembly=MSU7 -G $index $ref` unless (-e "$index.stidx");
`/opt/stampy/1.0.21-py2.7/stampy.py -g $index -H $index` unless (-e "$index.sthash");
#`/opt/stampy/1.0.21-py2.7/stampy.py -g $index -h $index -M $read1 $read2 > $prefix.stampy.sam` unless (-e "$prefix.stampy.sam");
#`/opt/stampy/1.0.21-py2.7/stampy.py -g $index -h $index -M $read3 $read4 > $prefix.stampy.clean_mPing.sam` unless (-e "$prefix.stampy.clean_mPing.sam");
`/opt/stampy/1.0.21-py2.7/stampy.py -g $index -h $index --xa-max=3 --xa-max-discordant=10 --bamkeepgoodreads -M $prefix.bwa.bam > $prefix.stampy.sam` unless (-e "$prefix.stampy.sam");
`/opt/stampy/1.0.21-py2.7/stampy.py -g $index -h $index --xa-max=3 --xa-max-discordant=10 --bamkeepgoodreads -M $prefix.bwa.clean_mPing.bam > $prefix.stampy.clean_mPing.sam` unless (-e "$prefix.stampy.clean_mPing.sam");

print "Clean temp in mapping";
if (1){
#`rm $read1.sai $read2.sai $read3.sai $read4.sai`;
#`rm $index.sthash $index.stidx`;
#`rm $read1 $read2`;
`rm $read3 $read4`;
#`rm $prefix.bwa.bam $prefix.bwa.sam $prefix.bwa.clean_mPing.sam $prefix.bwa.clean_mPing.bam`;
#`rm $prefix.stampy.sam $prefix.stampy.clean_mPing.sam`;
}
return $status;
}

sub remove
{
my ($read1,$read2,$mping)=@_;
my $prefix1=basename($read1,".fq");
my $prefix2=basename($read2,".fq");
my $prefix0=basename($read1,".p1.fq");

my %status; ## store information for pair, readname->read1or2->(status,five or three) status 0,1,2: 0 mean read in mping, 1 mean read cover partial mping, 2 mean read not cover mping. five or three mean five primer was mPing and three primer is mPing
print "Remove part in reads that match mPing:\n";
open OUT1, ">$prefix1.clean_mPing.fq" or die "%!";
open OUT2, ">$prefix2.clean_mPing.fq" or die "%!";
#open OUT0, ">$prefix0.Unpair.clean_mPing.fq" or die "%!";
open IN, "$read1" or die "$!";
open IN1, "$read2" or die "$!";
while(<IN>){
     my $head=$_;
     $head=$1 if $head=~/^@(.*)$/;
     my $seq =<IN>;
     my $qhead=<IN>;
     my $qual=<IN>;
     chomp $head;chomp $seq; chomp $qhead; chomp $qual;
     my $head1=<IN1>;my $seq1=<IN1>;$qhead1=<IN1>;$qual1=<IN1>;$head1=$1 if $head1=~/^@(.*)$/;
     chomp $head1; chomp $seq1; chomp $qhead1; chomp $qual1;
     ##subfunction choose a read in mping (1 or 2) and return the start of length to retrieve from the original reads
     my $flag1=1; $flag2=1;
     my $refhash=matchmping($mping->{$head},\%status);
     if ($refhash->{1}->[1] == 0){  ### read1 all in mping, now set this type as start=0, len=101 becasue the result will not influence analysis and keep pair-end are easy to handle
        $flag1=0;
     }else{
        $newseq=substr($seq,$refhash->{1}->[0],$refhash->{1}->[1]);
        $newqual=substr($qual,$refhash->{1}->[0],$refhash->{1}->[1]);
     }
     if ($refhash->{2}->[1] == 0){  ### read2 all in mping
        $flag2=0;
     }else{
        $newseq1=substr($seq1,$refhash->{2}->[0],$refhash->{2}->[1]);
        $newqual1=substr($qual1,$refhash->{2}->[0],$refhash->{2}->[1]);
     }
     if ($flag1 and $flag2){ ### read1 and read2 have sequence left
        print OUT1 "\@$head\n$newseq\n$qhead\n$newqual\n";
        print OUT2 "\@$head1\n$newseq1\n$qhead1\n$newqual1\n";
     }else{ ### read1 or/and read2 have no sequence left
        #print OUT0 "$head\n$newseq\n$qhead\n$newqual\n" if $flag1;
        #print OUT0 "$head1\n$newseq1\n$qhead1\n$newqual1\n" if $flag2;
     }
}
close IN;
close OUT1;
close OUT2;
#close OUT0;
return \%status;
}

sub matchmping
{
my ($align,$status)=@_;
print Dumper($align);
my %hash;
for(my $i=0;$i<2;$i++){
   my $start=0;
   my $len=0;
   my @unit=split("\t",$align->[$i]);
   print "$unit[0]\t$unit[1]\t$unit[2]\t$rank\n";
   my $flag=samflag($unit[1]);
   print "$flag\n";
   my $pair= $flag=~/\'64\'/ ? 1 : 2;    #1 for read 1, 2 for read2
   my $strand = $flag=~/\'16\'/ ? 1 : 0; #0 for forward, 1 for reverse
   print "$pair in Pair\nStrand: $strand\n";
   ##parse cigar to retrieve matches on mPing
   ## start and len: 0, 0 means entile reads match mPing, these reads are useless 
   ##                NA, NA means reads are unable to map to mPing or map is not reliable
   ##                start, len means read are maped to mPing, and start and len will get unmaped region on reads
   if ($unit[5]=~/^\d+M$/){ ### perfect match
      $start=0;
      $len  =length $unit[9];
      $status->{$unit[0]}->{$pair}=[0,0]; 
   }elsif($unit[5]=~/^(\d+.*?)(\d+)M$/){ ###match at right
      my $length=$2;
      my $subalign=$1;
      my $add; ### add length if subalign
      my $maxm=0; ### max Match length in subalign
      while($subalign=~/(\d+)(\D+)/g){
         print "Sub: $1\t$2\n";
         my $t=$1;
         #$add+=$1;
         if ($2=~/M/){
            $maxm= $t > $maxm ? $t : $maxm;
         }
      }
      $add=(length $unit[9])-$length;
      print "MAX match in sub: $maxm\tLength of match: $length\n";
      if ($maxm <= $length){ ### match is subalign smaller than match at the end
         if ($strand == 0){ ###forward
            $start=0;
            $len  =$add;
         }else{
            $start=$length;
            $len  =$add;
         }
         $direct= $strand == 1 ? 3 : 5;
         $status->{$unit[0]}->{$pair}=[1,$direct];
      }else{ ### match in subalign larger than match at the end, probably unreliable align
         $start=0;
         $len  =length $unit[9];
         $status->{$unit[0]}->{$pair}=[2,0];
         print "Strange subalign: $align->[$i]\n";
      }
   }elsif($unit[5]=~/^(\d+)M(\d+.*)$/){ ###match at left
      my $length=$1;
      my $subalign=$2;
      my $add=0;my $maxm=0;
      print "Subalign: $subalign\n";
      while($subalign=~/(\d+)(\D+)/g){
         #$add+=$1;
         my $t=$1;
         print "Add and \$1: $add\t$1\n";
         if ($2=~/M/){
            $maxm= $t > $maxm ? $t : $maxm;
         }
      } 
      $add=(length $unit[9])-$length; 
      if ($maxm < $length){ ### match is subalign smaller than match at the end
         if ($strand == 0){ ###forward strand
            $start=$length;
            $len  =$add;
         }else{
            $start=0;
            $len  =$add;
         }
         $direct= $strand == 1 ? 5 : 3;
         $status->{$unit[0]}->{$pair}=[1,$direct];
      }else{ ### match in subalign larger than match at the end, probably unreliable align
         $start=0;
         $len  =length $unit[9];
         $status->{$unit[0]}->{$pair}=[2,0];
         print "Strange subalign: $align->[$i]\n";
      }

   }elsif($unit[5]=~/\*/){ ### * no match
      $start=0;
      $len  =length $unit[9];
      $status->{$unit[0]}->{$pair}=[2,0];
   }else{ ### other unreliable match
      $start=0;
      $len  =length $unit[9];
      $status->{$unit[0]}->{$pair}=[2,0];
      print "Strange: $align->[$i]\n";
   }
   $hash{$pair}=[$start,$len];
   print "Read $pair: $start, $len\n";
}
return \%hash;
}

sub analyze
{
my ($align,$direct)=@_;
#print Dumper($align);
my @unit=split("\t",$align);
my $flag=samflag($unit[1]);
my $pair= $flag=~/\'64\'/ ? 1 : 2;    #1 for read 1, 2 for read2
my $strand = $flag=~/\'16\'/ ? 1 : 0; #0 for forward, 1 for reverse
#print "$pair in Pair\nStrand: $strand\n";

##parse cigar to retrieve matches on genome
##chr
##breakpoint
##strand
##startonread
##readsequence  
my ($chr,$breakpoint,$start,$read); 
if ($unit[5]=~/^(\d+)M$/){ ### perfect match
   $chr=$unit[2];
   if ($strand == 1){
      $breakpoint=$direct == 5 ? $unit[3] : $unit[3]+length $unit[9];
   }else{
      $breakpoint=$direct == 5 ? $unit[3]+length $unit[9] : $unit[3];
   }
   $start=0;
   $read=$unit[9];
}elsif($unit[5]=~/^(\d+.*?)(\d+)M$/){ ###match at right
      my $length=$2;
      my $subalign=$1;
      my $add; ### add length if subalign
      my $maxm=0; ### max Match length in subalign
      while($subalign=~/(\d+)(\D+)/g){
         print "Sub: $1\t$2\n";
         my $t=$1;
         if ($2=~/M/){
            $maxm= $t > $maxm ? $t : $maxm;
         }
      }
      $add=(length $unit[9])-$length;
      print "MAX match in sub: $maxm\tLength of match: $length\n";
      if ($maxm <= $length){ ### match is subalign smaller than match at the end
         $chr=$unit[2];
         if ($strand == 1){
            $breakpoint=$direct == 5 ? $unit[3] : $unit[3]+$length;
         }else{
            $breakpoint=$direct == 5 ? $unit[3]+$length : $unit[3];
         }
         $start=$add;
         $read=$unit[9];
      }else{ ### match in subalign larger than match at the end, probably unreliable align
         print "Strange subalign: $align\n";
         $chr="NA";        
         $breakpoint="NA";
         $start="NA";
         $read="NA";

      }
}elsif($unit[5]=~/^(\d+)M(\d+.*)$/){ ###match at left
      my $length=$1;
      my $subalign=$2;
      my $add=0;my $maxm=0;
      print "Subalign: $subalign\n";
      while($subalign=~/(\d+)(\D+)/g){
         my $t=$1;
         print "Add and \$1: $add\t$1\n";
         if ($2=~/M/){
            $maxm= $t > $maxm ? $t : $maxm;
         }
      } 
      $add=(length $unit[9])-$length; 
      if ($maxm < $length){ ### match is subalign smaller than match at the end
         $chr=$unit[2];  
         if ($strand == 1){
            $breakpoint=$direct == 5 ? $unit[3] : $unit[3]+$length;
         }else{
            $breakpoint=$direct == 5 ? $unit[3]+$length : $unit[3];
         }
         $start=0;
         $read=$unit[9];
      }else{ ### match in subalign larger than match at the end, probably unreliable align
         print "Strange subalign: $align\n";
         $chr="NA";       
         $breakpoint="NA";
         $start="NA";
         $read="NA";

      }
}elsif($unit[5]=~/\*/){ ### * no match
      $chr="NA";   
      $breakpoint="NA";
      $start="NA";
      $read="NA";
}else{ ### other unreliable match
      print "Strange: $align->[$i]\n";
      $chr="NA";               
      $breakpoint="NA";
      $start="NA";
      $read="NA";
}
return ($chr,$breakpoint,$start,$read,$strand);
}


sub readbam
{
my ($bam,$mping)=@_;
open OUT, ">temp.mping.ref.sam" or die "$!";
open IN, "/usr/local/bin/samtools view $bam |" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/ or $_=~/^@/);
   my $read1=$_;
   my @unit=split("\t",$read1);
   if (exists $mping->{$unit[0]}){
      print OUT "$read1\n";
   }
}
close IN;
close OUT;
`sort -k1,1 temp.mping.ref.sam > temp.mping.ref.sort.sam`;
}

sub readsam
{
my ($file)=@_;
my %hash;
my %count;
open OUT, ">HEG4_2.3.mPing.reads.list" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^@/);
    my $read=$_;
    my @unit=split("\t",$read);
    my $flag=samflag($unit[1]);
    my $pair= $flag=~/\'64\'/ ? 0 : 1;    #0 for read 1, 1 for read2
    $count{$unit[0]}++;
    $hash{$unit[0]}->[$pair]=$read;
    print OUT "$unit[0]\n" if ($count{$unit[0]}==1);
}
close IN;
close OUT;
return \%hash;
}

sub readsam2
{
my ($file)=@_;
my %hash;
my %count;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^@/);
    my $read=$_;
    my @unit=split("\t",$read);
    my $flag=samflag($unit[1]);
    my $pair= $flag=~/\'64\'/ ? 0 : 1;    #0 for read 1, 1 for read2
    $count{$unit[0]}++;
    $hash{$unit[0]}->[$pair]=$read;
}
close IN;
return \%hash;
}

sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}




sub samflag
{
my ($flag)=@_;
my %TABLE= qw(  0   '0'
1 '1'
2 '2'
3 '1''2'
4 '4'
5 '1''4'
6 '2''4'
7 '1''2''4'
8 '8'
9 '1''8'
10 '2''8'
11 '1''2''8'
12 '4''8'
13 '1''4''8'
14 '2''4''8'
15 '1''2''4''8'
16 '16'
17 '1''16'
18 '2''16'
19 '1''2''16'
20 '4''16'
21 '1''4''16'
22 '2''4''16'
23 '1''2''4''16'
24 '8''16'
25 '1''8''16'
26 '2''8''16'
27 '1''2''8''16'
28 '4''8''16'
29 '1''4''8''16'
30 '2''4''8''16'
31 '1''2''4''8''16'
32 '32'
33 '1''32'
34 '2''32'
35 '1''2''32'
36 '4''32'
37 '1''4''32'
38 '2''4''32'
39 '1''2''4''32'
40 '8''32'
41 '1''8''32'
42 '2''8''32'
43 '1''2''8''32'
44 '4''8''32'
45 '1''4''8''32'
46 '2''4''8''32'
47 '1''2''4''8''32'
48 '16''32'
49 '1''16''32'
50 '2''16''32'
51 '1''2''16''32'
52 '4''16''32'
53 '1''4''16''32'
54 '2''4''16''32'
55 '1''2''4''16''32'
56 '8''16''32'
57 '1''8''16''32'
58 '2''8''16''32'
59 '1''2''8''16''32'
60 '4''8''16''32'
61 '1''4''8''16''32'
62 '2''4''8''16''32'
63 '1''2''4''8''16''32'
64 '64'
65 '1''64'
66 '2''64'
67 '1''2''64'
68 '4''64'
69 '1''4''64'
70 '2''4''64'
71 '1''2''4''64'
72 '8''64'
73 '1''8''64'
74 '2''8''64'
75 '1''2''8''64'
76 '4''8''64'
77 '1''4''8''64'
78 '2''4''8''64'
79 '1''2''4''8''64'
80 '16''64'
81 '1''16''64'
82 '2''16''64'
83 '1''2''16''64'
84 '4''16''64'
85 '1''4''16''64'
86 '2''4''16''64'
87 '1''2''4''16''64'
88 '8''16''64'
89 '1''8''16''64'
90 '2''8''16''64'
91 '1''2''8''16''64'
92 '4''8''16''64'
93 '1''4''8''16''64'
94 '2''4''8''16''64'
95 '1''2''4''8''16''64'
96 '32''64'
97 '1''32''64'
98 '2''32''64'
99 '1''2''32''64'
100 '4''32''64'
101 '1''4''32''64'
102 '2''4''32''64'
103 '1''2''4''32''64'
104 '8''32''64'
105 '1''8''32''64'
106 '2''8''32''64'
107 '1''2''8''32''64'
108 '4''8''32''64'
109 '1''4''8''32''64'
110 '2''4''8''32''64'
111 '1''2''4''8''32''64'
112 '16''32''64'
113 '1''16''32''64'
114 '2''16''32''64'
115 '1''2''16''32''64'
116 '4''16''32''64'
117 '1''4''16''32''64'
118 '2''4''16''32''64'
119 '1''2''4''16''32''64'
120 '8''16''32''64'
121 '1''8''16''32''64'
122 '2''8''16''32''64'
123 '1''2''8''16''32''64'
124 '4''8''16''32''64'
125 '1''4''8''16''32''64'
126 '2''4''8''16''32''64'
127 '1''2''4''8''16''32''64'
128 '128'
129 '1''128'
130 '2''128'
131 '1''2''128'
132 '4''128'
133 '1''4''128'
134 '2''4''128'
135 '1''2''4''128'
136 '8''128'
137 '1''8''128'
138 '2''8''128'
139 '1''2''8''128'
140 '4''8''128'
141 '1''4''8''128'
142 '2''4''8''128'
143 '1''2''4''8''128'
144 '16''128'
145 '1''16''128'
146 '2''16''128'
147 '1''2''16''128'
148 '4''16''128'
149 '1''4''16''128'
150 '2''4''16''128'
151 '1''2''4''16''128'
152 '8''16''128'
153 '1''8''16''128'
154 '2''8''16''128'
155 '1''2''8''16''128'
156 '4''8''16''128'
157 '1''4''8''16''128'
158 '2''4''8''16''128'
159 '1''2''4''8''16''128'
160 '32''128'
161 '1''32''128'
162 '2''32''128'
163 '1''2''32''128'
164 '4''32''128'
165 '1''4''32''128'
166 '2''4''32''128'
167 '1''2''4''32''128'
168 '8''32''128'
169 '1''8''32''128'
170 '2''8''32''128'
171 '1''2''8''32''128'
172 '4''8''32''128'
173 '1''4''8''32''128'
174 '2''4''8''32''128'
175 '1''2''4''8''32''128'
176 '16''32''128'
177 '1''16''32''128'
178 '2''16''32''128'
179 '1''2''16''32''128'
180 '4''16''32''128'
181 '1''4''16''32''128'
182 '2''4''16''32''128'
183 '1''2''4''16''32''128'
184 '8''16''32''128'
185 '1''8''16''32''128'
186 '2''8''16''32''128'
187 '1''2''8''16''32''128'
188 '4''8''16''32''128'
189 '1''4''8''16''32''128'
190 '2''4''8''16''32''128'
191 '1''2''4''8''16''32''128'
192 '64''128'
193 '1''64''128'
194 '2''64''128'
195 '1''2''64''128'
196 '4''64''128'
197 '1''4''64''128'
198 '2''4''64''128'
199 '1''2''4''64''128'
200 '8''64''128'
201 '1''8''64''128'
202 '2''8''64''128'
203 '1''2''8''64''128'
204 '4''8''64''128'
205 '1''4''8''64''128'
206 '2''4''8''64''128'
207 '1''2''4''8''64''128'
208 '16''64''128'
209 '1''16''64''128'
210 '2''16''64''128'
211 '1''2''16''64''128'
212 '4''16''64''128'
213 '1''4''16''64''128'
214 '2''4''16''64''128'
215 '1''2''4''16''64''128'
216 '8''16''64''128'
217 '1''8''16''64''128'
218 '2''8''16''64''128'
219 '1''2''8''16''64''128'
220 '4''8''16''64''128'
221 '1''4''8''16''64''128'
222 '2''4''8''16''64''128'
223 '1''2''4''8''16''64''128'
224 '32''64''128'
225 '1''32''64''128'
226 '2''32''64''128'
227 '1''2''32''64''128'
228 '4''32''64''128'
229 '1''4''32''64''128'
230 '2''4''32''64''128'
231 '1''2''4''32''64''128'
232 '8''32''64''128'
233 '1''8''32''64''128'
234 '2''8''32''64''128'
235 '1''2''8''32''64''128'
236 '4''8''32''64''128'
237 '1''4''8''32''64''128'
238 '2''4''8''32''64''128'
239 '1''2''4''8''32''64''128'
240 '16''32''64''128'
241 '1''16''32''64''128'
242 '2''16''32''64''128'
243 '1''2''16''32''64''128'
244 '4''16''32''64''128'
245 '1''4''16''32''64''128'
246 '2''4''16''32''64''128'
247 '1''2''4''16''32''64''128'
248 '8''16''32''64''128'
249 '1''8''16''32''64''128'
250 '2''8''16''32''64''128'
251 '1''2''8''16''32''64''128'
252 '4''8''16''32''64''128'
253 '1''4''8''16''32''64''128'
254 '2''4''8''16''32''64''128'
255 '1''2''4''8''16''32''64''128');
return $TABLE{$flag};
}
