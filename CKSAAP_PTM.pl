#!/usr/bin/perl
die "usage: $0 <datadirname> <traindata> <testdata> <up> <down>\n" unless @ARGV == 5;
my ($dirname, $trainfile, $testfile, $up, $down) = @ARGV;
my ($ptrainnum, $ptestnum);
open(TR,"./$dirname/$trainfile") or die "can not open TR!\n$dirname\t$trainfile";
open(TE,"./$dirname/$testfile") or die "can not open TE!\n$dirname\t$testfile";
my @train = <TR>;
my @test = <TE>;

my $set = 5; # the value of k(i.e.kmax=4)

my $sidesum = 20;

my $start = $sidesum - $up;
my $length = $up + $down + 1;

$ptestnum= "";
$ptrainnum= "";
open(TRAIN,">./$dirname/$trainfile.$up.$down.encode")or die "can not open train file!";
open(TEST,">./$dirname/$testfile.$up.$down.encode")or die "can not open test file!";
for($i=0;$i<@test;$i++)
{
	if($test[$i]=~/\s+positive\s+/)
	{
		$ptestnum++;
	}
}
for($j=0;$j<@train;$j++)
{
	if($train[$j]=~/\s+positive\s+/)
	{
		$ptrainnum++;
	}
}

for($te=0;$te<@test;$te++)
{
	@l=();
	if($test[$te]=~/\S+\s+\d+\s+(\S+)/)
	{
		$rawseq=$1;
	}
	$seq=substr($rawseq,$start,$length);
	@l=kspace($seq,$length,$set);
	@etest=encoding_print($ptestnum,$te,@l);
	print TEST "@etest";
}

for($tr=0;$tr<@train;$tr++)
{
	@l=();
	if($train[$tr]=~/\S+\s+\d+\s+(\S+)/)
	{
		$rawtrseq=$1;
	}
	$trseq=substr($rawtrseq,$start,$length);
	@trl=kspace($trseq,$length,$set);
	@etrain=encoding_print($ptrainnum,$tr,@trl);
	print TRAIN "@etrain";
}

#to generate a SVM model and test with SVM_light software
#you have to change the path of SVM_light and the parameters can be changed if you want
print "./svm_learn -j 5 -t 2 -g 0.01 -c 0.8  ./$dirname/$trainfile.$up.$down.encode ./$dirname/$trainfile.$up.$down.model\n";
print "./svm_classify ./$dirname/$testfile.$up.$down.encode ./$dirname/$trainfile.$up.$down.model ./$dirname/$testfile.$up.$down.prediction\n";
`./svm_learn -j 5 -t 2 -g 0.01 -c 0.8  ./$dirname/$trainfile.$up.$down.encode ./$dirname/$trainfile.$up.$down.model`;
`./svm_classify ./$dirname/$testfile.$up.$down.encode ./$dirname/$trainfile.$up.$down.model ./$dirname/$testfile.$up.$down.prediction`;
calperformance($ptestnum, $dirname, $testfile, $up, $down);
####################################################################################################

#prediction peformance
sub calperformance{
	my($pnum, $dirname, $testfile, $up, $down);
	$pnum  = $_[0];
	$dirname = $_[1];
	$testfile = $_[2];
	$up = $_[3];
	$down = $_[4];
	open(PRE,"./$dirname/$testfile.$up.$down.prediction") or die "can not open prediction file!";
	open(RES,">./$dirname/$testfile.$up.$down.result") or die "can not open result file!";
	open(TE,"./$dirname/$testfile") or die "can not open test file!";
	my @testdata=<TE>;
	@file=<PRE>;
	#print RES "#the first column is the ID of sequence in Swissprot,\n#the second column is the position of glycosylation site,\n#the third column is the fragment with 40 amino acids around the glycosylation site,\n#the fourth column denotes whether the site in the middle is glycosylation site,\n#the fifth column is the predicted score,\n#the sixth column represents if the site was correctly predicted, T represents a correct prediction, F represents a wrong prediction\n"; 
	#print RES "#===============================================================================\n";
	for($f=0;$f<@file;$f++)
	{
		chomp($testdata[$f]);
		chomp($testdata[$f]);
		#print $testdata[$f];
		$score = sprintf("%.4f",$file[$f]);
		print RES "$testdata[$f]  $score\n";
	}
	print RES "#================================================================================\n";
	@result=cal_result($pnum,@file);
	$string=join" ",@result;
	@array=split/ /,$string;
	@meaturement=qw(Accuracy Sensitivity Specificity MCC);
	for($m=0;$m<@meaturement;$m++)
	{
		$value=sprintf("%.3f",$array[$m]);
		print RES "#$meaturement[$m] = $value\n";
	}
#`rm ./$dirname/*prediction`;
#END sub calperformance
}
####################################################################################################
sub kspace{
#need two parameters
#the first is the sequence 
#the second is the window size
#the third is the k-value
my($seq,@one,@two,$set,$kspace,$ti,$i,$j,$length,$w,$ab,$a,$b,@ab,$l,@l);
$seq=$_[0];
$length=$_[1];
$set=$_[2];
    for($ti=1;$ti<=$set;$ti++)
         {
       @one=qw(A C D E F G H I K L M N P Q R S T V W Y -);
        for($i=0;$i<@one;$i++)
          {
          @two=qw(A C D E F G H I K L M N P Q R S T V W Y -);
          for($j=0;$j<@two;$j++)
             {
               @ab=();
               for($w=0;$w<$length;$w++)
                  {
                   $ab="";
                   $a=substr($seq,$w,1);
                   $b=substr($seq,$w+$ti,1);
                   $ab=$a.$b;
                   push @ab,$ab;
                   }
                  $kspace="";
                  $kspace=$one[$i].$two[$j];
                  $abstring=join" ",@ab;
                  @km=$abstring=~/$kspace/g;
                  $l=@km;
                  push @l,$l;
             }
          }
      }
return @l;
}
####################################################################################################
sub encoding_print{
my($pos,$split,@test,$r,$ll,@l,@txt,@file,$i);
@txt=@_;
$split=$_[0];
$pos=$_[1];
@l=();
for($i=2;$i<@txt;$i++)
   {
    push @l,$txt[$i];
   }
         if($pos<$split)
           {
             push @test,"+1 ";
             for($ll=0;$ll<@l;$ll++)
                {
                $r=$ll+1;
                push @test,"$r:$l[$ll] ";
                }
                push @test,"\n";
              }
          else
            {
            push @test,"-1 ";
             for($ll=0;$ll<@l;$ll++)
                {
                $r=$ll+1;
                push @test,"$r:$l[$ll] ";
                }
                push @test,"\n";
             }
return @test;
}

####################################################################################################
sub cal_result{
my ($i,$j,$f,@txt,@file);
my $in=$_[0];
@txt=@_;
@file=();
for($f=1;$f<@txt;$f++)
   {
   push @file,$txt[$f];
   }
my $long=@file;
my $n=0;
my $m=0;
for($i=0;$i<$in;$i++)
   {
   if($file[$i]>=0)
     {
      $m++;
     }
   }
for($j=$in;$j<$long;$j++)
    {
   if($file[$j]<0)
     {
     $n=$n+1;
     }
    }
$tp=$m;
$tn=$n;
$fn=$in-$tp;
$fp=$long-$in-$tn;
$sn=$tp/($tp+$fn);
$sp=$tn/($tn+$fp);
$ac=($tp+$tn)/($tp+$fp+$tn+$fn);
$fp_rate=($long-$in-$tn)/($tn+$fp);
$ccpre = (($tn+$fn)*($tn+$fp)*($tp+$fn)*($tp+$fp))**(1/2);
if($ccpre == 0)
{
	$cc = "error_division_zero";
}
else
{
	$cc=($tp*$tn-$fp*$fn)/$ccpre;
}
return "$ac $sn $sp $cc";
}