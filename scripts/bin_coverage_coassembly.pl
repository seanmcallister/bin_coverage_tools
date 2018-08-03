#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Use the metabat depth output (x.depth.txt) and bins from DASTool to generate an average coverage per bin.

# - - - - - U S E R    V A R I A B L E S - - - - - - - -


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %options=();
getopts("b:c:o:x:n:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-b = path to bin folder\n";
        print "-c = depth coverage file from metabat\n";
        print "-o = output prefix\n";
        print "-x = bin file suffix (i.e. fasta, fa, fna, etc)\n";
        print "-n = Number of samples in coassembly.\n";
	print "-h = This help message\n\n";
	die;
    }

my %Sequences;
my %Bin_names;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
opendir(DIR, "$options{b}") or die "\n\nNada $options{b} you fool!!\n\n";
my $fastaread = 1;
while (my $file = readdir(DIR))
    {   my $binsuffix = qr/$options{x}/;
        if ($file =~ m/$binsuffix$/)
            {   $Bin_names{$file}{'Name'} = $file;
                $Bin_names{$file}{'FileNum'} = $fastaread;
                foreach my $b (1..$options{n})
                    {   $Bin_names{$file}{'total_cov_'.$b} = 0;
			$Bin_names{$file}{'total_len_'.$b} = 0;
                    }
                &FASTAread($file, $fastaread);
                $fastaread += 1;
            }
    }

#Process coverage file to split by bins.
open(IN2, "<$options{c}") or die "\n\nNADA $options{c} you FOOL!!!\n\n";
my @cov_data = <IN2>; close(IN2);
my $cov_head = shift(@cov_data); chomp($cov_head); my @cov_h = split('\t', $cov_head);
foreach my $line (@cov_data)
    {   chomp($line);
        my @data = split('\t', $line);
	foreach my $b (1..$options{n})
	     {	if (exists $Sequences{$data[0]})
		{ $Sequences{$data[0]}{'HeaderLib_'.$b} = $cov_h[(2*$b)+1];
		  $Sequences{$data[0]}{'COV_'.$b} = $data[(2*$b)+1];
		  $Sequences{$data[0]}{'COV_Len'} = $data[1];
		}
	     }
    } 
   
foreach my $k (sort keys %Sequences)
    {   foreach my $b (1..$options{n})
            {   my $math_cov = $Sequences{$k}{'COV_'.$b};
                my $math_len = $Sequences{$k}{'COV_Len'};
                my $inflated_cov = $math_cov * $math_len;
                #print "$Sequences{$k}{'HEAD'}\t$math_cov\t$math_len\t$inflated_cov\n";
                $Sequences{$k}{'sizeinflatedcov_'.$b} = $inflated_cov;
            }
    }

foreach my $i (sort keys %Bin_names)
    {   foreach my $j (sort keys %Sequences)
            {   if ($Sequences{$j}{'filename'} eq $Bin_names{$i}{'Name'})
                    {   foreach my $b (1..$options{n})
                            {   $Bin_names{$i}{'total_cov_'.$b} = $Bin_names{$i}{'total_cov_'.$b} + $Sequences{$j}{'sizeinflatedcov_'.$b};
                                $Bin_names{$i}{'total_len_'.$b} = $Bin_names{$i}{'total_len_'.$b} + $Sequences{$j}{'COV_Len'};
                                #print "$Bin_names{$i}{'Name'}\t$Bin_names{$i}{'total_cov'}\t$Bin_names{$i}{'total_len'}\n";
                            }
                    }                
            }
    }

open(OUT, ">".$options{o}."_lenweighted_avgbincov.txt");
print OUT "Bin_Name";
foreach my $b (1..$options{n})
    {   print OUT "\t$cov_h[(2*$b)+1]";}
print OUT "\n";
foreach my $m (sort keys %Bin_names)
    {	if ($Bin_names{$m}{'total_len_1'} == 0)
          {print "Warning: some bins appear to have a 0 length:\n$Bin_names{$m}{'Name'}\n";}
        else
	{   print OUT "$Bin_names{$m}{'Name'}";
            foreach my $b (1..$options{n})
            {   my $total_cov = $Bin_names{$m}{'total_cov_'.$b};
                my $total_len = $Bin_names{$m}{'total_len_'.$b};
                my $avg_cov = $total_cov / $total_len;
                print OUT "\t$avg_cov";
            }
            #print "$Bin_names{$m}{'Name'}\t$avg_cov\n";
            print OUT "\n";
	}
    }
close(OUT);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub FASTAread
{	#print "   Reading file . . . \n";
	# 1. Load FIlE . . . . . . . . . .
	$/=">";                                     # set input break string
	my $infile = $_[0];
        my $filenumber = $_[1];
        my $path = $options{b};
	open(IN, "<".$path."/".$infile) or die "\n\nNADA $path/$infile you FOOL!!!\n\n";
	my @DATA = <IN>; close(IN); shift(@DATA);	
	# 2. Parse sequence data . . . . . . . . . . . . .
	#my $unid = $filenumber.10000001;                           # string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		#my $seq = '';
		#foreach my $i (1..$#data)
		#{	$seq .= $data[$i];  }
		#$seq =~ s/>//;
		$Sequences{$data[0]}{'HEAD'}    = $data[0];       # store header
		my @shorthead = split(' ', $data[0]);
		$Sequences{$data[0]}{'SHORT_HEAD'} = $shorthead[0];
		#$Sequences{$unid}{'gappy-ntseq'}   = uc($seq);       # store aligned sequence
		#$Sequences{$unid}{'SIZE'}    = length($seq);   # store length
		#$seq =~ s/\.//;
                #$seq =~ s/\-//;
                #$Sequences{$unid}{'degapped-ntseq'} = uc($seq);     # store degapped sequence
                $Sequences{$data[0]}{'filenumber'} = $filenumber;
                $Sequences{$data[0]}{'filename'} = $infile;
              #  $unid += 1;
	}
	$/="\n";
}
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
