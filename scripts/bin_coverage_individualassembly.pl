#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Use the metabat depth output (x.depth.txt) and bins from DASTool to generate an average coverage per bin.

# - - - - - U S E R    V A R I A B L E S - - - - - - - -


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %options=();
getopts("b:c:o:x:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-b = path to bin folder\n";
        print "-c = depth coverage file from metabat\n";
        print "-o = output prefix\n";
        print "-x = bin file suffix (i.e. fasta, fa, fna, etc)\n";
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
                $Bin_names{$file}{'total_cov'} = 0;
                $Bin_names{$file}{'total_len'} = 0;
                &FASTAread($file, $fastaread);
                $fastaread += 1;
            }
    }

open(IN2, "<$options{c}") or die "\n\nNADA $options{c} you FOOL!!!\n\n";
my @cov_data = <IN2>; close(IN2); shift (@cov_data);
foreach my $line (@cov_data)
    {   chomp($line);
        my @data = split('\t', $line);
        foreach my $i (sort keys %Sequences)
            {   if ($Sequences{$i}{'HEAD'} eq $data[0])
                    {   $Sequences{$i}{'COV'} = $data[2];
                        $Sequences{$i}{'COV_Len'} = $data[1];
                    }
            }
    }

foreach my $k (sort keys %Sequences)
    {   my $math_cov = $Sequences{$k}{'COV'};
        my $math_len = $Sequences{$k}{'COV_Len'};
        my $inflated_cov = $math_cov * $math_len;
        #print "$Sequences{$k}{'HEAD'}\t$math_cov\t$math_len\t$inflated_cov\n";
        $Sequences{$k}{'sizeinflatedcov'} = $inflated_cov;
    }

foreach my $i (sort keys %Bin_names)
    {   foreach my $j (sort keys %Sequences)
            {   if ($Sequences{$j}{'filename'} eq $Bin_names{$i}{'Name'})
                    {   $Bin_names{$i}{'total_cov'} = $Bin_names{$i}{'total_cov'} + $Sequences{$j}{'sizeinflatedcov'};
                        $Bin_names{$i}{'total_len'} = $Bin_names{$i}{'total_len'} + $Sequences{$j}{'COV_Len'};
                        #print "$Bin_names{$i}{'Name'}\t$Bin_names{$i}{'total_cov'}\t$Bin_names{$i}{'total_len'}\n";
                    }                
            }
    }

open(OUT, ">".$options{o}."_lenweighted_avgbincov.txt");
foreach my $m (sort keys %Bin_names)
    {   my $total_cov = $Bin_names{$m}{'total_cov'};
        my $total_len = $Bin_names{$m}{'total_len'};
        if ($total_len == 0)
            {print "Warning: some bins appear to have a 0 length:\n$Bin_names{$m}{'Name'}\n";}
        else
        {   my $avg_cov = $total_cov / $total_len;
            #print "$Bin_names{$m}{'Name'}\t$avg_cov\n";
            print OUT "$Bin_names{$m}{'Name'}\t$avg_cov\n";
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
	my $unid = $filenumber.10000001;                           # string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		#my $seq = '';
		#foreach my $i (1..$#data)
		#{	$seq .= $data[$i];  }
		#$seq =~ s/>//;
		$Sequences{$unid}{'HEAD'}    = $data[0];       # store header
		my @shorthead = split(' ', $data[0]);
		$Sequences{$unid}{'SHORT_HEAD'} = $shorthead[0];
		#$Sequences{$unid}{'gappy-ntseq'}   = uc($seq);       # store aligned sequence
		#$Sequences{$unid}{'SIZE'}    = length($seq);   # store length
		#$seq =~ s/\.//;
                #$seq =~ s/\-//;
                #$Sequences{$unid}{'degapped-ntseq'} = uc($seq);     # store degapped sequence
                $Sequences{$unid}{'filenumber'} = $filenumber;
                $Sequences{$unid}{'filename'} = $infile;
                $unid += 1;
	}
	$/="\n";
}
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
