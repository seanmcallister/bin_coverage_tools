#!/usr/bin/perl -w
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Use idxstats and bins from DASTool to generate a total number of reads recruiting to each bin.

# - - - - - U S E R    V A R I A B L E S - - - - - - - -


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %options=();
getopts("b:i:x:o:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-b = path to bin folder\n";
        print "-x = bin file suffix (i.e. fasta, fa, fna, etc)\n";
	print "-i = idxstats file from samtools\n";
	print "-o = output prefix\n";
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
                $Bin_names{$file}{'total_reads'} = 0;
                &FASTAread($file, $fastaread);
                $fastaread += 1;
            }
    }

open(IN2, "<$options{i}") or die "\n\nNADA $options{i} you FOOL!!!\n\n";
my @idx_data = <IN2>; close(IN2);
foreach my $line (@idx_data)
    {   chomp($line);
        my @data = split('\t', $line);
        my $seq_header = $data[0];
	chomp($seq_header);
	if (exists $Sequences{$seq_header})
		{	$Sequences{$seq_header}{'READ_COUNT'} = $data[2];}
	else {print "This contig seems to be erroneous: $seq_header\n";}
	
    }

foreach my $i (sort keys %Bin_names)
    {   foreach my $j (sort keys %Sequences)
            {	my $test1 = $Sequences{$j}{'filename'};
		#print "Sequence file name: $test1\n";
		my $test2 = $Bin_names{$i}{'Name'};
		#print "Bin file name: $test2\n";
		if ($test1 eq $test2)
                    {   $Bin_names{$i}{'total_reads'} += $Sequences{$j}{'READ_COUNT'};
		    }                
            }
    }

open(OUT, ">".$options{o}."_readrecruitment.txt");
foreach my $m (sort keys %Bin_names)
    {   my $total_reads = $Bin_names{$m}{'total_reads'};
	print OUT "$Bin_names{$m}{'Name'}\t$total_reads\n";
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
	#my $unid = $filenumber.100000001;                            string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		#my $seq = '';
		#foreach my $i (1..$#data)
		#{	$seq .= $data[$i];  }
		#$seq =~ s/>//;
		$Sequences{$data[0]}{'HEAD'}    = $data[0];       # store header
		my @shorthead = split(' ', $data[0]);
		$Sequences{$data[0]}{'SHORT_HEAD'} = $shorthead[0];
		$Sequences{$data[0]}{'READ_COUNT'} = 0;
		#$Sequences{$unid}{'gappy-ntseq'}   = uc($seq);       # store aligned sequence
		#$Sequences{$unid}{'SIZE'}    = length($seq);   # store length
		#$seq =~ s/\.//;
                #$seq =~ s/\-//;
                #$Sequences{$unid}{'degapped-ntseq'} = uc($seq);     # store degapped sequence
                $Sequences{$data[0]}{'filenumber'} = $filenumber;
                $Sequences{$data[0]}{'filename'} = $infile;
                #$unid += 1;
	}
	$/="\n";
}
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
