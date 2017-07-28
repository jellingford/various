#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;

my $config;

($config) = &configure(scalar @ARGV);

main($config);

sub main{
    my ($config) = @_;
    my $chrOI = $config->{chr};
#    my $filename = `basename include/$config->{variants}`;
    open (OUT,">$config->{outputdir}/variants.vcf_$chrOI");

    open (FILE, $config->{variants}) or die "cannot open $config->{variants} : $!";
    while(<FILE>){
	chomp;
	my @d = split "\t";
	my ($chr) = $d[0];
	if ($_ =~/^#/){
	    print OUT "$d[0]\t$d[1]\t$d[2]\t$d[3]\t$d[4]\n";
	}
	if ($chr eq $chrOI){
	    print OUT "$_\n";
	}
    }
}

#######################
sub configure{
    my ($config) = @_;
    my $args = shift;
    my $config = {};
    my @samples;
    $config = {'samples' => \@samples};
    GetOptions($config, "variants=s",
               "outputdir=s",
               "log=s",
               "chr=s",
               "help|h!",)
	|| warn "error : $!\n";
    if (!defined $config->{variants}){
	print "\n ERROR:\n a varaints file needs to be provided with the --variants flag\n\n";
        usage();
    }
    if (!defined $config->{outputdir}){
	$config->{outputdir} = getcwd();
        print "using $config->{outputdir} as the output directory\n";
    }
    if (!defined $config->{chr}){
        print "\n ERROR:\n a chromosome of interest needs to be provided with the --chr flag [this should match the chromosome syntac in the variants file, e.g. chr1 or 1...etc]\n\n";
        usage();
    }
    return($config);
}
##################
sub usage{
    die "\n\n Not using script correctly \n
Usage:
--variants
--outputdir
--chr
\n\n";
}
