#!/usr/bin/perl
use strict;
use Getopt::Long;
use Cwd;

my $config;

($config) = &configure(scalar @ARGV);

main($config);

sub main{
    my $hash;
    ($hash)=make_gnomad_hash($config,$hash);
    ($hash)=make_vep_hash($config,$hash);
    compare_data($hash);
}
####################
sub compare_data{
    my ($h_compare) = @_;
    open (OUT, ">$config->{outputdir}/compared-gnomad-vep.txt");
    print OUT "Chr\tPos\tRef\tAlt\tID\tGene\tConsequence\tEffect\tTranscript\tVariant_Count\tVariant_Het_Count\tVariant_Hom_Count\tCount_Total_Gnomad\tCount_Hom_Gnomad\tCount_Hemi_Gnomad\tAF_Genomes_Gnomad\tAF_Exomes_Gnomad\tOld_ExAC_Freq\n";
    open (MISSED, ">$config->{outputdir}/compared-gnomad-vep-MISSED-VARIANTS.txt");
    open (MULTI, ">$config->{outputdir}/compared-gnomad-multi-annotations.txt");
    foreach my $variant (sort keys %{$h_compare->{gnomad}}){
	if ($h_compare->{vepNum}{$variant} > 1){
	    print MULTI "$variant\t$h_compare->{vepNum}{$variant}\n";
	}
	if (defined $h_compare->{vep}{$variant}){
	    foreach my $cDNA (keys %{$h_compare->{vep}{$variant}}){
		
		print OUT "$h_compare->{gnomadEssential}{$variant}\t$h_compare->{rs}{$variant}\t$h_compare->{gene}{$variant}\t$cDNA\t$h_compare->{effect}{$variant}\t$h_compare->{transcript}{$variant}\t$h_compare->{gnomadFrequencies}{$variant}\t$h_compare->{exac}{$variant}\n";
	    }
	    
	}
	else{
#	    print "ERROR: Some variants were not present in the filtered transcripts file - check these transcripts are present in the supplied filtered transcripts file\n\tmissed variants are present in >$config->{outputdir}/compared-gnomad-vep-MISSED-VARIANTS.txt\n\n";
	    print MISSED "$variant\n";
	}
    }

}
####################
sub make_gnomad_hash{
    my ($config,$h) = @_;
    open (GNO, $config->{gnomad}) or die "cannot open $config->{gnomad} : $!";
    while(<GNO>){
    chomp;
    my @d = split "\t";
    my ($c,$p,$ref,$alt)=($d[0],$d[1],$d[2],$d[3]);
    my $varGno = join("\t",$c,$p,$alt);
    my $varGnoPrint = join("\t",$c,$p,$ref,$alt);
    $h->{gnomad}{$varGno}++;
    $h->{gnomadEssential}{$varGno}=$varGnoPrint;
    my $freqs = join("\t",$d[4],$d[5],$d[6],$d[7],$d[8],$d[9],$d[10],$d[11]);
    $h->{gnomadFrequencies}{$varGno}=$freqs;
    }
    close GNO;
    return ($h);
}
###################
sub make_vep_hash{
    my ($config,$h) = @_;
    open (VEP, $config->{vep}) or die "cannot open $config->{vep} : $!";
    while(<VEP>){
    chomp;
    ## chr1:211654433 T NM_183059.2 intron_variant rs148178003 RD3 NM_183059.2:c.296+29G>A 0.003215
    my @dat = split "\t";
    my @CHROM=split ':', $dat[0];
    my ($c,$p,$alt,$trans,$eff,$id,$gene,$cDNA,$exac) = ($CHROM[0],$CHROM[1],$dat[1],$dat[2],$dat[3],$dat[4],$dat[5],$dat[6],$dat[7]);
    my $varVep = join("\t",$c,$p,$alt);
    $h->{vepNum}{$varVep}++;
    $h->{vep}{$varVep}{$cDNA}++;
    $h->{rs}{$varVep}=$id;
    $h->{gene}{$varVep}=$gene;
    $h->{exac}{$varVep}=$exac;
    $h->{transcript}{$varVep}=$trans;
    $h->{effect}{$varVep}=$eff;
    }
    close VEP;
    return ($h);
}
##################
sub configure{
    my ($config) = @_;
    my $args = shift;
    my $config = {};
    my @samples;
    $config = {'samples' => \@samples};
    GetOptions($config, "gnomad=s",
               "outputdir=s",
               "log=s",
               "vep=s",
               "help|h!",)
        || warn "error : $!\n";
    if (!defined $config->{gnomad}){
        print "\n ERROR:\n a varaints file with gnomad annotations needs to be provided with the --gnomad flag\n\n";
        usage();
    }
    if (!defined $config->{outputdir}){
        $config->{outputdir} = getcwd();
        print "using $config->{outputdir} as the output directory\n";
    }
    if (!defined $config->{vep}){
	print "\n ERROR:\n a varaints file with vep annotations needs to be provided with the --vep flag\n\n";
        usage();
    }
    return($config);
}
#####################
sub usage{
    die "\n\n Not using script correctly \n

Usage:
--gnomad
--outputdir
--vep
\n\n";

}
