#!/usr/bin/perl
use strict;

my $hash;

main();

sub main{
    ## define files
    my $org_file = "var-freq-RD.txt";
    my $vep_file = $ARGV[0];

    ## make hash from original uploaded variants
    ($hash)=original_variants($hash,$org_file);
    
    ## make hash from vep returned results
    ($hash)=vep_results($hash,$vep_file);
    
    ## compare hashes and get results
    comparison($hash);

}

#############
sub comparison{
    my ($compare)=@_;
    
### compare hashes at variant level (chr, position, alt)
    my ($count_match,$count_missing,$count_org) = (0,0,0);
    foreach my $variant (sort keys %{$compare->{original_variant}}){
	$count_org++;
	if (defined $compare->{vep_variant}{$variant}){
	    $count_match++;
	}
	else{
	    $count_missing++;
	}
    }
    print "\n\ntotal variants uploaded = $count_org\nvariants found = $count_match\nvariants missing = $count_missing\n\n";

### compare hashes at position level (chr, position, alt)
=head
    my ($count_match,$count_missing,$count_org) = (0,0,0);
    foreach my $position (sort keys %{$compare->{original_position}}){
        $count_org++;
        if (defined $compare->{vep_position}{$position}){
            $count_match++;
        }
        else{
            $count_missing++;
        }
    }
    print "\n\ntotal variants uploaded = $count_org\nvariants found = $count_match\nvariants missing = $count_missing\n\n";
=cut

}

#############
sub original_variants{
    my ($h,$file1) = @_;
    open (ORG, $file1) or die "cannot open $file1 : $!"; 
    while(<ORG>){
	chomp;
    # chr1    10032272        .       T       A       3 
	## could add to process vcf format automatically?
	## could add a check to make sure on fwd strand?

	my @d = split "\t";
	my ($c,$p,$id,$ref,$alt,$count) = ($d[0],$d[1],$d[2],$d[3],$d[4],$d[5]);
	my $variant = join("\t",$c,$p,$alt);
	$h->{original_variant}{$variant}++;
	my $position = join("\t",$c,$p);
	$h->{original_position}{$position}++;
    }
    close ORG;
    return($h);
}

#############

sub vep_results{
    my ($h2,$file2) = @_;
    open (VEP, $file2) or die "cannot open $file2 : $!";
    while(<VEP>){
	chomp;
#   .       chr1:10032272   A
	my @data = split "\t";
	my @pos = split ':', $data[1];
	my ($c1,$p1,$alt1)=($pos[0],$pos[1],$data[2]);
	my $variant_vep = join ("\t",$c1,$p1,$alt1);
	$h2->{vep_variant}{$variant_vep}++;
	my $position_vep = join ("\t",$c1,$p1);
	$h2->{vep_position}{$position_vep}++;
    }
    close VEP;
    return($h2);
}

############
