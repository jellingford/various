#!/usr/bin/perl
use strict;
use lib '/Users/mmmskje2/Documents/git/perl_modules/';
use tabix;
use Getopt::Long;
use Cwd;
#use diagnostics;


my $config;

($config) = &configure(scalar @ARGV);

main($config);
sub main{
    my ($config) = @_;
    my ($hash,$gnomad);
    
    ### declare and store fullpaths to ALL gnomad files
    ($hash)=find_gnomad_file($hash);
    
    ### collect variant details into hash
    ### needs a brancher module based on input type (bespoke, vcf, vep..etc).
    my $file_var = $config->{variants}; ## e.g. "/Users/mmmskje2/Desktop/intronics/var-freq-RD.txt";
    ($hash)=create_variant_hash($hash,$file_var);
    
    ### define gnomad file for use and perform gnomad comparison
    ### splits analysis by type of variant and then perform comparison through subroutine compare_gnomad_bulk
    my $wd = $config->{outputdir}; ## e.g. "/Users/mmmskje2/Desktop/intronics/gnomad-output";
    ($hash)=compare_gnomad($hash,$wd);
    

    ### check analysis complete for all variants
    ## link to check_variant_numbers script

#($gnomad)=tabix::query($file,$query);
}
############
sub compare_gnomad{
    my ($h_compare,$wd)=@_;
    ### declare output files ###
        ### overall ###
    open (REFMATCH, ">$wd/refmatch.txt");
    open (NONREFMATCH, ">$wd/nonrefmatch.txt");
    open (GNOMADFREQ, ">$wd/gnomadfreq.txt");
    open (ABSENTGNOMAD, ">$wd/gnomadabsent.txt");
    open (COMBINED, ">$wd/combined-Input-Gnomad.txt");
    print NONREFMATCH "Chr\tPosition_Variant\tReference_Variant\tAlt_Variant\tPosition_Gnomad\tReference_Gnomad\tAlt_Gnomad\tPosition_Difference\n";
    print GNOMADFREQ "Chr\tPosition\tAlt_Variant\tVariant_Count\tGnomad_Count\tGnomad_AF\n";
    print ABSENTGNOMAD "Chr\tPosition\tAlt_Variant\tVariant_Count\tGnomad_Count\tGnomad_AF\n";
    
        ### INDELS ###
    open (INDELSMATCH, ">$wd/indels-refmatch.txt");
    open (OVERLAP, ">$wd/overlaps.txt");
    print INDELSMATCH "Chr\tPosition\tAlt_Variant\tVariant_Count\tGnomad_Count\tGnomad_AF\n"; print OVERLAP "Chr\tPosition\tAlt_Variant\tVariant_Ref>Alt\tCohort_Count\tn_Gnomad_AltAlleles\tGnomad_AltAllele\tGnomad_AF\n";
    #chr1    110146523       T       C>T     1       1       TGCCAGTCTCTACTAAAAG>T   0.0322955690479266
    ### loop through each variant ###
    foreach my $variant(sort keys %{$h_compare->{variant}}){
	#print "$variant\n";
	my @splitVar = split "\t", $variant;
	my ($ref,$alt,$pos,$chrOI,$rsID) = ($h_compare->{variant_wRef}{$variant},$splitVar[2],$splitVar[1],$splitVar[0],$splitVar[3]);
	my $gnomadOI = $h_compare->{gnomadFile_Chr}{$chrOI};
	if (!-e $gnomadOI) { warn "gnomad file doesn't exist : $!"}; ## if gnomad file doesn't exist then kill analysis
	my ($file,$query) = ($gnomadOI,$h_compare->{gnomad_query}{$variant});
	my $gnomad=tabix::query($file,$query); ## gnomad record retrieved, now needs analysis #### COMMON SOURCE OF ERROR ### make sure tabix is loaded ###
	my $nHits = scalar@{$gnomad}; ## must declare scalar for array before subroutine, bug made it inaccessible to access after entered subroutine compare_gnomad_bulk
	print "hits=$nHits\n";
	#### the following section is prepared to enable splitting of SNVs and indels when compared to gnomad data ###
	#### this is done on the basis of reference and alternate allele length ####
	#### currently erronously misses 1/2 genotypes ###
#	foreach my $jay (@{$gnomad}){ ## added
#	next unless ($nHits > 1);

	my ($length_var,$length_alt) = (length($ref),length($alt));
	
	### LOGIC 1 - SNVs ###
	if ($length_var == 1 && $length_var eq $length_alt){
#	if ($length_var == 1){
	my($af,$hit,$overlap,$over,$altAllele)=compare_gnomad_bulk($gnomad,$nHits,$length_var,$length_alt,$ref,$alt);
	    print "$variant\t$hit\t$af\n";
	    foreach my $fullLine (keys %{$h_compare->{variant_line}{$variant}}){
		if ($hit eq 'Y'){
		    print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t$af\n";
		    print COMBINED "$fullLine\t$af\n";
		}
		elsif($hit eq 'N'){
		    print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
		    print ABSENTGNOMAD "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
		    print COMBINED "$fullLine\t$af\t$ref\n";
		}
		if ($overlap eq 'Y'){
		    foreach my $o (keys %{$over->{AltAllele}}){
			print OVERLAP "$variant\t$ref>$alt\t$h_compare->{variant_count}{$variant}\t$nHits\t$o\t$over->{AltAllele}{$o}\n";
			
		    }
		}
	    }
	}
	
	### LOGIC 2 - multi-nucleotide substitutions ###
	elsif ($length_var > 1 && $length_var eq $length_alt){ 
	    my($af,$hit,$overlap,$over,$altAllele)=compare_gnomad_bulk($gnomad,$nHits,$length_var,$length_alt,$ref,$alt);
#	    print "$variant\t$hit\t$af\n";
	    if ($hit eq 'Y'){
                print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t$af\t$ref\n";
            }
            elsif($hit eq 'N'){
                print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
                print ABSENTGNOMAD "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
            }
	    if ($overlap eq 'Y'){
                foreach my $o (keys %{$over->{AltAllele}}){
                    print OVERLAP "$variant\t$ref>$alt\t$h_compare->{variant_count}{$variant}\t$nHits\t$o\t$over->{AltAllele}{$o}\n";
                }
            }

	}
	
	### LOGIC 3 - deletions (i.e. overall loss) ###
	elsif ($length_var > $length_alt){ 
	    my($af,$hit,$overlap,$over,$altAllele)=compare_gnomad_bulk($gnomad,$nHits,$length_var,$length_alt,$ref,$alt);
 #           print "$variant\t$hit\t$af\n";
	    if ($hit eq 'Y'){
                print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t$af\t$ref\n";
            }
            elsif($hit eq 'N'){
                print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
                print ABSENTGNOMAD "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
            }
	    if ($overlap eq 'Y'){
		foreach my $o (keys %{$over->{AltAllele}}){
		    print OVERLAP "$variant\t$ref>$alt\t$h_compare->{variant_count}{$variant}\t$nHits\t$o\t$over->{AltAllele}{$o}\n";
		}
	    }

	}

	### LOGIC 4 - insertions (i.e. overall gain) ###
	elsif ($length_var < $length_alt){
	    my($af,$hit,$overlap,$over,$altAllele)=compare_gnomad_bulk($gnomad,$nHits,$length_var,$length_alt,$ref,$alt);
  #          print "$variant\t$hit\t$af\n";
	    if ($hit eq 'Y'){
                print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t$af\t$ref\n";
            }
            elsif($hit eq 'N'){
                print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
                print ABSENTGNOMAD "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
            }
	    if ($overlap eq 'Y'){
		foreach my $o (keys %{$over->{AltAllele}}){
		    print OVERLAP "$variant\t$ref>$alt\t$h_compare->{variant_count}{$variant}\t$nHits\t$o\t$over->{AltAllele}{$o}\n";
		}
	    }
	    
	}

	### ANYTHING missed by logics 1-4 ###
	else{
	    my($af,$hit,$overlap,$over,$altAllele)=compare_gnomad_bulk($gnomad,$nHits,$length_var,$length_alt,$ref,$alt);
   #         print "$variant\t$hit\t$af\n";
	    if ($hit eq 'Y'){
                print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t$af\t$ref\n";
	    }
            elsif($hit eq 'N'){
                print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
                print ABSENTGNOMAD "$variant\t$h_compare->{variant_count}{$variant}\t0\t$ref\n";
            }
	    if ($overlap eq 'Y'){
                foreach my $o (keys %{$over->{AltAllele}}){
                    print INDELSOVERLAP "$variant\t$ref>$alt\t$h_compare->{variant_count}{$variant}\t$nHits\t$o\t$over->{AltAllele}{$o}\n";
                }
            }
	}
    }
 #   }
    return($h_compare);
}

############
sub compare_gnomad_bulk{
    my ($gnomad,$nHits,$length_var,$length_ref,$ref,$alt) = @_;
    my ($af,$hit,$overlap,$af_overlap)=(0,'N','N',0);
    my ($altDets,$overL);
    my @AFoverlap;
#    print "scalar(@gnomad)\n";
#    print "$ref\t$hit\t$af\n";
    for (my $i = 0; $i < $nHits; $i++){
	($overlap) = 'Y';
	my (%match,%overlap);
	my @ln = split ('\t', ${$gnomad}[$i]);
	my @alt_col = split (',',$ln[4]);
	my ($gnomad_ref) = ($ln[3]);
   ### SECTION START - define the number of nucleotides in reference ###
	my ($length_var,$length_alt) = (length($ref),length($alt));
	my $length_gnomad = length($ln[3]);
	#next unless ($length_var == $length_gnomad);
	    ### SECTION END ###
   ### SECTION START - check if reference of variant matches with reference of alternate ###
	    if ($ln[3] eq $ref){
		#print REFMATCH "$variant\n";
	    }
	    else{
#		my $pos_diff = ($pos - $ln[1]);
		#print NONREFMATCH "$chrOI\t$pos\t$ref\t$alt\t$ln[1]\t$ln[3]\t$ln[4]\t$pos_diff\n";
	    }
   ### SECTION END ###
   ### SECTION START - compare variant to each possible alternate allele ###
	    for (my $j = 0; $j < scalar@alt_col; $j++){
		print "alt_gnomad=$alt_col[$j]\n";
		if ($ln[3] eq $ref && $alt_col[$j] eq $alt){# && $ln[6]){ #ne 'AC_Adj0_Filter'){                                                                                                        
		    print "comparison: alt=$alt alt_gnomad=$alt_col[$j]\n";
		  #  my $af;
		    my @info = split(';', $ln[7]);
		    foreach my $inf (@info){
			my ($tag,$val) = ($inf =~ /(.*?)=(.*)/);
                        my @vals = split(',',$val);
                        $match{$tag} = $vals[$j];
		    }
		    unless ($match{AC} == '0' || $match{AN} == '0'){
			($af) = ($match{AC}/$match{AN})*100;
			#print GNOMADFREQ "$variant\t$h_compare->{variant_count}{$variant}\t$data{AC}\t$af\n";
			($hit) = 'Y';
			#print "$af\t$hit\n";
		    }
		}
		else{
		    #print "#### overlap ###";
		    my $altVariant = $alt_col[$j];
		    ($altDets) = "$ln[3]>$altVariant";
		  #  push(@AltVar, $altDets);
		    my @info = split(';', $ln[7]);
		    foreach my $inf (@info){
                        my ($tag,$val) = ($inf =~ /(.*?)=(.*)/);
                        my @vals = split(',',$val);
                        $overlap{$tag} = $vals[$j];
                    }
                    unless ($overlap{AC} == '0' || $overlap{AN} == '0'){
                        ($af_overlap) = ($overlap{AC}/$overlap{AN})*100;
			($overlap) = 'Y';
			$overL->{AltAllele}{$altDets}=$af_overlap;
#			$overlap{AltFreq_Custom}=$af_overlap;
			push(@AFoverlap,$af_overlap);
		    }
		    
		}
	    }
    }
    return($af,$hit,$overlap,$overL,$altDets);
    #return($af,$hit,$overlap,$af_overlap,$altDets);
}

############
sub create_variant_hash{
    my ($h_var,$file1)=@_;
    open (ORG, $file1) or die "cannot open $file1 : $!";
	while(<ORG>){
	    #if ($config->{type} eq 'vcf'){
	    chomp;
	    # chr1    10032272        .       T       A       3
	    next if ($_ =~/^#/);
	    my @d = split "\t";
	    #my ($c,$p,$id,$ref,$alt,$count) = ($d[0],$d[1],$d[2],$d[3],$d[4],$d[5]);
	    my ($c,$p,$id,$ref,$alt,$count);
	    if ($config->{type} eq 'vcf'){
		($c,$p,$id,$ref,$alt,$count) = ($d[0],$d[1],$d[2],$d[3],$d[4],$d[5]);
	    }
	    elsif ($config->{type} eq 'vep'){
		my @pos = split ':', $d[1];
		($c,$p,$id,$ref,$alt,$count) = ($pos[0],$pos[1],'NA','NA',$d[2],'NA');
		#print "#$c#\t#$p#\t#$alt#\n";
	    }
	    my $variant = join("\t",$c,$p,$alt,$id);
	    $h_var->{variant}{$variant}++;
	    $h_var->{variant_wRef}{$variant}=$ref;
	    my $position = join("\t",$c,$p);
	    my $c_NoChr = $c; ($c_NoChr=~s/chr//);
	    my $gnomad_query = "$c_NoChr:$p-$p";
	    $h_var->{variant_position}{$position}++;
	    $h_var->{variant_chr}{$c}{$variant}++;
	    $h_var->{gnomad_query}{$variant}=$gnomad_query;
	    $h_var->{variant_count}{$variant}=$count;
	    $h_var->{variant_line}{$variant}{$_}++;
	}
	close ORG;
    return($h_var);
}

############
sub find_gnomad_file{
    my ($h_gnomad) = @_;
    my $gnomadPaths = "/Users/mmmskje2/Documents/gnomad/data/fullpaths-gnomad.txt";
    open (GPATH, $gnomadPaths) or die "cannot open $gnomadPaths : $!";
    while(<GPATH>){
	chomp;
	my @d = split "\t";
	my ($chr,$path) = ($d[0],$d[1]);
	#print "#$chr#\t#$path#\n";
	$h_gnomad->{gnomadFile_noChr}{$chr}=$path;
	$h_gnomad->{gnomadFile_Chr}{"chr$chr"}=$path;
    }
    close GPATH;
    return($h_gnomad);
}
##############
sub configure{
    my ($config) = @_;
    my $args = shift;
    my $config = {};
    my @samples;
    $config = {'samples' => \@samples};
    GetOptions($config, "variants=s",
	       "outputdir=s",
	       "log=s",
	       "type=s",
	       "help|h!",)
	|| warn "error : $!\n";
    if (!defined $config->{variants}){
	print "\n ERROR:\n a variants file needs to be provided with the --variants flag\n\n";
	usage();
    }
    if (!defined $config->{outputdir}){
        $config->{outputdir} = getcwd();
	print "using $config->{outputdir}\n";
    }
    if (!defined $config->{log}){
	$config->{log} = "$config->{outputdir}.log";
    }
    if (!defined $config->{type}){
        print "\n ERROR:\n a type of input in the variant file needs to be provided with the --type flag\n\n";
        usage();
    }
    return($config);
}
sub usage{
    die "\n\n Not using script correctly \n

Usage:
--variants
--outputdir
--log
--type    vcf/vep
\n\n";
}

