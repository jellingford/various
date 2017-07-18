package tabix;

sub query{
    use strict;
    
    my ($file, $query) = @_;

    my @hits = qx(tabix $file $query);
    #print "$hits[0]\n";
    return (\@hits);
}

1;
