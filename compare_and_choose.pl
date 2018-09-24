#!/usr/bin/env perl


# Feel free to contact me with questions or suggestions 
# Copyright Nicolas Magain njm18@duke.edu nicolas.magain@ulg.ac.be nicolas.magain@gmail.com

use Modern::Perl '2011';
use autodie;
use Smart::Comments '###';
use List::MoreUtils qw(uniq); 
use Getopt::Euclid;

my $threshold=$ARGV{'-t'};
my $outfile=$ARGV{'-o'};
my $format = $ARGV{'-f'};

my $loci0=$ARGV{'-l'};
my @loci;
for my $locus (@$loci0){
push (@loci, $locus)
}


my @nchars;
my @idlist;
my @locusnames;
my @selected;
my @storage;

my %main_hash_storage;


foreach my $locus (@loci){
	

	get_nchar_loop($locus);
	### @nchars

	store_sequences_loop($locus);
	#### %main_hash_storage
	
}	
update_ids();	
#### @idlist
unless($ARGV{'-c'} eq 'no'){compare()};
choose();
	






sub get_nchar_loop {
	say "I get the number of characters in ".  $_[0];
	open my $in, '<', $_[0];
	LINE:
	while (my $line = <$in>) {
		chomp $line;
		
		if ($line !~ /nchar/i ) {next LINE }
		
		else {
			my @arguments = split /nchar=/i, $line;
			if  ($arguments[1] =~ m/(?<nchar>[0-9]+)/ ) {
				push(@nchars, $+{nchar})
			}
		}
	}
	return(@nchars);
}


sub store_sequences_loop {
	my %main_hash;

	say "I am storing sequences for " . $_[0];
	my $inmatrix;
	open my $in, '<', $_[0];
	LINE:
	while (my $line = <$in>) {
		chomp $line;
		$inmatrix = $inmatrix . " " .  $line ;
		}
		
	$inmatrix =~ s/\[.+?\]//g;
	#### $inmatrix
	my ($temp1,$processed1) = split(/MATRIX/i, $inmatrix);
	#### $processed1
	my ($processed2, $temp2) = split(/;/,$processed1);

	$processed2 =~ s/^\s+|\s+$//g;

	%main_hash = split(/[\s]+/, $processed2);
	### %main_hash
	
	$main_hash_storage{$_[0]}=\%main_hash;
	
	#### @storage
	push(@storage, { %main_hash });
	#### @storage
	
	
	push(@locusnames,$_[0]);
	#### @locusnames
}

sub update_ids {
	say "I am updating the list of IDs";

	foreach my $hash (sort keys %main_hash_storage){	
		foreach my $secondlevel (sort keys $main_hash_storage{$hash})
			{push(@idlist, $secondlevel);
			@idlist = uniq @idlist;
			}
		}
}
	
sub sumchar {
	my $acc = 0;
	foreach (@_){
  	$acc += $_; }
	return($acc)
}
	
	
sub compare {
	say 'I compare the taxon IDs from ' . scalar @nchars . ' matrices.';
	say 'We have ' . scalar @nchars . " locus which contain @nchars characters, respectively";

	open(FILE, '>', $ARGV{'-c'}) || die "File not found";

	say FILE "ID \t" . join("\t",@locusnames) . "\t Number of loci";

	foreach my $seqid (sort @idlist) {
		my @presentabsent;
		my $presencenumber; 
		for my $locus_id (0 .. $#nchars) {
			if (exists ${ $storage[$locus_id] }{$seqid} ) 
				{ $presentabsent[$locus_id] = 'present';
				  $presencenumber ++;
				} else {$presentabsent[$locus_id] = 'absent'}
				
			}
		my @finalline = ($seqid, @presentabsent, $presencenumber);	
		say FILE join("\t", @finalline) ;
	}	
}
	
sub choose {
	say "I choose taxa that have at least $threshold loci out of the " .scalar(@loci). " loci ". join (" ", @loci);
	open(FILE, '>', $outfile) || die "File not found";
	say FILE '#NEXUS';
	say FILE 'BEGIN DATA;';

	foreach my $seqid (sort @idlist) {
		my @presentabsent;
		my $presencenumber; 
		for my $locus_id (0 .. $#nchars) {
			if (exists ${ $storage[$locus_id] }{$seqid} ) 
				{ 
				  $presencenumber ++;
				} 
			}
		if ($presencenumber >= $threshold)
		{push (@selected, $seqid);}
	}
		@selected = uniq @selected;

	say FILE 'DIMENSIONS NTAX='. scalar @selected . ' NCHAR=' . sumchar(@nchars) . ' ;';

	if ($format eq 'interleave') {
		say FILE 'FORMAT DATATYPE=DNA  MISSING=? GAP=-  interleave;';
		say FILE 'MATRIX';
		for my $locus_id (0 .. $#nchars) {
			foreach my $seqid (sort @selected) {
    		
        		if (exists ${ $storage[$locus_id] }{$seqid} ) 
    				{say FILE "$seqid \t ${ $storage[$locus_id] }{$seqid}";
    				}
    			else { say FILE $seqid . "\t " . 'N' x $nchars[$locus_id];
    				}
    		}
    	    	say FILE " " ;

    	}
	}
	 else {	say FILE 'FORMAT DATATYPE=DNA  MISSING=? GAP=-  interleave;';
		say FILE 'MATRIX';
		foreach my $seqid (sort @selected) {
			print FILE $seqid . "\t";
			for my $locus_id (0 .. $#nchars) {
				if (exists ${ $storage[$locus_id] }{$seqid} ) 
    				{print FILE "${ $storage[$locus_id] }{$seqid}";
    				}
    			else { print FILE 'N' x $nchars[$locus_id];
    				}
    			}
    		print FILE "\n";
    	}		
	}	
	say FILE ';';
	say FILE 'END;';
	close(FILE);	
}
	
	
=head1 NAME

Compare_and_choose.pl

=head1 USAGE


=head1 VERSION

Version 0.0
 
=head1 REQUIRED ARGUMENTS

=over

=item -l[oc[us|i]] <loci>...

provide the loci files in the order you want them to be combined.

=for Euclid:
	loci.type: readable

=item -t[hreshold] <threshold>

Minimum of loci of a taxa for it to be chosen.

=for Euclid:
	threshold.type: number



=head1 OPTIONS

=over 

=item -f[ormat] <format>

=for Euclid
	format.type: string
	format.default: 'interleave'

=over

=item -o[ut[put|file]] <outfile>

Provide the name of the outfile. 

=for Euclid
	outfile.type: writeable
	outfile.default: 'selectedtaxa.nex'

=item -c[omp[arison][file]] <comparisonfile>

Specify the name of the comparisonfile. If no, no file is produced.

=for Euclid
	comparisonfile.type: writeable
	comparisonfile.default: 'comparelocus.txt'

	
=head1 COPYRIGHT
Nicolas Magain.