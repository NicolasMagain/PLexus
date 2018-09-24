#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;
use Smart::Comments '###';
use Getopt::Euclid;
use List::MoreUtils qw(uniq); 

my $infile = $ARGV{'-m'};
my $outfile = $ARGV{'-o'};

my $threshold = $ARGV{'-t'};

my @sequences;
my @identitysets;
my $nchar;
my %main_hash;

my %number_n;

my %newhash;
my %to_collapse;
my %hashfinal;



store_sequences_one($infile);
get_nchar_one($infile);
collapse_multi();



sub collapse_multi {

# toutes les lettres des séquences sont séparées, et un numéro leur est attribué
foreach my $seqid (sort keys %main_hash) {
	my @letters = split('',$main_hash{$seqid});
	#### @letters
	my %subhash;
	
	foreach my $n (1..scalar(@letters)) {
		$subhash{$n} = $letters[$n-1];
	}
	#### %subhash
	$newhash{$seqid}= \%subhash
}
#### %newhash


# chaque séquence est comparée aux autres
my @temparray = sort keys %main_hash;

foreach my $temp1 (0..$#temparray) {
	
	my $seqid1 = $temparray[$temp1];
	
	#on compte les N dans chaque séquence
	my $control_n=0;
	
	
	foreach my $n (1..$nchar) {
				my $character = $newhash{$seqid1}{$n};
				if ($character !~ /[AaCcGgTt-]/) {$control_n++}
				}
	$number_n{$seqid1}=$control_n;
}


#on prend la première séquence
ROUND:
foreach my $temp1 (0..$#temparray) {
	
	my $seqid1 = $temparray[$temp1];
	

	
	if ($to_collapse{$seqid1}){next ROUND}
	#say "seqid1" . $seqid1;
	
	# on prend la deuxième séquence
	my $seqid2;
	my @identitygroup;
	LINE:
	foreach my $temp2 (0..$#temparray){
	 	$seqid2 = $temparray[$temp2];
		#say "seqid2" . $seqid2;
		
		#on ne compare pas une séquence qui a déjà été comparée
		if ($temp2 <= $temp1) {next LINE}
		 
		my $control_ident=0;
		my $control_notempty=0;
		
		
		my $control_n1_non2=0;
		my $control_n2_non1=0;
		
		#on ne compare pas une séquence avec elle-même
		if ($seqid1 eq $seqid2) {
			next LINE}
		
		
		else{
			my @consensus;
			#on compare la séquence character par character
			foreach my $n (1..$nchar) {
				my $character1 = $newhash{$seqid1}{$n};
				my $character2 = $newhash{$seqid2}{$n};
				
				
				
				#si les séquences ont un caractère ACGT- différent à une position elles ne sont pas identiques on passe à la suite
				if ( ($character1 =~ /[AaCcGgTt-]/) and ($character2 =~ /[AaCcGgTt-]/) and ($character1 ne $character2))
					{
					#say "le taxon $seqid1 n'est pas égal au taxon $seqid2 car $seqid1 a $newhash{$seqid1}{$n} en position $n tandis que $seqid2 a $newhash{$seqid2}{$n}";
					#say "$seqid1 is not equal to $seqid2";
					next LINE
					}
					
				#si les séquences ont au moins une position ACGT identique, on ne compare pas que avec du missing, le compteur qui vérifie ça passe à 1. On ajoute aussi 1 au compteur d'identité	
				elsif (	($character1 =~ /[AaCcGgTt]/) and ($character2 =~ /[AaCcGgTt]/) and ($character1 eq $character2))
					{$control_notempty = 1;
					$control_ident ++}
				#sinon, un des deux caractères est ambigus, on laisse passer le compteur	
				elsif (	($character1 =~ /[AaCcGgTt-]/) and ($character2 !~ /[AaCcGgTt-]/))
					{
					$control_ident ++;
					$consensus[$n]=$character1;
					$control_n2_non1 ++;
					}
				elsif (	($character1 !~ /[AaCcGgTt-]/) and ($character2 =~ /[AaCcGgTt-]/))
					{
					$control_ident ++;
					$consensus[$n]=$character2;
					$control_n1_non2 ++;
					}
				else { $control_ident++;
						} 	
				
				}
				
				
				#si le compteur de vérification vaut un, et le compteur d'identité = la longueur de la séquence, on décide de collapser
				if (($control_ident == $nchar) and ($control_notempty == 1)) 
					{
						say "$seqid1 is identical to $seqid2";
						unless ($control_n1_non2 > $threshold and $control_n2_non1 > $threshold) {
							$to_collapse{$seqid1}=1;
							$to_collapse{$seqid2}=1;
							push (@identitygroup, $seqid1);
							push (@identitygroup, $seqid2);
							#### @identitygroup
						
							#on égalise les 2 séquences
							foreach my $n (1..$nchar){
								if ($consensus[$n]) {
									$newhash{$seqid1}{$n}=$consensus[$n];
									$newhash{$seqid2}{$n}=$consensus[$n];
									}
								}
							}
						else {say "it seems that $seqid1 and $seqid2 sequences have mutually exclusive sets of locis";
							next LINE
							}	
						}	
				else{
				say "$seqid1 is not identical to $seqid2"
					}	
		}
	
	
	}
	#on groupe les séquences identiques dans un array, lui même dans un array
	@identitygroup=uniq(@identitygroup);
	if (@identitygroup){
		push(@identitysets,\@identitygroup);
		}

}
#### %to_collapse
#### @identitysets

#on commence à préparer l'output
foreach my $secondlevel (@identitysets){
	
	#le nouveau nom de la séquence est le groupement des anciens noms
	my @seqid_sorted = sort {$number_n{$a} <=> $number_n{$b}} @$secondlevel ;
	say join(' ', @seqid_sorted);
	my $newseqid =join('',@seqid_sorted);
	foreach my $seq (@seqid_sorted) {say "in the seqid_sorted $seq $number_n{$seq}"};
	say "newseqid" . $newseqid;
	#### $newseqid
	
	#on doit choisir quelle séquence garder
	my $seqtokeep;				
	#my $grandbout='';
	my $min_so_far=$nchar;
	
	FIRST:
	foreach my $seqid (@$secondlevel) {
		say "seqid is $seqid";
		say "number of n: $number_n{$seqid}";
		#si la séquence est complète sans ambigus, on la prend
		if ($main_hash{$seqid} =~ /[AaCcGgTt]{$nchar}/) {$seqtokeep = $main_hash{$seqid};
														say "seqtokeep 1st condition is $seqid";
																			}
						#ou sans ambigu mais avec gaps
		elsif ($main_hash{$seqid} =~ /[AaCcGgTt-]{$nchar}/) {$seqtokeep = $main_hash{$seqid};
															say "seqtokeep 2nd condition is $seqid";
																			}
		#sinon on prend celle qui a le moins de N							
		else { if ($number_n{$seqid} < $min_so_far) { 
			say "seqid : $seqid";
			say "number_n $number_n{$seqid}";
			say "minsofar before: $min_so_far";
			$min_so_far = $number_n{$seqid};
			say "minsofar after : $min_so_far";
			$seqtokeep = $main_hash{$seqid};
			say "seqtokeep 3rd condition in the loop is $seqtokeep";
				}
	 		
	 		}
	 		
		}
$hashfinal{$newseqid}=$seqtokeep;	
	
	
}
#### %hashfinal

# on recompte le nombre de séquences après le collapse
my $nseq;
foreach my $sequence (sort keys %main_hash) {
	unless ($to_collapse{$sequence}) {$nseq++}
	}
foreach my $sequence (sort keys %hashfinal) {$nseq++}

#on prépare l'output
open (OUT, '>', $outfile);
say OUT "#NEXUS";
say OUT "BEGIN DATA;";
say OUT "\t DIMENSIONS  NTAX=$nseq NCHAR=$nchar;";
say OUT "\t FORMAT DATATYPE=DNA  MISSING=? GAP=- ;";
say OUT "MATRIX";

#on prend la séquence originale sauf si ça a été collapsé
foreach my $sequence (sort keys %main_hash) {
	unless ($to_collapse{$sequence}) {say OUT $sequence . " " . $main_hash{$sequence}}
}
#on ajoute les collapsées
foreach my $sequence (sort keys %hashfinal) {
	say OUT $sequence . " " . $hashfinal{$sequence}
}
say OUT ";";
say OUT "END;";
}

# store_sequences_one is a function that store nexus sequences from one locus in a hash named %main_hash
sub store_sequences_one {
	say "I am storing sequences for " . $_[0];
	my $inmatrix='';
	open my $in, '<', $_[0];

	LINE:
	while (my $line = <$in>) {
		chomp $line;

		$inmatrix = $inmatrix . " " . $line ;
		}
		$inmatrix =~ s/\[.+?\]//g;
		#### $inmatrix
		my ($temp1,$processed1) = split(/MATRIX/i, $inmatrix);
		#### $processed1
		my ($processed2, $temp2) = split(/;/,$processed1);
		#### $processed2
		$processed2 =~ s/^\s+|\s+$//g;
		%main_hash = split(/[\s]+/, $processed2);
		#### %main_hash
	
	
		foreach my $sequence (keys %main_hash){push(@sequences,$main_hash{$sequence})};
		#### @sequences
	
		return %main_hash ;
}

sub get_nchar_one {
	#say "I get the number of characters in ".  $_[0];
	
	#ouverture du fichier nexus
	open my $in, '<', $_[0];
	LINE:
	while (my $line = <$in>) {
		chomp $line;
		
		# si la ligne contient ne contient pas nchar, on passe à la suivante
		if ($line !~ /nchar/i ) {next LINE }
		
		# si elle contient nchar, on capture les chiffres contenus après et on le store dans nchar
		else {
			my @arguments = split /nchar=/i, $line;
			if  ($arguments[1] =~ m/(?<nchar>[0-9]+)/ ) {
				$nchar= $+{nchar};
				}
			  
			  }
		
	}
	#### $nchar	
}



=head1 NAME

collapse_multi.pl

=head1 USAGE

collapse_multi.pl -m <matrix> -o <outfile>

=head1 VERSION

Version 0.1
 
=head1 REQUIRED ARGUMENTS

=over

=item -m <matrix>

provide the loci files in the order you want them to be combined.

=for Euclid:
	matrix.type: readable

=head1 OPTIONS

=over

=item -o[ut[put|file]] <outfile>

Provide the name of the outfile. 

=for Euclid
	outfile.type: writeable
	outfile.default: 'combinedpart.nex'

=over

=item -t[hreshold] <threshold>

Provide the number of ambiguous characters that each sequence should have not to be collapsed (should be an approximation of the shorter sequence for one locus in the dataset)
If 0, sequences with ambiguous characters are not collapsed against sequences without ambiguous characters.

=for Euclid
	threshold.type: integer
	threshold.default: 200

=head1 COPYRIGHT
Nicolas Magain.