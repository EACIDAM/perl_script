#! /usr/bin/perl

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Cds;






use vars qw($gb_file $region_file $help);
GetOptions("in=s" => \$gb_file ,				#input : genbank file
			"out=s" => \$region_file ,			#output : fasta file
			"h" => \$help) ;					# help

verif_options($gb_file, $region_file, $help);

				##############
				#### main ####
				##############

 
open(FILE, $gb_file) or die "Cannot open GenBank file";


my @tab_gene;


my $id;
my $acc;
my $desc;
my $debut;
my $fin;
my $seq="";
my $bool=0;
my @tab_gene;
my @tab_gene_trie;

while(<FILE>)
{

		
	chomp $_ ;

	my $ligne = $_ ;


	if($ligne =~ /^LOCUS\s+/)
	{
		my $tmp  = $' ;
		if($tmp =~ /(\w+)\s+/)
		{
			$id = $1 ;
			$desc = $' ;
		}
		next;
	}

	if($ligne =~/^ACCESSION\s+/)
	{
		$acc = $';
		next;
	}


	if($ligne =~ /^DEFINITION\s+/)
	{
		$desc .= "  $'";
		next;
	}


	if($ligne =~ /^\s+CDS\s+/)
	{
		my $tmp = $' ;

		if ($tmp =~ /\d+/)
		{

			$tmp =~ s/\>//g;
			$tmp =~ s/\<//g;


			if($tmp =~ /^complement\((\d+)\.\.(\d+)/)
			{
				$debut=$1;
				$fin=$2;
			}
			elsif($tmp =~ /(\d+)\.\.(\d+)/)
			{
				$debut=$1;
				$fin=$2;
			}


			my $description_cds="$id"."$desc"."| $debut..$fin";
			my $obj_cds = Cds -> new_1 ($description_cds, $debut, $fin);
			push(@tab_gene, $obj_cds);

		}
	}


	if($ligne =~ /\/\//)
	{
		my $taille = @tab_gene;
		my $obj = Bio::Seq -> new (
							- id => $acc ,
							- desc => $desc ,
							- seq => $seq ) ;


		if($taille>0)
		{

			@tab_gene_trie=tri_tableau(\@tab_gene);


			recuperation_region_intergenique(\@tab_gene_trie, $taille, $region_file, $obj);
		}
		


		$seq="";
		$bool=0;
		@tab_gene=();
		@tab_gene_trie=();
	}


	if($ligne =~ /^ORIGIN/)
	{
		$bool=1;
	}
	if($bool)
	{
		if($ligne =~ /\s+\d+/)
		{
			$seq .= $' ;
			$seq =~ s/ //g;
		}
	}
}

sub tri_tableau
{
	my ($tab_gene) = @_ ;


	my @tab = @$tab_gene;
	my @tab_trier;
	my $taille_tab= @tab; 
	
	my $i;
	my $u;
	my $v;


	for($i=0 ; $i<$taille_tab ; $i++)
	{
		my $bool=1;
		my $taille_tab_trier = @tab_trier;

		my $debut1 = $tab[$i] -> get_start();
		my $fin1 = $tab[$i] -> get_end();

		for($u=0 ; $u<$taille_tab_trier ; $u++)
		{
			my $debut2 = $tab_trier[$u] -> get_start();
			my $fin2 = $tab_trier[$u] -> get_end();
			if($debut1 <= $debut2)
			{
				for($v=$taille_tab_trier-1 ; $v>=$u ; $v--)
				{
					$tab_trier[$v+1] = $tab_trier[$v];
				}
				$tab_trier[$u] = $tab[$i];
				$bool=0;
			}
			if($bool==0)
			{
				last;
			}
		}


		if($bool)
		{
			my $test = $tab[$i]->get_start();
			push(@tab_trier, $tab[$i]);

		}
	}
	
	return @tab_trier;
}
sub recuperation_region_intergenique
{
	my($tab_cds, $taille, $output, $obj) = @_ ;


	my $desc= $obj -> desc();
	my $seq = $obj -> seq() ;
	my $id = $obj -> id() ;


	my @tab= @$tab_cds;


	my $taille_seq = length($seq);

	my $cds;
	my $start;
	my $end;
	my $start1;
	my $end1;
	my $start2;
	my $end2;


	my $out = Bio::SeqIO -> new (
							- file => ">>$output" ,
							- format => "fasta" ) ;


	if($taille==1)
	{

		$cds = $tab[0];
		$start = $cds -> get_start() ;
		$end = $cds -> get_end();
		


		my $region_amont = substr($seq, 0, $start);
		my $region_aval = substr($seq,$end);



		my $desc_amont = $desc." | 0..$start";
		my $desc_aval = $desc." | $end..$taille_seq";


		my $obj_amont = Bio::Seq -> new (
								- id => $id ,
								- desc => $desc_amont ,
								- seq => $region_amont ) ;
		my $obj_aval = Bio::Seq -> new (
								- id => $id, 
								- desc => $desc_aval ,
								- seq => $region_aval ) ;


		$out -> write_seq($obj_amont) if(length($region_amont)>20);
		$out -> write_seq($obj_aval) if(length($region_aval)>20);
	}
	else
	{

		$cds = $tab[0];
		$start1 = $cds -> get_start() ;
		$end1 = $cds -> get_end() ;


		my $region_amont = substr($seq, 0,$start1);


		my $desc_amont = $desc." | 0..$start1";


		my $obj_amont = Bio::Seq -> new (
								- id => $id ,
								- desc => $desc_amont ,
								- seq => $region_amont ) ;

		$out -> write_seq($obj_amont) if(length($region_amont)>20);



		for(my $i=1; $i<$taille; $i++)
		{

			$cds = $tab[$i];
			$start2 = $cds -> get_start() ;
			$end2 = $cds -> get_end() ;


			my $taille_region =	($start2 - $end1)+1;
			my $region_cds = substr($seq, $end1, $taille_region);


			my $desc_cds = $desc." | $end1..$start2";


			my $obj_cds = Bio::Seq -> new (
								- id => $id ,
								- desc => $desc_cds ,
								- seq => $region_cds ) ;



			if($start2<=$end1 && $end2<=$end1)
			{
				next;
			}





			if($start2>$end1)
			{
				$out -> write_seq($obj_cds) if(length($region_cds)>20);	
			}
			


			$start1 = $start2;
			$end1 = $end2 ;
		}


		my $region_aval = substr($seq, $end2);


		my $desc_aval = $desc." | $end2..$taille_seq";


		my $obj_aval = Bio::Seq -> new (
							- id => $id ,
							- desc => $desc_aval ,
							- seq => $region_aval ) ;


		$out -> write_seq($obj_aval) if(length($region_aval)>20);
	}
}

sub verif_options
{
	my ($gb_file, $region_file, $help) = @_ ;

	if($help)
	{
		help();
	}


	if(!$gb_file || !$region_file)
	{
		system('clear');
		print("Erreur. Syntaxe correct : region_intergenique -in genbank_file -out output\nregion_intergenique -h pour afficher l'aide.\n\n");
		exit;
	}

	if(-e $region_file)
	{
		system('clear');
		print("Attention le fichier $region_file fournis en argument existe dÃ©ja. Veuillez donner un autre nom.\nregion_intergenique -h pour afficher l'aide.\n\n");
		exit;
	}
}


sub help
{
	system('clear');

	print("
command line:  intergenic_extract.pl -in GenBank_file -out output_file
extraction of intergenic regions (in fasta format) from a genbank file with CDS features.
description line contain : accession number, nucleic sequence size
installation of the Cds.pm object: Copy the Cds.pm file in the directory containing perl classes (e.g. /usr/bin/perl/5.18) 
	
additional option 
	-h 		for help
");

	exit();
}
