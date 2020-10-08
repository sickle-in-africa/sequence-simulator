package Simulate;

use strict;
use warnings;

use Carp;
use Math::Random qw(random_normal random_exponential);
use Getopt::Long;
use List::Util qw(max min);
use lib qw(..);
use JSON::Parse 'parse_json';

#
#  functions
#
#	summary:
#		public:
#		+ initialize_parameters()
#		+ load_reference_sequence()
#		+ mutate_reference()
#		+ generate_reads()
#		private:
#		+ comp()
#		+ revcomp()
#		+ rand_dna()
#		+ rand_equals()
#		+ add_seq_errs()
#

sub initialize_parameters {

	# create reference to the input parameters hash 
	#	& specify default values

	my @fa_fn_ar = ();

	my %inputs = (
		fa_fn 		=> \@fa_fn_ar,
		ref_file 	=> "",		# input reference filename
		out_ref_file=> "",	# output reference filename
		rf 			=> "",		# reference sequence
		long 		=> 0,		# 1 -> generate long reads
		paired 		=> 1,		# 1 -> generate paired-end reads
		prefix 		=> "reads", # output files start with this string
		nreads 		=> 10000, 	# number of reads/end to generate
		rdlen_av 	=> 75,	# average to use when drawing from exponential\$inputs
		rdlen_exact	=> undef,	# exact length for all reads, overrides randomness
		rdlen_min 	=> 40,	# minimum read length (added to exponential draw)
		frag_av 	=> 250, 	# mean fragment len
		frag_sd 	=> 45,	# s.d. to use when drawing frag len from normal dist\$inputs
		verbose 	=> 0,		# be talkative

		input_json	=> 'none',	# input json filename for user-sepcified arguments

		has_snps   	=> 0,		# does the simulated sample genome have SNPS? true (1) /false (0)

		has_indels 	=> 0,		# does the simulated sample genome have indels? true (1) / false (0)
		has_rearr  	=> 0,		# does the simulated sample genome have rearrangements?

		psnp		=> 0.0012, 	# probability per base of a snp
		findel		=> 0.0005, 	# fraction of indels wrt original sequence length
		nrearr_av	=> 3,		# number of small dna rearrangements (an average value for a distribution)

		output_vcf	=> 'truth.g.vcf',	# output filename of truth vcf file
		cohort_id   => '_'
	);

	# update values with command line arguments
	GetOptions (
		"fasta|reference=s" => %inputs{fa_fn},
		"input_ref=s"		=> \$inputs{ref_file},
		"output_ref=s"		=> \$inputs{out_ref_file},
		"long"              => \$inputs{long},
		"verbose"           => \$inputs{verbose},
		"nreads=i"          => \$inputs{nreads},
		"read-avg=i"        => \$inputs{rdlen_av},
		"read-len=i"        => \$inputs{rdlen_exact},
		"read-min=i"        => \$inputs{rdlen_min},
		"frag-avg=i"        => \$inputs{frag_av},
		"frag-sd=i"         => \$inputs{frag_sd},
		"unpaired"          => sub { $inputs{paired} = 0; },
		"reads_prefix=s"    => \$inputs{prefix},
		"input_json=s" 		=> \$inputs{input_json},
		"output_vcf=s"		=> \$inputs{output_vcf},
		"cohort_id=s"		=> \$inputs{cohort_id}
	) || die "Bad option";

	# update values with input json file arguments
	my $from_json = parse_json(
		do {
		   open(my $json_fh, "<:encoding(UTF-8)", $inputs{input_json})
		      or die("Can't open \$inputs{input_json}\": $!\n");
		   local $/;
		   <$json_fh>
		}
	);
	$inputs{has_snps}   = $from_json->{has_snps}	if exists $from_json->{has_snps};
	$inputs{psnp}		= $from_json->{psnp}		if exists $from_json->{psnp};
	$inputs{has_indels} = $from_json->{has_indels}	if exists $from_json->{has_indels};
	$inputs{findel} 	= $from_json->{findel} 		if exists $from_json->{findel};
	$inputs{has_rearr}  = $from_json->{has_rearr} 	if exists $from_json->{has_rearr};
	$inputs{nrearr_av} 	= $from_json->{nrearr_av} 	if exists $from_json->{nrearr_av};
	$inputs{nreads}		= $from_json->{nreads}		if exists $from_json->{nreads};
	$inputs{rdlen_av}	= $from_json->{rdlen_av}	if exists $from_json->{rdlen_av};
	$inputs{rdlen_min}	= $from_json->{rdlen_min}	if exists $from_json->{rdlen_min};
	$inputs{cohort_id}	= $from_json->{cohort_id}	if exists $from_json->{cohort_id};

	if($inputs{long}) {
		$inputs{nreads}    = 6000;
		$inputs{rdlen_av}  = 300;
		$inputs{rdlen_min} = 40;
	}

	if ( scalar(@{$inputs{fa_fn}}) == 0 ) {
		push(@{$inputs{fa_fn}}, $inputs{ref_file});
	}

	# return reference to the input parameters hash
	return(\%inputs)
}


sub load_reference_sequence {
	my $inputs = $_[0];

	my $rf = "";

	scalar(@{$inputs->{fa_fn}}) > 0 || die "Must specify at least one reference FASTA file with --fasta";

	print STDERR "Loading reference files...\n";

	for my $fn (@{$inputs->{fa_fn}}) {
		if (! -e $fn) {
	    	print "simulate ERROR: the reference file $fn does not exist\n";
		}
		open(FN, $fn) || confess;
		my $name = "";
		while(<FN>) {
			chomp;
			$rf .= $_ if substr($_, 0, 1) ne ">";
		}
		close(FN);
	}

	return $rf;
}

sub print_reference_sequence {
	my $inputs = $_[0];
	my $rf    = $_[1];

	open(OREF, ">$inputs->{out_ref_file}") || die;

	print OREF ">NC_001416.1 $inputs->{cohort_id}, complete genome\n";

	print OREF "$_\n" for unpack '(A70)*', $$rf;

	close(OREF)

}

sub mutate_reference {
	my $inputs = $_[0];
	my $rf = $_[1];

	open(TVCF, ">$inputs->{output_vcf}") || die;

	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

	print TVCF 
		"##fileformat=VCFv4.3\n",
		"##fileDate=", sprintf("%02d%02d%02d\n", 1900+$year, $mon, $mday),
		"##source=simulate program\n",
		"##reference=$inputs->{ref_file}\n",
		'##FILTER=<ID=none,Description="No filters applied">',"\n",
		"##contig=<",
			"ID=NC_001416.1,",
			"length=", length($$rf), ",",
			"assembly=B36,",
			"md5=f126cdf8a6e0c7f379d618ff66beb2da,",
			"species=",'"Enterobacteria phage lambda"',",taxonomy=x",
		">\n",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

	my $nsnp = 0;
	if ($inputs->{has_snps}) {
		print STDERR "Adding single-base substitutions...\n";
		for(0..length($$rf)-1) {
			if(rand() < $inputs->{psnp}) {
				my $oldc = substr($$rf, $_, 1);
				substr($$rf, $_, 1) = substr("ACGT", int(rand(4)), 1);
				$nsnp++ if substr($$rf, $_, 1) ne $oldc;

				if (substr($$rf, $_, 1) ne $oldc) {
					my $chrom = "NC_001416.1";
					my $pos = $_ + 1;
					my $id  = ".";
					my $ref = $oldc;
					my $alt = substr($$rf, $_, 1);
					my $qual = "inf";
					my $filter = "none";
					my $info = ".";
					print TVCF "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\n";
				}

			}
		}
	} else {
		print STDERR "Skipping single-base substitutions...\n";
	}

	my $microgap = 0;
	if ($inputs->{has_indels}) {
		print STDERR "Adding microindels...\n";
		{
			my $newrf = "";
			my $nins = int(length($$rf) * $inputs->{findel} + 0.5);
			my $ndel = int(length($$rf) * $inputs->{findel} + 0.5);
			$microgap = $nins + $ndel;
			my %indel = ();
			for(1..$nins) {
				my $off = int(rand(length($$rf)));
				$indel{$off}{ty} = "ins";
				$indel{$off}{len} = int(random_exponential(1, 3))+1;
			}
			for(1..$ndel) {
				my $off = int(rand(length($$rf)));
				$indel{$off}{ty} = "del";
				$indel{$off}{len} = int(random_exponential(1, 3))+1;
			}
			my $lasti = 0;
			for my $rfi (sort {$a <=> $b} keys %indel) {
				if($rfi > $lasti) {
					$newrf .= substr($$rf, $lasti, $rfi - $lasti);
					$lasti = $rfi;
				}
				if($indel{$rfi}{ty} eq "ins") {
					# insert short DNA sequence
					my $chrom = "NC_001416.1";
					my $pos = $rfi;
					my $id  = ".";
					my $ref = substr($$rf, $pos-1, 1);
					my $alt = rand_dna($indel{$rfi}{len});
					my $qual = "inf";
					my $filter = "none";
					my $info = ".";
					$newrf .= $alt;
					$alt = $ref . $alt;
					print TVCF "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\n";
				} else {
					# make small deletion
					my $chrom = "NC_001416.1";
					my $pos = $rfi;
					my $id  = ".";
					my $ref = substr($$rf, $pos-1, $indel{$rfi}{len}+1);
					my $alt = substr($ref, 0, 1);
					my $qual = "inf";
					my $filter = "none";
					my $info = ".";
					$lasti += $indel{$rfi}{len};
					print TVCF "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\n"; 
				}
			}
			if($lasti < length($$rf)-1) {
				$newrf .= substr($$rf, $lasti, length($$rf) - $lasti - 1);
			}
			$$rf = $newrf;
		}
	} else {
		print STDERR "Skipping microindels...";
	}

	my $nrearr = 0;
	if ($inputs->{has_rearr}) {
		print STDERR "Adding larger rearrangements...\n";
		$nrearr = int(random_exponential(1, $inputs->{nrearr_av})+1);
		for(0..$nrearr) {
			my $break = int(rand(length($$rf)));
			my $before = substr($$rf, 0, $break);
			my $after = substr($$rf, $break);
			$after = revcomp($after) if int(rand()) == 0;
			$$rf = $after.$before;
		}
	} else {
		print STDERR "Skipping larger rearrangements...\n";
	}

	print STDERR "Added $nsnp SNPs\n";
	print STDERR "Added $microgap Microindels\n";
	print STDERR "Added $nrearr Rearrangements\n";

}


sub generate_reads {
	my $inputs = $_[0];
	my $rf     = $_[1];

	my @fraglens  = ();     # fragment lengths (for paired)
	my @readlens  = ();     # read/end lengths

	print STDERR "Picking read and fragment lengths...\n";
	# Pick random read lengths
	if(defined($inputs->{rdlen_exact})) {
		@readlens = ($inputs->{rdlen_exact}) x ($inputs->{nreads} * ($inputs->{paired} ? 2 : 1));
	} else {
		@readlens = random_exponential($inputs->{nreads} * ($inputs->{paired} ? 2 : 1), $inputs->{rdlen_av});
		@readlens = map int, @readlens;
		@readlens = map { int($_ + $inputs->{rdlen_min}) } @readlens;
	}
	if($inputs->{paired}) {
		# Pick random fragment and read lengths
		@fraglens = random_normal($inputs->{nreads}, $inputs->{frag_av}, $inputs->{frag_sd});
		@fraglens = map int, @fraglens;
		for(my $i = 0; $i < scalar(@readlens); $i += 2) {
			$fraglens[$i/2] = max($fraglens[$i/2], $readlens[$i] + $readlens[$i+1]);
		}
	}

	# Now simulate 
	print STDERR "Simulating reads...\n";
	my $rflen = length($$rf);
	if($inputs->{paired}) {
		open(RD1, ">$inputs->{prefix}_1P.fq") || die;
		open(RD2, ">$inputs->{prefix}_2P.fq") || die;
		for(my $i = 0; $i < scalar(@fraglens); $i++) {
			# Extract fragment
			my $flen = $fraglens[$i];
			my $off = int(rand($rflen - ($flen-1)));
			my $fstr = substr($$rf, $off, $flen);
			# Check if it has too many Ns
			my %ccnt = ();
			for my $j (1..$flen) {
				my $c = uc substr($fstr, $j, 1);
				$ccnt{tot}++;
				$ccnt{non_acgt}++ if ($c ne "A" && $c ne "C" && $c ne "G" && $c ne "T");
				$ccnt{$c}++;
			}
			# Skip if it has >10% Ns
			if(1.0 * $ccnt{non_acgt} / $ccnt{tot} > 0.10) {
				$i--;
				next;
			}
			# Possibly reverse complement
			$fstr = revcomp($fstr) if (int(rand(2)) == 0);
			# Get reads 1 and 2
			my $rdlen1 = min($readlens[2*$i], $flen);
			my $rdlen2 = min($readlens[2*$i+1], $flen);
			my $rd1 = substr($fstr, 0, $rdlen1);
			my $rd2 = substr($fstr, length($fstr)-$rdlen2);
			length($rd2) == $rdlen2 || die "Got ".length($rd2)." expected $rdlen2";
			# Reverse complement 2 to simulate --fr orientation
			$rd2 = revcomp($rd2);
			# Generate random quality values
			my $qu1 = rand_quals($rdlen1);
			$rd1 = add_seq_errs($rd1, $qu1);
			length($rd1) == length($qu1) || die;
			my $qu2 = rand_quals($rdlen2);
			$rd2 = add_seq_errs($rd2, $qu2);
			length($rd2) == length($qu2) || die;
			# Print
			print RD1 "\@r".($i+1)."\n$rd1\n+\n$qu1\n";
			print RD2 "\@r".($i+1)."\n$rd2\n+\n$qu2\n";
		}
		close(RD1);
		close(RD2);
		print STDERR "Made pairs: <base>_1P.fq and <base>_2P.fq where <base> is:\n $inputs->{prefix}\n";
	} else {
		open(RD1, ">$inputs->{prefix}.fq") || die;
		for(my $i = 0; $i < scalar(@readlens); $i++) {
			# Extract fragment
			my $rdlen = $readlens[$i];
			my $off = int(rand($$rflen - ($rdlen-1)));
			my $rd = substr($$rf, $off, $rdlen);
			# Check if it has too many Ns
			my %ccnt = ();
			for my $j (1..$rdlen) {
				my $c = uc substr($rd, $j, 1);
				$ccnt{tot}++;
				$ccnt{non_acgt}++ if ($c ne "A" && $c ne "C" && $c ne "G" && $c ne "T");
				$ccnt{$c}++;
			}
			# Skip if it has >10% Ns
			if(1.0 * $ccnt{non_acgt} / $ccnt{tot} > 0.10) {
				$i--;
				next;
			}
			length($rd) == $rdlen || die;
			# Possibly reverse complement
			$rd = revcomp($rd) if int(rand(2)) == 0;
			# Generate random quality values
			my $qu = rand_quals($rdlen);
			length($rd) == length($qu) || die "length(seq) = ".length($rd).", length(qual) = ".length($qu);
			$rd = add_seq_errs($rd, $qu);
			length($rd) == length($qu) || die "length(seq) = ".length($rd).", length(qual) = ".length($qu);
			# Print
			print RD1 "\@r".($i+1)."\n$rd\n+\n$qu\n";
		}
		close(RD1);
		print STDERR "Made unpaired reads: $inputs->{prefix}.fq\n";
	}

	print STDERR "DONE\n";
}

sub comp($) {
	my %revcompMap = (
	"A" => "T", "T" => "A", "a" => "t", "t" => "a",
	"C" => "G", "G" => "C", "c" => "g", "g" => "c",
	"R" => "Y", "Y" => "R", "r" => "y", "y" => "r",
	"M" => "K", "K" => "M", "m" => "k", "k" => "m",
	"S" => "S", "W" => "W", "s" => "s", "w" => "w",
	"B" => "V", "V" => "B", "b" => "v", "v" => "b",
	"H" => "D", "D" => "H", "h" => "d", "d" => "h",
	"N" => "N", "." => ".", "n" => "n" );

	my $ret = $revcompMap{$_[0]} || confess "Can't reverse-complement '$_[0]'";
	return $ret;
}

sub revcomp {
	my ($ret) = @_;
	$ret = reverse $ret;
	for(0..length($ret)-1) { substr($ret, $_, 1) = comp(substr($ret, $_, 1)); }
	return $ret;
}

sub rand_dna($) {
	my $ret = "";
	for(1..$_[0]) { $ret .= substr("ACGT", int(rand(4)), 1); }
	return $ret;
}

sub rand_quals($) {
	my $ret = "";
	my $upper = (rand() < 0.2 ? 11 : 40);
	$upper = 4 if rand() < 0.02;
	for(1..$_[0]) {
		$ret .= chr(33+int(rand($upper)));
	}
	return $ret;
}

sub add_seq_errs($$) {
	my($rd, $qu) = @_;
	my $origLen = length($rd);
	for(0..length($rd)-1) {
		my $c = substr($rd, $_, 1);
		my $q = substr($qu, $_, 1);
		$q = ord($q)-33;
		my $p = 10 ** (-0.1 * $q);
		if(rand() < $p) {
			$c = substr("ACGTNNNNNN", int(rand(10)), 1);
		}
		substr($rd, $_, 1) = $c;
		substr($qu, $_, 1) = $q;
	}
	length($rd) == $origLen || die;
	return $rd;
}

1;