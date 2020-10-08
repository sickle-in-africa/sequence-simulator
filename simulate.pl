#!/usr/bin/env perl

##
# Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
#
# This file is part of Bowtie 2.
#
# Bowtie 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Bowtie 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
#

#
# The paired-end data is made by (a) changing to the reads subdirectory and (b)
# running 'perl simulate.pl --ref=../reference/lambda_virus.fa'.
#

#
# The long-read data is made by (a) changing to the reads subdirectory and (b)
# running 'perl simulate.pl --ref=../reference/lambda_virus.fa --long
# --unpaired --prefix=longreads'.
#

use strict;
use warnings;
use FindBin 1.51 qw( $RealBin );
use lib $RealBin;
use Simulate;

sub main {

	my $mode = shift @ARGV or die "input ERROR: invalid simulate option";

	print "Running simulate in $mode mode", "\n";

	if( $mode eq "MutateReference" ){
		MutateReference();
	} elsif ( $mode eq "GenerateReads" ){
		GenerateReads();
	} else {
		die "input ERROR: invalid simulate option"
	}

}

sub MutateReference {

	# initialize input parameters
	my $inputs = Simulate::initialize_parameters();

	# load the reference sequence
	my $ref_string = Simulate::load_reference_sequence($inputs);

	# mutate the reference
	Simulate::mutate_reference($inputs, \$ref_string);

	# print the mutated reference
	Simulate::print_reference_sequence($inputs, \$ref_string);
}

sub GenerateReads {

	# initialize input parameters
	my $inputs = Simulate::initialize_parameters();

	# load the reference sequence
	my $ref_string = Simulate::load_reference_sequence($inputs);

	# generate the reads
	Simulate::generate_reads($inputs, \$ref_string);
}


#
#  run main
#
main();
