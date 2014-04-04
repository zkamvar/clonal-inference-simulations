 P#!/usr/bin/perl

# written by Zhian Kamvar 2012-02-20

#----------------------------------------------------------#
# This program allows the user to run multiple simuPOP batch 
# files at once without necessitating any new input from the 
# command line.
#----------------------------------------------------------#

use strict;
use warnings;
my ($config_file, @list);

#----------------------------------------------------------#
# Obtains a list of configuration files from the current dir
#----------------------------------------------------------#

@list = `ls | grep clone | grep \"\\.cfg\"`;
my @burns = `ls | grep BURNIN | grep \"\\.cfg\"`;
my @burnpops = `ls | grep BURNIN | grep \"\\.pop\"`;
system "clear";
my $configs = scalar(@burns);
my $pops = scalar(@burnpops);
print "\n\t$configs BURNIN config files\n\t$pops BURNIN populations\n\n";
print "Would you like to run all files as batch? (yes/no)\n";
my $yes_no = <STDIN>;
my $batch;
if($configs != 0 && $pops == 0){
    foreach my $burn (@burns){
	print "Burning in $burn\n";
    	system "multi_burnin.py --config=$burn";
    }
} elsif($pops == 0){
    die("There are no burnin configuration files\n");
}
#----------------------------------------------------------#
# Choosing no for this option will require the user to give
# a response for each configuration file.
#----------------------------------------------------------#

if ($yes_no =~ /^n|N.+?/) {
	foreach $config_file (@list){
		chomp $config_file;
		if($config_file =~ /"^BURNIN.+?cfg"/){
			next;
		}
		print "Would you like to run $config_file as batch?\n";
		my $yes_no = <STDIN>;
		if ($yes_no =~ /^y|Y.+?/) {
			$batch = "Sim_Rand_Seed.py --config=$config_file --gui=False";
		} else {
			$batch = "Sim_Rand_Seed.py --config=$config_file";
		}
		system "$batch";
	}

#----------------------------------------------------------#
# Choosing yes will run all the files as a batch.
# Note: This is only useful if you have no duplicate names.
#----------------------------------------------------------#

} else {
	foreach $config_file (@list){
		chomp $config_file;
		if($config_file =~ /"^BURNIN.+?cfg"/){
			next;
		}
		print "\n\nRunning $config_file...\n\n";
		$batch = "Sim_Rand_Seed.py --config=$config_file --gui=False";
		system "$batch";
	}
}
