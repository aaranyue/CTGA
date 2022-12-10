#!/usr/bin/perl
## perl ContigsToGenomeAssemble.pl high_quality_ref.fasta contigs_or_scaffolds.fasta polish_file_or_folder(or none)
## Author: Junyang Yue
use File::Copy;
use FileHandle;
use File::Basename;
use Cwd 'abs_path';
use threads;
use threads ('yield', 'stack_size' => 64*4096, 'exit' => 'threads_only', 'stringify');
use threads::shared;
use Getopt::Long;
use POSIX qw/mktime/;
use POSIX qw/strftime/;
use Time::HiRes qw(usleep);
use Time::HiRes qw(sleep time);
no strict 'refs';

my $localtime_start = `date "+%Y-%m-%d %H:%M:%S"`;
chomp ($localtime_start); #2022-01-20 14:58:58
my $job_start = &get_the_local_time();
my $worktime_start = $job_start -> {date_detail}; #20220120145858

our $thread = `cat /proc/cpuinfo | grep "processor" | wc -l`;
chomp ($thread);
our $memsize = `cat /proc/meminfo | grep MemTotal | grep -Eo [0-9]+`;
chomp ($memsize);
our $memory = int ($memsize/1200000);
chomp ($memory);

&checkInstalledSoft ("minimap2");
&checkInstalledSoft ("nucmer");

if (@ARGV != 3) {
	print "Usage: perl $0 high_quality_ref contigs_or_scaffolds polish_file_or_folder\n";
	exit(0);
}
my ($high_quality_ref, $contigs_or_scaffolds, $polish_file_or_folder) = @ARGV;

my $min_mummer_map = 50000; # default 50000
my $min_minimap_ratio = 0.4; # default 0.4
my $flanking_seq = 100000; # default 100000

our @suffix_list = qw(.gz .fastq .fq .fasta .fa .gff .gff3 .paf);
my $insert_n = "N"x100;
my $chr_numbers;

########## work_dir/ ##########
# 01_genome
#   |--[generated:] $prefix_high_quality_ref.fasta, $prefix_contigs_or_scaffolds.fasta
# 02_polish
#   |--[generated:] Hongyang.ont.fasta, Hongyang.hap1.fasta
# 03_split
#   |--[generated:] 
# 04_map
#   |--[generated:] 
# 05_mini
#   |--[generated:] 
# 06_concerned
#   |--[generated] 
# 07_adjust
#   |--[generated] 
# 08_temp
#   |--[temp] 
# 09_result
#   |--[generated]
###############################

my @provided_files = qw(01_genome 02_polish 03_split 04_map 05_mini 06_concerned 07_adjust 08_temp 09_result);

foreach (@provided_files) {
	chomp;
	if (! -e -d "$_") {
		system "mkdir $_";
	}
}

my $this_path = `pwd`;
chomp ($this_path);

my $abs_high_quality_ref = abs_path ($high_quality_ref);
my $base_high_quality_ref = basename ($high_quality_ref);
my $prefix_high_quality_ref = basename ($high_quality_ref, @suffix_list);
my $abs_contigs_or_scaffolds = abs_path ($contigs_or_scaffolds);
my $base_contigs_or_scaffolds = basename ($contigs_or_scaffolds);
my $prefix_contigs_or_scaffolds = basename ($contigs_or_scaffolds, @suffix_list);
my $abs_polish_file_or_folder = abs_path ($polish_file_or_folder);
my $base_polish_file_or_folder = basename ($polish_file_or_folder);

if ((-e -f "$abs_high_quality_ref") && (-e -f "$abs_contigs_or_scaffolds") && ($prefix_high_quality_ref ne $prefix_contigs_or_scaffolds)) {
	if (! -e -f "01_genome/$prefix_high_quality_ref.fasta") {
		system "ln -s $abs_high_quality_ref 01_genome/$prefix_high_quality_ref.fasta";
	}
	&convert_fasta_2_fasta("01_genome/$prefix_high_quality_ref.fasta");
	
	$chr_numbers = `cat 01_genome/$prefix_high_quality_ref.fasta | grep '^>' | wc -l`;
	chomp ($chr_numbers);
	
	if (! -e -f "01_genome/$prefix_contigs_or_scaffolds.fasta") {
		system "ln -s $abs_contigs_or_scaffolds 01_genome/$prefix_contigs_or_scaffolds.fasta";
	}
	&convert_fasta_2_fasta("01_genome/$prefix_contigs_or_scaffolds.fasta");
	
	if (-e -f "$abs_polish_file_or_folder") {
		print "The $abs_polish_file_or_folder is a file.\n";
		my $prefix_polish_file_or_folder = basename ($polish_file_or_folder, @suffix_list);
		
		if (! -e -f "02_polish/$prefix_polish_file_or_folder.fasta") {
			system "ln -s $abs_polish_file_or_folder 02_polish/$prefix_polish_file_or_folder.fasta";
		}
		&convert_fasta_2_fasta("02_polish/$prefix_polish_file_or_folder.fasta");
	} elsif (-e -d "$abs_polish_file_or_folder") {
		print "The $abs_polish_file_or_folder is a folder.\n";
		while (my $polish_seq_file = glob "$abs_polish_file_or_folder/*\.(fa|fasta)") {
			my $base_polish_seq_file = basename ($polish_seq_file);
			my $prefix_polish_seq_file = basename ($polish_seq_file, @suffix_list);
			if (! -e -f "02_polish/$prefix_polish_seq_file.fasta") {
				system "ln -s $abs_polish_file_or_folder/$base_polish_seq_file 02_polish/$prefix_polish_seq_file.fasta";
			}
			&convert_fasta_2_fasta("02_polish/$prefix_polish_seq_file.fasta");
		}
	}
	
	opendir (DIR, "02_polish") || die $!;
	my @files = readdir (DIR);
	closedir (DIR);
	my @files_select = grep /^\w/, @files;
	@files_select = sort @files_select;
	
	if (! -e -f "03_split/$prefix_contigs_or_scaffolds.split.fasta") {
		$/ = ">";
		open (IN, "01_genome/$prefix_contigs_or_scaffolds.fasta") || die $!;
		open (GAP, ">03_split/$prefix_contigs_or_scaffolds.split.gap.count");
		open (OUT, ">03_split/$prefix_contigs_or_scaffolds.split.fasta");
		open (LOG, ">03_split/$prefix_contigs_or_scaffolds.split.fasta.length");
		while (<IN>) {
			chomp;
			if ($_) {
				my @line = split m/\n/, $_;
				chomp ($line[0]);
				chomp ($line[1]);
				$line[1] =~ s/^N+//g;
				$line[1] =~ s/N+$//g;
				my @list = split m/N+/i, $line[1];
				print GAP $line[0]."\t".$#list."\n";
				for (my $i=0;$i<=$#list;$i++) {
					my $id = sprintf ("%03d", $i+1);
					print OUT ">".$line[0]."-".$id."\n".$list[$i]."\n";
					print LOG $line[0]."-".$id."\t".length($list[$i])."\n";
				}
			}
		}
		close (IN);
		close (GAP);
		close (OUT);
		close (LOG);
		$/ = "\n";
	}
	
	my $mummer1 = "mummer_".$prefix_high_quality_ref."_".$prefix_contigs_or_scaffolds.".split";
	if (! -e -f "03_split/$mummer1.pdf") {
		&run_mummer("$abs_high_quality_ref", "$this_path/03_split/$prefix_contigs_or_scaffolds.split.fasta");
		
		system "mv $mummer1.pdf 03_split/$mummer1.pdf"; # show
		system "mv $mummer1.filter 03_split/$mummer1.filter";
	}
	
	if (! -e -f "03_split/$mummer1.filter.mini") {
		my $title_line;
		open (IN, "03_split/$mummer1.filter") || die $!;
		open (OUT, ">03_split/$mummer1.filter.mini");
		while (<IN>) {
			chomp;
			if ($_ !~ m/(^\d+$)|(^\-\d+$)/) {
				if ($_ =~ m/^>/) {
					$title_line = $_;
				} else {
					my @line = split m/\s+/, $_;
					chomp ($line[0]);
					chomp ($line[1]);
					my $diff = $line[1] - $line[0];
					if ($diff > $min_mummer_map) {
						if ($title_line) {
							print OUT $title_line."\n";
							$title_line = "";
						}
						print OUT $_."\n";
					}
				}
			}
		}
		close (IN);
		close (OUT);
	}
	
	if (! -e -f "04_map/$prefix_contigs_or_scaffolds.map") {
		my %hash_count;
		my $label;
		my $i = 1;
		do {
			$label = 0;
			my $cut_value = $min_mummer_map * $i;
			
			$/ = ">";
			open (IN, "03_split/$mummer1.filter.mini") || die $!;
			open (OUT, ">04_map/$prefix_contigs_or_scaffolds.map.temp");
			while (<IN>) {
				chomp;
				if ($_) {
					my @line = split m/\n/, $_;
					chomp ($line[0]);
					my @id = split m/\s+/, $line[0];
					chomp ($id[0]);
					chomp ($id[1]);
					chomp ($id[2]);
					chomp ($id[3]);
					my ($sum, $sum1, $sum2);
					for (my $i=1;$i<=$#line;$i++) {
						chomp ($line[$i]);
						my @list = split m/\s+/, $line[$i];
						chomp ($list[0]);
						chomp ($list[1]);
						chomp ($list[2]);
						chomp ($list[3]);
						$sum += $list[0]; # for position
						my $diff = $list[1] - $list[0] + 1;
						if ($list[2] < $list[3]) {
							$sum1 += $diff;
						} else {
							$sum2 += $diff;
						}
					}
					if (($sum1 >= $cut_value) || ($sum2 >= $cut_value)) {
						my $position = int ($sum/$#line);
						if ($sum1 > $sum2) {
							print OUT $id[0]."\t".$id[1]."\t".$id[2]."\t".$id[3]."\t"."+"."\t".$position."\n";
							if (not exists $hash_count{$id[1]}) {
								$hash_count{$id[1]} = $id[0];
							} else {
								$label++;
							}
						} else {
							print OUT $id[0]."\t".$id[1]."\t".$id[2]."\t".$id[3]."\t"."-"."\t".$position."\n";
							if (not exists $hash_count{$id[1]}) {
								$hash_count{$id[1]} = $id[0];
							} else {
								$label++;
							}
						}
					}
				}
			}
			close (IN);
			close (OUT);
			$/ = "\n";
			
			$i++;
			undef %hash_count;
		} until ($label == 0);
		
		my $chr_map_numbers = `cut -f 1 04_map/$prefix_contigs_or_scaffolds.map.temp | sort | uniq | wc -l`;
		chomp ($chr_map_numbers);
		
		if ($chr_numbers != $chr_map_numbers) {
			print "Attention !!! Some chromosomes fail for mapping by any of the contigs. Please check the $mummer1.pdf in the folder \'03_split\'. Maybe a certain of contigs or the parameters \$min_mummer_map and/or \$flanking_seq need to be adjusted manually.\n";
			exit(0);
		}
		
		system "cd 04_map && cat $prefix_contigs_or_scaffolds.map.temp | sort -k1,1 -k6,6n > $prefix_contigs_or_scaffolds.map"; # revised
		system "rm -f 04_map/$prefix_contigs_or_scaffolds.map.temp";
	}
	
	## update the map configuration file
	if ((-e -f "04_map/$prefix_contigs_or_scaffolds.map") && (-e -f "04_map/$prefix_contigs_or_scaffolds.map.fasta")) {
		my $mtime_map = (stat "04_map/$prefix_contigs_or_scaffolds.map")[9];
		my $mdate_map = strftime "%Y%m%d%H%M%S", (localtime $mtime_map)[0..5];
		my $mdate_readable_map = strftime "%Y/%m/%d %H:%M:%S", (localtime $mtime_map)[0..5];
		
		my $mtime_fasta = (stat "04_map/$prefix_contigs_or_scaffolds.map.fasta")[9];
		my $mdate_fasta = strftime "%Y%m%d%H%M%S", (localtime $mtime_fasta)[0..5];
		my $mdate_readable_fasta = strftime "%Y/%m/%d %H:%M:%S", (localtime $mtime_fasta)[0..5];
		if ($mdate_map > $mdate_fasta + 10) {
			my $bak_file = "backup_".$worktime_start;
			system "mkdir -p backup/$bak_file";
			system "cd backup/$bak_file && mkdir 04_map 05_mini 06_concerned 07_adjust 09_result";
			system "mv 04_map/* backup/$bak_file/04_map/";
			system "mv 05_mini/* backup/$bak_file/05_mini/";
			system "mv 06_concerned/* backup/$bak_file/06_concerned/";
			system "mv 07_adjust/* backup/$bak_file/07_adjust/";
			system "mv 09_result/* backup/$bak_file/09_result/";
			system "cp backup/$bak_file/04_map/$prefix_contigs_or_scaffolds.map 04_map/";
			
			print "As $prefix_contigs_or_scaffolds.map.fasta \($mdate_readable_fasta\) was generated before $prefix_contigs_or_scaffolds.map \($mdate_readable_map\), we needs to rerun the pipeline.\n";
		}
	}
	
	if (! -e -f "04_map/$prefix_contigs_or_scaffolds.map.fasta") {
		my %hash_chr;
		my %hash_direaction1; # strand +
		my %hash_direaction2; # strand -
		my @all_new_id;
		my %hash_new_id;
		open (IN, "04_map/$prefix_contigs_or_scaffolds.map") || die $!;
		while (<IN>) {
			chomp;
			if ($_) {
				my @line = split m/\t/, $_;
				chomp ($line[0]);
				chomp ($line[1]);
				chomp ($line[4]);
				my $chr_id = $line[0];
				my ($hic_id, $hic_serial) = split m/-/, $line[1];
				
				my $new_id = $chr_id."#".$hic_id;
				if (not exists $hash_new_id{$new_id}) {
					push @all_new_id, $new_id;
					$hash_new_id{$new_id} = 1;
				}
				
				$hic_serial =~ s/^0//g;
				push @$new_id, $hic_serial;
				
				$hash_chr{$line[1]} = $line[0]; # no rep ids
				
				if ($line[4] eq "+") {
					if (exists $hash_direaction1{$new_id}) {
						$hash_direaction1{$new_id}++;
					} else {
						$hash_direaction1{$new_id} = 1;
					}
					if (not exists $hash_direaction2{$new_id}) {
						$hash_direaction2{$new_id} = 0;
					}
				} else {
					if (not exists $hash_direaction1{$new_id}) {
						$hash_direaction1{$new_id} = 0;
					}
					if (exists $hash_direaction2{$new_id}) {
						$hash_direaction2{$new_id}++;
					} else {
						$hash_direaction2{$new_id} = 1;
					}
				}
			}
		}
		close (IN);
		
		$/ = ">";
		my %hash_sequence; # split
		open (FA, "03_split/$prefix_contigs_or_scaffolds.split.fasta") || die $!;
		while (<FA>) {
			chomp;
			if ($_) {
				my @line = split m/\n/, $_;
				chomp ($line[0]);
				chomp ($line[1]);
				$hash_sequence{$line[0]} = $line[1];
			}
		}
		close (FA);
		$/ = "\n";
		
		my %hash_sequences; # map
		my %hash_count; # map
		my %hash_rep; # rep
		open (LOG, ">04_map/$prefix_contigs_or_scaffolds.map.log");
		foreach (@all_new_id) {
			my $this_id = $_;
			my ($chr_id, $hic_id) = split m/\#/, $this_id;
			my @each_serials = @$this_id;
			@each_serials = sort {$a <=> $b} @each_serials;
			my $start_serial = $each_serials[0] - 1;
			my $end_serial = $each_serials[$#each_serials] + 1;
			my $gap_count = 0;
			my $sequence;
			for ($start_serial .. $end_serial) {
				my $serial = sprintf ("%03d", $_);
				my $id = $hic_id."-".$serial;
				if (! $hash_sequence{$id}) {
					print LOG "The $id without sequences is not included\n";
				} elsif (($hash_chr{$id}) && ($hash_chr{$id} ne $chr_id)) {
					print LOG "The $id mapping to $hash_chr{$id} is not included\n";
				} elsif (! $sequence) {
					$sequence = $hash_sequence{$id};
					print LOG "The $id is included in $chr_id\n";
					if (not exists $hash_rep{$id}) {
						$hash_rep{$id} = $chr_id;
					} else {
						$hash_rep{$id} = $hash_rep{$id}."#".$chr_id;
					}
				} else {
					$sequence .= $insert_n.$hash_sequence{$id};
					print LOG "The $id is included in $chr_id\n";
					$gap_count++;
					if (not exists $hash_rep{$id}) {
						$hash_rep{$id} = $chr_id;
					} else {
						$hash_rep{$id} = $hash_rep{$id}."#".$chr_id;
					}
				}
			}
			
			if ($hash_direaction1{$this_id} < $hash_direaction2{$this_id}) { # strand -
				$sequence = reverse ($sequence);
				$sequence =~ tr/ACGTUacgtu/TGCAAtgcaa/;
			}
			
			if (not exists $hash_sequences{$chr_id}) {
				$hash_sequences{$chr_id} = $sequence;
				$hash_count{$chr_id} = $gap_count;
			} else {
				$hash_sequences{$chr_id} = $hash_sequences{$chr_id}.$insert_n.$sequence;
				$hash_count{$chr_id} = $hash_count{$chr_id} + 1 + $gap_count;
			}
			
		}
		
		foreach (sort keys %hash_rep) {
			if ($hash_rep{$_} =~ m/\#/) {
				print LOG "Attention : The $_ is included in $hash_rep{$_}\n";
			}
		}
		close (LOG);
		
		open (OUT, ">04_map/$prefix_contigs_or_scaffolds.map.fasta");
		open (GAP, ">04_map/$prefix_contigs_or_scaffolds.map.gap.count");
		foreach (sort keys %hash_sequences) {
			print OUT ">".$_."\n".$hash_sequences{$_}."\n";
			print GAP $_."\t".$hash_count{$_}."\n";
		}
		close (OUT);
		close (GAP);
		
		open (TMP, "03_split/$prefix_contigs_or_scaffolds.split.fasta.length") || die $!;
		open (LEN, ">04_map/$prefix_contigs_or_scaffolds.map.fasta.id");
		while (<TMP>) {
			chomp;
			if ($_) {
				my @line = split m/\t/, $_;
				chomp ($line[0]);
				chomp ($line[1]);
				if ($hash_chr{$line[0]}) {
					print LEN $line[0]."\t".$line[1]."\t".$hash_chr{$line[0]}."\n";
				} else {
					print LEN $line[0]."\t".$line[1]."\t"."Chr00"."\n";
				}
			}
		}
		close (TMP);
		close (LEN);
	}
	
	if (! -e -f "04_map/$prefix_contigs_or_scaffolds.map.gap.locus") {
		my %hash_map_gap;
		$/ = ">";
		open (FA, "04_map/$prefix_contigs_or_scaffolds.map.fasta") || die $!;
		while (<FA>) {
			chomp;
			if ($_) {
				my @line = split m/\n/, $_;
				chomp ($line[0]);
				chomp ($line[1]);
				$line[1] =~ s/^N+//g;
				$line[1] =~ s/N+$//g;
				if ($line[1] =~ m/N+/) {
					my @list = split m/N+/, $line[1];
					my $gaps_value;
					my $front_seq = 1;
					for (my $i=0;$i<$#list;$i++) {
						chomp ($list[$i]);
						my $gap_start = length ($list[$i]) + $front_seq;
						my $gap_end = $gap_start + 99;
						my $gap_value = $gap_start."-".$gap_end;
						if (! $gaps_value) {
							$gaps_value = $gap_value;
						} else {
							$gaps_value = $gaps_value."#".$gap_value;
						}
						$front_seq = $gap_end + 1;
					}
					$hash_map_gap{$line[0]} = $gaps_value;
				}
			}
		}
		close (FA);
		$/ = "\n";
		
		open (TMP, "04_map/$prefix_contigs_or_scaffolds.map.gap.count") || die $!;
		open (OUT, ">04_map/$prefix_contigs_or_scaffolds.map.gap.locus");
		while (<TMP>) {
			chomp;
			if ($_) {
				my @line = split m/\t/, $_;
				chomp ($line[0]);
				chomp ($line[1]);
				if ($hash_map_gap{$line[0]}) {
					print OUT $line[0]."\t".$line[1]."\t".$hash_map_gap{$line[0]}."\n";
				} else {
					print OUT $line[0]."\t".$line[1]."\t"."0"."\n";
				}
			}
		}
		close (TMP);
		close (OUT);
	}
	
	my $mummer2 = "mummer_".$prefix_high_quality_ref."_".$prefix_contigs_or_scaffolds.".map";
	if (! -e -f "04_map/$mummer2.pdf") {
		&run_mummer("$abs_high_quality_ref", "$this_path/04_map/$prefix_contigs_or_scaffolds.map.fasta");
		
		system "mv $mummer2.pdf 04_map/$mummer2.pdf"; # show
		system "mv $mummer2.filter 04_map/$mummer2.filter";
	}
	
	if (@files_select > 0) {
		if (! -e -f "05_mini/$prefix_contigs_or_scaffolds.mini.fasta") {
			my $insert_m = "M"x10;
			$/ = ">";
			open (IN, "04_map/$prefix_contigs_or_scaffolds.map.fasta") || die $!;
			open (OUT, ">05_mini/$prefix_contigs_or_scaffolds.mini.fasta");
			open (GAP, ">05_mini/$prefix_contigs_or_scaffolds.mini.gap.site");
			while (<IN>) {
				chomp;
				if ($_) {
					my @line = split m/\n/, $_;
					chomp ($line[0]);
					chomp ($line[1]);
					$line[1] =~ s/^N+//g;
					$line[1] =~ s/N+$//g;
					if ($line[1] =~ m/N+/i) {
						my $sequence;
						my @list = split m/N+/i, $line[1];
						for (my $i=0;$i<=$#list;$i++) {
							chomp ($list[$i]);
							my ($seq_up, $seq_down);
							if (length($list[$i]) >= $flanking_seq) {
								$seq_up = substr ($list[$i], 0, $flanking_seq);
								$seq_down = substr ($list[$i], 0-$flanking_seq, );
							} else {
								$seq_up = $list[$i];
								$seq_down = $list[$i];
							}
							if (! $sequence) {
								$sequence = $seq_up.$insert_m.$seq_down;
							} else {
								$sequence = $sequence.$insert_n.$seq_up.$insert_m.$seq_down;
							}
						}
						my @seqs = split m/M+/i, $sequence;
						#for (my $j=0;$j<=$#seqs;$j++) { # gaps + telomeres
						for (my $j=1;$j<$#seqs;$j++) { # only gaps
							chomp ($seqs[$j]);
							my ($seq_start, $seq_end) = split m/N+/i, $seqs[$j];
							my $site_start_n = length($seq_start);
							my $site_end_n = $site_start_n + 101; # Nx100
							my $serial = sprintf ("%02d", $j+1);
							print OUT ">".$line[0]."-".$serial."\n".$seqs[$j]."\n";
							print GAP $line[0]."-".$serial."\t".$site_start_n."#".$site_end_n."\n";
						}
					}
				}
			}
			close (IN);
			close (OUT);
			close (GAP);
			$/ = "\n";
		}
		
		for (my $i=0;$i<=$#files_select;$i++) {
			chomp ($files_select[$i]);
			my $prefix_query = basename ($files_select[$i], @suffix_list);
			my $prefix_paf = $prefix_contigs_or_scaffolds."-".$prefix_query;
			
			if (! -e -f "05_mini/$prefix_paf.paf") {
				system "minimap2 05_mini/$prefix_contigs_or_scaffolds.mini.fasta 02_polish/$prefix_query.fasta > 05_mini/$prefix_paf.paf";
			}
			
			if (! -e -f "05_mini/$prefix_paf.filter") {
				my %hash_site;
				open (SITE, "05_mini/$prefix_contigs_or_scaffolds.mini.gap.site") || die $!;
				while (<SITE>) {
					chomp;
					if ($_) {
						my @line = split m/\t/, $_;
						chomp ($line[0]);
						chomp ($line[1]);
						$hash_site{$line[0]} = $line[1];
					}
				}
				close (SITE);
				
				my %match_lable;
				my %match_seq;
				open (IN, "05_mini/$prefix_paf.paf") || die $!;
				open (OUT, ">05_mini/$prefix_paf.filter.temp");
				while (<IN>) {
					chomp;
					if ($_) {
						my @line = split m/\t/, $_;
						chomp ($line[0]); # query_id
						chomp ($line[1]);
						chomp ($line[2]); # query_start
						chomp ($line[3]); # query_end
						chomp ($line[4]);
						chomp ($line[5]); # chr_id
						chomp ($line[6]);
						chomp ($line[7]); # ref_start
						chomp ($line[8]); # ref_end
						chomp ($line[9]); # query_match_length
						chomp ($line[10]); # ref_match_length
						
						my ($nn_start, $nn_end) = split m/\#/, $hash_site{$line[5]};
						my $match_id = $line[0]."#".$line[5];
						my $ratio = $line[9]/$line[10];
						my $result_line = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[5]."\t".$line[6]."\t".$line[7]."\t".$line[8]."\t".$line[9]."\t".$line[10]."\t".$nn_start."\t".$nn_end;
						if (($ratio >= $min_minimap_ratio) && ($line[10] > 1000)) { # adjust
							if (($line[7] < $nn_start - 1000) && ($line[8] > $nn_end + 1000)) {
								my $query_matches = $line[3] - $line[2];
								my $ref_matches = $line[8] - $line[7];
								if ($query_matches > $ref_matches) {
									print OUT $result_line."\t"."Across"."\n";
								}
							} elsif (($line[7] > $nn_end) && ($line[7] - 1000 <= $nn_end) && ($ratio >= 0.8)) {
								$match_lable{$match_id} .= "D".$line[4]; # Downstream
								$match_seq{$match_id} .= $result_line."\t"."SegmentD"."\n";
							} elsif (($line[8] < $nn_start) && ($line[8] + 1000 >= $nn_start) && ($ratio >= 0.8)) {
								$match_lable{$match_id} .= "U".$line[4]; # Upstream
								$match_seq{$match_id} .= $result_line."\t"."SegmentU"."\n";
							}
						}
					}
				}
				close (IN);
				close (OUT);
				
				open (OUT, ">>05_mini/$prefix_paf.filter.temp");
				foreach (sort keys %match_lable) {
					if ((($match_lable{$_} =~ m/D\+/) && ($match_lable{$_} =~ m/U\+/)) || (($match_lable{$_} =~ m/D\-/) && ($match_lable{$_} =~ m/U\-/))) {
						print OUT $match_seq{$_};
					}
				}
				close (OUT);
				
				system "cd 05_mini && sort -k6,6 -k1,1 -k8,8n $prefix_paf.filter.temp > $prefix_paf.filter";
				system "rm -f 05_mini/$prefix_paf.filter.temp";
			}
			
			if (! -e -f "06_concerned/$prefix_paf.concerned.paf") {
				my %hash_ref_fasta;
				my %hash_query_fasta;
				$/ = ">";
				open (FA1, "04_map/$prefix_contigs_or_scaffolds.map.fasta") || die $!;
				while (<FA1>) {
					chomp;
					if ($_) {
						my @line = split m/\n/, $_;
						chomp ($line[0]);
						chomp ($line[1]);
						$hash_ref_fasta{$line[0]} = $line[1];
					}
				}
				close (FA1);
				
				open (FA2, "02_polish/$prefix_query.fasta") || die $!;
				while (<FA2>) {
					chomp;
					if ($_) {
						my @line = split m/\n/, $_;
						chomp ($line[0]);
						chomp ($line[1]);
						$hash_query_fasta{$line[0]} = $line[1];
					}
				}
				close (FA2);
				$/ = "\n";
				
				open (OUT1, ">06_concerned/$prefix_contigs_or_scaffolds.concerned.ref.fasta");
				my $ref_id = `cut -f 6 05_mini/$prefix_paf.filter | cut -d '-' -f 1 | sort | uniq`;
				chomp ($ref_id);
				my @ref_ids = split m/\n/, $ref_id;
				foreach (sort @ref_ids) {
					chomp;
					print OUT1 ">".$_."\n".$hash_ref_fasta{$_}."\n";
				}
				close (OUT1);
				
				open (OUT2, ">06_concerned/$prefix_query.concerned.query.fasta");
				my $query_id = `cut -f 1,5 05_mini/$prefix_paf.filter | sort | uniq`;
				chomp ($query_id);
				my @query_ids = split m/\n/, $query_id;
				foreach (sort @query_ids) {
					chomp;
					my ($id, $direction) = split m/\t/, $_;
					my $sequence = $hash_query_fasta{$id};
					if ($direction eq "-") {
						$sequence = reverse ($sequence);
						$sequence =~ tr/ACGTUacgtu/TGCAAtgcaa/;
					}
					print OUT2 ">".$id.$direction."\n".$sequence."\n";
				}
				close (OUT2);
				
				system "cd 06_concerned && minimap2 $prefix_contigs_or_scaffolds.concerned.ref.fasta $prefix_query.concerned.query.fasta > $prefix_paf.concerned.paf";
			}
			
			if (! -e -f "06_concerned/$prefix_paf.concerned.filter") {
				my %hash_id;
				my %hash_title;
				open (ID, "05_mini/$prefix_paf.filter") || die $!;
				while (<ID>) {
					chomp;
					if ($_) {
						my @line = split m/\t/, $_;
						chomp ($line[0]);
						chomp ($line[5]);
						chomp ($line[13]);
						my $id = $line[5];
						$id =~ s/\-\d+$//;
						if (not exists $hash_id{$line[0]}) {
							$hash_id{$line[0]} = $id;
						} else {
							$hash_id{$line[0]} .= "#".$id;
						}
						if (not exists $hash_title{$line[0]}) {
							$hash_title{$line[0]} = $line[13]
						} else {
							$hash_title{$line[0]} .= "#".$line[13];
						}
					}
				}
				close (ID);
				
				open (IN, "06_concerned/$prefix_paf.concerned.paf") || die $!;
				open (OUT, ">06_concerned/$prefix_paf.concerned.filter.temp");
				while (<IN>) {
					chomp;
					if ($_) {
						my @line = split m/\t/, $_;
						chomp ($line[0]);
						chomp ($line[5]);
						chomp ($line[9]);
						chomp ($line[10]);
						my $query_id = $line[0];
						$query_id =~ s/(\+|\-)$//;
						if (($hash_id{$query_id} =~ m/$line[5]/) && ($line[10] > $min_mummer_map) && ($line[9]/$line[10] >= $min_minimap_ratio)) {
							print OUT $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[5]."\t".$line[6]."\t".$line[7]."\t".$line[8]."\t".$line[9]."\t".$line[10]."\t".$hash_title{$query_id}."\n";
						}
					}
				}
				close (IN);
				close (OUT);
				
				system "cd 06_concerned && sort -k6,6 -k1,1 -k8,8n $prefix_paf.concerned.filter.temp > $prefix_paf.concerned.filter";
				system "rm -f 06_concerned/$prefix_paf.concerned.filter.temp";
			}
			
			my $mummer3 = "mummer_".$prefix_contigs_or_scaffolds.".concerned.ref_".$prefix_query.".concerned.query";
			if (! -e -f "06_concerned/$mummer3.pdf") {
				&run_mummer("$this_path/06_concerned/$prefix_contigs_or_scaffolds.concerned.ref.fasta", "$this_path/06_concerned/$prefix_query.concerned.query.fasta");
				
				system "mv $mummer3.pdf 06_concerned/$mummer3.pdf"; # show
				system "mv $mummer3.filter 06_concerned/$mummer3.filter";
				
				my $title_line;
				my $chr_id;
				my $i;
				open (IN, "06_concerned/$mummer3.filter") || die $!;
				open (OUT, ">06_concerned/$mummer3.filter.mini");
				while (<IN>) {
					chomp;
					if ($_ !~ m/(^\d+$)|(^\-\d+$)/) {
						if ($_ =~ m/^>/) {
							$title_line = $_;
							$title_line =~ s/^>//g;
							$i = 1; # need revised yuejy 
						} else {
							my @line = split m/\s+/, $_;
							chomp ($line[0]);
							chomp ($line[1]);
							my $diff = $line[1] - $line[0];
							if ($diff > $min_mummer_map) {
								if ($title_line) {
									my @title = split m/\s+/, $title_line;
									chomp ($title[0]);
									chomp ($title[1]);
									chomp ($title[2]);
									chomp ($title[3]);
									$chr_id = $title[0];
									my $this_site = $line[0] - 1; # for sorting
									print OUT ">"."\t".$chr_id."\t".$this_site."\t".$title[2]."\t".$title[3]."\t".$prefix_query."\#".$title[1]."\n";
									$title_line = "";
								}
								print OUT $i."\t".$chr_id."\t".$line[0]."\t".$line[1]."\t".$line[2]."\t".$line[3]."\n";
								$i++;
							}
						}
					}
				}
				close (IN);
				close (OUT);
				
				system "sort -k2,2 -k3,3n 06_concerned/$mummer3.filter.mini >> 07_adjust/$prefix_contigs_or_scaffolds.adjust.tab";
			}
		}
		
		if (! -e -f "07_adjust/$prefix_contigs_or_scaffolds.final.gap") {
			my %hash_gap;
			open (GAP, "04_map/$prefix_contigs_or_scaffolds.map.gap.locus") || die $!;
			while (<GAP>) {
				chomp;
				if ($_) {
					my @line = split m/\t/, $_;
					chomp ($line[0]);
					chomp ($line[2]);
					if ($line[2] != 0) {
						$hash_gap{$line[0]} = $line[2];
					}
				}
			}
			close (GAP);
			
			$/ = ">";
			open (IN, "07_adjust/$prefix_contigs_or_scaffolds.adjust.tab") || die $!;
			open (OUT, ">07_adjust/$prefix_contigs_or_scaffolds.final.gap.temp");
			while (<IN>) {
				chomp;
				if ($_) {
					my @line = split m/\n/, $_;
					chomp ($line[0]);
					my @value_id = split m/\t/, $line[0];
					chomp ($value_id[1]);
					chomp ($value_id[5]);
					my $chr_id = $value_id[1];
					my $query_id = $value_id[5];
					
					my $gaps = $hash_gap{$chr_id};
					my @sites = split m/\#/, $gaps;
					my $this_gaps = $gaps;
					$this_gaps =~ s/\-/\#/g;
					my @this_sites = split m/\#/, $this_gaps;
					
					my ($min_site, $max_site);
					chomp ($line[1]);
					my @value_1 = split m/\t/, $line[1];
					chomp ($value_1[2]);
					$min_site = $value_1[2];
					
					chomp ($line[$#line]);
					my @value_last = split m/\t/, $line[$#line];
					chomp ($value_last[3]);
					$max_site = $value_last[3];
					
					push @this_sites, $min_site;
					push @this_sites, $max_site;
					
					@this_sites = sort {$a <=> $b} @this_sites;
					my $new_sites = join ("-", @this_sites);
					
					$new_sites =~ s/^.+\-$min_site\-/\-/;
					$new_sites =~ s/^$min_site\-/\-/;
					$new_sites =~ s/\-$max_site\-.+$/\-/;
					$new_sites =~ s/\-$max_site/\-/;
					$new_sites =~ s/^\-//;
					$new_sites =~ s/\-$//;
					
					if ($new_sites =~ m/\-/) {
						$new_sites =~ s/\-.+\-/\-/;
						my @subs = split m/\-/, $new_sites;
						my $sub_start = $subs[0];
						my $sub_end = $subs[1];
						
						for (my $i=1;$i<=$#line;$i++) {
							chomp ($line[$i]);
							my @values = split m/\t/, $line[$i];
							chomp ($values[2]);
							chomp ($values[3]);
							push @subs, $values[2];
							push @subs, $values[3];
						}
						@subs = sort {$a <=> $b} @subs;
						my $new_subs = join ("-", @subs);
						
						my ($chr_start, $chr_end);
						if ($new_subs =~ m/(\d+)\-$sub_start\-/) {
							$chr_start = $1;
						}
						if ($new_subs =~ m/\-$sub_end\-(\d+)/) {
							$chr_end = $1;
						}
						
						my ($query_start, $query_end);
						for (my $i=1;$i<=$#line;$i++) {
							chomp ($line[$i]);
							my @values = split m/\t/, $line[$i];
							chomp ($values[2]);
							chomp ($values[3]);
							chomp ($values[4]);
							chomp ($values[5]);
							if ($chr_start == $values[2]) {
								$query_start = $values[4];
							} elsif ($chr_start == $values[3]) {
								$query_start = $values[5];
							}
							if ($chr_end == $values[2]) {
								$query_end = $values[4];
							} elsif ($chr_end == $values[3]) {
								$query_end = $values[5];
							}
						}
						print OUT $chr_id."\t".$chr_start."\t".$chr_end."\t".$query_id."\t".$query_start."\t".$query_end."\n";
					}
				}
			}
			close (IN);
			close (OUT);
			$/ = "\n";
			
			system "cd 07_adjust && sort -k1,1 -k2,2n $prefix_contigs_or_scaffolds.final.gap.temp > $prefix_contigs_or_scaffolds.final.gap";
			system "rm -f 07_adjust/$prefix_contigs_or_scaffolds.final.gap.temp";
		}
		
		if (! -e -f "07_adjust/$prefix_contigs_or_scaffolds.concerned.merge.query.fasta") {
			while (my $concerned_file = glob "06_concerned/*\.concerned\.query\.fasta") {
				my $base_concerned_file = basename ($concerned_file);
				my $this_file = $base_concerned_file;
				$this_file =~ s/\.concerned\.query\.fasta$//;
				
				open (IN, "06_concerned/$base_concerned_file") || die $!;
				open (OUT, ">>07_adjust/$prefix_contigs_or_scaffolds.concerned.merge.query.fasta");
				while (<IN>) {
					chomp;
					if ($_ =~ m/^>/) {
						my $id = $_;
						$id =~ s/^>//;
						my $new_id = $this_file."#".$id;
						print OUT ">".$new_id."\n";
					} else {
						print OUT $_."\n";
					}
				}
				close (IN);
				close (OUT);
			}
		}
		
		if (-e -f "07_adjust/$prefix_contigs_or_scaffolds.final.gap") { # revised
			my $head_line = `head -1 07_adjust/$prefix_contigs_or_scaffolds.final.gap | cut -f 4`;
			chomp ($head_line);
			if ($head_line !~ m/\#\#/) {
				my $raw_id_count = `cut -f 1 07_adjust/$prefix_contigs_or_scaffolds.final.gap | sort | uniq | wc -l`;
				chomp ($raw_id_count);
				my $new_id_count;
				do {
					system "cd 07_adjust && cp $prefix_contigs_or_scaffolds.final.gap $prefix_contigs_or_scaffolds.final.gap.$worktime_start";
					system "rm -f 07_adjust/$prefix_contigs_or_scaffolds.final.gap.temp";
					system "cd 07_adjust && mv $prefix_contigs_or_scaffolds.final.gap $prefix_contigs_or_scaffolds.final.gap.temp";
					
					my (%hash_start, %hash_end, %hash_merge);
					open (TEMP, "07_adjust/$prefix_contigs_or_scaffolds.final.gap.temp") || die $!;
					while (<TEMP>) {
						chomp;
						if ($_) {
							my @line = split m/\t/, $_;
							chomp ($line[0]);
							chomp ($line[1]);
							chomp ($line[2]);
							chomp ($line[3]);
							chomp ($line[4]);
							chomp ($line[5]);
							if (not exists $hash_start{$line[0]}) {
								my $chr_start = $line[1];
								my $chr_end = $line[2];
								my $query_start = $line[4];
								my $query_end = $line[5];
								
								$hash_start{$line[0]} = $chr_start;
								$hash_end{$line[0]} = $chr_end;
								
								$hash_merge{$line[0]} = $line[3]."##".$query_start."##".$query_end;
							} else {
								if ($line[2] <= $hash_end{$line[0]}) {
									#print "Nothing to do\n";
								} elsif (($line[1] <= $hash_end{$line[0]}) && ($line[2] > $hash_end{$line[0]})) {
									my $this_chr_start = $hash_end{$line[0]} + 1;
									my $chr_end = $line[2];
									my $this_query_start = $this_chr_start - $line[1] + $line[4];
									my $query_end = $line[5];
									
									$hash_end{$line[0]} = $chr_end;
									
									$hash_merge{$line[0]} = $hash_merge{$line[0]}."###".$line[3]."##".$this_query_start."##".$query_end;
								} elsif ($line[1] > $hash_end{$line[0]}) {
									my $this_chr_start = $hash_end{$line[0]} + 1; # cut and grab raw sequences from the $prefix_contigs_or_scaffolds.fasta
									my $this_chr_end = $line[1] - 1;
									my $this_query_id = $line[0];
									my $this_query_start = $this_chr_start;
									my $this_query_end = $this_chr_end;
									
									my $chr_start = $line[1];
									my $chr_end = $line[2];
									my $query_start = $line[4];
									my $query_end = $line[5];
									
									$hash_end{$line[0]} = $chr_end;
									
									$hash_merge{$line[0]} = $hash_merge{$line[0]}."###".$this_query_id."##".$this_query_start."##".$this_query_end."###".$line[3]."##".$query_start."##".$query_end;
								}
							}
						}
					}
					close (TEMP);
					
					open (OUT, ">07_adjust/$prefix_contigs_or_scaffolds.final.gap");
					foreach (sort keys %hash_merge) {
						print OUT $_."\t".$hash_start{$_}."\t".$hash_end{$_}."\t".$hash_merge{$_}."\n";
					}
					close (OUT);
					
					$new_id_count = `cut -f 1 07_adjust/$prefix_contigs_or_scaffolds.final.gap | sort | uniq | wc -l`;
					chomp ($new_id_count);
				} until ($raw_id_count == $new_id_count);
				
				system "rm -f 07_adjust/$prefix_contigs_or_scaffolds.final.gap.temp";
			}
		}
		
		if (! -e -f "09_result/$prefix_contigs_or_scaffolds.ctga.fasta") {
			$/ = ">";
			my %hash_raw_seq;
			open (FA1, "04_map/$prefix_contigs_or_scaffolds.map.fasta") || die $!;
			while (<FA1>) {
				chomp;
				if ($_) {
					my @line = split m/\n/, $_;
					chomp ($line[0]);
					chomp ($line[1]);
					$hash_raw_seq{$line[0]} = $line[1];
				}
			}
			close (FA1);
			
			my %hash_merge_seq;
			open (FA2, "07_adjust/$prefix_contigs_or_scaffolds.concerned.merge.query.fasta") || die $!;
			while (<FA2>) {
				chomp;
				if ($_) {
					my @line = split m/\n/, $_;
					chomp ($line[0]);
					chomp ($line[1]);
					$hash_merge_seq{$line[0]} = $line[1];
				}
			}
			close (FA2);
			$/ = "\n";
			
			my %hash_new_seq;
			open (TAB, "07_adjust/$prefix_contigs_or_scaffolds.final.gap") || die $!;
			while (<TAB>) {
				chomp;
				if ($_) {
					my @line = split m/\t/, $_;
					chomp ($line[0]);
					chomp ($line[1]);
					chomp ($line[2]);
					chomp ($line[3]);
					my @list = split m/###/, $line[3];
					my $seqs;
					for (my $i=0;$i<=$#list;$i++) {
						chomp ($list[$i]);
						my ($query_id, $query_start, $query_end) = split m/##/, $list[$i];
						if ($query_id =~ m/#/) {
							my $seq = substr ($hash_merge_seq{$query_id}, $query_start - 1, $query_end - $query_start + 1);
							$seqs .= $seq;
						} else {
							my $seq = substr ($hash_raw_seq{$query_id}, $query_start - 1, $query_end - $query_start + 1);
							$seqs .= $seq;
						}
					}
					my $chr_up = substr ($hash_raw_seq{$line[0]}, 0, $line[1] - 1);
					my $chr_down = substr ($hash_raw_seq{$line[0]}, $line[2], );
					my $sequence = $chr_up.$seqs.$chr_down;
					$hash_new_seq{$line[0]} = $sequence;
				}
			}
			close (TAB);
			
			open (OUT, ">09_result/$prefix_contigs_or_scaffolds.ctga.fasta");
			foreach (sort keys %hash_raw_seq) {
				chomp;
				if (! $hash_new_seq{$_}) {
					print OUT ">".$_."\n".$hash_raw_seq{$_}."\n";
				} else {
					print OUT ">".$_."\n".$hash_new_seq{$_}."\n";
				}
			}
			close (OUT);
		}
		
		my $mummer4 = "mummer_".$prefix_high_quality_ref."_".$prefix_contigs_or_scaffolds.".ctga";
		if (! -e -f "09_result/$mummer4.pdf") {
			&run_mummer("$abs_high_quality_ref", "$this_path/09_result/$prefix_contigs_or_scaffolds.ctga.fasta");
			
			system "mv $mummer4.pdf 09_result/$mummer4.pdf"; # show
			system "mv $mummer4.filter 09_result/$mummer4.filter";
		}
	} else {
		print "The polished files are not provided.\n";
	}
} else {
	print "The required files are not provided.\n";
	exit(0);
}

########## sub-program ##########
sub get_the_local_time {
	my $time = shift || time();
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
	$mon ++;
	$sec  = ($sec<10) ? "0".$sec:$sec;
	$min  = ($min<10) ? "0".$min:$min;
	$hour = ($hour<10) ? "0".$hour:$hour;
	$mday = ($mday<10) ? "0".$mday:$mday;
	$mon  = ($mon<10) ? "0".$mon:$mon;
	$year+=1900;
	
	my $weekday = ('Sun','Mon','Tue','Wed','Thu','Fri','Sat')[$wday];
	return {
		'second' => $sec,
		'minute' => $min,
		'hour'   => $hour,
		'day'    => $mday,
		'month'  => $mon,
		'year'   => $year,
		'weekNo' => $wday,
		'wday'   => $weekday,
		'yday'   => $yday,
		'date'   => "$year$mon$mday",
		'date_detail'   => "$year$mon$mday$hour$min$sec"
	};
}

sub checkInstalledSoft { ##(softname["name"]<string>)
	$/ = "\n";
	
	`$_[0] -version > soft_check.log 2>&1`;
	
	open (LOG, "soft_check.log");
	
	my @array;
	
	while (<LOG>) {
		chomp;
		if ($_) {
			push @array, $_;
		}
	}
	
	my $count = @array;
	my $softname = uc($_[0]);
	
	if (($count <= 2) && ("@array" =~ m/command not found/i)) {
		my $message = "\nError: The $softname software is not found. You should specify the path or install it firstly !\n\n";
		die $message;
	} else {
		print "The $softname software could be correctly called\. \n";
	}
	
	close (LOG);
	
	system "rm -f soft_check.log";
}

sub convert_fasta_2_fasta {
	if (@_ != 1) {
		print "Usage: perl $0 input_file\n";
		exit(0);
	}
	my ($input_file) = @_;
	
	my $input_file_seq_count = `cat $input_file | grep '>' | wc -l`;
	chomp ($input_file_seq_count);
	my $input_file_line_count = `cat $input_file | wc -l`;
	chomp ($input_file_line_count);
	
	if ($input_file_line_count > $input_file_seq_count * 2 + 1) {
		system "rm -f $input_file.temp";
		system "mv $input_file $input_file.temp";
		local $/ = "\n>";
		open (IN, "$input_file.temp") || die $!;
		open (OUT, ">$input_file");
		while (<IN>) {
			chomp;
			if ($_ =~ m/(\S+).*?\n(.+)/sx) {
				my $id = $1;
				my $seq = $2;
				$id =~ s/>//g;
				$seq =~ s/\r//g;
				$seq =~ s/\n//g;
				$seq =~ s/\*$//g;
				$seq =~ s/_$//g;
				my $sequence = uc ($seq);
				print OUT ">".$id."\n".$sequence."\n";
			}
		}
		close (IN);
		close (OUT);
		system "rm -f $input_file.temp";
	}
}

sub run_mummer {
	if (@_ != 2) {
		print "Usage: perl $0 species_a_name species_b_name\n";
		exit(0);
	}
	my ($species_a_name, $species_b_name) = @_;
	my $abs_species_a_name = abs_path ($species_a_name);
	my $base_species_a_name = basename ($species_a_name);
	my $prefix_species_a_name = basename ($species_a_name, @suffix_list);
	my $abs_species_b_name = abs_path ($species_b_name);
	my $base_species_b_name = basename ($species_b_name);
	my $prefix_species_b_name = basename ($species_b_name, @suffix_list);
	system "rm -rf mummer_temp";
	mkdir "mummer_temp";
	if (-e -f "$abs_species_a_name") {
		system ("ln -s $abs_species_a_name mummer_temp/$base_species_a_name");
	} else {
		die "Can not find $species_a_name: $!";
	}
	if (-e -f "$abs_species_b_name") {
		system ("ln -s $abs_species_b_name mummer_temp/$base_species_b_name");
	} else {
		die "Can not find $species_b_name: $!";
	}
	
	system "cd mummer_temp && nucmer -p nucmer -c 200 -g 200 $base_species_a_name $base_species_b_name"; #focus on similar
	system "cd mummer_temp && delta-filter -i 90 -l 2000 -m nucmer.delta > nucmer.filter";
	system "cd mummer_temp && mummerplot -p nplot nucmer.filter -t postscript";
	
	system "mv mummer_temp/nplot.gp mummer_temp/nplot.gp.temp";
	open (IN, "mummer_temp/nplot.gp.temp") || die "Can not open nplot.gp.temp: $!\n";
	open (OUT, ">mummer_temp/nplot.gp");
	while (<IN>) {
		chomp;
		if ($_) {
			my $line = $_;
			$line =~ s/w lp ls/w line ls/;
			print OUT $line."\n";
		}
	}
	close (IN);
	close (OUT);
	system "cd mummer_temp && gnuplot nplot.gp";
	system "cd mummer_temp && ps2pdf nplot.ps nplot.pdf";
	my $this_mummer = "mummer_".$prefix_species_a_name."_".$prefix_species_b_name;
	system "rm -f $this_mummer.filter $this_mummer.pdf";
	system "mv mummer_temp/nucmer.filter $this_mummer.filter";
	system "mv mummer_temp/nplot.pdf $this_mummer.pdf";
	system "rm -rf mummer_temp";
}
