#!/usr/bin/perl

use strict;
use warnings;

#use Genome;
use IO::File;
use Data::Dumper;

print_usage() and exit unless @ARGV;

#print "USAGE: $usage\n" and exit unless @ARGV;
`echo $0 @ARGV > lumpy_vcf_filter.log`;

# remove existing files from old runs to avoid potential confusion with earlier attempts with this script
&clean_up_old_output_files(); #exit;

my $input_params = parse_input_params();
#print Dumper $input_params; exit;

# TODO change to column_number_for_header
my %header_column_number; set_header_column_number();
#print Dumper \%header_column_number; exit;

# get individual names from #CHROM line in vcf
my( $ALL_INDIVIDUALS, $TARGET_INDIVIDUALS, $NON_TARGET_INDIVIDUALS ) = validated_individuals_to_process();
#print Dumper $ALL_INDIVIDUALS, $TARGET_INDIVIDUALS, $NON_TARGET_INDIVIDUALS; exit;

# chromosomes/refs to exclude
my $CHRS_TO_EXCLUDE = set_chromosomes_to_exclude();
# print Dumper $CHRS_TO_EXCLUDE; exit;

print "PARSING ".$input_params->{vcf_file}."\n";
my $vcf_fh = input_vcf_filehandle();
my $vcf_out_fh = output_vcf_filehandle();
while( my $line = $vcf_fh->getline ) {
    next if $line =~ /^#/; # skip headers

    my $sv_info = get_sv_info_from_line( $line );
    #print Dumper $sv_info; last;

    # skip all BND/INTERCHROM SV types
    next if record_is_BND_sv( $sv_info );

    add_to_later_summarize( $sv_info, 'PRE_FILTERED' );
    add_sv_for_stats( $sv_info, 'PRE_FILTERED' );

    # user input filtering
    # chrs
    next if sv_is_of_excluded_chromosome( $sv_info );
    # min/max length
    next unless sv_min_length_is_met( $sv_info );
    next unless sv_max_length_is_met( $sv_info );
    # sv type DEL, INV etc
    next unless sv_matches_input_sv_type( $sv_info );

    # sv length and type doesn't change when
    # buffered SV position is applied so do the
    # buffered filtering here
    add_to_buffered_pos_shared_summary( $sv_info ) if
	$input_params->{sv_position_buffer_length};

    # sr, pe, sr + pe read count filter
    next unless inds_pass_read_support_filter( $sv_info );
    # genotype 0/0, 1/1, 0/1 etc
    next unless inds_pass_genotype_filter( $sv_info );
    # filter min or specific number of inds that are 'shared' ie 1/1 or 0/1
    next unless inds_pass_shared_inds_count_filter( $sv_info );
    # filter by all homoz or all hetroz or 
    next unless inds_pass_fixed_genotype_filter( $sv_info );

    print_to_filtered_vcf( $line, $sv_info );

    add_to_later_summarize( $sv_info, 'FILTERED' );
    add_sv_for_stats( $sv_info, 'FILTERED' );
}

$vcf_fh->close;
$vcf_out_fh->close;
#print Dumper $_{PRE_FILTERED_SVS}; exit;
#print Dumper $_{FILTERED_SVS}; exit;
#print Dumper $_{SV_TYPE_STATS}; exit;
#print Dumper $_{INDIVIDUAL_STATS}; exit;
#print Dumper $_{BUFFERED_SVS}; exit;

# print pre-filtered SVs
{
    my @headers = (qw/ CHR START END VARIANT /);
    SV_TYPE: for my $sv_type(qw/ DEL INS INV DUP /) {
	# create empty file if sv_type is not found
	my $sv_type_full_name = full_name_for_sv_type( $sv_type );
	my $fh = IO::File->new( '>ALL.'.$sv_type_full_name.'.tsv' );
	my $svs = get_sorted_svs_for_type( $sv_type, 'PRE_FILTERED' );
	$fh->close and next SV_TYPE unless $svs;
	$fh->print( join("\t", @headers),"\n" );
	print "  Writing ALL $sv_type_full_name SVs\n";
	for my $sv_desc( @$svs ) {
	    my @tmp = split(/\s+/, $sv_desc);
	    $fh->print( join("\t", @tmp),"\n" );
	}
	$fh->close;
    }
}

# print filtered SVs
{
    my @headers = (qw/ CHR START END VARIANT HZ_RATE SHARED_INDS / );
    SV_TYPE: for my $sv_type(qw/ DEL INS INV DUP /) {
	# create empty file if no SVs exists for type
	my $sv_type_full_name = full_name_for_sv_type( $sv_type );

	my $fh = IO::File->new( '>FILTERED.'.$sv_type_full_name.'.tsv' );
	my $svs = get_sorted_svs_for_type( $sv_type, 'FILTERED' );
	# just create an empty file if no filtered SVs exists
	$fh->close and next SV_TYPE unless $svs;
	print "  writing FILTERED $sv_type_full_name SVs\n";
	$fh->print( join("\t", @headers),"\n" );
	for my $sv_desc( @$svs ) {
	    my @tmp = split(/\s+/, $sv_desc);
	    # tmp[0] = chr
	    # tmp[1] = start
	    # tmp[2] = end
	    # tmp[3] = sv summary
	    # tmp[4] = homozyg rate
	    # tmp[5] = shared individuals count
	    $fh->print( join("\t", @tmp),"\n" );
	}
	$fh->close;
    }
}

my $gene_positions; # global
if( $input_params->{genes_bed_file} ) {
    print "FILTERING GENE MATCHING SVs\n";
    my @headers = (qw/ CHR VAR_START VAR_END VARIANT GENE_START HZ_RATE SHARED_INDS GENE_END GENE_NAME GENE_LENGTH / );
    $gene_positions = load_gene_positions();
    SV_TYPE: for my $sv_type(qw/ DEL INS INV DUP /) {
	my @overlaps;
	# if no $sv_type exists, output an enpty file
	# unless input sv_types is specified and it's not one of them
	next SV_TYPE unless type_does_not_match_user_input( $sv_type );
	my $file_root = full_name_for_sv_type( $sv_type );
	my $file = 'FILTERED.'.$file_root.'.MATCHING_GENES.tsv';
	my $fh = IO::File->new( ">$file" );
	my $gene_list_file = 'FILTERED.'.$file_root.'.MATCHING_GENES_UNIQUE_LIST.txt';
	my $gene_list_fh = IO::File->new( ">$gene_list_file" );
	my @genes_seen;
	my $svs = get_sorted_svs_for_type( $sv_type, 'FILTERED' );
	print "   After filtering, no SVs are avilable for sv type: $sv_type\n" and next SV_TYPE unless $svs;
	$fh->print( join("\t", @headers),"\n" );
	print " Writing $file_root.".".MATCHING_GENES\n";
	for my $sv_desc( @$svs ) {
	    my @tmp = split(/\s+/, $sv_desc);
	    my $chr   = $tmp[0];
	    my $start = $tmp[1];
	    my $end   = $tmp[2];
	    my $sv_summary = $tmp[3];
	    my $homoz_rate = $tmp[4];
	    my $shared_inds = $tmp[5];
	    my @gene_pos = sv_maps_in_gene( $chr, $start, $end );
	    if( @gene_pos ) {
		for( @gene_pos ) {
		    my @ar = split(/\s+/, $_);
		    # ar[0] = gene_start
		    # ar[1] = gene_end
		    # ar[2] = gene_name
		    # for unique gene list
		    push @genes_seen, $ar[2] unless grep( /^$ar[2]$/, @genes_seen );
		    push @overlaps, $chr.' '.$start.' '.$end.' '.$sv_summary.' '.$homoz_rate.' '.$shared_inds.' '.$_;
		}
	    }
	}
	for( @overlaps ) {
	    my @tmp = split(/\s+/, $_);
	    $fh->print( join("\t", @tmp),"\n" );
	}
	# print unique gene names
	$gene_list_fh->print( join("\n", @genes_seen ),"\n" ) if @genes_seen;
	$fh->close;
	$gene_list_fh->close;
    }
}

# print INDIVIDUAL stats
{
    my @stat_types = (qw/ TOTAL HAS_CALL IS_HOMOZYGOUS SHARED_SVS_COUNT /);
    my $fh = IO::File->new( '>INDIVIDUAL_STATS.txt' );
    my $spacing = "%-20s%14s%14s%14s%14s%14s%14s%14s%14s%14s";
    print_individual_stats_file_headers( $fh, $spacing );

    for my $ind( @$TARGET_INDIVIDUALS ) {
	# add pre filtered stats
	my @prefiltered_stats;
	for my $stat_type( @stat_types ) {
	    my $count = ( $_{INDIVIDUAL_STATS}{PRE_FILTERED}{$ind}{TOTALS}{$stat_type} )
		? $_{INDIVIDUAL_STATS}{PRE_FILTERED}{$ind}{TOTALS}{$stat_type}
	        : 0 ;
	    push @prefiltered_stats, $count;
	}
	# add ratio
	for( 1 .. 3 ) {
	    my $ratio = get_ratio( $prefiltered_stats[0], $prefiltered_stats[$_] );
	    $prefiltered_stats[$_] = $prefiltered_stats[$_].'('.$ratio.')';
	}
	# add filtered stats
	my @filtered_stats;
	for my $stat_type( @stat_types ) {
	    my $count = ( $_{INDIVIDUAL_STATS}{FILTERED}{$ind}{TOTALS}{$stat_type} )
		? $_{INDIVIDUAL_STATS}{FILTERED}{$ind}{TOTALS}{$stat_type}
	        : 0 ;
	    push @filtered_stats, $count;
	}
	for( 1 .. 3 ) {
	    my $ratio = get_ratio( $filtered_stats[0], $filtered_stats[$_] );
	    $filtered_stats[$_] = $filtered_stats[$_].'('.$ratio.')';
	}
	my $filtered_ratio = get_ratio( $prefiltered_stats[0], $filtered_stats[0] );
	$fh->printf( "$spacing\n", $ind, @prefiltered_stats, $filtered_ratio, @filtered_stats );
    }
    $fh->close;
}

# print SV_TYPE stats
{
    my @stat_types = (qw/ COUNT HOMOZYGOUS SHARED / );
    my $stats_fh = IO::File->new( '>SV_TYPE_STATS.txt' );
    my $spacing = "%-14s%14s%14s%14s%14s%14s%14s%14s";
    print_sv_type_stats_headers( $stats_fh, $spacing );
    for my $sv_type( &sv_types_to_process ) {
	my $abbreviated_sv_name = abb_name_for_sv_type( $sv_type );
	next unless type_does_not_match_user_input( $abbreviated_sv_name );
	# prefiltered stats
	my $pre_filtered_stats = $_{SV_TYPE_STATS}{ PRE_FILTERED }{ SV_TYPES }{ $sv_type };
	my @stats_prefiltered;
	for my $stat_type( @stat_types ) {
	    my $value = ( $pre_filtered_stats->{$stat_type} )
		? $pre_filtered_stats->{$stat_type}
	        : 0 ;
	    push @stats_prefiltered, $value;
	}
	for( 1 .. 2 ) {
	    my $ratio = get_ratio( $stats_prefiltered[0], $stats_prefiltered[$_] );
	    $stats_prefiltered[$_] = $stats_prefiltered[$_].'('.$ratio.')';
	}
	# filtered stats
	my $filtered_stats = $_{SV_TYPE_STATS}{ FILTERED }{ SV_TYPES }{ $sv_type };
	my @stats_filtered;
	for my $stat_type( @stat_types ) {
	    my $value = ( $filtered_stats->{$stat_type} )
		? $filtered_stats->{$stat_type}
	        : 0 ;
	    push @stats_filtered, $value;
	}
	for( 1 .. 2 ) {
	    my $ratio = get_ratio( $stats_filtered[0], $stats_filtered[$_] );
	    $stats_filtered[$_] = $stats_filtered[0].'('.$ratio.')';
	}
	my $filter_ratio = get_ratio( $stats_prefiltered[0], $stats_filtered[0] );
	$stats_fh->printf( "$spacing\n", $sv_type, @stats_prefiltered, $filter_ratio, @stats_filtered );
    }
    $stats_fh->close;
    print "WROTE FILTERED STATS\n";
}

if( $input_params->{sv_position_buffer_length} ) {
    my @headers = (qw/ CHR START END INDIVIDUALS(COUNT) /);
    my $buf_len = $input_params->{sv_position_buffer_length};
    SV_TYPE: for my $sv_type(qw/ DEL INS INV DUP / ) {
	next SV_TYPE unless type_does_not_match_user_input( $sv_type );
	my $all_svs = $_{BUFFERED_SVS}{$sv_type};
	next SV_TYPE unless $all_svs;
	my $full_name_for_sv = full_name_for_sv_type( $sv_type );
	my $fh = IO::File->new( '>BUFFERED_AND_FILTERED.'.$full_name_for_sv.'.TXT' );
	$fh->print( join("\t", @headers),"\n" );
	for my $chr( sort {$a<=>$b} keys %$all_svs ) {
	    #print "CHR: $chr\n";
	    my $chr_svs = $all_svs->{$chr};
	    my %sorted_svs;
	    for my $svs( @$chr_svs ) {
		next if $svs->{is_processed};
		$svs->{is_processed} = 1;
		#print "  CHECKING EXPANDED POSITION ".$svs->{buffered_pos}." FOR OVERLAPS\n";
		my( $expanded_start, $expanded_end ) = split('-', $svs->{buffered_pos});
		my @expanded_inds = @{$svs->{individuals}};
		COMP: for my $compared_svs( @$chr_svs ) {
		    # prevent self-comparison
		    next COMP if $compared_svs->{is_processed};
		    my( $compared_start, $compared_end ) = split('-', $compared_svs->{buffered_pos});
		    if( overlaps( $expanded_start, $expanded_end, $compared_start, $compared_end ) ) {
			#print "    FOUND: $expanded_start $expanded_end <=> $compared_start $compared_end\n";
			$compared_svs->{is_processed} = 1;
			push @expanded_inds, @{$compared_svs->{individuals}};
		    }
		}
		#print "  POSITION ".$svs->{buffered_pos}." has ".scalar( @expanded_inds )." individuals\n";
		$sorted_svs{$expanded_start}{$expanded_end} = \@expanded_inds;
	    }
	    # print RESULT
	    for my $start( sort {$a<=>$b} keys %sorted_svs ) {
		for my $end( sort {$a<=>$b} keys %{$sorted_svs{$start}} ) {
		    #print "RESULT: $sv_type $chr $start $end ".scalar( @{$sorted_svs{$start}{$end}} )."\n";
		    my $inds = $sorted_svs{$start}{$end};
		    next unless buffered_inds_pass_read_support_filter( $inds );
		    next unless buffered_inds_pass_genotype_filter( $inds );
		    next unless buffered_inds_pass_shared_inds_count_filter( $inds );
		    next unless buffered_inds_pass_fixed_genotype_filter( $inds );
		    #print Dumper $inds;
		    my $target_inds_and_counts = target_individuals_found_and_counts( $inds );
		    $fh->print( join("\t", $chr, $start, $end, $target_inds_and_counts ),"\n" );
		}
	    }
	}
	$fh->close;
    }
}

exit;

sub target_individuals_found_and_counts {
    my $inds = shift;
    my %matches;
    for my $ind( @$inds ) {
	my $ind_name = $ind->{name};
	$matches{ $ind_name }++ if grep(/^$ind_name$/, @$TARGET_INDIVIDUALS);
    }
    my @counts;
    for my $target_ind( @$TARGET_INDIVIDUALS ) {
	my $ct = ( exists $matches{ $target_ind } ) ? $matches{ $target_ind } : 0 ;
	push @counts, $target_ind.'('.$ct.')';
    }
    return join(',', @counts);
}

sub buffered_inds_pass_fixed_genotype_filter {
    my $inds = shift;
    return 1 unless $input_params->{fixed_homo_var} || $input_params->{fixed_hetro} || $input_params->{fixed_homo_ref};
    # fixed_homo_var = 1/1
    # fixed_homo_ref = 0/0
    # fixed_hetro    = 0/1 
    my %non_target_inds_genotypes;
    if( $input_params->{output_svs_only_unique_to_input_individuals} ) {
	for my $ind( @$inds ) {
	    my $ind_name = $ind->{name};
	    next unless grep(/^$ind_name$/, @$NON_TARGET_INDIVIDUALS);
	    my $ind_genotype = $ind->{GT};
	    $non_target_inds_genotypes{ $ind_genotype }++;
	}
    }
    my %target_inds_genotypes;
    for my $ind( @$inds ) {
	my $ind_name = $ind->{name};
	next unless grep(/^$ind_name$/, @$TARGET_INDIVIDUALS);
	my $ind_genotype = $ind->{GT};
	$target_inds_genotypes{ $ind_genotype }++;
    }
    my $target_inds_genotypes_count = scalar( keys %target_inds_genotypes );
    if( $input_params->{fixed_homo_var} ) {
        # expect only 1/1
	return if $non_target_inds_genotypes{'1/1'};
	return unless $target_inds_genotypes_count == 1 && exists $target_inds_genotypes{'1/1'};
    }
    elsif( $input_params->{fixed_hetro} ) {
        # expect only 1/0
        return if $non_target_inds_genotypes{'1/0'};
	return unless $target_inds_genotypes_count == 1 && exists $target_inds_genotypes{'1/0'};
    }
    elsif( $input_params->{fixed_homo_ref} ) {
	# expect only 0/0
        return if $non_target_inds_genotypes{'0/0'};
        return unless $target_inds_genotypes_count == 1 && exists $target_inds_genotypes{'0/0'};
    }
    return 1;
}

sub buffered_inds_pass_shared_inds_count_filter {
    my $inds = shift;
    return 1 unless $input_params->{target_shared_individuals} || $input_params->{min_shared_individuals};
    # SHARED = GT 1/1 OR 0/1
    if( $input_params->{output_svs_only_unique_to_input_individuals} ) {
	for my $ind( @$inds ) {
	    my $ind_name = $ind->{name};
	    next unless grep(/^$ind_name$/, @$NON_TARGET_INDIVIDUALS);
	    return if $ind->{GT} eq '1/1' || $ind->{GT} eq '0/1';
	}
    }
    my %shared_inds;
    for my $ind( @$inds ) {
	my $ind_name = $ind->{name};
	next unless grep( /^$ind_name$/, @$TARGET_INDIVIDUALS );
	$shared_inds{$ind_name}++ if $ind->{GT} eq '1/1' || $ind->{GT} eq '0/1';
    }
    my $shared_inds_count = scalar( keys %shared_inds );
    if( $input_params->{target_shared_individuals} ) {
	return unless $shared_inds_count == $input_params->{target_shared_individuals};
    }
    if( $input_params->{min_shared_individuals} ) {
	return unless $shared_inds_count >= $input_params->{min_shared_individuals};
    }
    return 1;
}

sub buffered_inds_pass_genotype_filter {
    my $inds = shift;
    return 1 unless $input_params->{genotype};
    if( $input_params->{output_svs_only_unique_to_input_individuals} ) {
	for my $ind( @$inds ) {
	    my $ind_name = $ind->{name};
	    next unless grep(/^$ind_name$/, @$NON_TARGET_INDIVIDUALS);
	    return if $input_params->{genotype} eq 'heterozygous' && $ind->{GT} eq '1/1';
	    return if $input_params->{genotype} eq 'homozygous'   && $ind->{GT} eq '0/1';
	}
    }
    for my $ind( @$inds ) {
	my $ind_name = $ind->{name};
	next unless grep( /^$ind_name$/, @$TARGET_INDIVIDUALS );
	return 1 if $input_params->{genotype} eq 'heterozygous' && $ind->{GT} eq '1/1';
	return 1 if $input_params->{genotype} eq 'homozygous'   && $ind->{GT} eq '0/1';
    }
    return 1;
}

sub buffered_inds_pass_read_support_filter {
    my $inds = shift;
    return 1 unless $input_params->{min_support} || $input_params->{min_sr_support} || $input_params->{min_pe_support};
    if( $input_params->{output_svs_only_unique_to_input_individuals} ) {
	# if any of the non-user specified ind passes filter, reject the whole record
	for my $ind( @$inds ) {
	    my $ind_name = $ind->{name};
	    next unless grep(/^$ind_name$/, @$NON_TARGET_INDIVIDUALS);
	    return if ind_passes_min_support_filter( $ind );
	    return if ind_passes_min_pe_support_filter( $ind );
	    return if ind_passes_min_sr_support_filter( $ind );
	}
    }
    for my $ind( @$inds ) {
	my $ind_name = $ind->{name};
	next unless grep( /^$ind_name$/, @$TARGET_INDIVIDUALS );
	return unless ind_passes_min_support_filter( $ind );
	return unless ind_passes_min_pe_support_filter( $ind );
	return unless ind_passes_min_sr_support_filter( $ind );
    }
    return 1;
}

sub input_vcf_filehandle {
    if( $input_params->{vcf_file} =~ /\.gz$/ ) {
	return IO::File->new( 'zcat '.$input_params->{vcf_file}.' |' );
    }
    else {
	return IO::File->new( $input_params->{vcf_file} );
    }
    return;
}

sub output_vcf_filehandle {
    my $fh = IO::File->new( '>'.$input_params->{vcf_file}.'.FILTERED.VCF' );
    add_headers_to_vcf_out( $fh );
    return $fh;
}

sub inds_pass_fixed_genotype_filter {
    my $sv_info = shift;
    return 1 unless $input_params->{fixed_homo_var} || $input_params->{fixed_hetro} || $input_params->{fixed_homo_ref};
    # fixed_homo_var = 1/1
    # fixed_homo_ref = 0/0
    # fixed_hetro    = 0/1
    my %non_target_inds_genotypes;
    if( $input_params->{output_svs_only_unique_to_input_individuals} ) {
	for my $ind( @$NON_TARGET_INDIVIDUALS ) {
	    my $sv = $sv_info->{individuals}->{$ind};
	    my $ind_genotype = $sv->{GT};
	    $non_target_inds_genotypes{ $ind_genotype }++;
	}
    }
    my %target_inds_genotypes;
    for my $ind( @$TARGET_INDIVIDUALS ) {
	my $sv = $sv_info->{individuals}->{$ind};
	my $ind_genotype = $sv->{GT};
	$target_inds_genotypes{ $ind_genotype }++;
    }
    my $target_inds_genotypes_count = scalar( keys %target_inds_genotypes );

    if( $input_params->{fixed_homo_var} ) {
	# expect only 1/1
	return if $non_target_inds_genotypes{'1/1'};
	return unless $target_inds_genotypes_count == 1 && exists $target_inds_genotypes{'1/1'};
    }
    elsif( $input_params->{fixed_hetro} ) {
	# expect only 1/0
	return if $non_target_inds_genotypes{'1/0'};
	return unless $target_inds_genotypes_count == 1 && exists $target_inds_genotypes{'1/0'};
    }
    elsif( $input_params->{fixed_homo_ref} ) {
	# expect only 0/0
	return if $non_target_inds_genotypes{'0/0'};
	return unless $target_inds_genotypes_count == 1 && exists $target_inds_genotypes{'0/0'};
    }
    return 1;
}

sub get_ratio {
    my $value1 = shift;
    my $value2 = shift;
    return '0.0' if $value1 == 0 || $value2 == 0;
    my $ratio = sprintf( "%.1f", $value2 / $value1 * 100 );
    return $ratio;
}

sub print_sv_type_stats_headers {
    my $fh = shift;
    my $spacing = shift;
    
    my $filtered_shared_header = 'SHARED(>=2)';
    $filtered_shared_header = 'SHARED(=='.$input_params->{target_shared_individuals}.')' if $input_params->{target_shared_individuals};
    $filtered_shared_header = 'SHARED(>='.$input_params->{min_shared_individuals}.')' if $input_params->{min_shared_individuals};
    my @headers = ( 'SV_COUNT', 'HOMOZYGOUS', 'SHARED(>=2)', 'RATE', 'SV_COUNT', 'HOMOZYGOUS', $filtered_shared_header );
    $fh->printf( "$spacing\n", '', 'PRE_FILTERED', 'PRE_FILTERED','PRE_FILTERED', 'FILTER', 'FILTERED', 'FILTERED','FILTERED' );
    $fh->printf( "$spacing\n", 'SV_TYPES', @headers );
    $fh->print( '=' x 112,"\n" );
    return;
}

sub print_individual_stats_file_headers {
    my $fh = shift;
    my $spacing = shift;
    my $filtered_shared_header = 'SHARED(>=2)';
    $filtered_shared_header = 'SHARED(=='.$input_params->{target_shared_individuals}.')' if $input_params->{target_shared_individuals};
    $filtered_shared_header = 'SHARED(>='.$input_params->{min_shared_individuals}.')' if $input_params->{min_shared_individuals};
    my @headers1 = (qw/ INDIVIDUAL PRE_FILTERED PRE_FILTERED PRE_FILTERED PRE_FILTERED FILTER FILTERED FILTERED FILTERED FILTERED /);
    my @headers2 = ('NAMES', 'TOTAL', 'WITH_CALL', 'HOMOZYGOUS', 'SHARED(>=2)', 'RATE', 'TOTAL', 'WITH_CALL', 'HOMOZYGOUS', $filtered_shared_header );
    $fh->printf( "$spacing\n", @headers1 );
    $fh->printf( "$spacing\n", @headers2 );
    $fh->print( '=' x 146,"\n" );
    return;
}

sub add_sv_for_stats {
    my $sv = shift;
    my $type = shift; # PRE OR POST FILTERED

    # these hash key differ depending on the FILTERED or PRE_FILTERED
    # individuals also differ
    my( $shared_inds_key, $homozygous_rate_key, @inds_to_process );
    if( $type eq 'FILTERED' ) {
	@inds_to_process     = @$TARGET_INDIVIDUALS;
	$shared_inds_key     = 'filtered_shared_individuals';
	$homozygous_rate_key = 'filtered_homozygous_rate';
    }
    elsif( $type eq 'PRE_FILTERED' ) {
	@inds_to_process     = @$ALL_INDIVIDUALS;
	$shared_inds_key     = 'pre_filtered_shared_individuals';
	$homozygous_rate_key = 'pre_filtered_homozygous_rate';
    }
    else {
	die "Invalid processed type, $type, expected either FILTERED OR PRE_FILTERED\n";
    }

    # SV TYPE STATS
    my $sv_type = full_name_for_sv_type( $sv->{SVTYPE} );
    $_{SV_TYPE_STATS}{ $type }{SV_TYPES}{ $sv_type }{ COUNT }++;
    $_{SV_TYPE_STATS}{ $type }{SV_TYPES}{ $sv_type }{ HOMOZYGOUS }++ if $sv->{ $homozygous_rate_key } > 0;
    $_{SV_TYPE_STATS}{ $type }{SV_TYPES}{ $sv_type }{ SHARED }++ if $sv->{ $shared_inds_key } >= 2;

    # IND STATS
    # ./. = missing call
    # 0/0 = homozygous reference
    # 0/1 = heterozygous
    # 1/1 = homozygous variant
    for my $ind_id( @inds_to_process ) {
	my $ind_info = $sv->{individuals}->{$ind_id};
	$_{INDIVIDUAL_STATS}{ $type }{ $ind_id }{ TOTALS }{ TOTAL }++;
	$_{INDIVIDUAL_STATS}{ $type }{ $ind_id }{ SV_TYPES }{ $sv_type }{ TOTAL }++;
	# NOT missing call
	unless( $ind_info->{GT} eq './.' ) {
	    $_{INDIVIDUAL_STATS}{ $type }{ $ind_id }{ TOTALS }{ HAS_CALL }++;
	    $_{INDIVIDUAL_STATS}{ $type }{ $ind_id }{ SV_TYPES }{ $sv_type }{ HAS_CALL }++;
	}
	# individual's sv is homozygous
	if( $ind_info->{GT} eq '1/1' ) {
	    $_{INDIVIDUAL_STATS}{ $type }{ $ind_id }{ TOTALS }{ IS_HOMOZYGOUS }++;
	    $_{INDIVIDUAL_STATS}{ $type }{ $ind_id }{ SV_TYPES }{ $sv_type }{ IS_HOMOZYGOUS }++;
	}
	# shared
	if( $sv->{ $shared_inds_key } >= 2 ) {
	    $_{INDIVIDUAL_STATS}{ $type }{ $ind_id }{ TOTALS }{ SHARED_SVS_COUNT }++;
	    $_{INDIVIDUAL_STATS}{ $type }{ $ind_id }{ SV_TYPES }{ $sv_type }{ SHARED_SVS_COUNT }++;
	}
    }
    return;
}

sub sv_types_to_process {
    return (qw/ DELETIONS INSERTIONS INVERSIONS DUPLICATIONS /);
}

sub clean_up_old_output_files {
    # remove any existing results files from previous attempts to avoid any confusion
    for my $sv_type( &sv_types_to_process ) {
	unlink 'ALL.'.$sv_type.'.tsv';
	unlink 'FILTERED.'.$sv_type.'.tsv';
    }
    return;
}

sub sv_is_of_excluded_chromosome {
    my $sv_info = shift;
    my $sv_chr = $sv_info->{chrom};
    return 1 if $CHRS_TO_EXCLUDE->{ $sv_chr };
    return;
}

sub has_failed_inds_among_user_input_inds {
    my $sv_info = shift;
    for my $ind( @$TARGET_INDIVIDUALS ) {
	my $sv = $sv_info->{individuals}->{$ind};
	return if $sv->{param_filter_status} && $sv->{param_filter_status} eq 'failed';
    }
    return 1;
}

sub has_failed_inds_among_all_inds {
    my $sv_info = shift;
    for my $ind( keys %{$sv_info->{individuals}} ) {
	my $sv = $sv_info->{individuals}->{$ind};
	return if $sv->{param_filter_status} && $sv->{param_filter_status} eq 'failed';
    }
    return 1;
}

sub print_to_filtered_vcf {
    my $old_line = shift;
    my $sv_info = shift;

    chomp $old_line;
    my @processed_ind_svs;
    for my $ind( @$TARGET_INDIVIDUALS ) {
	my $sv = $sv_info->{individuals}->{$ind}->{vcf_entry};
	die "Can't find sv info for individual $ind\n" if not defined $sv;
	push @processed_ind_svs, $sv;
    }

    my @tmp = split(/\s+/, $old_line);
    my $new_line = join("\t", @tmp[0 .. 8], @processed_ind_svs);
    $vcf_out_fh->print( "$new_line\n" );
}

sub get_sorted_svs_for_type {
    my $sv_type = shift;
    my $filter_status = shift; # FILTERED OR PRE_FILTERED
    my $svs = $_{$filter_status.'_SVS'}{$sv_type};
    return if not $svs;
    my @sorted_svs;
    for my $chr( sort {$a<=>$b} keys %$svs ) {
	for my $start( sort {$a<=>$b} keys %{$svs->{$chr}} ) {
	    for my $end( sort {$a<=>$b} keys %{$svs->{$chr}->{$start}} ) {
		for my $sv_summary( @{$svs->{$chr}->{$start}->{$end}} ) {
		    push @sorted_svs, $chr.' '.$start.' '.$end.' '.$sv_summary;
		}
	    }
	}
    }
    return \@sorted_svs;
}

sub type_does_not_match_user_input {
    my $sv_type = shift;
    my $input_sv_types = $input_params->{sv_types};
    return 1 if not $input_sv_types; # no user input so all types are okay
    my @input_types = split(',', $input_sv_types);
    return 1 if grep(/^$sv_type$/, @input_types);
    return;
}

sub sv_maps_in_gene {
    my $chr = shift;
    my $sv_start = shift;
    my $sv_end = shift;
    return unless $gene_positions->{$chr};
    my @overlaps;
    for my $gene( @{$gene_positions->{$chr}} ) {
	# look for overlap
	if( $sv_start >= $gene->{start} && $sv_end <= $gene->{end} ||  $sv_end >= $gene->{start} && $sv_end <= $gene->{end} || $gene->{start} >= $sv_start && $gene->{start} <= $sv_end || $gene->{end} >= $sv_start && $gene->{end} <= $sv_start ) {
	    my $gene_length = abs( $gene->{end} - $gene->{start} );
	    push @overlaps, $gene->{start}.' '.$gene->{end}.' '.$gene->{name}.' '.$gene_length;
	}
    }

    return @overlaps if @overlaps;
    return;
}

sub load_gene_positions {
    my $fh = IO::File->new( $input_params->{genes_bed_file} ) ;
    my %pos;
    while( my $line = $fh->getline ) {
	chomp $line;
	next if $line =~ /^\s+$/;
	my @tmp = split(/\s+/, $line);
	next unless $tmp[1] =~ /\d+/;  # header??
	# chr      start    end        gene id
	my %h = (
	    'chr'   => $tmp[0],
	    'start' => $tmp[1],
            'end'   => $tmp[2],
            'name'  => $tmp[3],
            );
	push @{$pos{$tmp[0]}}, \%h;
    }
    $fh->close;
    return \%pos;
}

sub common_column_id_headers {
    return(qw/ CHROM POS ID REF ALT QUAL FILTER INFO FORMAT /);
}

sub add_headers_to_vcf_out {
    my $vcf_out_fh = shift;
    my $fh = IO::File->new( $input_params->{vcf_file} );
    while( my $line = $fh->getline ) {
	last unless $line =~ /^#/;
	if( $line =~ /^##/ ) {
	    $vcf_out_fh->print( $line ); # generic headers
	}
	elsif( $line =~ /^#CHROM/ ) {
	    my @common_headers = common_column_id_headers();
	    my $new_header = join("\t", @common_headers, @$TARGET_INDIVIDUALS );
	    $vcf_out_fh->print( '#'.$new_header."\n" );
	}
	else {
	    # reached end of headers
	    last;
	}
    }
    $fh->close;
    return;
}

sub abb_name_for_sv_type {
    my $sv_type = shift;
       my %names = (
	   'DELETIONS'    => 'DEL',
	   'INSERTIONS'   => 'INS',
	   'INVERSIONS'   => 'INV',
	   'DUPLICATIONS' => 'DUP',
	   'BNDS'         => 'BND',
       );
    my $abb = $names{$sv_type};
    die "Can't determine abbreviated sv name for $sv_type\n" unless defined $abb;
    return $abb;
}

sub full_name_for_sv_type {
    my $sv_type = shift;
    my %names = (
	'DEL' => 'DELETIONS',
	'INS' => 'INSERTIONS',
	'INV' => 'INVERSIONS',
	'DUP' => 'DUPLICATIONS',
	'BND' => 'BNDS',
    );
    my $name = $names{$sv_type};
    die "Can't determine file name for sv type: $sv_type\n" if not defined $name;
    return $name;
}

sub add_to_later_summarize {
    my $sv_info = shift;
    my $filter_status = shift; # FILTERED OR PRE_FILTERED
    # summary eg +-sD(100)
    my $sv_summary = create_sv_summary( $sv_info ); 
    my $type = $sv_info->{SVTYPE};
    my $chr  = $sv_info->{chrom};
    my $start= $sv_info->{start};
    my $len  = $sv_info->{SVLEN};
    my $end  = $start + $len;
    my @summary = ( $sv_summary );
    if( $filter_status eq 'FILTERED' ) {
	push @summary, $sv_info->{filtered_homozygous_rate}    if exists $sv_info->{filtered_homozygous_rate};
	push @summary, $sv_info->{filtered_shared_individuals} if exists $sv_info->{filtered_shared_individuals};
    }
    push @{$_{$filter_status.'_SVS'}{ $type }{ $chr }{ $start }{ $end }}, join(' ', @summary);
    return;
}

sub add_to_buffered_pos_shared_summary {
    my $sv_info = shift;

    my $type = $sv_info->{SVTYPE};
    my $chr  = $sv_info->{chrom};
    my $start = $sv_info->{start};
    my $end = $start + $sv_info->{SVLEN};
    my $buf_len = $input_params->{sv_position_buffer_length};
    $start = ( ($start - $buf_len) > 0 ) ? $start - $buf_len : 1 ;
    $end += $buf_len;
    my $inds = get_individuals_in_record( $sv_info );
    my %h = (
	chr => $chr,
	buffered_pos => $start.'-'.$end,
	individuals => $inds,
    );
    push @{$_{BUFFERED_SVS}{$type}{$chr}}, \%h;
    return;
}

sub get_individuals_in_record {
    my $sv_info = shift;
    my $inds = $sv_info->{individuals};
    my @individuals;
    for my $ind_name( keys %{$sv_info->{individuals}} ) {
	my %h = %{$sv_info->{individuals}->{$ind_name}};
	$h{name} = $ind_name;
	push @individuals, \%h;
    }
    return \@individuals;
}

sub overlaps {
    my( $new_start, $new_end, $start, $end ) = @_;
    if( $new_start >= $start && $new_start <= $end || $new_end >= $start && $new_end <= $end || $start >= $new_start && $start <= $new_end || $end >= $new_start && $end <= $new_end ) {
	return 1;
    }
    return;
}

sub create_sv_summary {
    my $sv_info = shift;
    my %abbrevs = (
	'DEL' => 'D',
	'INS' => 'I',
	'INV' => 'V',
	'DUP' => 'U',  
	'BND' => 'T',
    );
    my $abbrev = $abbrevs{ $sv_info->{SVTYPE} };
    die "Can't determine SVTYPE abbreviation for ".$sv_info->{SVTYPE}."\n" if not defined $abbrev;
    return $sv_info->{SVSTR}.'s'.$abbrev.'('.$sv_info->{SVLEN}.')';
}

sub inds_pass_genotype_filter {
    my $sv_info = shift;
    return 1 unless $input_params->{genotype};
    my $inds = $sv_info->{individuals};
    if( $input_params->{output_svs_only_unique_to_input_individuals} ) {
	# none of the non-target individuals must pass this filter
	# if output_svs_only_unique_to_input_individual is set
	for my $ind( @$NON_TARGET_INDIVIDUALS ) {
	    return if $input_params->{genotype} eq 'heterozygous' && $inds->{$ind}->{GT} eq '1/1';
	    return if $input_params->{genotype} eq 'homozygous'  && $inds->{$ind}->{GT} eq '0/1';
	}
    }
    for my $ind( @$TARGET_INDIVIDUALS ) {
	# passes if at least one individual is homoz or hetroz
	return 1 if $input_params->{genotype} eq 'heterozygous' && $inds->{$ind}->{GT} eq '1/1';
	return 1 if $input_params->{genotype} eq 'homozygous'  && $inds->{$ind}->{GT} eq '0/1';
    }
    return;
}

sub inds_pass_read_support_filter {
    my $sv_info = shift;
    # check for NON-target inds that pass this
    return 1 unless $input_params->{min_support} || $input_params->{min_sr_support} || $input_params->{min_pe_support};
    {
	next;
	# not sure if user wants to actually do this since this seems to
	# filter everything out
	if( $input_params->{output_svs_only_unique_to_input_individuals} ) {
	    # this fails if any of the non-target individuas pass this filter
	    # since we're only looking for SVs where only target inds pass
	    for my $ind( @$NON_TARGET_INDIVIDUALS ) {
		my $ind_info = $sv_info->{individuals}->{$ind};
		return if ind_passes_min_support_filter( $ind_info );
		return if ind_passes_min_pe_support_filter( $ind_info );
		return if ind_passes_min_sr_support_filter( $ind_info );
	    }
	}
    }
    for my $ind( @$TARGET_INDIVIDUALS ) {
	my $ind_info = $sv_info->{individuals}->{$ind};
	return unless ind_passes_min_support_filter( $ind_info );
	return unless ind_passes_min_pe_support_filter( $ind_info );
	return unless ind_passes_min_sr_support_filter( $ind_info );
    }
    
    return 1;
}

sub inds_pass_shared_inds_count_filter {
    my $sv_info = shift;
    return 1 unless $input_params->{target_shared_individuals} || $input_params->{min_shared_individuals};

    if( $input_params->{output_svs_only_unique_to_input_individuals} ) {
	my $non_target_inds_count = 0;
	for my $ind( @$NON_TARGET_INDIVIDUALS ) {
	    my $ind_info = $sv_info->{individuals}->{$ind};
	    $non_target_inds_count++ if $ind_info->{GT} eq '1/1' || $ind_info->{GT} eq '0/1';
	}
	return if $non_target_inds_count > 0;
    }

    my $read_support_passing_inds_count = 0;
    for my $ind( @$TARGET_INDIVIDUALS ) {
	my $ind_info = $sv_info->{individuals}->{$ind};
	$read_support_passing_inds_count++ if $ind_info->{GT} eq '1/1' || $ind_info->{GT} eq '0/1';
    }
    if( $input_params->{target_shared_individuals} ) {
	return unless $read_support_passing_inds_count == $input_params->{target_shared_individuals};
    }
    if( $input_params->{min_shared_individuals} ) {
	return unless $read_support_passing_inds_count >= $input_params->{min_shared_individuals};
    }
    return 1;
}

sub ind_passes_min_sr_support_filter {
    my $ind_info = shift;
    my $min_sr = $input_params->{min_sr_support};
    return 1 if not $min_sr;
    return if $ind_info->{SR} < $min_sr;
    return 1;
}

sub ind_passes_min_pe_support_filter {
    my $ind_info = shift;
    my $min_pe = $input_params->{min_pe_support};
    return 1 if not $min_pe;
    return if $ind_info->{PE} < $min_pe;
    return 1;
}

sub ind_passes_min_support_filter {
    my $ind_info = shift;
    my $min_support = $input_params->{min_support};
    return 1 if not $min_support;
    return if $ind_info->{SUP} < $min_support;
    return 1;
}

sub sv_matches_input_sv_type {
    my $sv_info = shift;
    my $sv_type = $sv_info->{SVTYPE};
    # track unique sv types
    my $input_sv_types = $input_params->{sv_types};
    return 1 unless defined $input_sv_types;
    my @types_to_process = split(',', $input_sv_types);
    return 1 if grep( /^$sv_type$/, @types_to_process );
    return;
}

sub exclude_svs_not_unique_to_input_inds {
    my $sv_info = shift;
    return unless $input_params->{output_svs_only_unique_to_input_individuals};
    my @desired_inds = split(',', $input_params->{individuals});
    for my $ind( keys %{$sv_info->{individuals}} ) {
	return unless grep( /^$ind$/, @desired_inds );
    }
    return 1;
}

sub sv_max_length_is_met {
    my $sv_info = shift;
    return 1 unless $input_params->{max_sv_length};
    return if $sv_info->{SVLEN} > $input_params->{max_sv_length};
    return 1;
}

sub sv_min_length_is_met {
    my $sv_info = shift;
    return 1 unless $input_params->{min_sv_length};
    return if $sv_info->{SVLEN} < $input_params->{min_sv_length};
    return 1;
}

sub record_is_BND_sv {
    my $sv_info = shift;
    return 1 if $sv_info->{SVTYPE} eq 'BND';
    return;
}

sub get_line_sv_info {
    my $line = shift;
    chomp $line;
    my @tmp = split(/\s+/, $line);
    my @tmp_svinfo = split(';', $tmp[7]);
    my %info;
    for( @tmp_svinfo ) {
	if( $_ =~ /^(\S+)=(\S+)$/ ) {
	    $info{$1} = $2;
	}
	else {
	    $info{$_} = 1;
	}
    }
    die "Can't deteremine sv type because SVTYPE is not found in line: $line\n" if not $info{SVTYPE};
    return \%info;
}

sub get_sv_info_from_line {
    my $line = shift;
    chomp $line;
    my @tmp = split(/\s+/, $line);
    my %info;
    $info{chrom} = $tmp[0];
    $info{start} = $tmp[1];

    my( $svtype, $svlen, $svstr ) = parse_sv_type_and_length( $line );
    $info{SVLEN} = $svlen;
    $info{SVTYPE} = $svtype;
    $info{SVSTR} = $svstr;
    $info{pre_filtered_shared_individuals} = 0;

    # load ALL individuals
    for my $ind_name( @$ALL_INDIVIDUALS ) {
	my $column_num = $header_column_number{ $ind_name };
	my $sv_for_ind = $tmp[ $column_num - 1 ];
	my $ind_sv_info = parse_ind_sv_info( $sv_for_ind );
	# eg:
	#'vcf_entry' => '0/1:10:5:5',
	#'SR'  => '5',
	#'PE'  => '5',
	#'SUP' => '10',
	#'GT'  => '0/1'
	$info{individuals}{$ind_name} = $ind_sv_info;
	$info{pre_filtered_shared_individuals}++ if $ind_sv_info->{GT} eq '1/1' || $ind_sv_info->{GT} eq '0/1'
    }
    my $pre_filtered_homozygous_rate = sprintf( "%.2f", $info{pre_filtered_shared_individuals} / scalar( @$ALL_INDIVIDUALS) );
    $info{pre_filtered_homozygous_rate} = $pre_filtered_homozygous_rate;

    # add filtered homozygous count
    $info{filtered_shared_individuals} = 0;
    for my $ind_name( @$TARGET_INDIVIDUALS ) {
	my $ind_sv_info = $info{individuals}{$ind_name};
	$info{filtered_shared_individuals}++ if$ind_sv_info->{GT} eq '1/1' || $ind_sv_info->{GT} eq '0/1'
    }
    my $filtered_homozygous_rate = sprintf( "%.2f", $info{filtered_shared_individuals} / scalar( @$TARGET_INDIVIDUALS) );
    $info{filtered_homozygous_rate} = $filtered_homozygous_rate;

    return \%info;
}

sub parse_ind_sv_info {
    my $vcf_entry = shift;
    my %info;
    my @tmp = split(':', $vcf_entry);
    for(qw/ GT SUP PE SR /) {
	$info{$_} = shift @tmp;
    }
    $info{vcf_entry} = $vcf_entry;
    return \%info;
}

sub parse_sv_type_and_length {
    my $line = shift;
    my @tmp_line = split(/\s+/, $line);
    my @tmp_info = split(';', $tmp_line[7]);
    my( $svtype, $svlen, $svstr );
    for my $info( @tmp_info ) {
	my( $name, $value ) = split('=', $info);
	next unless defined $name && defined $value;
	$svtype = $value if $name eq 'SVTYPE';
	$svlen = abs( $value ) if $name eq 'SVLEN';
	if( $name eq 'STR' || $name eq 'STRANDS' ) {
	    $svstr = $value;
	    $svstr =~ s/:\d+//g;
	}
    }
    # expect these to exist in all records
    die "Can not parse SVTYPE from line:\n\t$line\n" if not defined $svtype;
    die "Can not parse STR from line:\n\t$line\n" if not defined $svstr;

    return $svtype, $svlen, $svstr;
}

sub validated_individuals_to_process {
    my $available_inds = get_individuals_in_vcf();
    my( @target_individuals, @non_target_individuals );
    if( $input_params->{individuals} ) {
	for my $ind( split(',', $input_params->{individuals}) ) {
	    die "Invalid individual name input .. $ind is not found in vcf file.  Valid names are: ".join(' ', @$available_inds)."\n" unless grep( /^$ind$/, @$available_inds );
	    push @target_individuals, $ind;
	}
    }
    #elsif( $input_params->{genome_specie_names} ) {
	#my @species_names = split(',', $input_params->{genome_specie_names});
	#for my $specie_name( @species_names ) {
	    #print "Finding samples for specie: $specie_name\n";
	    #my $taxon = Genome::Taxon->get( name => $specie_name );
	    #print "No genome specie found for name, $specie_name\n" and exit unless $taxon;
	    #my @samples = $taxon->samples;
	    #print "No genome samples found for species name: $specie_name\n" and exit unless @samples;
	    #for my $sample( @samples ) {
		#my $name = $sample->name;
		#print " skipping sample, $name, not found in vcf file\n" and next unless grep( /^$name$/i, @$available_inds );
		#push @target_individuals, $name;
	    #}
	#}
    #}
    else {
	@target_individuals = @$available_inds;
    }
    for my $ind( @$available_inds ) {
	push @non_target_individuals, $ind unless grep(/^$ind$/, @target_individuals );
    }
    return $available_inds, \@target_individuals, \@non_target_individuals;
}

sub get_column_id_header_in_vcf {
    my $fh = IO::File->new( $input_params->{vcf_file} );
    my @ids;
    while( my $line = $fh->getline ) {
	if( $line =~ /^#CHROM/ ) {
	    chomp $line;
	    $line =~ s/^#//;
	    @ids = split(/\s+/, $line);
	    last;
	}
    }
    $fh->close;
    die "Couldn't find column ids header in vcf file\n" if not @ids;;
    return @ids;
}

sub get_individuals_in_vcf {
    my @expected_headers = common_column_id_headers();
    my @column_headers = get_column_id_header_in_vcf();
    for my $expected_header( @expected_headers ) {
	my $header = shift @column_headers;
	die "Expected header to be $expected_header but is $header\n" unless $header eq $expected_header;
    }
    # what's left should be individual just names
    die "Did not find any individual names in column id header in vcf file\n" if not @column_headers;
    return \@column_headers;
}

sub valid_input_params {
    return( qw/
exclude_chromosomes
fixed_hetero
fixed_homo_ref
fixed_homo_var
genes_bed_file
genome_specie_names
genotype
individuals
max_sv_length
min_pe_support
min_shared_individuals
min_sr_support
min_support
min_sv_length
sv_types
target_shared_individuals
target_support_passing_individuals
vcf_file
output_svs_only_unique_to_input_individuals
sv_position_buffer_length
/);
}

sub parse_input_params {
    my %params;
    my @valid_inputs = valid_input_params();
    for( @ARGV ) {
	my( $name, $value ) = split('=', $_);
	die "Invalid input param name: $name\n" unless grep( /^$name$/, @valid_inputs );
	die "Can't get name and value from param: $_ .. expected format is 'name=value'\n" unless defined $name && defined $value;
	if( $name =~ /file/ ) {
	    die "\nCan't find $value or $value is empty: $value\n" unless -s $value;
	}
	if( $name =~ /^min_/ || $name =~ /^max_/ ) {
	    die "\nExpected value for $name to be a number, not $value\n" unless $value =~ /^\d+$/;
	}
	if( $name =~ /^sv_position_buffer_length/ ) {
	    print "\nExpected value for $name to be a number > 0, not $value\n" and exit unless $value =~ /^\d+$/ && $value > 0;
	}
	if( $name eq 'genotype' ) {
	    die "\nExpected value for $name to be either heterozygous or homozygous not $value\n" unless lc $value eq 'heterozygous' || lc $value eq 'homozygous';
	}
	print "Param $name has been entered multiple times\n" and exit if $params{ $name };
	$params{ $name } = $value;
    }
    die "Need vcf_file=?? input\n" if not $params{vcf_file};

    if( $params{min_shared_individuals} && $params{target_shared_individuals} ) {
	print "\nCan't input both min_shared_individuals and target_shared_individuals .. choose only one of the two\n";
	exit;
    }
    if( $params{min_shared_individuals} || $params{target_shared_individuals} ) {
	print "\nmin_shared_individuals or target_shared_individuals require specifying of genotype=homozygous or genotype=heterozygou\n" and exit unless $params{genotype};
    }
    if( $params{individuals} && $params{genome_specie_names} ) {
	print "\nCan't input both individuals and genome_specie_names to parse individuals\n";
	exit;
    }
    if( $params{output_svs_only_unique_to_input_individuals} ) {
	print "\noutput_svs_only_unique_to_input_individuals param must be accompanied by either individuals or genome_specie_names to filter SVs that only match these individuals and no other individuals\n" and exit unless
	    $params{individuals} || $params{genome_specie_names};
    }
    # can't input more than 1 of the following: fixed_homo_var, fixed_hetro, fixed_homo_ref since they are mutually exclusive
    my $fixed_genotype_input_counts = 0;
    my @fixed_genotype_inputs = (qw/ fixed_homo_var fixed_hetro fixed_homo_ref / );
    for( @fixed_genotype_inputs ) { 
	$fixed_genotype_input_counts++ if $params{ $_ };
    }
    if( $fixed_genotype_input_counts > 1 ) {
	print "Can't input more than one of the following: ". join(',', @fixed_genotype_inputs)." they are mutually exclusive\n";
	exit;
    }
    if( $fixed_genotype_input_counts > 0 && $params{genotype} ) {
	print "Can't specify ".join(',', @fixed_genotype_inputs)," and genotype since they are mutually exclusive\n";
	exit;
    }
    return \%params;
}

sub set_chromosomes_to_exclude {
    my %chrs;
    if( $input_params->{exclude_chromosomes} ) {
	%chrs = map{ $_=> 1} ( split(',', $input_params->{exclude_chromosomes}) );
    }
    return \%chrs;
}

sub set_header_column_number {
    my $fh = IO::File->new( $input_params->{vcf_file} );
    while( my $line = $fh->getline ) {
	next if $line =~ /^##/;
	if( $line =~ /^#CHROM/ ) {
	    chomp $line;
	    $line =~ s/^#//;
	    my @tmp = split(/\t+/, $line);
	    my $column_number = 0;
	    for( @tmp ) {
		$header_column_number{ $_ } = ++$column_number;
	    }
	    last;
	}
    }
    $fh->close;
    return;
}

sub print_usage {
    print "USAGE: $0\n".
	"\tINPUT LUMPY VCF FILE\n".
	"\t vcf_file=<INPUT_VCF_FILE>            #input vcf file\n".
	"\tFILTER FOR SAMPLES/INDIVIDUALS\n".
	"\t individuals=<?NAME1,NAME2..>         #consider only these individuals\n".
	#"\t genome_specie_names=<?NAME>          #consider samples for this specie\n".
	"\tFILTER FOR VARIANT TYPES\n".
	"\t sv_types=<?TYPE1,TYPE2>              #consider only these types of SVs (valid inputs: DEL INS DUP INV)\n".
	"\tFILTER FOR VARIANT LENGTHS\n".
	"\t min_sv_length=<?LENGTH>              #exclude SVs less than this length\n".
	"\t max_sv_length=<?LENGTH>              #exclude SVs greater than this length\n".
	"\tFILTER FOR CHROMOSOMES/REFERENCE SEQUENCES\n".
	"\t exclude_chromosomes=<?,CHR?>         #exclude these chromosomes\n".
	"\tFILTER FOR VARIANT SUPPORTING READS\n".
	"\t min_support=<?COUNT>                 #exclude SVs with less than this number of pe + sr read support\n".
	"\t min_pe_support=<?COUNT>              #exclude SVs with less than this number of pe read support\n".
	"\t min_sr_support=<?COUNT>              #exclude SVs with less than this number of sr read support\n".
	"\tFILTER FOR GENOTYPES\n".
	"\t genotype=<?homozygous/heterozygous>  #homozygous = at least one 1/1 individual; heterozygous = at least one 0/1 individual\n".
	"\t target_shared_individuals=<?COUNT>   #exclude unless SV has this exact number of shared individuals, ie, 1/1 or 0/1\n".
	"\t min_shared_individuals=<?COUNT>      #exclude unless SV has this min number of shared individuals, ie, 1/1 or 0/1\n".
	"\t fixed_homo_var=1                     #exclude unless all individuals are homozygous variant(1/1)\n".
	"\t fixed_hetero=1                       #exclude unless all individuals are heterozygous(0/1)\n".
	"\t fixed_homo_ref=1                     #exclude unless all individuals are homozygous reference(0/0)\n".
	"\tBUFFER SV COORDINATES TO INCREASE LIKELY HOOD OF FINDING SHARED SVS\n".
	"\t sv_position_buffer_length=<?NUMBER>  #extend sv position by this length up and down stream to increase the likelyhood of finding shared SVS\n".
	"\tFILTER FOR FIXED INDIVIDUALS\n".
	"\t output_svs_only_unique_to_input_individuals=<?1> #output only SVs that are unique to input individuals\n".
	"\tOUTPUT SVS THAT OVERLAP GENE REGIONS IN BED FILE\n".
	"\t genes_bed_file=<?FILE>               #output variants that overlap gene coordinates in this bed file\n";
}
