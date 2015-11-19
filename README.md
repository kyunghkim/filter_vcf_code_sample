# filter_vcf_code_sample

filter_lumpy_vcf.pl allows filtering of lumpy derived VCF file is variety of ways.  Input options are listed below:

USAGE:  filter_lumpy_vcf.pl
       INPUT LUMPY VCF FILE
        vcf_file=<INPUT_VCF_FILE>            #input vcf file
       FILTER FOR SAMPLES/INDIVIDUALS
        individuals=<?NAME1,NAME2..>         #consider only these individuals
       FILTER FOR VARIANT TYPES
        sv_types=<?TYPE1,TYPE2>              #consider only these types of SVs (valid inputs: DEL INS DUP INV)
       FILTER FOR VARIANT LENGTHS
        min_sv_length=<?LENGTH>              #exclude SVs less than this length
        max_sv_length=<?LENGTH>              #exclude SVs greater than this length
       FILTER FOR CHROMOSOMES/REFERENCE SEQUENCES
        exclude_chromosomes=<?,CHR?>         #exclude these chromosomes
       FILTER FOR VARIANT SUPPORTING READS
        min_support=<?COUNT>                 #exclude SVs with less than this number of pe + sr read support
        min_pe_support=<?COUNT>              #exclude SVs with less than this number of pe read support
        min_sr_support=<?COUNT>              #exclude SVs with less than this number of sr read support
       FILTER FOR GENOTYPES
        genotype=<?homozygous/heterozygous>  #homozygous = at least one 1/1 individual; heterozygous = at least one 0/1 individual
        target_shared_individuals=<?COUNT>   #exclude unless SV has this exact number of shared individuals, ie, 1/1 or 0/1
        min_shared_individuals=<?COUNT>      #exclude unless SV has this min number of shared individuals, ie, 1/1 or 0/1
        fixed_homo_var=1                     #exclude unless all individuals are homozygous variant(1/1)
        fixed_hetero=1                       #exclude unless all individuals are heterozygous(0/1)
        fixed_homo_ref=1                     #exclude unless all individuals are homozygous reference(0/0)
       BUFFER SV COORDINATES TO INCREASE LIKELY HOOD OF FINDING SHARED SVS
        sv_position_buffer_length=<?NUMBER>  #extend sv position by this length up and down stream to increase the likelyhood of finding shared SVS
       FILTER FOR FIXED INDIVIDUALS
        output_svs_only_unique_to_input_individuals=<?1> #output only SVs that are unique to input individuals
       OUTPUT SVS THAT OVERLAP GENE REGIONS IN BED FILE
        genes_bed_file=<?FILE>               #output variants that overlap gene coordinates in this bed file

Examples:
 filter_lumpy_vcf.pl vcf_file=file.vcf individuals=IND1,IND3 sv_types=DEL,INS min_sv_length=1000 max_sv_length=1000000 min_support=3
 #   process only IND1 and IND3 individuals
 #   process only deletions and insertions
 #   ignore SVs less than 1000 bases or greater than 1000000 bases
 #   igmore unless SV has pe_support + sr_support 3 or greater

 filter_lumpy_vcf.pl vcf_file=file.vcf individuals=IND1,IND3 sv_types=DEL,INS min_sv_length=1000 max_sv_length=1000000 min_support=3 genes_bed_file=genes.bed
 #   process only IND1 and IND3 individuals
 #   process only deletions and insertions
 #   ignore SVs less than 1000 bases or greater than 1000000 bases
 #   igmore unless SV has pe_support + sr_support 3 or greater
 #   output SVs that overlap genes in genes.bed file

 filter_lumpy_vcf.pl vcf_file=file.vcf individuals=IND1,IND3 sv_types=DEL,INS min_sv_length=1000 max_sv_length=1000000 min_support=3 genes_bed_file=genes.bed sv_position_buffer_length=500
 #   process only IND1 and IND3 individuals
 #   process only deletions and insertions
 #   ignore SVs less than 1000 bases or greater than 1000000 bases
 #   igmore unless SV has pe_support + sr_support 3 or greater
 #   output SVs that overlap genes in genes.bed file
 #   expand SV's position by 500 on each side to increase the likelyhood of finding SVs that overlap genes in genes.bed file