version 1.0
## Copyright Broad Institute, 2020
## This script takes as input a set of raw de novo SNV calls and runs them through a filtering pipeline
## Filter codes include:
##  MAF = population minor allele frequency (gnomAD)
##  COD = protein-coding
##  MUC-HLA = MUC- and HLA- genes; tend to be very noisy
##  CAF = cohort allele frequency
##  IGV = IGV manual review result
##  dbSNP = fail if variant present in dbSNP
##  OUT = outlier_flag result
##  
## Variants passing all filters are then passed to EM-mosaic for mosaic variant discovery
##  
## TESTED: 
## Versions of other tools on this image at the time of testing:
##
## LICENSING : This script is released under the WDL source code license (BSD-3) (see LICENSE in https://github.com/broadinstitute/wdl). 
## Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script. 
## Please see the docker for detailed licensing information pertaining to the included programs.
##

###########################################################################
#WORKFLOW DEFINITION
###########################################################################
workflow ApplyFilters {

	input {
		File anno 
		String output_prefix

		String CAF_outprefix

		File estimation_script
		Int cohort_size

		File script_CAF

		File outlier_script

		Int expected_dnsnvs 
		Int case_cutoff 

		File script_update_filter_col

		File script_filtct 

		File script_printpass 
	}

	#########################################################
	## generate cohort allele frequency and outlier flags
	#########################################################
	# generate allele counts (# carriers) file
	call estimate_cohort_AF {
		input:
			script = estimation_script,
			infile = anno,
			cohort_size = cohort_size
	}

	# run CAF (cohort allele frequency)
	call filter_CAF {
		input:
			infile = anno,
			script = script_CAF,
			caf_file = estimate_cohort_AF.out
	}

	#run outlier filter
	call filter_outlier{
		input:
			infile = filter_CAF.out,
			script = outlier_script,
			cohort_size = cohort_size,
			exp = expected_dnsnvs,
			cutoff = case_cutoff
	}

	#########################################################
	## parse filter flags, summarize filtering, output variants passing all filters
	#########################################################
	#run update_filter_column script to combine filter flags into single column
	call update_filt_col {
		input:
			infile = filter_outlier.out,
			script = script_update_filter_col
	}

	#run script to summarize counts of variants flagged by each filter
	call summarize_counts{
		input:
			infile = update_filt_col.outfile,
			script = script_filtct
	}

	#run script to write out variants passing all filters, to be used as input to EM-mosaic
	call print_pass_vars{
		input:
			infile = update_filt_col.outfile,
			script = script_printpass
	}


	output {
		File denovos_all = update_filt_col.outfile # full denovos table with all annotations
		File denovos_PASS = print_pass_vars.outfile # only denovos passing all filters, to be passed to EM-mosaic scoring step
		File filter_summary = summarize_counts.outfile # summary table detailing how many variants failed each filter step
	}


}

###########################################################################
#Task Definitions
###########################################################################

# estimates cohort AF
task estimate_cohort_AF{

	input {
		File script
		File infile
		Int cohort_size
	}

	String outprefix = basename(infile, '.denovo.VEP.PV4.SB.FDR.RR.VC.txt')
	String outfname = "AC.~{outprefix}.txt"

	command {
		python ~{script} -i ~{infile} -n ~{cohort_size} -o ~{outfname}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outfname}"
	}
}


#Adds cohort allele frequency filter-related columns to variants file for downstream filtering
task filter_CAF {

	input {
		File infile
		File script
		File caf_file
	}

	String outprefix = basename(infile, '.txt')

	command {
		python ~{script} -i ~{infile} -c ~{caf_file} -o "~{outprefix}.CAF.txt"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.CAF.txt"
	}

}

#Adds Outlier filter-related columns to variants file for downstream filtering
task filter_outlier {

	input {
		File infile
		File script
		Int cohort_size
		Int exp
		Int cutoff
	}

	String outprefix = basename(infile, '.txt')

	command {
		Rscript ~{script} ~{infile} ~{outprefix} ~{cohort_size} ~{exp} ~{cutoff}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.OUT.txt"
	}

}

#Parses filter information from upstream steps and combines into an updated filter column
#Note: assumes 'filter' col already exists (from upstream reformat_vcf.py step)
#Note: also applies filters for 
#   (1) popfreq (MAF)
#   (2) protein-coding (COD)
#   (3) MUC/HLA genes (MUC-HLA)
#   (4) dbSNP (dbSNP)
task update_filt_col {

	input {
		File infile 
		File script 
	}

	String outsuffix = basename(infile, '.VEP.PV4.SB.FDR.RR.VC.CAF.OUT.txt')
	String outfname = "ADfile.~{outsuffix}.txt"

	command {
		python ~{script} -i ~{infile} -o ~{outfname}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File outfile = "~{outfname}"
	}

}


#Parses filter column and summarize how many variants are flagged by each filter
task summarize_counts {

	input {
		File infile 
		File script
	}

	String outprefix = basename(infile, '.txt')

	command {
		python ~{script} ~{infile} "~{outprefix}.SUMMARY_COUNTS.txt"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File outfile = "~{outprefix}.SUMMARY_COUNTS.txt"
	}
}

#Prints out all variants that pass all filters and that belong to samples that are not flagged as outliers
task print_pass_vars {

	input {
		File infile 
		File script 
	}
	
	String outprefix = basename(infile, '.txt')

	command {
		python ~{script} ~{infile} "~{outprefix}.PASS.txt"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File outfile = "~{outprefix}.PASS.txt"
	}
}

