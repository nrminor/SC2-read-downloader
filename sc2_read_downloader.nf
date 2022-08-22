#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// NOTE that helpful runtime parameters are in the file `nextflow.config`
// We recommend that you modify any file paths, inputs, or outputs there.


// Defining workflow for generating non-supplementary figures for Halfmann et al. 2022
workflow {


	// INPUT CHANNELS
	ch_ont_reads = Channel // this channel streams in data on the SRA runs sequenced on ONT instruments
		.fromPath ( params.samplesheet )
		.splitCsv ( sep: ",", header: true )
		.map      { row -> tuple(row.sra_id, row.file_basename, row.platform, file(row.primer_set)) }
		.filter   { it[2] == "ont" }

	ch_illumina_reads = Channel // this channel streams in data on the SRA runs sequenced on Illumina instruments
		.fromPath ( params.samplesheet )
		.splitCsv ( sep: ",", header: true )
		.map      { row -> tuple(row.sra_id, row.file_basename, row.platform, file(row.primer_set)) }
		.filter   { it[2] == "illumina" }


	// WORKFLOW STEPS
	// ONT read portion
	GET_ONT_READS (
		ch_ont_reads
	)

	// Illumina read portion
	GET_ILL_READS (
		ch_illumina_reads
	)

}


process GET_ONT_READS {

	// Using sra-tools fasterq-dump to full these files now, though note:
	// https://github.com/ncbi/sra-tools/issues/463

	tag "${filename}"
	publishDir params.results, pattern: '*.fastq.gz', mode: 'move'

	input:
	tuple val(sra_id), val(filename), val(platform), file(primers)

	output:
	tuple val(filename), file("*.fastq.gz"), file(primers)

	script:
	"""

	prefetch ${sra_id}
	fasterq-dump ${sra_id}/${sra_id}.sra \
	--concatenate-reads --skip-technical
	gzip ${sra_id}.sra.fastq
	mv ${sra_id}.sra.fastq.gz ${filename}.fastq.gz

	"""

}


process GET_ILL_READS {

	tag "${filename}"
	publishDir params.results, pattern: '*.fastq.gz', mode: 'move'

	input:
	tuple val(sra_id), val(filename), val(platform), file(primers)

	output:
	tuple val(filename), file("*_R1.fastq.gz"), file("*_R2.fastq.gz"), file(primers)

	script:
	"""

	prefetch ${sra_id}
	fasterq-dump ${sra_id}/${sra_id}.sra \
	--split-files --skip-technical
	gzip ${sra_id}.sra_1.fastq ; mv ${sra_id}.sra_1.fastq.gz ${filename}_R1.fastq.gz
	gzip ${sra_id}.sra_2.fastq ; mv ${sra_id}.sra_2.fastq.gz ${filename}_R2.fastq.gz

	"""

}

