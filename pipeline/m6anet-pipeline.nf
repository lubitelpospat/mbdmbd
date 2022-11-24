#!/usr/bin/env nextflow

params.sequencingSummary = 'NO_FILE'
params.genome = ''
params.fastqDirectory = ''
params.fast5Directory = ''

process mapGenome{
    cpus params.t

    input:
    path genome
    path joinedFastq
    
    output:
    path "${genome.baseName}.sam", emit: mappedGenome
    path "${genome.baseName}.mmi", emit: genomeIndex

    script:
    if (params.mode == 'ncrna')
	    """
	    minimap2 -ax sr -d ${genome.baseName}.mmi $genome
	    minimap2 -ax sr -uf -t ${params.t} --secondary=no ${genome.baseName}.mmi ${joinedFastq} > ${genome.baseName}.sam
	    """
	else if (params.mode == 'cdna')
		"""
	    minimap2 -ax map-ont -d ${genome.baseName}.mmi $genome
	    minimap2 -ax map-ont -uf -t ${params.t} --secondary=no ${genome.baseName}.mmi ${joinedFastq} > ${genome.baseName}.sam
	    """
}

process joinFastq{
    input:
    path fastqDirectory
    
    output:
    path "joinedFastq.fastq", emit: fastq

    script:
    """
    cat $fastqDirectory/*.fastq > joinedFastq.fastq
    """

}

process sam2bam{
    input:
    path sam

    output:
    path "${sam.baseName}.bam", emit: bam
    path "${sam.baseName}.bam.bai", emit: bai

    script:
    """
    samtools view -Sb ${sam} | samtools sort -o ${sam.baseName}.bam -
    samtools index ${sam.baseName}.bam
    """
}

process nanopolishEventAlign {
    cpus params.t

    input:
    path fast5
    path joinedFastq
    path mappedBam
    path indexedBai //Just needs to be in the working directory
    path genome

    output:
    path 'summary.txt', emit: summary
    path 'eventalign.txt', emit: eventAlign


    """
    nanopolish index -d $fast5 $joinedFastq
    nanopolish eventalign --scale-events --signal-index --summary=summary.txt --threads=${params.t} --reads=$joinedFastq --bam=$mappedBam --genome=$genome  > eventalign.txt
    """
}

process m6anetDataprep{
    cpus params.t
    conda params.condaenv

    input:
    path nanopolishEventalign

    output:
    path "./dataprep_output/data.index", emit: index
    path "./dataprep_output/data.json", emit: json
    path "./dataprep_output/data.log", emit: log
    path "./dataprep_output/data.readcount", emit: readcount
    path "./dataprep_output/eventalign.index", emit: eventalignIndex
    path "./dataprep_output/", emit: outDir
    
    script:
    """
    m6anet-dataprep --eventalign $nanopolishEventalign --out_dir ./dataprep_output/ --n_processes ${params.t}
    """
}

process m6anetInference{
    cpus params.t
    conda params.condaenv

    input:
    path dataDir

    output:
    path './inference_output/*', emit: inference 
    
    """
    m6anet-run_inference --input_dir $dataDir --out_dir ./inference_output/ --infer_mod_rate --n_processes 96
    """
}

workflow{
    joinFastq( Channel.fromPath( params.fastqDirectory ) )
    mapGenome( Channel.fromPath( params.genome ), joinFastq.out.fastq )
    sam2bam( mapGenome.out.mappedGenome )
    nanopolishEventAlign( Channel.fromPath( params.fast5Directory ), joinFastq.out.fastq, sam2bam.out.bam, sam2bam.out.bai, Channel.fromPath( params.genome ) )
    m6anetDataprep( nanopolishEventAlign.out.eventAlign )
    m6anetInference( m6anetDataprep.out.outDir )
}
