#!/usr/bin/env nextflow

/////////////////////////////////////////////////////////////////
///////////////////////// setup /////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////// setup the sample data ////////////////////////
//// MUST BE nested array where the 1st element is codename and second is path to the parent folder 

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

fh = new File(params.inputFile)
def csv_content = fh.getText('utf-8')
 
def data_iterator = parseCsv(csv_content, separator: '\t', readFirstLine: false)
// println data_iterator.getClass()  // class com.xlson.groovycsv.CsvIterator
 
def sampleArray = []
def sampleArrayPairs = []
i = 0
for (line in data_iterator) {
    sampleArray[i] = [line[0],line[1]]
    sampleArrayPairs[i] = [line[0],line[2],line[3]]
    i += 1;
}

rundata = Channel.from(sampleArray).unique()

////////////////////////////////////////////////////////////
//////////////////////// dataprep //////////////////////////
////////////////////////////////////////////////////////////

//// MAIN CHANNELS ////
rundata.into { runDataChannel; runDataChanMirror }
genomeChannel = Channel.fromPath( params.genome )

// once for N samples (they use the same index)
if (params.mode == "cdna") {
  process genomeIndexCdna {
    publishDir "./output/genomeIndexed"
    memory params.mem
    cpus params.t

    output:
    file('*.mmi') into minimapIndex

    script:
    """
    minimap2 -ax map-ont -d genome_index.mmi ${params.genome}
    """
  }
}
if (params.mode == "ncrna") {
  process genomeIndexNcrna {
    publishDir "./output/genomeIndexed"
    memory params.mem
    cpus params.t

    output:
    file('*.mmi') into minimapIndex

    script:
    """
    minimap2 -ax sr -d genome_index.mmi ${params.genome}
    """
  }
}

// once per sample
process catfq {
  publishDir "./output/$runname/catFq"

  memory params.mem
  cpus params.t

  input:
  tuple runname, filepath from runDataChannel

  output:
  tuple runname, '*.fastq' into joinFq
  tuple runname, '*.fastq' into joinFqMirror

  script:
  """
  cat ${filepath}/fastq_pass/*.fastq > ${runname}.fastq
  """
}

minimapChan = joinFq.combine(minimapIndex)

if (params.mode == "cdna") {
  process miniMappingCdna {
    publishDir "./output/$runname/mapped"
    memory params.mem
    cpus params.t

    input:
    tuple runname, catFq, indexFile from minimapChan

    output:
    tuple runname, '*.sam' into samFiles

    script:
    """
    minimap2 -ax map-ont -uf -t ${params.t} --secondary=no ${indexFile} ${catFq} > ${indexFile.baseName}.sam
    """
  }
}
if (params.mode == "ncrna") {
  process miniMappingNcrna {
    publishDir "./output/$runname/mapped"
    memory params.mem
    cpus params.t

    input:
    tuple runname, catFq, indexFile from minimapChan

    output:
    tuple runname, '*.sam' into samFiles

    script:
    """
    minimap2 -ax sr -uf -t ${params.t} --secondary=no ${indexFile} ${catFq} > ${indexFile.baseName}.sam
    """
  }
}

process sam2bam {
  publishDir "./output/$runname/bam"
  memory params.mem
  cpus params.t

  input:
  tuple runname, samFile from samFiles

  output:
  tuple runname, '*.bam' into bamFiles

  script:
  """
  samtools view -Sb ${samFile} | samtools sort -o ${samFile.baseName}.bam -
  samtools index ${samFile.baseName}.bam
  """
}

///////// combine using shared keys (runname) ////////////////////
npInputChan = runDataChanMirror.join(joinFqMirror).join(bamFiles)

process nanopolish {
  publishDir "./output/$runname/np"
  memory params.mem
  cpus params.t

  input:
  tuple runname, filepath, catFq, bamFile from npInputChan

  output:
  tuple runname, 'eventalign.txt' into eventAlign
  tuple runname, 'summary.txt' into summary

  script:
  """
  nanopolish index -d "${filepath}/fast5_pass" ${catFq}
  nanopolish eventalign --reads ${catFq} --bam ${bamFile} --genome ${params.genome} --signal-index --scale-events --summary summary.txt --threads ${params.t} > eventalign.txt
  """
}

process dataprep {
  conda './m6aenv.yml'
  publishDir "./output/$runname/${runname}_dataprep", mode: 'copy'
  memory params.mem
  cpus params.t

  input:
  tuple val(runname), path(eventalign) from eventAlign

  output: 
  file('*') into xporeOut
  tuple val(runname), val("${params.homeDir}/output/${runname}/${runname}_dataprep") into dataprep

  shell:
  """
  m6anet-dataprep --eventalign $eventalign --n_processes ${params.t} --out_dir "."
  """
}


