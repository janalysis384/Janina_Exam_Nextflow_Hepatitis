// to run the nextflow scirpt properly please run this with the singualrity container. nextflow.config is provided
// nextflow run hepatitis_delta_pipeline_compare_to_ref.nf -profile singularity
// with --accession an accession can be given
// with --downloadurl provide a different url
params.accession = "M21012" // this is the default setting, it is the hepatitis delta complete genome
params.downloadurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=M21012&rettype=fasta&retmode=text" // this is the link to M21012

params.out = "${projectDir}/output"
params.store = "${projectDir}/downloads" // for the sequences, also combined 
params.storeDir = "${projectDir}/cache" // for the "worked" on data after aligne trimal etc
params.input = "${projectDir}/input" // store the FASTA files you like to analyse here, script eill look for them here
// params.run_fastqc = false // need to check if I need that

// this process downdloads the reference genome
// I have no clue why I had to give the "" again in the wget, I had scirpts 
// where it worked without the "", if it is part of the 
// dowloadurl, but it works now
process downloadRefGenome {
	storeDir params.store

	publishDir params.out, mode: 'copy', overwrite: true

	output:
		path "${params.accession}.fasta"

	""" 
	wget "${params.downloadurl}" -O ${params.accession}.fasta
	"""
}

// combine the RevGenome with the TestGenome from my imaginary collegues  
process combineFASTA {
    storeDir params.store

	publishDir params.out, mode: 'copy', overwrite: true

    input:
		path infile
	output:
		path "RefGenome_${params.accession}_and_TestGenomes.fasta"

    """
    cat ${infile} ${params.input}/*.fasta > RefGenome_${params.accession}_and_TestGenomes.fasta
    """


}

// aligne sequences with mafft in a singularity container 
// mafft settings are default (gap opening penalty 1.53, gap extension 0.0)
// scirpt could also adjusted for global and local alignment
process mafftAligner {
    storeDir params.storeDir
	publishDir params.out, mode: 'copy', overwrite: true
	container "https://depot.galaxyproject.org/singularity/mafft%3A7.525--h031d066_1"
	
	input: 
		path combinedFASTA

	output:
		path "RefGenome_${params.accession}_and_TestGenomes_mafft_align.fasta"
	"""
    mafft ${combinedFASTA} > RefGenome_${params.accession}_and_TestGenomes_mafft_align.fasta
	"""
}

// trimAl is a tool for the automated removal of spurious sequences 
// or poorly aligned regions from a multiple sequence alignment
process trimal {
    storeDir params.storeDir
	publishDir params.out, mode: 'copy', overwrite: true
	container "https://depot.galaxyproject.org/singularity/trimal%3A1.5.0--h9948957_2"

    input:
        path mafftAligned
    
    output:
    path "RefGenome_${params.accession}_and_TestGenomes_mafft_align_trimal.fasta"
    path "RefGenome_${params.accession}_and_TestGenomes_mafft_align_trimal.html"

    """
    trimal -in ${mafftAligned} -htmlout RefGenome_${params.accession}_and_TestGenomes_mafft_align_trimal.html -out RefGenome_${params.accession}_and_TestGenomes_mafft_align_trimal.fasta -automated1
    """



}

workflow {
	RefGenome = downloadRefGenome()
    allFASTAinOne = combineFASTA(RefGenome)
    mafftAligned = mafftAligner(allFASTAinOne)
    trimal(mafftAligned)
}

