require(sevenbridges)
###############################################################################
# FASTQ Splitter
###############################################################################
in.lst<-list(
		
			input(id="file", 
				  type="File",
				  description="the input fastq file",
				 # inputBinding = CommandLineBinding(
						  position = 3L
				 # )
		  	),
			input(id="size",
					type="int",
					#inputBinding = CommandLineBinding(
							prefix = "--lines=",
							position = 1L,
							separate = FALSE
					#)
			),
			input(id="suffix",
					type="string",
					#inputBinding = CommandLineBinding(
							prefix = "--additional-suffix=",
							position = 2L,
							separate = FALSE
					#)
			)
	)
		
out.lst<-list(
		output(id="output", type="array", glob="*$(inputs.suffix)*")
		)

		
splitter<-Tool(
				id="fastq_splitter",
				baseCommand="split",
				inputs=in.lst,
				outputs=out.lst
)

###############################################################################
# Trimmer
###############################################################################
library(yaml)
yaml.file<-yaml::yaml.load_file("/ngs_share/tools/BS-seq-pipelines/CWL/tools/trimmomatic.cwl")
proc.yaml.file<-yaml.file
proc.yaml.file$id<-'Trimmomatic'
proc.yaml.file$inputs<-lapply(yaml.file$inputs, function(inp) {for(field in c("position","prefix","separate")) inp[[field]]<-inp$inputBinding[[field]]; inp$inputBinding<-NULL ;do.call("input", inp)})
proc.yaml.file$outputs<-lapply(yaml.file$outputs, function(outp) {for(field in c("glob")) outp[[field]]<-outp$outputBinding[[field]]; outp$outputBinding<-NULL ;do.call("output", outp)})
proc.yaml.file$requirements<-NULL
proc.yaml.file$cwlVersion<-NULL
trimmer<-do.call("Tool", proc.yaml.file)


###############################################################################
# Aligner
###############################################################################

### bwa-meth

in.lst<-list(
		
		input(id="threads",
				required = FALSE, 
				default = 1L,
				type="integer",
				prefix = "--threads",
				position = 1L,
				separate = TRUE
		),
		
#		input(id="prefix",
#				type="string",
#				prefix = "--prefix",
#				position = 2L,
#				separate = FALSE
#		),
		
		input(id="reference", 
				required = TRUE, 
				type="file",
				description="the reference fasta file",
				position = 2L,
				prefix = "--reference",
				secondaryFiles=c(
						".bwameth.c2t",
						".bwameth.c2t.amb",
						".bwameth.c2t.ann",
						".bwameth.c2t.bwt",
						".bwameth.c2t.pac",
						".bwameth.c2t.sa")
		),
		
		input(id="read1", 
				required = TRUE, 
				type="file",
				description="the input fastq file with the first mate",
				position = 3L,
		),
		
		input(id="read2", 
				required = TRUE, 
				type="file",
				description="the input fastq file with the second mate",
				position = 4L,
		)
		
)

out.lst<-list(
		output(id="alignment", type="file")
)


bwameth<-Tool(
		id="bwameth",
		baseCommand="bwameth.py",
		inputs=in.lst,
		outputs=out.lst
)
cat(bwameth$toYAML(), file="/ngs_share/tools/BS-seq-pipelines/CWL/tools/bwameth.cwl")

###############################################################################
# Methylation caller
###############################################################################

### PileOMeth

#https://github.com/dpryan79/PileOMeth




###############################################################################
###############################################################################
# Workflow
###############################################################################
###############################################################################

flow<- splitter %>>% trimmer
