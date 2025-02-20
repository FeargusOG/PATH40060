#! nextflow

params.help = null
def helpMessage() {
  log.info"""
  ------------------------------------------------------------------------------
                            MANEline PIPELINE
  ------------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/tutorial_25

  Description:

  Get fastq from Google Drive, unsplit, prep for Sarek.

  Mandatory arguments:

    --runID         [str]       Identifier for the run

    --gdPath        [str]       Path to Google Drive files (all files therein taken!)

    --email         [str]       Email address to send reports

    --outDir        [str]       Output directory (top level)

  Optional:

    --help          [str]       Display help
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()

//requires
if(params.runID == null){
  exit 1, "Please specify --runID"
}

if(params.gdPath == null){
  exit 1, "Please specify --gdPath, path of Google Drive dir in which are fastq files"
}

if(params.email == null){
  exit 1, "Please specify --email"
}

if(params.outDir == null){
  exit 1, "Please specify --outDir"
}

//whole workflow
workflow {
  get_fq()
  set_fq(get_fq.out.collect())
  fin_fq(set_fq.out)
}

//process to get data
process get_fq {  
  publishDir 'data/fastq', mode: 'copy', overwrite: false
  
  output:
  path 'split/*.fq'

  script:
  """
  gdown ${params.gdPath} -O split --folder
  """
}

//Get the fastq bits, and operate on them
process set_fq {
    publishDir 'data/fastq/unsplit', mode: 'copy', overwrite: false

    input:
    path fqs
    
    output:
    path "${fastq_gz}"

    script:
    fastq_gz = params.runID + ".fastq.gz"
    """
    cat *.fq | gzip > ${fastq_gz}
    """
}

process fin_fq {

    input:
    path fqs

    script:
    """
    echo "Workflow complete"
    """
}

//when completed...
workflow.onComplete {
    sleep(100)
    def subject = """\
        [brucemoran/tutorial_25] SUCCESS [$workflow.runName]
        """
        .stripIndent()
    if (!workflow.success) {
        subject = """\
            [brucemoran/tutorial_25] FAILURE [$workflow.runName]
            """
            .stripIndent()
    }

    def msg = """\
        Pipeline execution summary
        ---------------------------
        RunName     : ${workflow.runName}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: "${params.email}",
                subject: subject,
                body: msg)
}
