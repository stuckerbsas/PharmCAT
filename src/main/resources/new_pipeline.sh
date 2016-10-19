#!/usr/bin/env bash

FASTA1='fasta1.gz'
FASTA2='fasta2.gz'
RG='@RG\tID:SRR504516_small_1.fq.gz\tPL:ILLUMINA\tSM:SRR504516'
WORKING_DIR='$WORKING_DIR'
REF='/usr/share/docker_data/root/work/utils/bwakit-0.7.12/hs37d5.fa'
GATK='/usr/share/docker_data/root/work/utils/gatk-3.4-0/GenomeAnalysisTK.jar'
PICARD='/usr/share/docker_data/root/work/utils/picard-1.130/picard.jar'
BWA='/usr/share/docker_data/root/work/utils/bwakit-0.7.12/run-bwamem'
JAVA_TEMP='-Djava.io.tmpdir=$WORKING_DIRtmp'
JAVA_SETTINGS="'-Xmx32768m'  '-XX:+UseParallelOldGC'  '-XX:ParallelGCThreads=4'  '-XX:GCTimeLimit=50'  '-XX:GCHeapFreeLimit=10' "
KNOWN_INDELS_MILLS='/usr/share/docker_data/root/work/ext_data/gatk_bundle/mills_and_1000G/Mills_and_1000G_gold_standard.indels.b37.vcf'
KNOWN_INDELS='/usr/share/docker_data/root/work/ext_data/gatk_bundle/1000G_indels/1000G_phase1.indels.b37.vcf'
KNOWN_DBSNP='/usr/share/docker_data/root/work/ext_data/gatk_bundle/pre-129_dbsnp/dbsnp_138.b37.excluding_sites_after_129.vcf'

'$BWA'  '-o' '$WORKING_DIRoutput'  '-R' '@RG\tID:SRR504516_small_1.fq.gz\tPL:ILLUMINA\tSM:SRR504516'  '-t' '36'  '-a'  '-d'  '-s'  '-S'  $REF  '$WORKING_DIR$FASTA1'  '$WORKING_DIR$FASTA2'  |sh
'java'  $JAVA_SETTINGS '$JAVA_TEMP'  '-jar' '$PICARD'  'BuildBamIndex'  'TMP_DIR=$WORKING_DIRtmp'  'VALIDATION_STRINGENCY=STRICT'  'MAX_RECORDS_IN_RAM=8000000'  INPUT='$WORKING_DIRoutput.aln.bam'  OUTPUT='$WORKING_DIRoutput.aln.bam.bai'
'mv'  '$WORKING_DIRoutput.aln.bam'  '$WORKING_DIRoutput.moved.bam'
'mv'  '$WORKING_DIRoutput.aln.bam.bai'  '$WORKING_DIRoutput.moved.bam.bai'
'java'  $JAVA_SETTINGS '$JAVA_TEMP'  '-jar' '$GATK'  '-T' 'RealignerTargetCreator'  '-I' '$WORKING_DIRoutput.moved.bam'  '-ip' '100'  '-R' '$REF'  '-nt' '36'  '-o' '$WORKING_DIRoutput.realigner.intervals'  '-known' '$KNOWN_INDELS_MILLS' '-known' '$KNOWN_INDELS'
'java'  $JAVA_SETTINGS '$JAVA_TEMP'  '-jar' '$GATK'  '-T' 'IndelRealigner'  '-I' '$WORKING_DIRoutput.moved.bam'  '-ip' '100'  '-R' '$REF'  '-known' '$KNOWN_INDELS_MILLS' '-known' '$KNOWN_INDELS'  '-targetIntervals' '$WORKING_DIRoutput.realigner.intervals'  '-o' '$WORKING_DIRoutput.realigned.bam'  '-filterNoBases'
'java'  $JAVA_SETTINGS '$JAVA_TEMP'  '-jar' '$GATK'  '-T' 'BaseRecalibrator'  '-I' '$WORKING_DIRoutput.realigned.bam'  '-ip' '100'  '-R' '$REF'  '-nct' '36'  '-knownSites' '$KNOWN_INDELS_MILLS' '-knownSites' '$KNOWN_INDELS' '-knownSites' '$KNOWN_DBSNP'  '-o' '$WORKING_DIRoutput.recal.before.table'
'java'  $JAVA_SETTINGS '$JAVA_TEMP'  '-jar' '$GATK'  '-T' 'PrintReads'  '-I' '$WORKING_DIRoutput.realigned.bam'  '-ip' '100'  '-R' '$REF'  '-nct' '36'  '-o' '$WORKING_DIRoutput.bam'
'java'  $JAVA_SETTINGS '$JAVA_TEMP'  '-jar' '$GATK'  '-T' 'BaseRecalibrator'  '-I' '$WORKING_DIRoutput.realigned.bam'  '-ip' '100'  '-R' '$REF'  '-BQSR' '$WORKING_DIRoutput.recal.before.table'  '-nct' '36'  '-knownSites' '$KNOWN_INDELS_MILLS' '-knownSites' '$KNOWN_INDELS' '-knownSites' '$KNOWN_DBSNP'  '-o' '$WORKING_DIRoutput.recal.after.table'
'java'  $JAVA_SETTINGS '$JAVA_TEMP'  '-jar' '/usr/share/docker_data/root/work/utils/picard-1.130/picard.jar'  'CollectMultipleMetrics'  'TMP_DIR=$WORKING_DIRtmp'  'VALIDATION_STRINGENCY=STRICT'  'MAX_RECORDS_IN_RAM=8000000'  'INPUT=$WORKING_DIRoutput.bam'  'OUTPUT=$WORKING_DIRstats/output.metrics'  'PROGRAM=' 'QualityScoreDistribution' 'PROGRAM=' 'MeanQualityByCycle' 'PROGRAM=' 'CollectAlignmentSummaryMetrics' 'PROGRAM=' 'CollectInsertSizeMetrics' 'PROGRAM=' 'CollectBaseDistributionByCycle'
'java'  $JAVA_SETTINGS '$JAVA_TEMP'  '-jar' '$GATK'  '-T' 'AnalyzeCovariates'  '-ip' '100'  '-R' '$REF'  '-before' '$WORKING_DIRoutput.recal.before.table'  '-after' '$WORKING_DIRoutput.recal.after.table'  '-plots' '$WORKING_DIRstats/output.recal.pdf'  '-csv' '$WORKING_DIRstats/output.recal.csv'

