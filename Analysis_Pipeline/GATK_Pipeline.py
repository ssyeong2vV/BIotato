import pandas as pd
import os
from datetime import datetime
#-----------------------------------------------------------------------------#
Sample = pd.read_csv("SampleSheet.txt", sep="\t", header=None)
Name = Sample.iloc[0, 0]
R1 = Sample.iloc[0, 1]
R2 = Sample.iloc[0, 2]
barcode = f'{datetime.now().strftime("%Y%m%d_%H%M%S")}_{Name}'
if os.path.isfile("ControlSheet.txt") :
    Control = pd.read_csv("ControlSheet.txt", sep="\t", header=None)
    CName = Control.iloc[0,0]
    CR1 = Control.iloc[0,1]
    CR2 = Control.iloc[0,2]
else :
    pass
#-----------------------------------------------------------------------------#
BATCH = {}
with open(f'{Name}_batchconfig.txt','r') as batch :
        for info in batch :
                info = info.strip()
                splitted = info.split('=')
                BATCH[splitted[0]] = splitted[1]
#-----------------------------------------------------------------------------#
def PreQC(r1, r2):
    if os.path.isdir("00.PreQC"):
        pass
    else:
        command = "mkdir 00.PreQC"
        os.system(command)

    command =f"fastqc -o 00.PreQC \
            -t {BATCH['CPU']} \
            {r1}\
            {r2}"
    os.system(command)
#-----------------------------------------------------------------------------#
def Trimming(r1, r2) :
    if os.path.isdir("01.Trimmed"):
        pass
    else:
        command = "mkdir 01.Trimmed"
        os.system(command)
    command = f"trimmomatic PE -phred33 -threads {BATCH['CPU']} \
        {r1} {r2} \
        ./01.Trimmed/{Name}_trim_paired_R1.fastq ./01.Trimmed/{Name}_trim_unpaired_R1.fastq \
        ./01.Trimmed/{Name}_trim_paired_R2.fastq ./01.Trimmed/{Name}_trim_unpaired_R2.fastq \
        ILLUMINACLIP:/home/lsy/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:{BATCH['maximum_mismatch_counts']}:{BATCH['PE_adapter_reads_similarity']}:{BATCH['SE_adapter_reads_similarity']} \
        LEADING:{BATCH['Leading']} TRAILING:{BATCH['Trailing']} SLIDINGWINDOW:{BATCH['WindowSize']}:{BATCH['meanQ']} MINLEN:{BATCH['minLen']}"
    os.system(command)
#-----------------------------------------------------------------------------#
def PostQC(name) :
    if os.path.isdir("02.PostQC"):
        pass
    else:
        command = "mkdir 02.PostQC"
        os.system(command)
    command =f"fastqc -o 02.PostQC \
            -t {BATCH['CPU']} \
            ./01.Trimmed/{name}_trim_paired_R1.fastq \
            ./01.Trimmed/{name}_trim_paired_R2.fastq"
    os.system(command)
#-----------------------------------------------------------------------------#
def Mapping(name) : #/Bioinformatics/01.Reference/b37/b37.fa
    if os.path.isdir("03.Mapped"):
        pass
    else:
        command = "mkdir 03.Mapped"
        os.system(command)
    command = f"/Bioinformatics/00.Tools/bwa-0.7.17/bwa mem -t {BATCH['CPU']} \
        {BATCH['Ref_path']} \
        ./01.Trimmed/{name}_trim_paired_R1.fastq \
        ./01.Trimmed/{name}_trim_paired_R2.fastq | samtools view -Sb > ./03.Mapped/{Name}_trimmed_bwa.bam"
    os.system(command)
#-----------------------------------------------------------------------------#
def AReadGroup(name) : #/media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar
    command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar AddOrReplaceReadGroups \
        -I ./03.Mapped/{name}_trimmed_bwa.bam \
        -O ./03.Mapped/{name}_trimmed_bwa_rg.bam \
        -LB {BATCH['LB']} \
        -PL {BATCH['PL']} \
        -PU {barcode} \
        -SM {name}"
    os.system(command)
#-----------------------------------------------------------------------------#
def Sorting(name) :
    command = f"samtools sort -@ {BATCH['CPU']} ./03.Mapped/{name}_trimmed_bwa_rg.bam \
        -o ./03.Mapped/{name}_trimmed_bwa_rg_sorted.bam"
    os.system(command)
#-----------------------------------------------------------------------------#
def Deduplication(name) :
    command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar MarkDuplicates \
        -I ./03.Mapped/{name}_trimmed_bwa_rg_sorted.bam \
        -O ./03.Mapped/{name}_trimmed_bwa_rg_sorted_dedup.bam \
        -M ./03.Mapped/{name}_trimmed_bwa_rg_sorted_dedup.metrics"
    os.system(command)
#-----------------------------------------------------------------------------#
def BaseRecalibration(name) :
    command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar BaseRecalibrator \
        -I ./03.Mapped/{name}_trimmed_bwa_rg_sorted_dedup.bam \
        --known-sites /Bioinformatics/01.Reference/hg19/Homo_sapiens_assembly19.dbsnp138.vcf \
        --known-sites /media/src/hg19/03.db/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
        -R {BATCH['Ref_path']} \
        -O ./03.Mapped/{name}_trimmed_bwa_rg_sorted_dedup_recal.table"
    os.system(command)

    command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar ApplyBQSR \
        -I ./03.Mapped/{name}_trimmed_bwa_rg_sorted_dedup.bam \
        -R {BATCH['Ref_path']} \
        --bqsr-recal-file ./03.Mapped/{name}_trimmed_bwa_rg_sorted_dedup_recal.table \
        -O ./03.Mapped/{name}_trimmed_bwa_rg_sorted_dedup_recal.bam"
    os.system(command)
#-----------------------------------------------------------------------------#
def Vcalling() :
    if os.path.isdir("04.Variants"):
        pass
    else:
        command = "mkdir 04.Variants"
        os.system(command)

    if BATCH['Type'] == 'Somatic' :
        if BATCH['Run_type'] =='WGS' :
            if os.path.isfile('ControlSheet.txt') :
                command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar Mutect2 \
                    -R {BATCH['Ref_path']} \
                    -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal.bam \
                    -I ./03.Mapped/{CName}_trimmed_bwa_rg_sorted_dedup_recal.bam \
                    -normal {CName} \
                    --germline-resource /media/src/hg19/03.db/af-only-gnomad.raw.sites.vcf.gz \
                    --panel-of-normals /Bioinformatics/01.Reference/hg19/Mutect2-WGS-panel-b37.vcf\
                    -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect.vcf.gz"
                os.system(command)
            else :
                command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar Mutect2 \
                    -R {BATCH['Ref_path']} \
                    -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal.bam \
                    --germline-resource /media/src/hg19/03.db/af-only-gnomad.raw.sites.vcf.gz \
                    --panel-of-normals /Bioinformatics/01.Reference/hg19/Mutect2-WGS-panel-b37.vcf\
                    -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect.vcf.gz"
                os.system(command)
        elif BATCH['Run_type'] == 'TARGET' :
            if os.path.isfile('ControlSheet.txt') :
                command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar Mutect2 \
                    -R {BATCH['Ref_path']} \
                    -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal.bam \
                    -I ./03.Mapped/{CName}_trimmed_bwa_rg_sorted_dedup_recal.bam \
                    -normal {CName} \
                    --germline-resource /media/src/hg19/03.db/af-only-gnomad.raw.sites.vcf.gz \
                    --panel-of-normals /Bioinformatics/01.Reference/hg19/Mutect2-WGS-panel-b37.vcf\
                    -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect.vcf.gz \
                    -L {BATCH['Panel']}"
                os.system(command)
            else :
                command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar Mutect2 \
                    -R {BATCH['Ref_path']} \
                    -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal.bam \
                    --germline-resource /media/src/hg19/03.db/af-only-gnomad.raw.sites.vcf.gz \
                    --panel-of-normals /Bioinformatics/01.Reference/hg19/Mutect2-WGS-panel-b37.vcf\
                    -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect.vcf.gz \
                    -L {BATCH['Panel']}"
                os.system(command)
    elif BATCH['Type'] =='Germline' :
        if BATCH['Run_type'] =='WGS' :
            command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller \
                -R {BATCH['Ref_path']} \
                -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal.bam \
                -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_g.vcf.gz -ERC GVCF"
            os.system(command)
            # CombineGVCFs 생략 -> 용량이 더 큰 샘플에서 앞의 모든 과정이 종료될 때까지 다른 샘플들도 대기해야함 ?
            # 추후에 Trio 분석 등을 위해 gvcf 생성
            # GenotypeGVCFs
            command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar GenotypeGVCFs \
                -R {BATCH['Ref_path']} \
                -V ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_g.vcf.gz \
                -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller.vcf.gz"
            os.system(command)
        elif BATCH['Run_type'] == 'TARGET' :
            command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar HaplotypeCaller \
                -R {BATCH['Ref_path']} \
                -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal.bam \
                -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_g.vcf.gz -ERC GVCF \
                -L {BATCH['Panel']}"
            os.system(command)

            command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar GenotypeGVCFs \
                -R {BATCH['Ref_path']} \
                -V ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_g.vcf.gz \
                -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller.vcf.gz \
                -L {BATCH['Panel']}"
            os.system(command)
#-----------------------------------------------------------------------------#
def Vfiltering() :
    if BATCH['Type'] == 'Somatic' :
        # Estimation contamination
        command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar GetPileupSummaries \
            -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal.bam \
            -L /media/src/hg19/03.db/af-only-gnomad.raw.sites.vcf.gz \
            -O ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_pileup.table \
            -V /media/src/hg19/03.db/af-only-gnomad.raw.sites.vcf.gz"
        os.system(command)
        command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar CalculateContamination \
            -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_pileup.table \
            -O ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_pileup_contam.table"
        os.system(command)
        # F1R2 model (Estimation strand bias)
        command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar CollectF1R2Counts \
            -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal.bam \
            -O ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_F1R2.tar.gz \
            -R {BATCH['Ref_path']}"
        os.system(command)
        command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar LearnReadOrientationModel \
            -I ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_F1R2.tar.gz \
            -O ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_F1R2_model.tar.gz"
        os.system(command)
        # Filter Mutect2 Call
        if BATCH['Run_type'] == 'WGS' :
            command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar FilterMutectCalls \
                -R {BATCH['Ref_path']} \
                -V ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect.vcf.gz \
                -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered.vcf.gz \
                --contamination-table ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_pileup_contam.table \
                -ob-priors ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_F1R2_model.tar.gz"
            os.system(command)
        elif BATCH['Run_type'] == 'TARGET' :
            command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar FilterMutectCalls \
                -R {BATCH['Ref_path']} \
                -V ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect.vcf.gz \
                -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered.vcf.gz \
                --contamination-table ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_pileup_contam.table \
                -ob-priors ./03.Mapped/{Name}_trimmed_bwa_rg_sorted_dedup_recal_F1R2_model.tar.gz"
            os.system(command)

    elif BATCH['Type'] =='Germline' :
        if BATCH['Run_type'] == 'WGS' :
            command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar VariantFiltration \
                -R {BATCH['Ref_path']} \
                -V ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller.vcf.gz \
                -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered.vcf.gz \
                --filter-name 'LowQD' --filter-expression 'QD < {BATCH['QD']}' \
                --filter-name 'LowMQ' --filter-expression 'MQ < {BATCH['MQ']}' \
                --filter-name 'LowMQRankSum' --filter-expression 'MQRankSum < {BATCH['MQRankSum']}' \
                --filter-name 'LowReadPosRankSum' --filter-expression 'ReadPosRankSum < {BATCH['ReadPosRankSum']}' \
                --filter-name 'HighFS' --filter-expression 'FS > {BATCH['HighFS']}' \
                --filter-name 'HighSOR' --filter-expression 'SOR > {BATCH['HighSOR']}'"
            os.system(command)
        elif BATCH['Run_type'] == 'TARGET' :
            command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar VariantFiltration \
                -R {BATCH['Ref_path']} \
                -V ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller.vcf.gz \
                -O ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered.vcf.gz \
                --filter-name 'LowQD' --filter-expression 'QD < {BATCH['QD']}' \
                --filter-name 'LowMQ' --filter-expression 'MQ < {BATCH['MQ']}' \
                --filter-name 'LowMQRankSum' --filter-expression 'MQRankSum < {BATCH['MQRankSum']}' \
                --filter-name 'LowReadPosRankSum' --filter-expression 'ReadPosRankSum < {BATCH['ReadPosRankSum']}' \
                --filter-name 'HighFS' --filter-expression 'FS > {BATCH['HighFS']}' \
                --filter-name 'HighSOR' --filter-expression 'SOR > {BATCH['HighSOR']}' \
                -L {BATCH['Panel']}"
            os.system(command)
#-----------------------------------------------------------------------------#
def Vannotation() :
    # Funcotator
    # if BATCH['Type'] == 'Somatic' :
    #     command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar Funcotator \
    #     --variant ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered.vcf.gz \
    #     --reference {BATCH['Ref_path']} \
    #     --ref-version {BATCH['Ref_ver']} \
    #     --data-sources-path /media/src/hg19/06.Annotation/funcotator_dataSources.v1.7.20200521s \
    #     --output ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered_funcotated.vcf \
    #     --output-file-format VCF"
    #     os.system(command)
    # elif BATCH['Type'] == 'Germline' :
    #     command = f"java -jar /media/src/Tools/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar Funcotator \
    #     --variant ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered.vcf.gz \
    #     --reference {BATCH['Ref_path']} \
    #     --ref-version {BATCH['Ref_ver']} \
    #     --data-sources-path /media/src/hg19/06.Annotation/funcotator_dataSources.v1.7.20200521g \
    #     --output ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered_funcotated.vcf \
    #     --output-file-format VCF"
    #     os.system(command)
    if BATCH['Type'] == 'Somatic' :
        command = f"java -Xmx4g -jar /media/src/Tools/snpEff/snpEff.jar -v hg19 \
        ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered.vcf.gz \
        > ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered_snpeff.vcf"
        os.system(command)

        command = f"perl /Bioinformatics/00.Tools/annovar/convert2annovar.pl -format vcf4 \
        ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered_snpeff.vcf \
        -includeinfo \
        -outfile ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered_snpeff.avinput"
        os.system(command)

        command = f"perl /Bioinformatics/00.Tools/annovar/table_annovar.pl \
        ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered_snpeff.avinput \
        /Bioinformatics/00.Tools/annovar/humandb/ \
        -buildver {BATCH['Ref_ver']} \
        -out ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered_snpeff_annovar \
        -remove -protocol refGeneWithVer,clinvar_20240917,gnomad211_exome,avsnp147,dbnsfp47a,cytoBand \
        -operation g,f,f,f,f,r \
        -otherinfo \
        -nastring . -polish"
        os.system(command)

        command = f"python /labmed/01.ALL/03.python/lsy/00.CODE/VCF_result.py ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered_snpeff.vcf ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_mutect_filtered_snpeff_annovar.{BATCH['Ref_ver']}_multianno.txt ./{Name}_Mutect2"
        os.system(command)

    elif BATCH['Type'] == 'Germline' :
        command = f"java -Xmx4g -jar /media/src/Tools/snpEff/snpEff.jar -v hg19 \
        ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered.vcf.gz \
        > ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered_snpeff.vcf"
        os.system(command)

        command = f"perl /Bioinformatics/00.Tools/annovar/convert2annovar.pl -format vcf4 \
        ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered_snpeff.vcf \
        -includeinfo \
        -outfile ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered_snpeff.avinput"
        os.system(command)
        
        command = f"perl /Bioinformatics/00.Tools/annovar/table_annovar.pl \
        ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered_snpeff.avinput \
        /Bioinformatics/00.Tools/annovar/humandb/ \
        -buildver {BATCH['Ref_ver']} \
        -out ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered_snpeff_annovar \
        -remove -protocol refGeneWithVer,clinvar_20240917,gnomad211_exome,avsnp147,dbnsfp47a,cytoBand \
        -operation g,f,f,f,f,r \
        -otherinfo \
        -nastring . -polish"
        os.system(command)

        command = f"python /labmed/01.ALL/03.python/lsy/00.CODE/VCF_result.py ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered_snpeff.vcf ./04.Variants/{Name}_trimmed_bwa_rg_sorted_dedup_recal_HCaller_filtered_snpeff_annovar.{BATCH['Ref_ver']}_multianno.txt ./{Name}_HCaller"
        os.system(command)
#-----------------------------------------------------------------------------#
if BATCH['Type'] == 'Somatic' :
    if os.path.isfile('ControlSheet.txt') :
        if BATCH['Step'] == 'ALL' :
            PreQC(R1,R2)
            PreQC(CR1,CR2)
            Trimming(R1,R2)
            Trimming(CR1,CR2)
            PostQC(Name)
            PostQC(CName)
            Mapping(Name)
            Mapping(CName)
            AReadGroup(Name)
            AReadGroup(CName)
            Sorting(Name)
            Sorting(CName)
            Deduplication(Name)
            Deduplication(CName)
            BaseRecalibration(Name)
            BaseRecalibration(CName)
            Vcalling()
            Vfiltering()
            Vannotation()

        elif BATCH['Step'] =='PreQC' :
            PreQC(R1,R2)
            PreQC(CR1,CR2)

        elif BATCH['Step'] =='Trimming' :
            Trimming(R1,R2)
            Trimming(CR1,CR2)

        elif BATCH['Step'] == 'PostQC' :
            PostQC(Name)
            PostQC(CName)

        elif BATCH['Step'] == 'Mapping' :
            Mapping(Name)
            Mapping(CName)

        elif BATCH['Step'] == 'AReadGroup' :
            AReadGroup(Name)
            AReadGroup(CName)

        elif BATCH['Step'] == 'Sorting' :
            Sorting(Name)
            Sorting(CName)

        elif BATCH['Step'] == 'Deduplication' :
            Deduplication(Name)
            Deduplication(CName)

        elif BATCH['Step'] == 'BaseRecalibration' :
            BaseRecalibration(Name)
            BaseRecalibration(CName)

        elif BATCH['Step'] == 'Vcalling' :
            Vcalling()
            Vfiltering()
            Vannotation()

        elif BATCH['Step'] == 'Vfiltering' :
            Vfiltering()
            Vannotation()

        elif BATCH['Step'] == 'Vannotation' :
            Vannotation()

    else :
        if BATCH['Step'] == 'ALL' :
            PreQC(R1,R2)
            Trimming(R1,R2)
            PostQC(Name)
            Mapping(Name)
            AReadGroup(Name)
            Sorting(Name)
            Deduplication(Name)
            BaseRecalibration(Name)
            Vcalling()
            Vfiltering()
            Vannotation()

        elif BATCH['Step'] =='PreQC' :
            PreQC(R1,R2)

        elif BATCH['Step'] =='Trimming' :
            Trimming(R1,R2)

        elif BATCH['Step'] == 'PostQC' :
            PostQC(Name)

        elif BATCH['Step'] == 'Mapping' :
            Mapping(Name)

        elif BATCH['Step'] == 'AReadGroup' :
            AReadGroup(Name)

        elif BATCH['Step'] == 'Sorting' :
            Sorting(Name)

        elif BATCH['Step'] == 'Deduplication' :
            Deduplication(Name)

        elif BATCH['Step'] == 'BaseRecalibration' :
            BaseRecalibration(Name)

        elif BATCH['Step'] == 'Vcalling' :
            Vcalling()
            Vfiltering()
            Vannotation()

        elif BATCH['Step'] == 'Vfiltering' :
            Vfiltering()
            Vannotation()

        elif BATCH['Step'] == 'Vannotation' :
            Vannotation()
#-----------------------------------------------------------------------------#
elif BATCH['Type'] =='Germline' :
    if BATCH['Step'] == 'ALL' :
        PreQC(R1,R2)
        Trimming(R1,R2)
        PostQC(Name)
        Mapping(Name)
        AReadGroup(Name)
        Sorting(Name)
        Deduplication(Name)
        BaseRecalibration(Name)
        Vcalling()
        Vfiltering()
        Vannotation()

    elif BATCH['Step'] =='PreQC' :
        PreQC(R1,R2)

    elif BATCH['Step'] =='Trimming' :
        Trimming(R1,R2)

    elif BATCH['Step'] == 'PostQC' :
        PostQC(Name)

    elif BATCH['Step'] == 'Mapping' :
        Mapping(Name)

    elif BATCH['Step'] == 'AReadGroup' :
        AReadGroup(Name)

    elif BATCH['Step'] == 'Sorting' :
        Sorting(Name)

    elif BATCH['Step'] == 'Deduplication' :
        Deduplication(Name)

    elif BATCH['Step'] == 'BaseRecalibration' :
        BaseRecalibration(Name)

    elif BATCH['Step'] == 'Vcalling' :
        Vcalling()
        Vfiltering()
        Vannotation()

    elif BATCH['Step'] == 'Vfiltering' :
        Vfiltering()
        Vannotation()

    elif BATCH['Step'] == 'Vannotation' :
        Vannotation()
#-----------------------------------------------------------------------------#
