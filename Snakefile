import configparser
import yaml
import re

fastq_dir = config['all']['fastq_dir']
anvio_output_dir = config['all']['root'] + '/anvio_output'
coassembly_yaml = config['all']['coassembly_yaml']
workdir: config['all']['root']

#######
### list of samples
### obtained from coassembly yaml file
#######

with open(coassembly_yaml, 'r') as file:
    coassemble_group = yaml.load(file, Loader = yaml.FullLoader)
    samples = list(coassemble_group.values())[0]    

rule all:
    input:
        anvio_output_dir + '/merged_profile/PROFILE.db'

rule merge_profile:
    input: 
        contig_db = ancient(anvio_output_dir + '/contigs.db')
    output:
        anvio_output_dir + '/merged_profile/PROFILE.db'
    params:
        profile_dir = anvio_output_dir + '/profile',
        outdir = anvio_output_dir + '/merged_profile'
    shell:
        """
            profiles=$(find {params.profile_dir} -name *PROFILE.db)
            anvi-merge ${{profiles}} -o {params.outdir} -c {input.contig_db} 
        """

rule profile:
    input:
        contig_db = ancient(anvio_output_dir + '/contigs.db'),
        bam = anvio_output_dir + '/{sample}/sorted.bam',
        sentinel = anvio_output_dir + '/.DONEcentrifuge'
    output:
        anvio_output_dir + '/profile/{sample}/PROFILE.db'
    params:
        anvio_output_dir + '/profile/{sample}'
    threads:
        config['profile']['threads'] 
    shell:
        """
            anvi-profile -i {input.bam} -o {params} -c {input.contig_db} -T {threads} -W
        """

rule import_taxonomy:
    input:
        contig_db = ancient(anvio_output_dir + '/contigs.db'),
        hits = anvio_output_dir + '/centrifuge_hits.tsv',
        report = anvio_output_dir + '/centrifuge_report.tsv'
    output:
        anvio_output_dir + '/.DONEcentrifuge'
    shell:
        """
            anvi-import-taxonomy-for-genes \
                -c {input.contig_db} -i {input.report} {input.hits} -p centrifuge && \
            touch {output}
        """

rule centrifuge:
    input:
        anvio_output_dir + '/gene-calls.fa'
    output:
        hits = anvio_output_dir + '/centrifuge_hits.tsv',
        report = anvio_output_dir + '/centrifuge_report.tsv'
    params:
        config['centrifuge']['centrifuge_db']
    threads:
        config['centrifuge']['threads']
    shell:
        """
            centrifuge -f -x {params} {input} -S {output.hits} --report-file {output.report} -p {threads}
        """

rule get_seq_for_gene_calls:
    input:
        contig_db = ancient(anvio_output_dir + '/contigs.db'),
        sentinel = anvio_output_dir + '/.DONEncbi_cog'
    output:
        anvio_output_dir + '/gene-calls.fa'
    shell:
        """
            anvi-get-sequences-for-gene-calls \
                -c {input.contig_db} -o {output}
        """

########
### If you are running COGs for the first time, 
### you will need to set them up on your computer using 
### anvi-setup-ncbi-cogs
### Refer to http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-ncbi-cogs
######## 

rule ncbi_cogs:
    input:
        contig_db = ancient(anvio_output_dir + '/contigs.db'),
        sentinel = anvio_output_dir + '/.DONEhmm'
    output:
        anvio_output_dir + '/.DONEncbi_cog'
    threads:
        config['ncbi_cog']['threads']
    shell:
        """
            anvi-run-ncbi-cogs --contigs-db {input.contig_db} \
                --num-threads {threads} --search-with blastp && \
            touch {output}
        """

rule run_hmm:
    input:
        ancient(anvio_output_dir + '/contigs.db')
    output:
        anvio_output_dir + '/.DONEhmm'
    threads:
        config['hmm']['threads']
    shell:
        """
            anvi-run-hmms --contigs-db {input} \
                --num-threads {threads} && \
            touch {output}
        """

rule get_contig_db:
    input:
        anvio_output_dir + '/reformatted_contigs.fa'
    output:
        anvio_output_dir + '/contigs.db'
    params:
        config['all']['project_name']
    shell:
        """
            anvi-gen-contigs-database --contigs-fasta {input} \
                --output-db-path {output} \
                --project-name {params}
        """

rule contigs_mapping:
    input:
        contig = anvio_output_dir + '/reformatted_contigs.fa',
        reads = expand(fastq_dir + '/{{sample}}_{rp}.fastq.gz', rp = [1, 2])
    output:
        bam = anvio_output_dir + '/{sample}/sorted.bam',
        bai = anvio_output_dir + '/{sample}/sorted.bam.bai'
    params:
        temp = temp(anvio_output_dir + '/{sample}/sorted.tmp')
    threads: 
        config['mapping']['threads']
    shell:
        """
        ~/miniconda3/envs/sunbeam/bin/minimap2 -ax sr -t {threads} {input.contig} {input.reads} | \
            ~/miniconda3/envs/sunbeam/bin/samtools sort -@ {threads} -o {output.bam} -T {params.temp} - 
        ~/miniconda3/envs/sunbeam/bin/samtools index {output.bam} {output.bai}
        """

rule reformat_fasta:
    input: 
        config['all']['contig']
    output: 
        anvio_output_dir + '/reformatted_contigs.fa'
    params:
        min_contig_length = config['reformat']['min_contig_length'],
        report_file = anvio_output_dir + '/simplify-names.txt'
    shell:
        """
            anvi-script-reformat-fasta {input} \
                --output-file {output} \
                --min-len {params.min_contig_length} \
                --simplify-names \
                --report-file {params.report_file}
        """

onsuccess:
    print('Workflow finished, no error')
    shell('mail -s "Workflow finished successfully" ' + config['all']['admin_email'] + ' < {log}')

onerror:
    print('An error occurred')
    shell('mail -s "An error occurred" ' + config['all']['admin_email'] + ' < {log}')

