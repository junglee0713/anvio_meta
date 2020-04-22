import configparser
import yaml
import re

contig_dir = config['all']['contig_dir']
anvio_output_dir = config['all']['root'] + '/anvio_output'
genome_storage_fp = anvio_output_dir + '/' + config['all']['project_name'] + '-GENOMES.db'

files, = glob_wildcards(contig_dir + '/{file}.fa')

workdir: config['all']['root']

########
### function to create the external genome path tsv
########

def get_external_genome_path(contig_db_list, output_fp):
    with open(output_fp, 'w') as f:
        f.write('name\tcontigs_db_path\n')
        for db in contig_db_list:
            name = re.sub(anvio_output_dir + '/', '', db)
            name = re.sub('/contigs.db', '', name)
            name = re.sub('-contigs', '', name)
            name = re.sub('-|\.', '_', name) 
            f.write(name + '\t' + db + '\n')

rule all:
    input: 
        genome_storage_fp

rule build_genome_storage:
    input:
        anvio_output_dir + '/external_genome_path.tsv'
    output:
        genome_storage_fp
    shell:
        """
            anvi-gen-genomes-storage \
                --external-genomes {input} \
                --output-file {output}
        """

rule external_genome_path:
    input:
        contig_db = expand(anvio_output_dir + '/{file}/contigs.db', file = files),
        sentinel = expand(anvio_output_dir + '/{file}/.DONEncbi_cog', file = files)
    output:
        anvio_output_dir + '/external_genome_path.tsv'
    run:
        get_external_genome_path(input.contig_db, output[0])    

########
### If you are running COGs for the first time, 
### you will need to set them up on your computer using 
### anvi-setup-ncbi-cogs
### Refer to http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-ncbi-cogs
######## 

rule ncbi_cogs:
    input:
        contig_db = ancient(anvio_output_dir + '/{file}/contigs.db'),
        sentinel = anvio_output_dir + '/{file}/.DONEhmm'
    output:
        anvio_output_dir + '/{file}/.DONEncbi_cog'
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
        ancient(anvio_output_dir + '/{file}/contigs.db')
    output:
        anvio_output_dir + '/{file}/.DONEhmm'
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
        anvio_output_dir + '/{file}/reformatted_contigs.fa'
    output:
        anvio_output_dir + '/{file}/contigs.db'
    params:
        file = '{file}'
    shell:
        """
            anvi-gen-contigs-database --contigs-fasta {input} \
                --output-db-path {output} \
                --project-name {params.file}
        """

rule reformat_fasta:
    input: 
        contig_dir + '/{file}.fa'
    output: 
        anvio_output_dir + '/{file}/reformatted_contigs.fa'
    params:
        min_contig_length = config['reformat']['min_contig_length'],
        report_file = anvio_output_dir + '/{file}/simplify-names.txt'
    shell:
        """
            anvi-script-reformat-fasta {input} \
                --output-file {output} \
                --min-len {params.min_contig_length} \
                --simplify-name \
                --report-file {params.report_file}
        """

onsuccess:
    print('Workflow finished, no error')
    shell('mail -s "Workflow finished successfully" ' + config['all']['admin_email'] + ' < {log}')

onerror:
    print('An error occurred')
    shell('mail -s "An error occurred" ' + config['all']['admin_email'] + ' < {log}')

