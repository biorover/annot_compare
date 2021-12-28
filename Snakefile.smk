import os,sys

configfile: 'config.yml'
genome_dir = config['genome_dir']
species_dict = config['species']
outdir = config.get('outdir','results')
happy_hmms = config['happy_hmms']
max_threads = config.get('max_threads',36)
exonerate_queries = config['exonerate_queries']

#rules
rule all:
    input:
        expand(outdir + '/HAPpy-ABCENTH/{species}/ABCENTH.gtf',species = species_dict.keys() )
        #expand(outdir + '/exonerate/{species}/exonerate.txt',species = species_dict.keys() )

rule happy_abcenth:
    output:
        gtf = '{outdir}/HAPpy-ABCENTH/{species}/ABCENTH.gtf',
        pep = '{outdir}/HAPpy-ABCENTH/{species}/ABCENTH.pep',
    conda: 'envs/abcenth.yml'
    threads: max_threads
    params:
        target = lambda w: genome_dir + '/' + species_dict[w.species],
    shell:
        """
       HAPpy \
            --genome {params.target} \
            --hmm_dir {happy_hmms} \
            --threads {threads} \
            --annotator ABCENTH \
            --output_dir $(dirname {output.gtf}) \
            --overwrite True \
            --ABCENTH_options '--orf_finder_E 0.01 --full_pseudoexon_search True'
        
        MAGOT gff2fasta {output.gtf} {params.target} --seq-type lorfaa | sed -e 's/>/>{wildcards.species}_/' > {output.pep}
        """

rule cat_exonerate:
    input:
        expand('{{outdir}}/exonerate/{{species}}/{qsp}.exonerate.txt',qsp = exonerate_queries.keys() )
    output:
        '{outdir}/exonerate/{species}/exonerate.txt'
    shell:
        """
        cat {input} > {output}
        """

rule exonerate:
    output:
        outfile = '{outdir}/exonerate/{species}/{qsp}.exonerate.txt'
    threads: min((12,max_threads))
    params: 
        query = lambda w: exonerate_queries[w.qsp],
        target = lambda w: genome_dir + '/' + species_dict[w.species],
    conda: 'envs/py37.yml'
    shell:
        """
        exonerate --model protein2genome \
            -c {threads} \
            --maxintron 3600 \
            --showtargetgff T \
            {params.query} {params.target} > {output}}
        """