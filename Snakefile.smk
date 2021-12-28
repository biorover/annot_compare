import os,sys

config: 'config.yml'
genome_dir = config['genome_dir']
species_dict = config['species']
outdir = config.get('outdir','results')
happy_hmms = config['happy_hmms']

#rules
rule all:
    input:
        expand(outdir + '/HAPpy-ABCENTH/{species}/ABCENTH.gtf',species = species_dict.keys() )
        #expand(outdir + '/exonerate/{species}/exonerate.txt',species = species_dict.keys() )

rule happy_abcenth:
    output:
        gtf = outdir + '/HAPpy-ABCENTH/{species}/ABCENTH.gtf',
        pep = outdir + '/HAPpy-ABCENTH/{species}/ABCENTH.pep',
    conda: 'envs/abcenth.yml'
    threads: 36
    params:
        target = lambda w: species_dict[w.species],
    shell:
        """
        python ~/tools/HAPpy-ABCENTH/HAP.py \
            --genome {params.target} 
            --hmm_dir {happy_hmms} \
            --threads {threads} \
            --annotator ABCENTH \
            --output_dir $(dirname {output.gtf}) \
            --overwrite True \
            --ABCENTH_options '--orf_finder_E 0.01 --full_pseudoexon_search True'
        
        MAGOT gff2fasta {output.gtf} {params.target} --seq_type aa | sed -e 's/>/>{species}_/' > {output.pep}
        """