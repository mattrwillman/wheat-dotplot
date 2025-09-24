#This workflow aligns a query genome assembly against a reference assembly.
#To manage the large size of the wheat genome, the full query assembly is aligned to a single chromosome of the reference assembly.

import os

QUERY_NAME = config['QUERY_NAME']
REF_NAME = config['REF_NAME']
QUERY_FASTA = config['QUERY_FASTA']
REF_FASTA = config['REF_FASTA']
CHRS = config['CHRS']
BREAKLEN = config['BREAKLEN']

rule all:
    input:
        expand("results/ref_chrs/{chr}.fa", chr = CHRS), 
        expand("results/plots/20kb/chr_full/{ref}_{query}_{chr}.png", ref = REF_NAME, query = QUERY_NAME, chr = CHRS), 
        expand("results/filter20kb_coords/{ref}_{query}_{chr}.coords", ref = REF_NAME, query = QUERY_NAME, chr = CHRS), 
        expand("results/Rplots/pairwise_{ref}_{query}_{chr}.coords.jpeg", ref = REF_NAME, query = QUERY_NAME, chr = CHRS)

rule change_gaps:
    input:
        ref = REF_FASTA, 
        query = QUERY_FASTA
    output:
        ref = expand("results/ref_adjusted_gaps/{ref}.fasta", ref = REF_NAME), 
        query = expand("results/query_adjusted_gaps/{query}.fasta", query = QUERY_NAME)
    params:
        breaklen = config['BREAKLEN']
    shell:
        "bash workflow/scripts/standardize_assembly_breaks.sh {input.ref} {output.ref} {params.breaklen} ; "
        "bash workflow/scripts/standardize_assembly_breaks.sh {input.query} {output.query} {params.breaklen}"

rule subset_ref_chrs:
    input:
        expand("results/ref_adjusted_gaps/{ref}.fasta", ref = REF_NAME)
    output:
        "results/ref_chrs/{chr}.fa"
    log:
        "logs/subset_ref_{chr}.log"
    shell:
        "seqkit grep -rp '{wildcards.chr}' {input} > {output} 2> {log}"

rule nucmer:
    input:
        ref = "results/ref_chrs/{chr}.fa", 
        sample = expand("results/query_adjusted_gaps/{query}.fasta", query = QUERY_NAME)
    output:
        expand("results/nucmer/{ref}_{query}_{{chr}}.delta", ref = REF_NAME, query = QUERY_NAME)
    params:
        prefix = expand("results/nucmer/{ref}_{query}_{{chr}}", ref = REF_NAME, query = QUERY_NAME), 
        minmatch = "100", 
        mincluster = "1000", 
        breaklen = config['BREAKLEN']
    log:
        "logs/nucmer_{chr}.log"
    threads: 12
    shell:
        "nucmer --prefix {params.prefix} "
        "-l {params.minmatch} "
        "-c {params.mincluster} "
        "-b {params.breaklen} "
        "--threads={threads} "
        "--mum "
        "{input.ref} "
        "{input.sample} "
        "2> {log}"

rule filter:
    input:
        expand("results/nucmer/{ref}_{query}_{{chr}}.delta", ref = REF_NAME, query = QUERY_NAME)
    output:
        expand("results/filter20kb/{ref}_{query}_{{chr}}.delta", ref = REF_NAME, query = QUERY_NAME)
    log:
        "logs/filter20kb_{chr}.log"
    threads: 12
    shell:
        "delta-filter -r -q -l 20000 {input} > {output} 2> {log}"

rule extract_coords:
    input:
        expand("results/filter20kb/{ref}_{query}_{{chr}}.delta", ref = REF_NAME, query = QUERY_NAME)
    output:
        coords = expand("results/filter20kb_coords/{ref}_{query}_{{chr}}.coords", ref = REF_NAME, query = QUERY_NAME), 
    log:
        "logs/extract20kb_{chr}.log"
    threads: 12
    shell:
        "show-coords -r -c -l {input} > {output.coords} 2> {log} ; "

rule extract_snps:
    input:
        expand("results/filter20kb/{ref}_{query}_{{chr}}.delta", ref = REF_NAME, query = QUERY_NAME)
    output:
        snps = expand("results/filter20kb_snps/{ref}_{query}_{{chr}}.snps", ref = REF_NAME, query = QUERY_NAME)
    log:
        "logs/extract20kb_{chr}.log"
    threads: 12
    shell:
        "show-snps -C {input} > {output.snps} 2>> {log}"

rule plotscript:
    input:
        expand("results/filter20kb/{ref}_{query}_{{chr}}.delta", ref = REF_NAME, query = QUERY_NAME)
    output:
        gp_full = expand("results/plots/20kb/chr_full/{ref}_{query}_{{chr}}.gp", ref = REF_NAME, query = QUERY_NAME), 
        gp_chr = expand("results/plots/20kb/chr_chr/{ref}_{query}_{{chr}}.gp", ref = REF_NAME, query = QUERY_NAME)
    params:
        prefix_full = expand("results/plots/20kb/chr_full/{ref}_{query}_{{chr}}", ref = REF_NAME, query = QUERY_NAME), 
        prefix_chr = expand("results/plots/20kb/chr_chr/{ref}_{query}_{{chr}}", ref = REF_NAME, query = QUERY_NAME), 
        chr = "{chr}", 
        title = expand("{query} (y) aligned to {ref} (x)", query = QUERY_NAME, ref = REF_NAME)
    threads: 12
    log: 
    shell:
        "mummerplot --png -large -title '{params.title}' -p {params.prefix_full} {input} ; "
        "mummerplot --png -large -title '{params.title}' -p {params.prefix_chr} -q {params.chr} {input}"

rule gnuplot:
    input:
        gp_full = expand("results/plots/20kb/chr_full/{ref}_{query}_{{chr}}.gp", ref = REF_NAME, query = QUERY_NAME), 
        gp_chr = expand("results/plots/20kb/chr_chr/{ref}_{query}_{{chr}}.gp", ref = REF_NAME, query = QUERY_NAME)
    output:
        gp_full = expand("results/plots/20kb/chr_full/{ref}_{query}_{{chr}}_edit.gp", ref = REF_NAME, query = QUERY_NAME), 
        gp_chr = expand("results/plots/20kb/chr_chr/{ref}_{query}_{{chr}}_edit.gp", ref = REF_NAME, query = QUERY_NAME), 
        plot_full = expand("results/plots/20kb/chr_full/{ref}_{query}_{{chr}}.png", ref = REF_NAME, query = QUERY_NAME), 
        plot_chr = expand("results/plots/20kb/chr_chr/{ref}_{query}_{{chr}}.png", ref = REF_NAME, query = QUERY_NAME)
    shell:
        "sed 's/set terminal png tiny size 1400,1400/set terminal png large size 1000,1000/' {input.gp_full} > {output.gp_full} ; "
        "sed 's/set terminal png tiny size 1400,1400/set terminal png large size 1000,1000/' {input.gp_chr} > {output.gp_chr} ; "
        "gnuplot {output.gp_full} > {output.plot_full} ; "
        "gnuplot {output.gp_chr} > {output.plot_chr}"

rule Rplot:
    input:
        coords = expand("results/filter20kb_coords/{ref}_{query}_{chr}.coords", ref = REF_NAME, query = QUERY_NAME, chr = CHRS)
    output:
        expand("results/Rplots/pairwise_{ref}_{query}_{chr}.coords.jpeg", ref = REF_NAME, query = QUERY_NAME, chr = CHRS)
    params:
        in_dir = "results/filter20kb_coords", 
        out_dir = "results/Rplots", 
        ref = REF_NAME, 
        query = QUERY_NAME
    shell:
        "Rscript workflow/scripts/plot_coords.R {params.in_dir} {params.out_dir} {params.ref} {params.query}"
