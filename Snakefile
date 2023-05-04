configfile: 'config.yaml'

from configparser import ConfigParser
from itertools import combinations

DATASET = ConfigParser()
DATASET.readfp(open(config['dataset_file']))

G = 'global'

BAM_DIR = config['bam_dir']
COOLER_BIN = config['cooler_bin']
GENOME = DATASET.get(G, 'genome_name')
HIC_MAPS = [x for x in DATASET.sections() if x != G]
HIC_MAPS_DIR = config['hic_maps_dir']
MIN_INTS = config['min_interactions']
NORM_ALGS = config['normalization_algorithms']
PRJ_DIR = config['hifive_project_dir']
REPLICATES = sorted(set(map(lambda x: x.split(':', 1)[0], HIC_MAPS)))
RESOLUTIONS = config['resolutions']
RESTRICTION_DIGEST_FILE = DATASET.get(G, 'fragments_file')
RESTRICTION_ENZYME = DATASET.get(G, 'restriction_enzyme')
RESTRICT_CHROMOSOMES = DATASET.get(G, 'restrict_chromosomes', fallback='None')
SAMTOOLS_BIN = config['samtools_bin']
SCRIPT_DIR = config['script_dir']
SUPER_RES = config['super_resolution_factors']
STATS_DIR = config['stats_dir']
TISSUES = sorted(set(map(lambda x: x.split(':', 1)[1], HIC_MAPS)))

FEND_FILE = '%s_%s.fend' %(GENOME, RESTRICTION_ENZYME)

rule all:
    input:
#        expand('%s/m{min_int}_r{res}_s{superres}:{alg}/{hic_map}.trv' %HIC_MAPS_DIR,
#                min_int = MIN_INTS, res = RESOLUTIONS, superres = SUPER_RES,
#                alg = NORM_ALGS, hic_map = HIC_MAPS),
        expand('%s/{sample}:{tissue}/m{min_int}_r{res}_s{superres}:{alg}.csv' %STATS_DIR,
                sample = map(lambda x: '-'.join(x), combinations(sorted(
                set(map(lambda x: x.split(':', 1)[0], HIC_MAPS))), 2)), 
                tissue = set(map(lambda x: x.split(':', 1)[1], HIC_MAPS)),
                min_int = MIN_INTS, res = RESOLUTIONS, superres = SUPER_RES,
                alg = NORM_ALGS),
        expand('%s/{alg}/m{min_int}_r{res}_s{superres}/{hic_map}.pdf' %HIC_MAPS_DIR,
                min_int = MIN_INTS, res = RESOLUTIONS[:2], superres = SUPER_RES,
                alg = NORM_ALGS, hic_map = HIC_MAPS + TISSUES),
        expand('%s/{alg}/m{min_int}_r{res}_s{superres}/{hic_map}_{chr}.pdf' %HIC_MAPS_DIR,
                min_int = MIN_INTS, res = RESOLUTIONS[:2], superres = SUPER_RES,
                alg = NORM_ALGS, hic_map = HIC_MAPS + TISSUES, chr=[1,3,4,5]),
        expand('%s/{alg}/{sample}:{tissue}.pdf' %STATS_DIR, 
                sample = map(lambda x: '-'.join(x), combinations(REPLICATES,
                2)), tissue = TISSUES, alg=NORM_ALGS),
        expand('%s/{alg}/{tissue}_m{min_int}_s{super_res}.cool' %HIC_MAPS_DIR,
                alg=NORM_ALGS, tissue=TISSUES, min_int=MIN_INTS,
                super_res=SUPER_RES),
        expand('%s/{alg}/{sample_tissue}_m{min_int}_s{super_res}.cool' %HIC_MAPS_DIR, 
                alg=NORM_ALGS, sample_tissue=HIC_MAPS, min_int=MIN_INTS,
                super_res=SUPER_RES)

rule digest_to_fend:
    input:
        RESTRICTION_DIGEST_FILE
    output:
        txt = temp('%s_%s.txt' %(GENOME, RESTRICTION_ENZYME)),
        fend = FEND_FILE
    shell:
        '%s/digest_to_fend.py {input} > {output.txt};' %SCRIPT_DIR +
        'hifive fends -g %s -r %s -F {output.txt} {output.fend}' %(GENOME, RESTRICTION_ENZYME)

rule link_bam_files:
    input:
        hic_maps = chain(*[map(lambda x: x[1], DATASET.items(s)) for s in
            HIC_MAPS])
    output:
        expand('%s/{hic_map}_{part}.orig.bam' %BAM_DIR, hic_map=HIC_MAPS,
                part = ('part1', 'part2'))
    run:
        from os import symlink, path
        for m in HIC_MAPS:
            symlink(path.relpath(DATASET.get(m, '1'), BAM_DIR), \
                    path.join(BAM_DIR, '%s_part1.orig.bam' %m))
            symlink(path.relpath(DATASET.get(m, '2'), BAM_DIR), \
                    path.join(BAM_DIR, '%s_part2.orig.bam' %m))


rule sort_bam_files:
    input:
        '%s/{hic_map}_{part}.orig.bam' %BAM_DIR
    output:
        '%s/{hic_map,[^_.]+}_{part,[^_.]+}.bam' %BAM_DIR
    shell:
        '%s sort {input} > {output}' %SAMTOOLS_BIN

rule index_bam_files:
    input:
        '%s/{hic_map}_{part}.bam' %BAM_DIR
    output:
        '%s/{hic_map,[^_.]+}_{part,[^_.]+}.bam.bai' %BAM_DIR
    shell:
        '%s index {input}' %SAMTOOLS_BIN
        
rule hic_data:
    input:
        part1 = '%s/{hic_map}_part1.bam' %BAM_DIR,
        part2 = '%s/{hic_map}_part2.bam' %BAM_DIR,
        idx1 = '%s/{hic_map}_part1.bam.bai' %BAM_DIR,
        idx2 = '%s/{hic_map}_part2.bam.bai' %BAM_DIR,
        fend = FEND_FILE
    output:
        '%s/{hic_map,.*:.*}.data' %PRJ_DIR
    log:
        'logs/hic-data_{hic_map}.log'
    shell:
        'hifive hic-data --bam {input.part1} {input.part2} {input.fend} '
        '{output} 2> {log}'


rule hic_data_full:
    input:
        part1 = lambda wildcards: list(map(lambda x: '%s/%s:%s_part1.bam' %(
                BAM_DIR, x, wildcards.tissue), REPLICATES)),
        part2 = lambda wildcards: list(map(lambda x: '%s/%s:%s_part2.bam' %(
                BAM_DIR, x, wildcards.tissue), REPLICATES)),
        idx1 = lambda wildcards: list(map(lambda x: '%s/%s:%s_part1.bam.bai' %(
                BAM_DIR, x, wildcards.tissue), REPLICATES)),
        idx2 = lambda wildcards: list(map(lambda x: '%s/%s:%s_part2.bam.bai' %(
                BAM_DIR, x, wildcards.tissue), REPLICATES)),
        fend = FEND_FILE
    output:
        '%s/{tissue}.data' %PRJ_DIR
    log:
        'logs/hic-data_{tissue}.log' 
    shell:
        'hifive hic-data --bam %s {input.fend} ' %(' --bam '.join(map(
                lambda x: '{input.part1[%s]} {input.part2[%s]}' %(x, x),
                range(len(REPLICATES))))) + 
        '{output} 2> {log}'


rule hic_project:
    input:
        '%s/{hic_map}.data' %PRJ_DIR
    params:
        chrom = RESTRICT_CHROMOSOMES
    output:
        '%s/m{min_int}/{hic_map}.prj' %PRJ_DIR
    log:
        'logs/hic-project_{hic_map}_m{min_int}.log'
    shell:
        'hifive hic-project -j 200 -n 250 -c {params.chrom} -f {wildcards.min_int} {input} '
        '{output} 2> {log}'


rule hic_normalize:
    input:
        '%s/m{min_int}/{hic_map}.prj' %PRJ_DIR
    params:
        chrom = RESTRICT_CHROMOSOMES
    output:
        '%s/m{min_int}/{hic_map}::{algorithm}.prj' %PRJ_DIR
    log:
        'logs/hic-normalize_{hic_map}_m{min_int}_{algorithm}.log'
    shell:
#        'hifive hic-normalize {wildcards.algorithm} -c {params.chrom} -r 100 '
#        '-y all -o {output} {input} 2> {log}'
# parameter settings for binning-probability algorithm
        'hifive hic-normalize {wildcards.algorithm} -c {params.chrom} -b 100 '
        '-y all -r 400 -o {output} {input} 2> {log}'
# parameter settings for probability normalization
#        'hifive hic-normalize {wildcards.algorithm} -c {params.chrom} -b 100 '
#        '-o {output} {input} 2> {log}'
# parameter settings for express normalization
#        'hifive hic-normalize {wildcards.algorithm} -f {wildcards.min_int} '
#        '-c {params.chrom} -e 400 -w all -o {output} {input} 2> {log}'


rule hic_heatmap:
    input:
        '%s/m{min_int}/{hic_map}::{alg}.prj' %PRJ_DIR
    params:
        chrom = RESTRICT_CHROMOSOMES
    output:
        map = '%s/{alg}/m{min_int}_r{res}_s{super_res}/{hic_map,[^/]+}.hdf5' %HIC_MAPS_DIR,
        img = '%s/{alg}/m{min_int}_r{res}_s{super_res}/{hic_map}.png' %HIC_MAPS_DIR
    log:
        'logs/hic-heatmap_{alg}_{hic_map}_m{min_int}_r{res}_s{super_res}.log'
    shell:
        'SUPER_RES=$(bc <<< \'{wildcards.res} * {wildcards.super_res}\');'
        'hifive hic-heatmap -t -d enrichment -f {wildcards.min_int} '
        '-b {wildcards.res} -c {params.chrom} -a $SUPER_RES -y -i {output.img} '
        '{input} {output.map} 2> {log}'
        #' -k min_color=\#450256 mid_color=\#21908D max_color=\#F9E721'


rule hdf_to_bed_bins:
    input:
        '%s/{alg}/{hic_map_dir}/{sample}:{tissue}.hdf5' %HIC_MAPS_DIR,
    output:
        pixels = temp('%s/{alg}/{hic_map_dir}/{sample,[^_:]+}:{tissue,[^_:]+}.bg2' %HIC_MAPS_DIR),
        bins = temp('%s/{alg}/{hic_map_dir}/{sample}_{tissue}.bed' %HIC_MAPS_DIR)
    log:
        'logs/hdf5_to_flatfile_{alg,[^/]+}_{sample,[^_/:]+}:{tissue,[^/:]+}_{hic_map_dir,[^_/]+}.log'
    shell:
        '%s/hdf5_to_flatfile.py -n -f bg2 -b {output.bins} {input} ' %SCRIPT_DIR +
        '> {output.pixels} 2> {log}' 


rule bed_to_cool:
    input:
        pixels = '%s/{alg}/{hic_map_dir}/{sample}:{tissue}.bg2' %HIC_MAPS_DIR,
        bins = '%s/{alg}/{hic_map_dir}/{sample}_{tissue}.bed' %HIC_MAPS_DIR
    output:
        '%s/{alg,[^/]+}/{hic_map_dir,[^/]+}/{sample,[^/:_]+}:{tissue,[^/:_]+}.cool' %HIC_MAPS_DIR,
    log:
        'logs/cooler_load_{alg}_{sample}:{tissue}_{hic_map_dir}.log'
    shell:
        '%s load -f bg2 --field count:dtype=float {input.bins} {input.pixels} ' %COOLER_BIN + 
        '{output} 2> {log}' 


rule hdf_to_bed_bins_full:
    input:
        '%s/{alg}/{hic_map_dir}/{tissue}.hdf5' %HIC_MAPS_DIR,
    output:
        pixels = '%s/{alg,[^/]+}/{hic_map_dir,[^/]+}/{tissue,[^/:_]+}.bg2' %HIC_MAPS_DIR,
        bins = '%s/{alg}/{hic_map_dir}/{tissue}.bed' %HIC_MAPS_DIR
    log:
        'logs/hdf5_to_flatfile_{alg}_{tissue}_{hic_map_dir}.log'
    shell:
        '%s/hdf5_to_flatfile.py -n -f bg2 -b {output.bins} {input} ' %SCRIPT_DIR +
        '> {output.pixels} 2> {log}' 


rule bed_to_cool_full:
    input:
        pixels = '%s/{alg}/{hic_map_dir}/{tissue}.bg2' %HIC_MAPS_DIR,
        bins   = '%s/{alg}/{hic_map_dir}/{tissue}.bed' %HIC_MAPS_DIR
    output:
        '%s/{alg,[^/]+}/{hic_map_dir,[^/]+}/{tissue,[^/:_]+}.cool' %HIC_MAPS_DIR,
    log:
        'logs/cooler_load_{alg}_{tissue}_{hic_map_dir}.log'
    shell:
        '%s load -f bg2 --field count:dtype=float {input.bins} {input.pixels} ' %COOLER_BIN + 
        '{output} 2> {log}' 


rule balance_matrix:
    input:
        '%s/{alg}/{hic_map_dir}/{hic_map}.cool' %HIC_MAPS_DIR,
    output:
        temp('%s/{alg}/{hic_map_dir}/{hic_map}.cool.dummy' %HIC_MAPS_DIR),
    log:
        'logs/cooler_balance_{alg,[^/]+}_{hic_map,[^/]+}_{hic_map_dir,[^_/]+}.log'
    shell:
        '%s balance -f --tol 1e-09 {input} 2> {log};' %COOLER_BIN + 
        'echo 1 > {output}' 


rule plot_matrix:
    input:
        mtrx = '%s/{alg}/{hic_map_dir}/{tissue}.cool' %HIC_MAPS_DIR,
        dummy = '%s/{alg}/{hic_map_dir}/{tissue}.cool.dummy' %HIC_MAPS_DIR,
    output:
        '%s/{alg,[^/]+}/{hic_map_dir,[^/]+}/{tissue,[^/]+}.pdf' %HIC_MAPS_DIR,
    log:
        'logs/plotMatrix_{alg}_{tissue}_{hic_map_dir}.log'
    shell:
        '%s/plotMatrix.py {input.mtrx} > {output} 2> {log}' %SCRIPT_DIR


rule show_chromosome:
    input:
        mtrx = '%s/{alg}/{hic_map_dir}/{tissue}.cool' %HIC_MAPS_DIR,
        dummy = '%s/{alg}/{hic_map_dir}/{tissue}.cool.dummy' %HIC_MAPS_DIR,
    output:
        '%s/{alg,[^/]+}/{hic_map_dir,[^/]+}/{tissue,[^/_]+}_{chr,[^/_]+}.pdf' %HIC_MAPS_DIR,
    shell:
        '%s show -b -s log10 --cmap viridis --dpi 600 ' %COOLER_BIN +
        '-o {output} {input.mtrx} {wildcards.chr}'


rule cool_to_multi:
    input:
        maps = expand('%s/{{alg}}/m{{min_int}}_r{res}_s{{super_res}}/' %HIC_MAPS_DIR +
                '{{sample}}:{{tissue}}.cool', res = RESOLUTIONS),
        dummy = expand('%s/{{alg}}/m{{min_int}}_r{res}_s{{super_res}}/' %HIC_MAPS_DIR +
                '{{sample}}:{{tissue}}.cool.dummy', res = RESOLUTIONS)
    output:
        '%s/{alg,[^/]+}/{sample,[^_/:]+}:{tissue,[^_:]+}' %HIC_MAPS_DIR +
                '_m{min_int,[0-9.]+}_s{super_res,[0-9.]+}.cool'
    log:
        'logs/cool2multi_{alg}_{sample}:{tissue}_m{min_int}_s{super_res}.cool'
    shell:
        '%s/cool2multi.py -o {output} {input.maps} 2> {log}' %SCRIPT_DIR


rule cool_to_multi_full:
    input:
        maps = expand('%s/{{alg}}/m{{min_int}}_r{res}_s{{super_res}}/' %HIC_MAPS_DIR +
                '{{tissue}}.cool', res = RESOLUTIONS),
        dummy = expand('%s/{{alg}}/m{{min_int}}_r{res}_s{{super_res}}/' %HIC_MAPS_DIR +
                '{{tissue}}.cool.dummy', res = RESOLUTIONS)
    output:
        '%s/{alg,[^/]+}/{tissue,[^/:_]+}' %HIC_MAPS_DIR +
                '_m{min_int,[0-9.]+}_s{super_res,[0-9.]+}.cool'
    log:
        'logs/cool2multi_{alg}_{tissue}_m{min_int}_s{super_res}.cool'
    shell:
        '%s/cool2multi.py -o {output} {input.maps} 2> {log}' %SCRIPT_DIR


rule dump_cool_to_bed:
    input:
        cool = '%s/{alg}/{hic_map_dir}/{hic_map}.cool' %HIC_MAPS_DIR,
        dummy = '%s/{alg}/{hic_map_dir}/{hic_map}.cool.dummy' %HIC_MAPS_DIR,
    output:
        temp('%s/{alg}/{hic_map_dir}/{hic_map}.cool.bg2' %HIC_MAPS_DIR),
    log:
        'logs/cooler_dump_{alg,[^/]+}_{hic_map,[^/]+}_{hic_map_dir,[^_/]+}.log'
    shell:
        '%s dump -t pixels --join -b {input.cool} > {output} ' %COOLER_BIN +
        '2> {log}'


rule quantify_difference:
    input:
        map1 = '%s/{alg}/{hic_map_dir}/{sample1}:{tissue}.cool' %HIC_MAPS_DIR,
        map2 = '%s/{alg}/{hic_map_dir}/{sample2}:{tissue}.cool' %HIC_MAPS_DIR,
        dummy1 = '%s/{alg}/{hic_map_dir}/{sample1}:{tissue}.cool.dummy' %HIC_MAPS_DIR,
        dummy2 = '%s/{alg}/{hic_map_dir}/{sample2}:{tissue}.cool.dummy' %HIC_MAPS_DIR,
    output:
        '%s/{alg,[^/]+}/{sample1,[^_/]+}-{sample2,[^_/]+}:{tissue,[^/]+}/{hic_map_dir,[^/]+}.csv' %STATS_DIR,
    shell:
        '%s/quantifyDifference.py {input.map1} {input.map2} > {output}' %SCRIPT_DIR


rule summarize_fullstats:
    input:
        expand('%s/{{alg}}/{{tissue}}/m{min_int}_r{resolution}_s{superres}.csv' %STATS_DIR,
                min_int = MIN_INTS, resolution = RESOLUTIONS, superres =
                SUPER_RES)
    output:
        '%s/{alg,[^/]+}/{tissue,[^/]+}.csv' %STATS_DIR
    shell:
        'IFS=$\'\\t\';'
        'for m in %s; do' %' '.join(map(str, MIN_INTS)) +
        '   for r in %s; do ' %' '.join(map(str, RESOLUTIONS)) +
        '       for s in %s; do ' %' '.join(map(str, SUPER_RES)) +
        '           DATA=$(cat %s/{wildcards.alg}/{wildcards.tissue}/m${{m}}_r${{r}}_s${{s}}.csv);' %STATS_DIR +
        '           echo -e "$m\\t$r\\t$s\\t$DATA";'
        '       done;'
        '   done;'
        'done > {output} || true'

rule visualize_stats:
    input:
        '%s/{stat_file}.csv' %STATS_DIR
    output:
        '%s/{stat_file}.pdf' %STATS_DIR
    shell:
        '%s/plotDifferences.py -p {input} > {output}' %SCRIPT_DIR

