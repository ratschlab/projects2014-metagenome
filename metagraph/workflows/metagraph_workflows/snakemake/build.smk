from metagraph_workflows.resource_management import BuildGraphResources, ResourceConfig
from metagraph_workflows import cfg_utils, constants, rule_utils

if build_primary:
    ruleorder: build_joint_primary > build
else:
    ruleorder: build > build_joint_primary


BUILD_RULE="build"
rule build:
    input: seqs_file_list_path
    output: graph_path
    threads: max_threads
    resources:
        mem_mb=BuildGraphResources(BUILD_RULE).get_mem(config),
        disk_mb=BuildGraphResources(BUILD_RULE).get_disk(config),
    params:
        k=config['k'],
        tempdir_opt=cfg_utils.temp_dir_config(config),
        mem_cap=BuildGraphResources(BUILD_RULE).get_mem_cap(config),
        disk_cap=BuildGraphResources(BUILD_RULE).get_disk_cap(config),
    log: cfg_utils.get_log_path(BUILD_RULE, config)
    shell:
        """
        cat {input} | {exec_cmd} build {verbose_opt} \
        --parallel {threads} \
        -k {params.k} \
        -o {output} \
        --mem-cap-gb {params.mem_cap} \
        --disk-cap-gb {params.disk_cap} \
        {params.tempdir_opt} > {log} 2>&1
        """


### Build Primary

canonical_graphs_dir=wdir/f'canonical_graphs'
canonical_graph_path=wdir/f'{graph}_canonical.dbg'

joint_contigs_path=wdir/f'{graph}_primary.fasta.gz'


sample_ids_spec = False
orig_samples_path=wdir/'orig_samples'



STAGE_SAMPLES_RULE="stage_samples"
rule stage_samples:
    output: temp(orig_samples_path/f"{{sample_id}}{config[constants.SAMPLE_STAGING_FILE_ENDING]}")
    params:
        staging_script_path=config[constants.SAMPLE_STAGING_SCRIPT_PATH],
        additional_options=config[constants.SAMPLE_STAGING_SCRIPT_ADDITIONAL_OPTIONS],
    log: cfg_utils.get_log_path(STAGE_SAMPLES_RULE, config, ['sample_id'])
    shell:
        """
        bash {params.staging_script_path} {wildcards.sample_id} {output} {params.additional_options} > {log} 2>&1
        """


BUILD_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE="build_canonical_graph_single_sample"
rule build_canonical_graph_single_sample:
    input: rule_utils.get_build_single_sample_input(config, orig_samples_path, seq_ids_dict)
    output: temp(canonical_graphs_dir/"{sample_id}.dbg")
    threads: max_threads
    resources:
        mem_mb=BuildGraphResources(BUILD_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE).get_mem(config),
        disk_mb=BuildGraphResources(BUILD_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE).get_disk(config),
    params:
        k=config['k'],
        tempdir_opt=cfg_utils.temp_dir_config(config),
        mem_cap=BuildGraphResources(BUILD_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE).get_mem_cap(config),
        disk_cap=BuildGraphResources(BUILD_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE).get_disk_cap(config),
    log: cfg_utils.get_log_path(BUILD_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE, config, ['sample_id'])
    shell:
        """
        echo "{input}" | {exec_cmd} build {verbose_opt} \
        --parallel {threads} \
        --mode canonical \
        -k {params.k} \
        -o {output} \
        --mem-cap-gb {params.mem_cap} \
        --disk-cap-gb {params.disk_cap} \
        {params.tempdir_opt} > {log} 2>&1
        """


PRIMARIZE_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE="primarize_canonical_graph_single_sample"
rule primarize_canonical_graph_single_sample:
    input: canonical_graphs_dir/"{sample_id}.dbg"
    output: temp(contigs_dir/"{sample_id}_primary.fasta.gz")
    threads: max_threads
    resources:
        mem_mb=ResourceConfig(PRIMARIZE_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE).get_mem(config),
    log: cfg_utils.get_log_path(PRIMARIZE_CANONICAL_GRAPH_SINGLE_SAMPLE_RULE, config, ['sample_id'])
    shell:
        """
        echo "{input}" | {exec_cmd} transform {verbose_opt} \
        --to-fasta \
        --primary-kmers \
        --parallel {threads} \
        -o {output} > {log} 2>&1
        """


BUILD_JOINT_GRAPH_RULE="build_joint_graph"
rule build_joint_graph:
    input: rule_utils.get_build_joint_input(config, contigs_dir, seq_ids_dict, seqs_file_list_path)
    output: temp(canonical_graph_path)
    threads: max_threads
    resources:
        mem_mb=BuildGraphResources(BUILD_JOINT_GRAPH_RULE).get_mem(config),
        disk_mb=BuildGraphResources(BUILD_JOINT_GRAPH_RULE).get_disk(config),
    params:
        k=config['k'],
        separate_build=str(bool(config[constants.PRIMARIZE_SAMPLES_SEPARATELY])).lower(),
        tempdir_opt=cfg_utils.temp_dir_config(config),
        mem_cap=BuildGraphResources(BUILD_JOINT_GRAPH_RULE).get_mem_cap(config),
        disk_cap=BuildGraphResources(BUILD_JOINT_GRAPH_RULE).get_disk_cap(config),
    log: cfg_utils.get_log_path(BUILD_JOINT_GRAPH_RULE, config)
    shell:
        """
        if {params.separate_build}; then
            SEQ_PATHS={wdir}/seqs_paths.txt
            echo "{input}" | tr ' ' '\n' > $SEQ_PATHS
        else
            SEQ_PATHS="{input}"
        fi

        cat $SEQ_PATHS | {exec_cmd} build {verbose_opt} \
        --parallel {threads} \
        --mode canonical \
        -k {params.k} \
        -o {output} \
        --mem-cap-gb {params.mem_cap} \
        --disk-cap-gb {params.disk_cap} \
        {params.tempdir_opt} > {log} 2>&1
        
        """

PRIMARIZE_JOINT_GRAPH_RULE="primarize_joint_graph"
rule primarize_joint_graph:
    input: canonical_graph_path
    output: temp(joint_contigs_path)
    threads: max_threads
    resources:
        mem_mb=ResourceConfig(PRIMARIZE_JOINT_GRAPH_RULE).get_mem(config),
    log: cfg_utils.get_log_path(PRIMARIZE_JOINT_GRAPH_RULE, config)
    shell:
        """
        echo "{input}" | {exec_cmd} transform {verbose_opt} \
        --to-fasta \
        --primary-kmers \
        --parallel {threads} \
        -o {output} > {log} 2>&1
        """


BUILD_JOINT_PRIMARY_RULE="build_joint_primary"
rule build_joint_primary:
    input: joint_contigs_path
    output: graph_path
    threads: max_threads
    resources:
        mem_mb=BuildGraphResources(BUILD_JOINT_PRIMARY_RULE).get_mem(config),
        disk_mb=BuildGraphResources(BUILD_JOINT_PRIMARY_RULE).get_disk(config),
    params:
        k=config['k'],
        tempdir_opt=cfg_utils.temp_dir_config(config),
        mem_cap=BuildGraphResources(BUILD_JOINT_PRIMARY_RULE).get_mem_cap(config),
        disk_cap=BuildGraphResources(BUILD_JOINT_PRIMARY_RULE).get_disk_cap(config)
    log: cfg_utils.get_log_path(BUILD_JOINT_PRIMARY_RULE, config)
    shell:
        """
        {exec_cmd} build {verbose_opt} \
        --parallel {threads} \
        --mode primary \
        -k {params.k} \
        -o {output} \
        --mem-cap-gb {params.mem_cap} \
        --disk-cap-gb {params.disk_cap} \
        {input} \
        {params.tempdir_opt} > {log} 2>&1
        """
