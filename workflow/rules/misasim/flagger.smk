import os
import json


all_sm_data = {}
output_dir = config["output_dir"]
os.makedirs(output_dir, exist_ok=True)

for cfg_sm in config["samples"]:
    sm_id = cfg_sm["sm_id"]
    asm_1 = cfg_sm["asm_1"]
    asm_2 = cfg_sm["asm_2"]
    reads = cfg_sm["reads"]
    alpha = cfg_sm["alpha"]
    threads = cfg_sm.get("threads", 24)
    aligner = cfg_sm.get("aligner", "minimap2")
    aligner_preset = cfg_sm.get("aligner_preset", "lr:hqae")
    kmer_size = cfg_sm.get("kmer_size", 19)

    cfg_file = os.path.join(output_dir, f"{sm_id}.json")
    # https://github.com/mobinasri/flagger/tree/v1.1.0/test_wdls/toil_on_slurm#steps-for-executing-workflows-on-slurm-with-toil
    params = {
        "HMMFlaggerEndToEndWithMapping.sampleName": sm_id,
        "HMMFlaggerEndToEndWithMapping.suffixForFlagger": f"HiFi_{sm_id}_flagger",
        "HMMFlaggerEndToEndWithMapping.suffixForMapping": f"HiFi_{sm_id}_map",
        "HMMFlaggerEndToEndWithMapping.hap1AssemblyFasta": asm_1,
        "HMMFlaggerEndToEndWithMapping.hap2AssemblyFasta": asm_2,
        "HMMFlaggerEndToEndWithMapping.readFiles": reads,
        "HMMFlaggerEndToEndWithMapping.aligner": aligner,
        "HMMFlaggerEndToEndWithMapping.presetForMapping": aligner_preset,
        "HMMFlaggerEndToEndWithMapping.kmerSize": kmer_size,
        "HMMFlaggerEndToEndWithMapping.alphaTsv": alpha,
    }

    with open(cfg_file, "wt") as fh:
        json.dump(params, fh)

    sm_data = {}
    sm_data["config"] = cfg_file
    sm_data["threads"] = threads
    sm_data["inputs"] = [asm_1, asm_2, *reads, alpha]

    all_sm_data[sm_id] = sm_data


wildcard_constraints:
    sm="|".join(all_sm_data.keys()),


rule run_flagger:
    input:
        wdl="flagger/wdls/workflows/hmm_flagger_end_to_end_with_mapping.wdl",
        json=lambda wc: all_sm_data[wc.sm]["config"],
        inputs=lambda wc: all_sm_data[wc.sm]["inputs"],
    output:
        jobs_store=directory(os.path.join(output_dir, "{sm}", "jobStore")),
        json=os.path.join(output_dir, "{sm}_outputs.json"),
        output_dir=directory(os.path.join(output_dir, "{sm}")),
    params:
        working_dir=lambda wc, output: os.path.dirname(output.jobs_store),
        logs_dir=directory(os.path.join(output_dir, "logs", "{sm}")),
    threads: lambda wc: all_sm_data[wc.sm]["threads"]
    log:
        "logs/flagger/run_flagger_{sm}.log",
    resources:
        mem="100GB",
    conda:
        "../../envs/misasim.yaml"
    shell:
        """
        toil-wdl-runner \
        --jobStore {output.jobs_store} \
        --stats \
        --clean=never \
        --workDir {params.working_dir} \
        --batchSystem single_machine \
        --maxCores "{threads}" \
        --batchLogsDir "{params.logs_dir}" \
        "{input.wdl}" \
        "{input.json}" \
        --outputDirectory {output.output_dir} \
        --outputFile "{output.json}" \
        --runLocalJobsOnWorkers \
        --retryCount 1 \
        --container singularity \
        --disableProgress 2>&1 | tee {log}
        """


rule flagger:
    input:
        expand(rules.run_flagger.output, sm=all_sm_data.keys()),
    default_target: True
