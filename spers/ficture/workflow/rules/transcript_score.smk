rule score_overlapped_hex_bins:
    input:
        tsv = targets["transcripts"],
        mdl = targets["coarse_mdl"]
    output:
        targets["overlapped_hex_scores"]
    params:
        params = config["transcript_score"]
    log:
        os.path.join(dirs["logs"], "score_overlapped_hex_bins.txt")
    benchmark:
        os.path.join(dirs["bench"], "score_overlapped_hex_bins.txt")
    conda:
        os.path.join(dirs["envs"], "pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"], "transcript_rescore.py")

