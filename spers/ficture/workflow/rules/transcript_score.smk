rule score_overlapped_hex_bins:
    input:
        tsv = targets["transcripts"],
        mdl = targets["coarse_mdl"]
    output:
        targets["overlapped_hex_scores"]
    params:
        params = config["transcript_score"]
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem=config["resources"]["big"]["mem"],
        time=config["resources"]["big"]["time"]
    log:
        os.path.join(dirs["logs"], "score_overlapped_hex_bins.txt")
    benchmark:
        os.path.join(dirs["bench"], "score_overlapped_hex_bins.txt")
    conda:
        os.path.join(dirs["envs"], "pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"], "transcript_rescore.py")


rule plot_scored_transcripts:
    input:
        scr = targets["overlapped_hex_scores"],
        trn = targets["transcripts"]
    output:
        png = targets["scored_png"],
        # svg = targets["scored_svg"]
    params:
        params=config["transcript_plot"]
    log:
        os.path.join(dirs["logs"], "plot_scored_transcripts.txt")
    benchmark:
        os.path.join(dirs["bench"], "plot_scored_transcripts.txt")
    threads:
        config["resources"]["ram"]["cpu"]
    resources:
        mem=config["resources"]["ram"]["mem"],
        time=config["resources"]["ram"]["time"]
    conda:
        os.path.join(dirs["envs"], "pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"], "scored_transcript_plot.py")