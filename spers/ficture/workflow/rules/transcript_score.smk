rule transcript_grid_score:
    input:
        tsv = targets["transcripts"],
        mdl = targets["model_mdl"]
    output:
        targets["overlapped_hex_scores"]
    params:
        params = config["transcript_score"]["grid_score"]
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem=config["resources"]["big"]["mem"],
        time=config["resources"]["big"]["time"]
    log:
        os.path.join(dirs["logs"], "transcript_grid_score.txt")
    benchmark:
        os.path.join(dirs["bench"], "transcript_grid_score.txt")
    conda:
        os.path.join(dirs["envs"], "pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"], "transcript_grid_score.py")


rule plot_scored_transcripts:
    input:
        scr = targets["overlapped_hex_scores"],
        trn = targets["transcripts"]
    output:
        png = targets["scored_png"],
    params:
        params=config["transcript_score"],
        plot=config["plot"]
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
        os.path.join(dirs["scripts"], "plot_scored_transcripts.py")