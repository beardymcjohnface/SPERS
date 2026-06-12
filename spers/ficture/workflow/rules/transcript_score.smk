# Score each sample with the shared joint model (one job per sample)
rule transcript_grid_score:
    input:
        tsv = patterns["transcripts"],
        mdl = model_files["model_mdl"]
    output:
        patterns["rescore"]
    params:
        params = config["transcript_score"]["grid_score"],
        platform = config["args"]["platform"]
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem=config["resources"]["big"]["mem"],
        time=config["resources"]["big"]["time"]
    log:
        os.path.join(dirs["logs"], "{sample}_transcript_grid_score.txt")
    benchmark:
        os.path.join(dirs["bench"], "{sample}_transcript_grid_score.txt")
    conda:
        os.path.join(dirs["envs"], "pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"], "transcript_grid_score.py")


rule plot_scored_transcripts:
    input:
        scr = patterns["rescore"],
        trn = patterns["transcripts"]
    output:
        top = patterns["scored_top_png"],
        factors = expand(patterns["scored_factor_png"], k=range(n_factors), allow_missing=True)
    params:
        params=config["transcript_score"],
        plot=config["plot"],
        platform=config["args"]["platform"]
    log:
        os.path.join(dirs["logs"], "{sample}_plot_scored_transcripts.txt")
    benchmark:
        os.path.join(dirs["bench"], "{sample}_plot_scored_transcripts.txt")
    threads:
        config["resources"]["ram"]["cpu"]
    resources:
        mem=config["resources"]["ram"]["mem"],
        time=config["resources"]["ram"]["time"]
    conda:
        os.path.join(dirs["envs"], "pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"], "plot_scored_transcripts.py")
