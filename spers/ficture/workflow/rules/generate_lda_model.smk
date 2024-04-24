rule generate_lda_model:
    input:
        tsv = targets["transcripts"],
    output:
        fit = targets["model_fit"],
        res = targets["model_res"],
        coh = targets["model_coh"],
        pos = targets["model_pos"],
        mtx = targets["model_mtx"],
        mdl = targets["model_mdl"]
    params:
        params=config["lda_model"]
    log:
        os.path.join(dirs["logs"],"generate_lda_model.txt")
    benchmark:
        os.path.join(dirs["bench"],"generate_lda_model.txt")
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem=config["resources"]["big"]["mem"],
        time=config["resources"]["big"]["time"]
    conda:
        os.path.join(dirs["envs"],"pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"],"generate_lda_model.py")


rule model_hex_bin_plot:
    input:
        fit = targets["model_fit"],
    output:
        png = targets["model_png"],
    params:
        hex_width=config["lda_model"]["bin"]["hex_width"],
        plot=config["plot"]
    log:
        os.path.join(dirs["logs"],"model_hex_bin_plot.txt")
    benchmark:
        os.path.join(dirs["bench"],"model_hex_bin_plot.txt")
    threads:
        config["resources"]["ram"]["cpu"]
    resources:
        mem=config["resources"]["ram"]["mem"],
        time=config["resources"]["ram"]["time"]
    conda:
        os.path.join(dirs["envs"],"pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"],"model_hex_bin_plot.py")
