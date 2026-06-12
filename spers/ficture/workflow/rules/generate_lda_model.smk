# Joint model: trained once on all samples together (shared factors)
rule generate_lda_model:
    input:
        tsv = targets["transcripts"],
    output:
        fit = model_files["model_fit"],
        res = model_files["model_res"],
        coh = model_files["model_coh"],
        pos = model_files["model_pos"],
        mtx = model_files["model_mtx"],
        mdl = model_files["model_mdl"]
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


# One coarse model plot per sample (samples share coordinates, so plot separately)
rule model_hex_bin_plot:
    input:
        fit = model_files["model_fit"],
    output:
        png = patterns["model_png"],
    params:
        hex_width=config["lda_model"]["bin"]["hex_width"],
        plot=config["plot"],
        sample=lambda w: w.sample
    log:
        os.path.join(dirs["logs"],"{sample}_model_hex_bin_plot.txt")
    benchmark:
        os.path.join(dirs["bench"],"{sample}_model_hex_bin_plot.txt")
    threads:
        config["resources"]["ram"]["cpu"]
    resources:
        mem=config["resources"]["ram"]["mem"],
        time=config["resources"]["ram"]["time"]
    conda:
        os.path.join(dirs["envs"],"pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"],"model_hex_bin_plot.py")
