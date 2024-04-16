rule coarse_model_fit:
    input:
        hex = targets["coarse_hex"],
        bch = targets["batched_matrix"]
    output:
        fit = targets["coarse_fit"],
        res = targets["coarse_res"],
        coh = targets["coarse_coh"],
        pos = targets["coarse_pos"],
        mtx = targets["coarse_mtx"],
        mdl = targets["coarse_mdl"]
    params:
        params=config["coarse_model"]
    log:
        os.path.join(dirs["logs"],"coarse_model_fit.txt")
    benchmark:
        os.path.join(dirs["bench"],"coarse_model_fit.txt")
    threads:
        config["resources"]["big"]["cpu"]
    resources:
        mem=config["resources"]["big"]["mem"],
        time=config["resources"]["big"]["time"]
    conda:
        os.path.join(dirs["envs"],"pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"],"coarse_model_fit.py")


rule coarse_hex_plot:
    input:
        fit = targets["coarse_fit"],
        hex = targets["coarse_hex"],
        trn = targets["transcripts"]
    output:
        png = targets["coarse_png"],
        svg = targets["coarse_svg"]
    params:
        params=config["coarse_plot"]
    log:
        os.path.join(dirs["logs"],"coarse_hex_plot.txt")
    benchmark:
        os.path.join(dirs["bench"],"coarse_hex_plot.txt")
    threads:
        config["resources"]["ram"]["cpu"]
    resources:
        mem=config["resources"]["ram"]["mem"],
        time=config["resources"]["ram"]["time"]
    conda:
        os.path.join(dirs["envs"],"pyscripts.yaml")
    script:
        os.path.join(dirs["scripts"],"coarse_hex_plot.py")
