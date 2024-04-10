rule coarse_model_fit:
    input:
        targets["coarse_hex"]
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
        os.path.join(dirs["envs"],"coarse.yaml")
    script:
        os.path.join(dirs["scripts"],"coarse_model_fit.py")