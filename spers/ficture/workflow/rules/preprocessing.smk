if config["args"]["platform"] == "visiumhd":

    # Visium HD ships binned outputs rather than a per-transcript CSV.
    # --input should point at the Space Ranger "outs/" directory.
    visiumhd_params = config["visiumhd"]
    visiumhd_bin_dir = os.path.join(
        config["args"]["input"], "binned_outputs", visiumhd_params["bin_size"]
    )

    rule visiumhd_to_tsv:
        input:
            mtx = os.path.join(visiumhd_bin_dir, "filtered_feature_bc_matrix", "matrix.mtx.gz"),
            barcodes = os.path.join(visiumhd_bin_dir, "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
            features = os.path.join(visiumhd_bin_dir, "filtered_feature_bc_matrix", "features.tsv.gz"),
            positions = os.path.join(visiumhd_bin_dir, "spatial", "tissue_positions.parquet"),
            scalefactors = os.path.join(visiumhd_bin_dir, "spatial", "scalefactors_json.json"),
        output:
            targets["transcripts"],
        params:
            params = visiumhd_params
        log:
            os.path.join(dirs["logs"], "visiumhd_to_tsv.txt")
        benchmark:
            os.path.join(dirs["bench"], "visiumhd_to_tsv.txt")
        threads:
            config["resources"]["big"]["cpu"]
        resources:
            mem=config["resources"]["big"]["mem"],
            time=config["resources"]["big"]["time"]
        conda:
            os.path.join(dirs["envs"], "pyscripts.yaml")
        script:
            os.path.join(dirs["scripts"], "visiumhd_to_tsv.py")

else:

    rule transcript_csv_to_tsv:
        input:
            config["args"]["input"]
        output:
            targets["transcripts"],
        params:
            params = config[config["args"]["platform"]]
        log:
            os.path.join(dirs["logs"], "transcript_csv_to_tsv.txt")
        benchmark:
            os.path.join(dirs["bench"], "transcript_csv_to_tsv.txt")
        conda:
            os.path.join(dirs["envs"], "pyscripts.yaml")
        script:
            os.path.join(dirs["scripts"], "transcript_csv_to_tsv.py")
