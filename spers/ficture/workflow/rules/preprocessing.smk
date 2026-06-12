# Resolve a sample's input path (one converter job per sample)
def sample_input(wildcards):
    return config["args"]["samples"][wildcards.sample]


if config["args"]["platform"] == "visiumhd":

    # Visium HD ships binned outputs rather than a per-transcript CSV.
    # Each --input should point at a Space Ranger "outs/" directory.
    def visiumhd_file(*parts):
        def _resolve(wildcards):
            return os.path.join(
                config["args"]["samples"][wildcards.sample],
                "binned_outputs", config["visiumhd"]["bin_size"], *parts)
        return _resolve

    rule visiumhd_to_tsv:
        input:
            mtx = visiumhd_file("filtered_feature_bc_matrix", "matrix.mtx.gz"),
            barcodes = visiumhd_file("filtered_feature_bc_matrix", "barcodes.tsv.gz"),
            features = visiumhd_file("filtered_feature_bc_matrix", "features.tsv.gz"),
            positions = visiumhd_file("spatial", "tissue_positions.parquet"),
            scalefactors = visiumhd_file("spatial", "scalefactors_json.json"),
        output:
            patterns["transcripts"],
        params:
            params = config["visiumhd"]
        log:
            os.path.join(dirs["logs"], "{sample}_visiumhd_to_tsv.txt")
        benchmark:
            os.path.join(dirs["bench"], "{sample}_visiumhd_to_tsv.txt")
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
            sample_input
        output:
            patterns["transcripts"],
        params:
            params = config[config["args"]["platform"]]
        log:
            os.path.join(dirs["logs"], "{sample}_transcript_csv_to_tsv.txt")
        benchmark:
            os.path.join(dirs["bench"], "{sample}_transcript_csv_to_tsv.txt")
        conda:
            os.path.join(dirs["envs"], "pyscripts.yaml")
        script:
            os.path.join(dirs["scripts"], "transcript_csv_to_tsv.py")
