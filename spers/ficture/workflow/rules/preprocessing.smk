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
