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


# rule coarse_hex:
#     input:
#         targets["transcripts"]
#     output:
#         targets["coarse_hex"]
#     params:
#         params = config["coarse_hex"]
#     log:
#         os.path.join(dirs["logs"], "coarse_hex.txt")
#     benchmark:
#         os.path.join(dirs["bench"], "coarse_hex.txt")
#     conda:
#         os.path.join(dirs["envs"], "pyscripts.yaml")
#     script:
#         os.path.join(dirs["scripts"], "hex_bin.py")


# rule spatial_minibatch:
#     input:
#         targets["transcripts"]
#     output:
#         targets["batched_matrix"]
#     params:
#         params = config["minibatch"]
#     log:
#         os.path.join(dirs["logs"], "spatial_minibatch.txt")
#     benchmark:
#         os.path.join(dirs["bench"], "spatial_minibatch.txt")
#     conda:
#         os.path.join(dirs["envs"], "pyscripts.yaml")
#     script:
#         os.path.join(dirs["scripts"], "spatial_minibatch.py")