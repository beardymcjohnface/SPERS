
# Directories
dirs = {
    "logs": os.path.join(config["args"]["output"], "ficture", "logs"),
    "bench": os.path.join(config["args"]["output"], "ficture", "bench"),
    "temp": os.path.join(config["args"]["output"], "ficture", "temp"),
    "results": os.path.join(config["args"]["output"], "ficture", "results"),
    "model": os.path.join(config["args"]["output"], "ficture", "results", "lda_model"),
    "rescore": os.path.join(config["args"]["output"], "ficture", "results", "scored_transcripts"),
    "envs": os.path.join(workflow.basedir, "envs"),
    "scripts": os.path.join(workflow.basedir, "scripts")
}


# Targets
targets = {
    # Convert to ficture TSV format
    "transcripts": os.path.join(dirs["results"], "transcripts.tsv.gz"),

    # lda model files
    "model_fit": os.path.join(dirs["model"], "fit.tsv.gz"),
    "model_res": os.path.join(dirs["model"], "results.pkl"),
    "model_coh": os.path.join(dirs["model"], "coherence.tsv.gz"),
    "model_pos": os.path.join(dirs["model"], "posterior_counts.tsv.gz"),
    "model_mtx": os.path.join(dirs["model"], "matrix.tsv.gz"),
    "model_mdl": os.path.join(dirs["model"], "model.pkl"),
    "model_png": os.path.join(dirs["model"], "plot.png"),

    # Overlap rescore files
    "overlapped_hex_scores": os.path.join(dirs["rescore"], "transcripts.rescored.tsv.gz"),
    "scored_png": os.path.join(dirs["rescore"], "plot.png"),
}


# Misc
target_rules = []


def targetRule(fn):
    """Mark rules as target rules for rule print_targets"""
    assert fn.__name__.startswith("__")
    target_rules.append(fn.__name__[2:])
    return fn


def copy_log_file():
    """Concatenate Snakemake log to output log file"""
    import glob

    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if files:
        current_log = max(files, key=os.path.getmtime)
        shell("cat " + current_log + " >> " + config["args"]["log"])


onsuccess:
    copy_log_file()

onerror:
    copy_log_file()
