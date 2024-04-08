
# Directories
dirs = {
    "logs": os.path.join(config["args"]["output"], "ficture", "logs"),
    "bench": os.path.join(config["args"]["output"], "ficture", "bench"),
    "temp": os.path.join(config["args"]["output"], "ficture", "temp"),
    "results": os.path.join(config["args"]["output"], "ficture", "results"),
    "envs": os.path.join(workflow.basedir, "envs"),
    "scripts": os.path.join(workflow.basedir, "scripts")
}


# Targets
targets = {
    # Convert to ficture TSV format
    "transcripts": os.path.join(dirs["results"], "transcripts.tsv.gz"),
    "transcript_minmax": os.path.join(dirs["temp"], "transcript_minmax.tsv"),
    "transcript_counts": os.path.join(dirs["temp"], "transcript_counts.tsv.gz"),

    # Preprocessing
    "batched_matrix": os.path.join(dirs["temp"], "transcripts.batched.matrix.tsv.gz")
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
