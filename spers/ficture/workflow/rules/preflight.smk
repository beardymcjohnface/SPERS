
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


n_factors = config["lda_model"]["lda"]["n_components"]


# Per-sample path patterns (carry the {sample} wildcard; used by rules)
patterns = {
    "transcripts": os.path.join(dirs["results"], "{sample}", "transcripts.pkl"),
    "hexbins": os.path.join(dirs["results"], "{sample}", "hexbins.pkl"),
    "rescore": os.path.join(dirs["results"], "{sample}", "scored_transcripts", "transcripts.rescored.pkl"),
    "scored_top_png": os.path.join(dirs["results"], "{sample}", "scored_transcripts", "top_factor.png"),
    "scored_factor_png": os.path.join(dirs["results"], "{sample}", "scored_transcripts", "factor_{k}.png"),
    "model_png": os.path.join(dirs["model"], "{sample}", "plot.png"),
}


# Shared (single) joint-model outputs
model_files = {
    "model_fit": os.path.join(dirs["model"], "fit.tsv.gz"),   # combined, has a 'sample' column
    "model_res": os.path.join(dirs["model"], "results.pkl"),
    "model_coh": os.path.join(dirs["model"], "coherence.tsv.gz"),
    "model_pos": os.path.join(dirs["model"], "posterior_counts.tsv.gz"),
    "model_mtx": os.path.join(dirs["model"], "matrix.tsv.gz"),
    "model_mdl": os.path.join(dirs["model"], "model.pkl"),
}


# Targets for `rule all` (expanded over samples / factors)
targets = {
    "transcripts": expand(patterns["transcripts"], sample=SAMPLES),
    **model_files,
    "model_png": expand(patterns["model_png"], sample=SAMPLES),
    "rescore": expand(patterns["rescore"], sample=SAMPLES),
    "scored_top_png": expand(patterns["scored_top_png"], sample=SAMPLES),
    "scored_factor_png": expand(patterns["scored_factor_png"], sample=SAMPLES, k=range(n_factors)),
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
