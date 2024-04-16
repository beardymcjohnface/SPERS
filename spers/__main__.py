"""
Entrypoint for SPERS

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from snaketool_utils.cli_utils import OrderedCommands, run_snakemake, copy_config, echo_click


def snake_base(rel_path):
    """Get the filepath to a Snaketool system file (relative to __main__.py)"""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    """Read and print the version from the version file"""
    with open(snake_base("spers.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("spers.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="spers.out",
            show_default=True,
        ),
        click.option(
            "--threads", help="Number of threads to use", default=16, show_default=True
        ),
        click.option(
            "--profile",
            default=None,
            help="Snakemake profile to use",
            show_default=False,
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default="conda",
            help="Custom conda env directory [default: (outputDir)/conda]",
            callback=default_to_output,
            show_default=False,
        ),
        click.option(
            "--workflow-profile",
            default="spers.profile",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/spers.profile/]",
        ),
        click.option(
            "--system-workflow-profile",
            default=snake_base(os.path.join("profile", "config.yaml")),
            help="Default workflow profile",
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """SPERS: SPatialomics Enhanced Research Suite

    \b
    For more options, run:
    spers command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
spers [pipeline] ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           spers [pipeline] --input [file]
Specify threads:    spers [pipeline] ... --threads [threads]
Disable conda:      spers [pipeline] ... --no-use-conda 
Add Snakemake args: spers [pipeline] ... --dry-run --keep-going --touch
Specify targets:    spers [pipeline] ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", help="Input file/directory", type=str, required=True)
@click.option("--platform", help="Spatialomics platform", default="xenium", show_default=True, type=click.Choice(["xenium", "cosmx"]))
@common_options
@click.option("--configfile",
              default=os.path.join("ficture", "config.yaml"),
              show_default=False,
              callback=default_to_output,
              help="Custom config file [default: (outputDir)/ficture/config.yaml]",)
@click.option("--system-config", default=snake_base(os.path.join("ficture", "config", "config.yaml")), hidden=True, )
@click.option("--log",
              default=os.path.join("ficture", "spers.ficture.log"),
              callback=default_to_output,
              hidden=True)
def ficture(**kwargs):
    """Run ficture pipeline"""
    merge_config = {
        "spers": {
            "args": kwargs
        }
    }

    run_snakemake(
        snakefile_path=snake_base(os.path.join("ficture", "workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs
    )


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@common_options
@click.option("--configfile",
              default=os.path.join("ficture", "config.yaml"),
              show_default=False,
              callback=default_to_output,
              help="Custom config file [default: (outputDir)/ficture/config.yaml]",)
@click.option("--system-config", default=snake_base(os.path.join("ficture", "config", "config.yaml")), hidden=True, )
@click.option("--log",
              default=os.path.join("ficture", "spers.ficture.log"),
              callback=default_to_output,
              hidden=True)
def ficture_test(**kwargs):
    """Run ficture test on xenium slice"""
    merge_config = {
        "spers": {
            "args": {
                "input": snake_base(os.path.join("ficture", "test_data", "xenium.smol.csv.gz")),
                "platform": "xenium"
            }
        }
    }

    run_snakemake(
        snakefile_path=snake_base(os.path.join("ficture", "workflow", "Snakefile")),
        merge_config=merge_config,
        **kwargs
    )



@click.command()
def citation(**kwargs):
    """Print the citation(s)"""
    print_citation()


cli.add_command(ficture)
cli.add_command(ficture_test)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
