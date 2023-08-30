"""
Entrypoint for TCGA_CPTAC_Taxonomic_Profiler

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
    with open(snake_base("tcga_cptac_taxonomic_profiler.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    """Read and print the Citation information from the citation file"""
    with open(snake_base("tcga_cptac_taxonomic_profiler.CITATION"), "r") as f:
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
            default="tcga_cptac_taxonomic_profiler.out",
            show_default=True,
        ),
        click.option(
            "--tmpdir",
            help="Temporary directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="tmp",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=1, show_default=True
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--rerun-incomplete",
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="tcga_cptac_taxonomic_profiler.log",
            callback=default_to_output,
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
    """Snakemake Pipeline to Mine TCGA and CPTAC WGS data for bacterial reads
    \b
    For more options, run:
    tcga_cptac_taxonomic_profiler command --help"""
    pass

help_msg_extra_extract = """
\b
CLUSTER EXECUTION:
tcga_cptac_taxonomic_profiler kraken ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           tcga_cptac_taxonomic_profiler extract --input [file]
Specify threads:    tcga_cptac_taxonomic_profiler extract ... --threads [threads]
Disable conda:      tcga_cptac_taxonomic_profiler extract ... --no-use-conda 
Change defaults:    tcga_cptac_taxonomic_profiler extract ... --snake-default="-k --nolock"
Add Snakemake args: tcga_cptac_taxonomic_profiler extract ... --dry-run --keep-going --touch
Specify targets:    tcga_cptac_taxonomic_profiler extract ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra_extract,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input file/directory", type=str, required=True)
@common_options
def extract(_input, output, log, **kwargs):
    """Run TCGA_CPTAC_Taxonomic_Profiler extract"""
    # Config to add or update in configfile
    merge_config = {
        "input": _input,
        "output": output,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_extract.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )



help_msg_extra_kraken = """
\b
CLUSTER EXECUTION:
tcga_cptac_taxonomic_profiler install_host ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           tcga_cptac_taxonomic_profiler install_host --database [path]
Specify threads:    tcga_cptac_taxonomic_profiler install_host ... --threads [threads]
Disable conda:      tcga_cptac_taxonomic_profiler install_host ... --no-use-conda 
Change defaults:    tcga_cptac_taxonomic_profiler install_host ... --snake-default="-k --nolock"
Add Snakemake args: tcga_cptac_taxonomic_profiler install_host ... --dry-run --keep-going --touch
Specify targets:    tcga_cptac_taxonomic_profiler install_host ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra_kraken,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
            "--database",
            help="Database directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="host_db",
            show_default=True,
        )
@common_options
def install_host(database, log, **kwargs):
    """Install host DB"""
    # Config to add or update in configfile
    merge_config = {
        "database": database,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_instlal_host.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


help_msg_extra_kraken = """
\b
CLUSTER EXECUTION:
tcga_cptac_taxonomic_profiler kraken ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           tcga_cptac_taxonomic_profiler kraken --input [file]
Specify threads:    tcga_cptac_taxonomic_profiler kraken ... --threads [threads]
Disable conda:      tcga_cptac_taxonomic_profiler kraken ... --no-use-conda 
Change defaults:    tcga_cptac_taxonomic_profiler kraken ... --snake-default="-k --nolock"
Add Snakemake args: tcga_cptac_taxonomic_profiler kraken ... --dry-run --keep-going --touch
Specify targets:    tcga_cptac_taxonomic_profiler kraken ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra_kraken,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input directory", type=str, required=True)
@common_options
def kraken(_input, output, log, **kwargs):
    """Run TCGA_CPTAC_Taxonomic_Profiler with kraken & bracken"""
    # Config to add or update in configfile
    merge_config = {
        "input": _input,
        "output": output,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_kraken.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


help_msg_extra_mmseqs = """
\b
CLUSTER EXECUTION:
tcga_cptac_taxonomic_profiler mmseqs ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           tcga_cptac_taxonomic_profiler mmseqs --input [file]
Specify threads:    tcga_cptac_taxonomic_profiler mmseqs ... --threads [threads]
Disable conda:      tcga_cptac_taxonomic_profiler mmseqs ... --no-use-conda 
Change defaults:    tcga_cptac_taxonomic_profiler mmseqs ... --snake-default="-k --nolock"
Add Snakemake args: tcga_cptac_taxonomic_profiler mmseqs ... --dry-run --keep-going --touch
Specify targets:    tcga_cptac_taxonomic_profiler mmseqs ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra_mmseqs,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input directory", type=str, required=True)
@common_options
def mmseqs(_input, output, log, **kwargs):
    """Run TCGA_CPTAC_Taxonomic_Profiler with mmseqs"""
    # Config to add or update in configfile
    merge_config = {
        "input": _input,
        "output": output,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_mmseqs.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


help_msg_extra_assembly = """
\b
CLUSTER EXECUTION:
tcga_cptac_taxonomic_profiler assembly ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           tcga_cptac_taxonomic_profiler assembly --input [file]
Specify threads:    tcga_cptac_taxonomic_profiler assembly ... --threads [threads]
Disable conda:      tcga_cptac_taxonomic_profiler assembly ... --no-use-conda 
Change defaults:    tcga_cptac_taxonomic_profiler assembly ... --snake-default="-k --nolock"
Add Snakemake args: tcga_cptac_taxonomic_profiler assembly ... --dry-run --keep-going --touch
Specify targets:    tcga_cptac_taxonomic_profiler assembly ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra_mmseqs,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input directory", type=str, required=True)
@common_options
def assembly(_input, output, log, **kwargs):
    """Run TCGA_CPTAC_Taxonomic_Profiler to assemble MAGs """
    # Config to add or update in configfile
    merge_config = {
        "input": _input,
        "output": output,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_assembly.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )

help_msg_extra_binning = """
\b
CLUSTER EXECUTION:
tcga_cptac_taxonomic_profiler binning ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           tcga_cptac_taxonomic_profiler binning --input [file]
Specify threads:    tcga_cptac_taxonomic_profiler binning ... --threads [threads]
Disable conda:      tcga_cptac_taxonomic_profiler binning ... --no-use-conda 
Change defaults:    tcga_cptac_taxonomic_profiler binning ... --snake-default="-k --nolock"
Add Snakemake args: tcga_cptac_taxonomic_profiler binning ... --dry-run --keep-going --touch
Specify targets:    tcga_cptac_taxonomic_profiler binning ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra_binning,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input directory", type=str, required=True)
@common_options
def binning(_input, output, log, **kwargs):
    """Run TCGA_CPTAC_Taxonomic_Profiler to bin MAGs """
    # Config to add or update in configfile
    merge_config = {
        "input": _input,
        "output": output,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_binning.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )

help_msg_extra_mge = """
\b
CLUSTER EXECUTION:
tcga_cptac_taxonomic_profiler mge ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           tcga_cptac_taxonomic_profiler mge --input [file]
Specify threads:    tcga_cptac_taxonomic_profiler mge ... --threads [threads]
Disable conda:      tcga_cptac_taxonomic_profiler mge ... --no-use-conda 
Change defaults:    tcga_cptac_taxonomic_profiler mge ... --snake-default="-k --nolock"
Add Snakemake args: tcga_cptac_taxonomic_profiler mge ... --dry-run --keep-going --touch
Specify targets:    tcga_cptac_taxonomic_profiler mge ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra_mge,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input directory", type=str, required=True)
@common_options
def mge(_input, output, log, **kwargs):
    """Run TCGA_CPTAC_Taxonomic_Profiler to find mges """
    # Config to add or update in configfile
    merge_config = {
        "input": _input,
        "output": output,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_mge.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )










help_msg_extra_annotate = """
\b
CLUSTER EXECUTION:
tcga_cptac_taxonomic_profiler mge ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           tcga_cptac_taxonomic_profiler annotate --input [file]
Specify threads:    tcga_cptac_taxonomic_profiler annotate ... --threads [threads]
Disable conda:      tcga_cptac_taxonomic_profiler annotate ... --no-use-conda 
Change defaults:    tcga_cptac_taxonomic_profiler annotate ... --snake-default="-k --nolock"
Add Snakemake args: tcga_cptac_taxonomic_profiler annotate ... --dry-run --keep-going --touch
Specify targets:    tcga_cptac_taxonomic_profiler annotate ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""


@click.command(
    epilog=help_msg_extra_annotate,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input directory", type=str, required=True)
@common_options
def annotate(_input, output, log, **kwargs):
    """Run TCGA_CPTAC_Taxonomic_Profiler to annotate MAGs """
    # Config to add or update in configfile
    merge_config = {
        "input": _input,
        "output": output,
        "log": log
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "run_annotate.smk")),
        system_config=snake_base(os.path.join("config", "config.yaml")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )




@click.command()
@common_options
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=snake_base(os.path.join("config", "config.yaml")))


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(kraken)
cli.add_command(install_host)
cli.add_command(mmseqs)
cli.add_command(assembly)
cli.add_command(binning)
cli.add_command(mge)
cli.add_command(annotate)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
