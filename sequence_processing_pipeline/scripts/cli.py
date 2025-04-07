import click
from sequence_processing_pipeline.Commands import demux_cmd


@click.group()
def cli():
    pass


@cli.command()
@click.option('--id-map', type=click.Path(exists=True), required=True)
@click.option('--infile', type=click.Path(exists=True), required=True)
@click.option('--output', type=click.Path(exists=True), required=True)
@click.option('--task', type=int, required=True)
@click.option('--maxtask', type=int, required=True)
def demux(id_map, infile, output, task, maxtask):
    demux_cmd(id_map, infile, output, task, maxtask)


if __name__ == '__main__':
    cli()
