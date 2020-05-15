import click
import os
from plex_stats.primer_stats import get_primer_stats




@click.command()
@click.argument('primer_fasta', type=click.Path(exists=True, dir_okay=False))
@click.option('--outpath', '-o', 
    type=click.Path(dir_okay=True, file_okay=False), default=os.getcwd(),
    help="Output directory.")
@click.option('--prefix', '-p', type=str, default="primer_stats",
    help="Output file name prefix")
@click.option('--Na', '-Na', default=50, type=click.IntRange(min=0),
    help="Millimolar concentration of Na")
@click.option('-K', '--K', default=0, type=click.IntRange(min=0),
    help="Millimolar concentration of K")
@click.option('-Tris', '--Tris', default=0, type=click.IntRange(min=0),
    help="Millimolar concentration of Tris")
@click.option('-Mg', '--Mg', default=0, type=click.IntRange(min=0),
    help="Millimolar concentration of Mg")
@click.option('-dNTPs', '--dNTPs', default=0, type=click.IntRange(min=0),
    help="Millimolar concentration of dNTPs")
def cli(primer_fasta, outpath, prefix, na, k, tris, mg, dntps):
    tm_params = {'Na': na, 'K': k, 'Tris': tris, 'Mg': mg, 'dNTPs': dntps}
    get_primer_stats(primer_fasta, outpath, prefix, tm_params)
    
if __name__ == '__main__':
    cli()

