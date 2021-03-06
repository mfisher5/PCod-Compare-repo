"""
Blast output conversion in GFF requires a BLAST+ tabular format which can be
obtained by using the `--outfmt 6` option with the default columns, as
specified in :func:`mgkit.io.blast.parse_blast_tab`. The script can get data
from the standard in and ouputs GFF lines on the standard output by default.

Uniprot
*******

The Function :func:`mgkit.io.blast.parse_uniprot_blast` is used, which filters
BLAST hits based on bitscore and adds by default a *db* attribute to the
annotation with the value `UNIPROT-SP`, indicating that the SwissProt db is
used and a *dbq* attribute with the value 10. The feature type used in the GFF
is CDS.

.. blockdiag::

    {
        "BLAST+" [color = "#377eb8" , textcolor = 'white', shape = flowchart.input];
        "parse_uniprot_blast" [color = "#e41a1c" , textcolor = 'white', width=200, fontsize=16];
        "GFF" [color = "#4daf4a" , textcolor = 'white'];
        "BLAST+"  -> "parse_uniprot_blast" -> GFF;
    }

BlastDB
*******

If a BlastDB, such as *nt* or *nr* was used, the **blastdb** command offers
some quick defaults to parse BLAST results.

Changes
*******

.. versionadded:: 0.2.3
    added *--fasta-file* option, added more data from a blsat hit

.. versionadded:: 0.2.2
    added *blastdb* command

.. versionchanged:: 0.2.1

    added *-ft* option

.. versionchanged:: 0.1.13

* added *-n* parameter to *uniprot* command
* added *-k* option to *uniprot* command

.. versionadded:: 0.1.12

"""
import sys
import argparse
import logging
from .. import logger
from . import utils
from ..io import blast, fasta


LOG = logging.getLogger(__name__)


def parse_attr_arg(value):
    values = value.split(':')
    if len(values) != 2:
        raise argparse.ArgumentTypeError(
            "Wrong filter format, must be 'key:value' 'key:value'"
        )

    return values[0], values[1]


def set_common_options(parser):
    parser.add_argument(
        '-dbq',
        '--db-quality',
        action='store',
        type=int,
        help='Quality of the DB used',
        default=10
    )
    parser.add_argument(
        '-b',
        '--bitscore',
        action='store',
        type=float,
        help='Minimum bitscore to keep the annotation',
        default=0.0
    )
    parser.add_argument(
        '-k',
        '--attr-value',
        action='append',
        type=parse_attr_arg,
        help='''Additional attribute and value to add to each annotation,
                in the form attr:value''',
        default=None
    )
    parser.add_argument(
        '-ft',
        '--feat-type',
        action='store',
        default='CDS',
        help='Feature type to use in the GFF'
    )
    parser.add_argument(
        '-a',
        '--fasta-file',
        type=argparse.FileType('r'),
        default=None,
        help="""
        Fasta file with nucleotide sequences, used to calculate the frame, if
        not used, the frame on the '-' strand will always be 0
        """
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=argparse.FileType('r'),
        default='-',
        help='BLAST+ output file in tabular format, defaults to stdin'
    )
    parser.add_argument(
        'output_file',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Output GFF file, defaults to stdout'
    )


def set_uniprot_parser(parser):
    parser.add_argument(
        '-db',
        '--db-used',
        action='store',
        type=str,
        default='UNIPROT-SP',
        help='Uniprot database used with BLAST'
    )
    parser.add_argument(
        '-n',
        '--no-split',
        action='store_true',
        default=False,
        help='''if used, the script assumes that the sequence header contains
                only the gene id'''
    )

    parser.set_defaults(func=convert_from_uniprot)


def set_blastdb_parser(parser):
    """
    .. versionadded:: 0.2.2
    """
    parser.add_argument(
        '-db',
        '--db-used',
        action='store',
        type=str,
        default='NCBI-NT',
        help='blastdb used'
    )

    parser.set_defaults(func=convert_from_blastdb)


def load_fasta_file(file_name):
    if file_name is None:
        return None

    return dict(
        (name, len(seq))
        for name, seq in fasta.load_fasta(file_name)
    )


def convert_from_blastdb(options):
    """
    .. versionadded:: 0.2.2
    """

    seqs = load_fasta_file(options.fasta_file)

    iterator = blast.parse_uniprot_blast(
        options.input_file,
        bitscore=options.bitscore,
        db=options.db_used,
        dbq=options.db_quality,
        name_func=None,
        feat_type=options.feat_type,
        seq_lengths=seqs
    )

    for annotation in iterator:
        if options.attr_value is not None:
            for key, value in options.attr_value:
                annotation.set_attr(key, value)

        annotation.to_file(options.output_file)


def convert_from_uniprot(options):

    seqs = load_fasta_file(options.fasta_file)

    if options.no_split is True:
        name_func = lambda x: x
    else:
        name_func = None

    iterator = blast.parse_uniprot_blast(
        options.input_file,
        bitscore=options.bitscore,
        db=options.db_used,
        dbq=options.db_quality,
        name_func=name_func,
        feat_type=options.feat_type,
        seq_lengths=seqs
    )

    for annotation in iterator:
        if options.attr_value is not None:
            for key, value in options.attr_value:
                annotation.set_attr(key, value)

        annotation.to_file(options.output_file)


def set_parser():
    """
    .. versionchanged:: 0.2.2
        added *blastdb* command

    Sets command line arguments parser
    """
    parser = argparse.ArgumentParser(
        description='Convert BLAST output to a GFF file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers()
    parser_u = subparsers.add_parser(
        'uniprot',
        help='Blast results from a Uniprot database, by default SwissProt'
    )

    set_uniprot_parser(parser_u)
    set_common_options(parser_u)

    parser_blastdb = subparsers.add_parser(
        'blastdb',
        help='Blast results from a NCBI database, like *nt*'
    )

    set_blastdb_parser(parser_blastdb)
    set_common_options(parser_blastdb)

    utils.add_basic_options(parser)

    return parser


def main():
    "Main function"

    options = set_parser().parse_args()

    logger.config_log(options.verbose)
    options.func(options)