#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Cataloga as sequências do arquivo fasta informado com sua data e local de coleta.

DESCRIÇÃO:
    Esse script suporta metadados oriundos do banco BV-BRC.
    Ele permite a identificação e remoção de sequências duplicadas,
    além de comparar se a data de coleta das mesmas é coerente e extrair a data em comum.
    Opcionalmente, pode-se extrair também o local de coleta, usando diferentes critérios de desempate.

EXEMPLO:
    Adiciona a data aos ids das sequências sem gerar uma tabela `.tsv`:

        $ python add_metadata.py -i sequences.fasta -m metadata.csv --clean-name --deduplicate -f output.fasta
"""

import os
import csv
import datetime
import argparse
import logging
from collections import Counter, defaultdict

try:
    from argcomplete import autocomplete
except ImportError:
    autocomplete = lambda *_: None


def read_metadata(file: str, sep: str = ","):
    DATE_FORMATS: list[tuple[str, str]] = [
        ("%Y-%m-%d", "%Y-%m-%d"),
        ("%Y-%m", "%Y-%m"),
        ("%Y", "%Y"),
        ("%b-%Y", "%Y-%m"),
        ("%d-%b-%Y", "%Y-%m-%d"),
    ]

    def to_iso_date(date: str):
        if not date:
            return ""

        for in_fmt, out_fmt in DATE_FORMATS:
            try:
                dt = datetime.datetime.strptime(date, in_fmt)
                return dt.strftime(out_fmt)
            except ValueError:
                pass

        raise ValueError(f"Unable to parse date: {date!r}")

    with open(file, "r", newline="") as handle:
        metadata_lines = [*csv.DictReader(handle, delimiter=sep)]

    metadata: dict[str, dict[str, str]] = {}
    for line in metadata_lines:
        metadata[line["GenBank Accessions"]] = {
            # YYYY-MM-DD, YYYY-MM, YYYY
            "date": to_iso_date(line["Collection Date"]),
            # replace spaces with underlines
            "local": "_".join(line["Isolation Country"].strip().split()),
        }

    return metadata


def read_fasta(file: str):
    with open(file, "r") as handle:
        seqs: list[str] = []
        seq_ids: list[str] = []

        for line in handle.read().splitlines():
            if line.startswith(">"):
                seqs.append("")
                seq_ids.append(line[1:].strip())
            else:
                seqs[-1] += line

    return [*zip(seq_ids, seqs)]


def deduplicate_fasta(seqs: list[tuple[str, str]]):
    data: defaultdict[str, list[str]] = defaultdict(list)
    [
        data[seq].append(seq_id)
        for multi_id, seq in seqs
        for seq_id in multi_id.split("__")
    ]
    return [("__".join(ids), seq) for seq, ids in data.items()]


def unify_metadata(
    seqs: list[tuple[str, str]],
    metadata: dict[str, dict[str, str]],
    *,
    most_common_local: bool,
):
    new_metadata: dict[str, dict[str, str]] = {}

    for multi_id, _ in seqs:
        seq_ids = [seq_id.removeprefix("accn|").removesuffix("(translated)") for seq_id in multi_id.split("__")]

        dates = Counter(metadata[seq_id]["date"] for seq_id in seq_ids)
        locs = Counter(metadata[seq_id]["local"] for seq_id in seq_ids)

        # discard missing data
        dates.pop("-", None)
        dates.pop("", None)
        locs.pop("-", None)
        locs.pop("", None)

        # picks the common date between sequences
        date = "-".join(os.path.commonprefix([dt.split("-") for dt in dates]))
        local = locs.most_common(1)[0][0] if locs else ""

        if not dates:
            logger.info(f"missing date: {seq_ids}")

        if dates and not date:
            logger.info(f"no common prefix between dates ({len(seq_ids)}): {seq_ids} {dates}")

        if not locs:
            logger.info(f"missing local: {seq_ids}")

        if len(locs) > 1 and not most_common_local:
            logger.info(f"many locals: {seq_ids} {locs}")
            local = ""

        new_metadata[multi_id] = {"date": date, "local": local}

    return new_metadata


def write_fasta(seqs: list[tuple[str, str]], file: str):
    with open(file, "w+") as handle:
        handle.writelines(f">{seq_id}\n{seq}\n" for seq_id, seq in seqs)


def write_tsv(data: list[tuple[str, ...]], file: str):
    with open(file, "w+") as handle:
        lines = ["\t".join(line) for line in data]
        handle.writelines(f"{line}\n" for line in lines)


class NameSpace(argparse.Namespace):
    input_fasta: str
    metadata_csv: str

    output_fasta: str | None
    output_tsv: str | None
    log: str | None

    discard_first: bool
    deduplicate: bool
    use_most_common_local: bool
    save_local: bool
    keep_missing: bool
    clean_name: bool


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-fasta", required=True)
    parser.add_argument("-m", "--metadata-csv", required=True)

    parser.add_argument("-f", "--output-fasta", required=False)
    parser.add_argument("-t", "--output-tsv", required=False)

    parser.add_argument(
        "--discard-first",
        action="store_true",
        help="Descarta a primeira sequência (RefSeq)",
    )
    parser.add_argument(
        "--deduplicate",
        action="store_true",
        help="Concatena os IDs de sequências duplicadas",
    )
    parser.add_argument(
        "--use-most-common-local",
        action="store_true",
        help="Se uma sequência duplicada tiver mais de um local de coleta, usa o mais frequente",
    )
    parser.add_argument(
        "--save-local",
        action="store_true",
        help="Anota o local de coleta das sequências",
    )
    parser.add_argument(
        "--keep-missing",
        action="store_true",
        help="Mantém sequências sem metadados"
    )
    parser.add_argument(
        "--clean-name",
        action="store_true",
        help="Remove `accn|` e `(translated)` dos IDs das sequências"
    )
    parser.add_argument("-l", "--log", help="Arquivo de log (default: STDOUT)")

    autocomplete(parser)
    return parser.parse_args(namespace=NameSpace())


def get_sequences(fasta_file: str, discard_first: bool, deduplicate: bool, clean_name: bool):
    seqs = read_fasta(fasta_file)

    # refseq ignorada de proposito
    if discard_first:
        seqs.pop(0)

    if deduplicate:
        seqs = deduplicate_fasta(seqs)

    if clean_name:
        seqs = [
            ("__".join(seq_id.removeprefix("accn|").removesuffix("(translated)") for seq_id in multi_id.split("__")), seq)
            for multi_id, seq in seqs
            ]
            
    return seqs


def main(args: NameSpace):
    seqs = get_sequences(args.input_fasta, args.discard_first, args.deduplicate, args.clean_name)

    metadata = read_metadata(args.metadata_csv)

    # remove metadado de sequencias duplicadas com vários locais
    metadata = unify_metadata(seqs, metadata, most_common_local=False)

    tsv, fasta = [], []
    for seq_id, seq in seqs:
        meta = metadata[seq_id]

        if not meta["date"] and not args.keep_missing:
            continue

        if args.save_local:
            if not meta["local"] and not args.keep_missing:
                continue
            tsv.append((seq_id, meta["date"], meta["local"]))
            fasta.append((f"{seq_id}__{meta['date']}__{meta['local']}", seq))
        else:
            tsv.append((seq_id, meta["date"]))
            fasta.append((f"{seq_id}__{meta['date']}", seq))

    if args.output_fasta:
        write_fasta(fasta, args.output_fasta)

    if args.output_tsv:
        write_tsv(tsv, args.output_tsv)


if __name__ == "__main__":
    args = parse_args()
    LOGGER_FORMAT = "%(message)s"

    logging.basicConfig(format=LOGGER_FORMAT, level=logging.INFO, filename=args.log)

    logger = logging.getLogger(__name__)

    main(args)
