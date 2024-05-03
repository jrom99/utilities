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

"""Identifica as regiões constantes ao longo do tempo em uma árvore MCC.

DESCRIÇÃO:
    Esse script suporta árvores geradas pelo TreeAnnotator e requere o Biopython.
    Ele identifica o consenso entre as sequências em cada nó da árvore, podendo ser filtrado por altura.

EXEMPLO:
        $ python constant_in_type.py mcc_tree.nex --consensus highest_posterior
"""

import argparse
import textwrap
import re
import sys
import time
import typing
from collections import defaultdict
from datetime import timedelta
from typing import Literal, NamedTuple

try:
    from argcomplete import autocomplete
except ImportError:
    def autocomplete(*args):
        return None

from Bio.Phylo.BaseTree import Clade, Tree


def parse_comment(comment: str) -> dict[str, list[str]]:
    comment = comment.removeprefix("[").removesuffix("]")
    data = [(k[::-1], v[::-1]) for v, k in re.findall(r",?(.+?)=([\d\w\.%]+)", comment[::-1])][::-1]
    data_dct = defaultdict(list)
    for k,v in data:
        data_dct[k].append(v)
    return data_dct


def get_consensus(sequences: dict[str, float] | list[str], weighted: bool = False, min_value: float = 100, no_cons: str = "*") -> str:
    assert all(len(next(iter(sequences))) == len(seq) for seq in sequences)
    assert 0 <= min_value <= 100

    weights: list[dict[str, float]] = [defaultdict(int) for _ in next(iter(sequences))]

    if weighted:
        assert isinstance(sequences, dict)
        for seq, prob in sequences.items():
            for idx, c in enumerate(seq):
                weights[idx][c] += prob
    else:
        for seq in sequences:
            for idx, c in enumerate(seq):
                weights[idx][c] += 1

    consensus = []
    for data in weights:
        total = sum(data.values())
        c = max(data, key=data.__getitem__)

        consensus.append(c if 100*data[c]/total >= min_value else no_cons)

    return "".join(consensus)


def get_data(node: "Clade", height_property: str):
    if isinstance(node, Newick.Clade):
        if node.comment is None:
            raise ValueError("Newick node missing comment")

        data = parse_comment(node.comment)
    else:
        raise TypeError("Unable to parse node")

    key1 = [k for k in data if k.endswith(".set.prob")][0]  # "OROV_L_sequences.set.prob" "OROV_M_seqs.set.prob"
    key2 = [k for k in data if k.endswith(".set")][0]  # "OROV_L_sequences.set" "OROV_M_seqs.set"

    seqs_prob = [*map(float, data[key1][0].strip("'{}").split(","))]
    seqs = data[key2][0].strip("'{}").replace('\\"', "").split(",")
    posteriors = {seq: prob for seq, prob in zip(seqs, seqs_prob)}

    height = float(data[height_property][0])

    return height, posteriors


def get_single_consensus(posteriors: dict[str, float], consensus: str, min_value: float):
    if consensus in ("weighted", "normal"):
        consensus = get_consensus(posteriors, weighted=(consensus == "weighted"), min_value=min_value)
    else:
        consensus = max(posteriors, key=posteriors.__getitem__)
    return consensus


def parse_tree(tree: "Tree", height_property: str, consensus_type: str, min_consensus: float):
    branches: list[tuple[Clade, Clade]] = []
    tips_consensus: list[str] = []
    heights: dict[Clade, float] = {}
    consensus: dict[Clade, str] = {}

    queue = [tree.root]
    while queue:
        node = queue.pop()
        heights[node], node_data = get_data(node, height_property)
        consensus[node] = get_single_consensus(node_data, consensus_type, min_consensus)

        if node.is_terminal():
            tips_consensus.append(consensus[node])

        queue.extend(node.clades)
        branches.extend((node, child) for child in node.clades)

    edges: list[tuple[NodeData, NodeData]] = []
    for n1, n2 in branches:
        if heights[n1] > heights[n2]:
            n1, n2 = n2, n1
        edges.append((
            NodeData(heights[n1], consensus[n1]),
            NodeData(heights[n2], consensus[n2])
        ))

    _, root_data = get_data(tree.root, height_property)
    root_consensus = get_single_consensus(root_data, consensus_type, min_consensus)

    return tips_consensus, edges, root_consensus


def parse_nucleotide_tree(tree: "Tree", height_property: str, consensus_type: str, min_consensus: float):
    from Bio.Seq import translate

    branches: list[tuple[Clade, Clade]] = []
    tips_consensus: list[str] = []
    heights: dict[Clade, float] = {}
    consensus: dict[Clade, str] = {}

    queue = [tree.root]
    while queue:
        node = queue.pop()
        _, node_data = get_data(node, height_property)

        heights[node], node_data = get_data(node, height_property)

        prot_data = defaultdict(float)
        for seq, post in node_data.items():
            prot_data[str(translate(seq))] += post

        consensus[node] = get_single_consensus(prot_data, consensus_type, min_consensus)

        if node.is_terminal():
            tips_consensus.append(consensus[node])

        queue.extend(node.clades)
        branches.extend((node, child) for child in node.clades)

    edges: list[tuple[NodeData, NodeData]] = []
    for n1, n2 in branches:
        if heights[n1] > heights[n2]:
            n1, n2 = n2, n1
        edges.append((
            NodeData(heights[n1], consensus[n1]),
            NodeData(heights[n2], consensus[n2])
        ))

    _, root_data = get_data(tree.root, height_property)

    prot_data: dict[str, float] = defaultdict(float)
    for seq, post in root_data.items():
        prot_data[str(translate(seq))] += post

    root_consensus = get_single_consensus(prot_data, consensus_type, min_consensus)

    return tips_consensus, edges, root_consensus


class NameSpace(argparse.Namespace):
    tree_file: str
    tree_format: Literal["newick", "nexus"]

    start_time: float | None
    end_time: float | None
    at_time: float | None

    height_property: Literal["height", "height_median"]
    min_consensus: float
    consensus: Literal["normal", "weighted", "highest_posterior"]


def _add_mutually_exclusive_groups(parser: argparse.ArgumentParser, args1: list[tuple[str, dict]], args2: tuple[str, dict], **kwargs):
    grp1 = parser.add_mutually_exclusive_group()
    grp2 = parser.add_mutually_exclusive_group()

    _args1 = [arg[0] for arg in args1]
    kwargs1 = [arg[1] for arg in args1]

    arg2, kwarg2 = args2

    grp1.add_argument(_args1[0], **kwargs1[0], **kwargs)
    act = grp1.add_argument(arg2, **kwarg2, **kwargs)

    # See: https://bugs.python.org/issue10984#msg219660
    grp2._group_actions.append(act)
    grp2.add_argument(_args1[1], **kwargs1[1], **kwargs)

    def _replace(func: typing.Callable[[], str]):
        def inner():
            return func().replace(
                f"[{_args1[0]} {kwargs1[0]['metavar']} | [{arg2} {kwarg2['metavar']} | {_args1[1]} {kwargs1[1]['metavar']}]",
                f"[{arg2} {kwarg2['metavar']} | {_args1[0]} {kwargs1[0]['metavar']} {_args1[1]} {kwargs1[1]['metavar']}]"
            )
        return inner
    
    parser.format_usage = _replace(parser.format_usage)
    parser.format_help = _replace(parser.format_help)


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("tree_file", metavar="FILE")

    _add_mutually_exclusive_groups(
        parser,
        [("--start-time", dict(metavar="TIME")), ("--end-time", dict(metavar="TIME"))],
        ("--at-time", dict(metavar="TIME")),
        type=float,
        default=None
    )

    parser.add_argument("-f", "--tree-format", choices=["newick", "nexus"], default="nexus", help="Tree file format (default: %(default)s)")
    parser.add_argument("--height", "--height-property", choices=["height", "height_median"], default="height", dest="height_property")
    parser.add_argument("--min-consensus", type=float, default=100)
    parser.add_argument("--consensus", choices=["normal", "weighted", "highest_posterior"], default="normal",
                        help=textwrap.dedent("""\
                        Node consensus algorithm (default: %(default)s)
                            - `highest_posterior` will look at the sequence with the highest posterior value at each node
                            - `normal` will use a simple majority rule
                            - `weighted` will rank sequences by their posterior values
                        """))
    parser.add_argument("--nucleotide", action="store_true", dest="is_nucleotide", help="If the tree data should be translated before analysis (default: %(default)s)")

    autocomplete(parser)
    return parser.parse_args(namespace=NameSpace())


class Colors:
    """ ANSI color codes """
    BLACK = "\033[0;30m"
    RED = "\033[0;31m"
    GREEN = "\033[0;32m"
    BROWN = "\033[0;33m"
    BLUE = "\033[0;34m"
    PURPLE = "\033[0;35m"
    CYAN = "\033[0;36m"
    LIGHT_GRAY = "\033[0;37m"
    DARK_GRAY = "\033[1;30m"
    LIGHT_RED = "\033[1;31m"
    LIGHT_GREEN = "\033[1;32m"
    YELLOW = "\033[1;33m"
    LIGHT_BLUE = "\033[1;34m"
    LIGHT_PURPLE = "\033[1;35m"
    LIGHT_CYAN = "\033[1;36m"
    LIGHT_WHITE = "\033[1;37m"
    BOLD = "\033[1m"
    FAINT = "\033[2m"
    ITALIC = "\033[3m"
    UNDERLINE = "\033[4m"
    BLINK = "\033[5m"
    NEGATIVE = "\033[7m"
    CROSSED = "\033[9m"
    END = "\033[0m"

    # cancel SGR codes if we don't write to a terminal
    if not sys.stdout.isatty():
        for _ in dir():
            if isinstance(_, str) and _[0] != "_":
                locals()[_] = ""
    else:
        # set Windows console in VT mode
        if __import__("platform").system() == "Windows":
            kernel32 = __import__("ctypes").windll.kernel32
            kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)
            del kernel32


class NodeData(NamedTuple):
    height: float
    consensus: str


if __name__ == "__main__":
    args = parse_args()

    from Bio import Phylo
    from Bio.Phylo import Newick
    from Bio.Phylo.BaseTree import Tree, Clade

    print("Reading tree...")

    # MCC nexus format may require conversion nexus -> figtree -> nexus
    # biopython doesn't understand "translate" block
    t0 = time.perf_counter()
    tree: Tree = Phylo.read(args.tree_file, format=args.tree_format)  # type: ignore

    t1 = time.perf_counter()
    print(f"Tree read in {timedelta(seconds=t1-t0)}")

    if args.is_nucleotide:
        tips_consensus, edges, root_consensus = parse_nucleotide_tree(tree, args.height_property, args.consensus, args.min_consensus)
    else:
        tips_consensus, edges, root_consensus = parse_tree(tree, args.height_property, args.consensus, args.min_consensus)

    START, END = "ROOT", "TIP"

    if args.at_time is not None:
        edges_at = [edge for edge in edges if edge[0].height <= args.at_time <= edge[1].height]
        START, END = "BEFORE", "AFTER"
        start_consensus = get_consensus([st.consensus for st, _ in edges_at])
        end_consensus = get_consensus([ed.consensus for _, ed in edges_at])
        full_consensus = get_consensus([start_consensus, end_consensus])
    else:
        if args.start_time is not None and args.end_time is not None:
            args.start_time, args.end_time = sorted([args.start_time, args.end_time])

        start_consensus = root_consensus
        end_consensus = get_consensus(tips_consensus)

        if args.start_time is not None:
            START = "START"
            edges = [edge for edge in edges if args.start_time <= edge[1].height]
            start_consensus = get_consensus([edge[0].consensus for edge in edges if edge[0].height <= args.start_time])
        if args.end_time is not None:
            END = "END"
            edges = [edge for edge in edges if edge[0].height <= args.end_time]
            end_consensus = get_consensus([edge[1].consensus for edge in edges if args.end_time <= edge[1].height])

        full_consensus = get_consensus([node.consensus for edge in edges for node in edge])

    print("FULL CONSENSUS")
    print(full_consensus.replace("*", f"{Colors.RED}*{Colors.END}"))
    print(f"\n{START} CONSENSUS")
    print(start_consensus.replace("*", f"{Colors.RED}*{Colors.END}"))
    print(f"\n{END} CONSENSUS")
    print(end_consensus.replace("*", f"{Colors.RED}*{Colors.END}"))
