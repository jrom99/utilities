# %% IMPORTS AND DEFS
from math import dist

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure


def get_contact_map(protein: "Chain", max_dist: float = 6):
    def get_distance_map(protein: "Chain"):
        matrix = np.zeros((len(protein), len(protein)))
        for x, res1 in enumerate(protein):
            for y, res2 in enumerate(protein):
                matrix[x, y] = dist(res1["CA"].coord, res2["CA"].coord)

        return matrix

    contact_map = get_distance_map(protein) < max_dist

    diag = np.arange(len(contact_map))
    np.fill_diagonal(contact_map, False)  # residue with itself
    contact_map[(diag[:-1], diag[1:])] = False  # residue with next (upper triangle)
    contact_map[(diag[1:], diag[:-1])] = False  # residue with next (lower triangle)

    return contact_map


def get_distance_map(protein: "Chain"):
    matrix = np.zeros((len(protein), len(protein)))
    for x, res1 in enumerate(protein):
        for y, res2 in enumerate(protein):
            matrix[x, y] = dist(res1["CA"].coord, res2["CA"].coord)

    return matrix


def read_protein(filename: str) -> Model:
    struc: "Structure" = PDBParser(QUIET=True).get_structure("", filename)  # type: ignore
    assert struc is not None
    return struc[0]


def get_sequence(protein: Chain):
    return [protein_letters_3to1[res.resname] for res in protein]


# %% LOAD DATA
protein = read_protein("native.pdb")["A"]

sequence = get_sequence(protein)


# %% FULL SQUARE DISTANCE MAP
sns.set_context("talk")

plt.figure(figsize=(10, 8))
ax = sns.heatmap(get_distance_map(protein), cmap="viridis", vmax=50, cbar_kws={"label": "Distance (angstroms)"})
ax.set(title="Cɑ-Cɑ distance map")
ax.xaxis.tick_top()

rr = range(0, len(sequence), 50)  # ticks every 50 steps
ax.set_xticks(rr, map(str, rr))
ax.set_yticks(rr, map(str, rr))

plt.savefig("native_protein_full_square.png")

# %% FULL SQUARE CONTACT MAP
sns.set_context("talk")

plt.figure(figsize=(8, 8))
ax = sns.heatmap(get_contact_map(protein, 4), cmap=["white", "black"], square=True, cbar=False)
ax.set(title="Contact map")
ax.xaxis.tick_top()

rr = range(0, len(sequence), 50)  # ticks every 50 steps
ax.set_xticks(rr, map(str, rr))
ax.set_yticks(rr, map(str, rr))

sns.despine(top=False, right=False)

plt.savefig("native_protein_full_square.png")
