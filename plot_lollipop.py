from math import dist
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure


def read_protein(filename: str) -> Model:
    struc: "Structure" = PDBParser(QUIET=True).get_structure("", filename)  # type: ignore
    assert struc is not None
    return struc[0]

def get_contact_map2(protein: "Chain", max_dist: float = 4):
    """How many atoms from the remaining of the protein are near each residue, ignoring adjacent residues?"""
    matrix = np.zeros(len(protein))
    for x, res1 in enumerate(protein):
        matrix[x] = sum(
            [
                # how many atoms from res2 are near res1?
                sum(any(dist(at1.coord, at2.coord) < max_dist for at1 in res1) for at2 in res2)
                for y, res2 in enumerate(protein)
                if abs(x - y) >= 2
            ]
        )
    return matrix


def get_distance_map(protein: "Chain"):
    matrix = np.zeros((len(protein), len(protein)))
    for x, res1 in enumerate(protein):
        for y, res2 in enumerate(protein):
            matrix[x, y] = dist(res1["CA"].coord, res2["CA"].coord)

    return matrix


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


def get_sequence(protein: Chain):
    return [protein_letters_3to1[res.resname] for res in protein]


def paired_property(data_dict: dict[str, pd.DataFrame], pep: str, properties: list[str]):
    """Convert two or more dataframes in a dict into a new paired dataframe.
    The dataframe must have a RESIDUE1 column

    Columns of returned dataframe:
        name: str, key from data_dict
        position: int
        label: str
    """
    # position aligns the sequences so that zero is equal to the start of the peptide
    # label masks the sequence before and after the peptide as -2, -1, +1, +2
    dfs: list[pd.DataFrame] = []
    for name, data in data_dict.items():
        seq = "".join(data["RESIDUE1"])
        a = seq.index(pep)
        label = [f"{i:+}" for i in range(-a, -a + len(seq) - len(pep) + 1)]
        label[a : a + 1] = pep
        position = [*range(-a, -a + len(seq))]

        dfs.append(
            pd.DataFrame(
                {
                    "name": name,
                    "position":position,
                    "label":label,
                    **{p: data[p].tolist() for p in properties}
                }
            )
        )

    return pd.concat(dfs, ignore_index=True)


RgbColor = tuple[float, float, float]

design = read_protein("design.pdb")["A"]
protein = read_protein("native.pdb")["A"]

sequence = get_sequence(protein)
contact_map = get_contact_map2(protein, 4)

pep = sequence[:]


# %% LOLLIPOP MULTIPLE PROPERTIES

df = paired_property({"native": df_native, "design": df_design}, pep, ["SASA_RELATIVE", "CONTACTS2"])
grouped = df.pivot(index=("position", "label"), columns="name", values=["SASA_RELATIVE", "CONTACTS2"]).reset_index()

# use sequence or R1, R2, R3... to label peptide
grouped.loc[grouped["position"].between(0, len(pep)-1), "label"] = [f"$R_{{{n}}}$" for n in range(1, len(pep)+1)]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5))

for (ax, p, ylabel) in [(ax1, "SASA_RELATIVE", "Relative SASA (%)"), (ax2, "CONTACTS2", r"Number of backbone contacts (≤4Å)")]:
    ax.axvline(-.5, color="black")
    ax.axvline(len(pep) - .5, color="black")

    ax.vlines(x=grouped["position"], ymin=grouped[(p, "native")], ymax=grouped[(p, "design")], color="#9A9996", lw=1.5)
    ax.scatter(x=grouped["position"], y=grouped[(p, "native")], alpha=1, s=30, label="Protein 1", color="#0072B2", zorder=2)
    ax.scatter(x=grouped["position"], y=grouped[(p, "design")], alpha=1, s=30, label="Protein 2", color="#D55E00", zorder=2)

    ax.set(xlabel="Relative position", ylabel=ylabel)
    ax.set_xticks(grouped["position"], grouped["label"])
    ax.legend(title="Protein", loc="upper left", bbox_to_anchor=(1, 1))

    ax.set_xlim(-2.5, len(pep) + 1.5)

fig.tight_layout()
fig.savefig("contact_map/lollipop.svg")
