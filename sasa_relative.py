#!/usr/bin/env python3
import sys
from pathlib import Path
from pymol import cmd

args = sys.argv[1:]

cmd.load(args[0])

# Execute o comando get_sasa_relative para calcular a SASA relativa da cadeia polimérica

sasa_data = cmd.get_sasa_relative("polymer")
fasta = {
    head.split("_")[-1]: "".join(data)
    for head, *data in map(str.splitlines, cmd.get_fastastr("polymer").split(">")[1:])
}

# Nome do arquivo de saída
output_file = args[1] if len(args) > 1 else f"sasa_relative_{Path(args[0]).stem}.txt"
dsasa_data = {}

for k, v in sasa_data.items():
    *_, ch, i = k
    dsasa_data.setdefault(ch, []).append(
       (int(i), v)
    )

nsasa_data = []
for ch, seq in sorted(fasta.items()):
    data = sorted(dsasa_data[ch])
    if len(seq) == len(data):
        nsasa_data.extend(
        (f"{ch}{i: 3d}: {r}", s) for (i, s), r in zip(data, seq)
        )
    else:
        nsasa_data.extend(
        (f"{ch}{i: 3d}", s) for i, s in data
        )

# Abre o arquivo para escrita
with open(output_file, "w") as f:
    # Escreve os dados no arquivo
    for residue_num, sasa_value in nsasa_data:
        f.write(f"Resíduo {residue_num}: SASA relativa = {sasa_value:.2f} Å^2\n")

print(f"Arquivo de saída '{output_file}' gerado com sucesso.")
