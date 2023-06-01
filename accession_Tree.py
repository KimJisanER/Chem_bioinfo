from Bio import Phylo
from tkinter import filedialog
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import matplotlib.pyplot as plt
import pandas as pd


fasta_file = filedialog.askopenfilename(title='Select the fastafile')

aln = AlignIO.read(open(f'{fasta_file}'), 'fasta')
constructor = DistanceTreeConstructor()
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)
njtree = constructor.nj(dm)
print(njtree)

fig, ax = plt.subplots(figsize=(24, 15))

print(dm)
dm_df = pd.DataFrame(list(dm))

dm_df.to_csv(fasta_file + '_nj_identitiy_distance_map' + '.csv')
Phylo.draw(njtree, do_show=False, axes=ax, branch_labels=lambda c: round(c.branch_length, 4))
plt.savefig(fasta_file + '_with_score_nj_identity62' + '.png', dpi=300, bbox_inches='tight')

