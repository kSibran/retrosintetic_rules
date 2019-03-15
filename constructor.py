from CGRtools.files import RDFread, RDFwrite, MRVread, SDFwrite
from collections import defaultdict
import copy
from itertools import product
from CGRtools.containers import ReactionContainer, MoleculeContainer, CGRContainer

def constructor(center_atoms, new_reaction, fg_fg, all_prot, all_coming, n):
    new_cgr = ~new_reaction
    c_a = new_cgr.center_atoms

    prot_gr = new_cgr.substructure(all_prot).split()
    coming_gr = new_cgr.substructure(all_coming).split()

    all_re = MoleculeContainer()
    all_pr =MoleculeContainer()
    for x in prot_gr:
        if center_atoms.intersection(x):
            center_atoms.update(x)
    for x in coming_gr:
        if center_atoms.intersection(x):
            center_atoms.update(x)

    for mol in new_reaction.reactants:
        all_re._node.update(mol._node)
        all_re._adj.update(mol._adj)
    for mol in new_reaction.products:
        all_pr._node.update(mol._node)
        all_pr._adj.update(mol._adj)

    for_dell_re = set(all_re).difference(center_atoms)
    for_dell_pr =set(all_pr).difference(center_atoms)
    [all_re.delete_atom(x) for x in for_dell_re]
    [all_pr.delete_atom(x) for x in  for_dell_pr]

    flag = 0
    for x in all_re.split():
        if set(c_a).intersection(x):
            flag += 1

    fg = ReactionContainer(reactants=(all_re,) , products=(all_pr,))

    if flag == 2:
        fg.meta['composition'] = '2'
    elif flag == 1:
        fg.meta['composition'] = '1'
    elif flag > 2:
        fg.meta['composition'] = 'multi'
    str_graph = str(fg)

    if str_graph not in fg_fg:
        fg_fg[str_graph] = fg
        fg_fg[str_graph].meta['id'] = {n}
    else:
        fg_fg[str_graph].meta['id'].add(n)
    return fg_fg
