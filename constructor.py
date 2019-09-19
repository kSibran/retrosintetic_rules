from collections import defaultdict
from CGRtools.containers import ReactionContainer, QueryContainer
# def diff_atoms(mols_1, mols_2):
#     set_1 = set()
#     set_2 = set()
#     for mols in mols_1:
#         for at in mols:
#             set_1.add(at)
#     for mols in mols_2:
#         for at in mols:
#             set_2.add(at)
#
#     return set_1.difference(set_2)

def diff_atoms(mols_1, mols_2):
    set_1 = set()
    set_2 = set()
    [set_1.update(set(mol)) for mol in mols_1]
    [set_2.update(set(mol)) for mol in mols_2]
    return set_1.difference(set_2)

def constructor(center_atoms, new_reaction, fg_fg,n):

    new_cgr = ~new_reaction
    c_a = new_cgr.center_atoms

    all_clevege = diff_atoms(new_reaction.reactants, new_reaction.products)
    all_coming = diff_atoms(new_reaction.products, new_reaction.reactants)

    if all_clevege:
        clevege_gr = new_cgr.substructure(all_clevege).copy().split()
        for x in clevege_gr:
            if center_atoms.intersection(x):
                center_atoms.update(x)
    if all_coming:
        coming_gr = new_cgr.substructure(all_coming).split()
        for x in coming_gr:
            if center_atoms.intersection(x):
                center_atoms.update(x)

    new_re = QueryContainer()
    new_pr = QueryContainer()

    new_center = new_cgr.substructure(center_atoms)
    new_order = {x[0]: n for n, x in enumerate(sorted(new_center.atoms_order.items(), key= lambda x: x[1]), start=1)}
    c_a = [new_order[a] for a in c_a]
    for mol in new_reaction.reactants:
        for x in new_order:
            if x in mol.atoms_numbers:
                new_re.add_atom(mol.atom(x),new_order[x], charge=mol.atom(x).charge, neighbors = mol.atom(x).neighbors, hybridization= mol.atom(x).hybridization)

        for x in mol.bonds():
            if new_order.get(x[0]) and new_order.get(x[1]):
                if new_order[x[0]] in new_re.atoms_numbers and new_order[x[1]] in new_re.atoms_numbers:
                    new_re.add_bond(new_order[x[0]],new_order[x[1]],x[2])

    for mol in new_reaction.products:
        for x in new_order:
            if x in mol.atoms_numbers:
                new_pr.add_atom(mol.atom(x), new_order[x], charge=mol.atom(x).charge,neighbors = mol.atom(x).neighbors, hybridization= mol.atom(x).hybridization)

        for x in mol.bonds():
            if new_order.get(x[0]) and new_order.get(x[1]):
                if new_order[x[0]] in new_pr.atoms_numbers and new_order[x[1]] in new_pr.atoms_numbers:
                    new_pr.add_bond(new_order[x[0]],new_order[x[1]],x[2])

    flag = 0
    for x in new_re.split():
        if set(c_a).intersection(x):
            flag += 1

    # for at_t in c_a:
    #     if at_t in new_re.atoms_numbers:
    #         new_re.atom(at_t).neighbors = None
    #     if at_t in new_pr.atoms_numbers:
    #         new_pr.atom(at_t).neighbors = None
    #
    # for atom_co_pr in all_clevege.union(all_coming).difference(c_a):
    #     if atom_co_pr in new_re.atoms_numbers:
    #         new_re.atom(atom_co_pr).hybridization = None
    #     if atom_co_pr in new_pr.atoms_numbers:
    #         new_pr.atom(atom_co_pr).hybridization = None

    fg = ReactionContainer(reactants=(new_re,) , products=(new_pr,))
    all_clevege = diff_atoms(fg.reactants, fg.products)
    all_coming = diff_atoms(fg.products, fg.reactants)
    center_atoms = [new_order[a] for a in center_atoms]
    for at_env in set(center_atoms).difference(all_clevege.union(all_coming).union(c_a)):

        if at_env in new_re.atoms_numbers:
            fg.reactants[0]._neighbors[at_env] = ()
            #fg.reactants[0]._hybridizations[at_env] = ()
            # new_re.atom(at_env).neighbors = None
        if at_env in new_pr.atoms_numbers:
            fg.products[0]._neighbors[at_env] = ()
           # fg.products[0]._hybridizations[at_env] = ()

        # if at_env in new_re.atoms_numbers:
        #     vzv = new_re.atom(at_env)
        #     new_re._neighbors[at_env] = ()
        #     #new_re.atom(at_env).neighbors = None
        # if at_env in new_pr.atoms_numbers:
        #     new_pr._neighbors[at_env] = ()
            #new_pr.atom(at_env).neighbors = None

    if flag == 2:
        fg.meta['composition'] = '2'

    elif flag == 1:
        fg.meta['composition'] = '1'
    elif flag > 2:
        fg.meta['composition'] = 'multi'
    str_graph = str(fg)

    if not fg.meta.get('composition'):
        return False
    else:
        if str_graph not in fg_fg:
            fg_fg[str_graph] = fg
            fg_fg[str_graph].meta['id'] = {n}
            fg_fg[str_graph].meta['years']= {int(new_reaction.meta['year'])}
            fg_fg[str_graph].meta['d_a'] = defaultdict(set)
            fg_fg[str_graph].meta['d_a'][new_reaction.meta['year']].add(n)
        else:
            fg_fg[str_graph].meta['id'].add(n)
            fg_fg[str_graph].meta['years'].add(int(new_reaction.meta['year']))
            fg_fg[str_graph].meta['d_a'][new_reaction.meta['year']].add(n)
        return fg_fg
