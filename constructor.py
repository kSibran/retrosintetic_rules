from collections import defaultdict
from CGRtools.containers import ReactionContainer, QueryContainer

def constructor(center_atoms, new_reaction, fg_fg, all_clevege, all_coming, n):
    new_reaction.reset_query_marks()

    new_cgr = ~new_reaction
    c_a = new_cgr.center_atoms

    clevege_gr = new_cgr.substructure(all_clevege).split()
    coming_gr = new_cgr.substructure(all_coming).split()

    new_re = QueryContainer()
    new_pr = QueryContainer()

    for x in clevege_gr:
        if center_atoms.intersection(x):
            center_atoms.update(x)
    for x in coming_gr:
        if center_atoms.intersection(x):
            center_atoms.update(x)


    # for mol in new_reaction.reactants:
    #     for x in mol.atoms():
    #         pre_re.add_atom(x[1], x[0])
    #     for x in mol.bonds():
    #         pre_re.add_bond(*x)
        #all_re._node.update(mol._node) ну очень грязный хак
        #all_re._adj.update(mol._adj)

    # for mol in new_reaction.products:
    #     for x in mol.atoms():
    #         pre_pr.add_atom(x[1], x[0])
    #     for x in mol.bonds():
    #         pre_pr.add_bond(*x)
        # all_pr._node.update(mol._node) ну очень грязный хак
        # all_pr._adj.update(mol._adj)

    # re_order = pre_re.atoms_order
    # pr_order = pre_pr.atoms_order

    new_center = new_cgr.substructure(center_atoms, as_view= False)

    new_order = {x[0]: n for n, x in enumerate(sorted(new_center.atoms_order.items(), key= lambda x: x[1]), start=1)}

    for mol in new_reaction.reactants:
        for x in new_order:
            if x in mol.atoms_numbers:
                new_re.add_atom(mol.atom(x),new_order[x])

        for x in mol.bonds():
            if new_order.get(x[0]) and new_order.get(x[1]):
                if new_order[x[0]] in new_re.atoms_numbers and new_order[x[1]] in new_re.atoms_numbers:
                    new_re.add_bond(new_order[x[0]],new_order[x[1]],x[2])

    for mol in new_reaction.products:
        for x in new_order:
            if x in mol.atoms_numbers:
                new_pr.add_atom(mol.atom(x), new_order[x])

        for x in mol.bonds():
            if new_order.get(x[0]) and new_order.get(x[1]):
                if new_order[x[0]] in new_pr.atoms_numbers and new_order[x[1]] in new_pr.atoms_numbers:
                    new_pr.add_bond(new_order[x[0]],new_order[x[1]],x[2])

    # for_dell_re = set(all_re).difference(center_atoms)
    # for_dell_pr =set(all_pr).difference(center_atoms)
    # [all_re.delete_atom(x) for x in for_dell_re]
    # [all_pr.delete_atom(x) for x in for_dell_pr]

    flag = 0
    for x in new_re.split():
        if set(c_a).intersection(x):
            flag += 1
    for at_env in center_atoms.difference(all_clevege.union(all_coming).union(c_a)):
        if at_env in new_re.atoms_numbers:
            new_re.atom(at_env).neighbors = None
        if at_env in new_pr.atoms_numbers:
            new_pr.atom(at_env).neighbors = None
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
        fg_fg[str_graph].meta['years']=[int(new_reaction.meta['year'])]
        fg_fg[str_graph].meta['d_a'] = defaultdict(list)
        fg_fg[str_graph].meta['d_a'][new_reaction.meta['year']].append(n)
    else:
        fg_fg[str_graph].meta['id'].add(n)
        fg_fg[str_graph].meta['years'].append(int(new_reaction.meta['year']))
        fg_fg[str_graph].meta['d_a'][new_reaction.meta['year']].append(n)
    return fg_fg
