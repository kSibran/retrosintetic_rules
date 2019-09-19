from collections import defaultdict
from itertools import product
from functools import reduce
from operator import   or_
from CGRtools.containers import ReactionContainer, MoleculeContainer

def diff_atoms(mols_1, mols_2):
    set_1 = set()
    set_2 = set()
    [set_1.update(set(mol)) for mol in mols_1]
    [set_2.update(set(mol)) for mol in mols_2]
    return set_1.difference(set_2)

def prot_come(rc, prot_or_come, gr):
    out_set = set()
    if set(rc).intersection(prot_or_come):
        for g in gr:
            if set(rc).intersection(g):
                out_set.update(set(g).difference(rc))
    return out_set

def atom_re_pr(set_atom, mol_source_re, mol_source_pr):
    out = defaultdict(dict)
    for atom in set_atom:
        re = None
        pr = None

        if atom in mol_source_re.atoms_numbers:
            re = mol_source_re.atom(atom).copy() #.copy  нужна для того чтобы при изменении меток гибридизации и меток соседних атомов, не изменялись те же самые метки в других реакциях( но возможно это бесполезно так как метки удаляються совсем в другом месте)

        if atom in mol_source_pr.atoms_numbers:
            pr = mol_source_pr.atom(atom).copy()

        out[atom] ={'re':re, 'pr' : pr}
    return out

def bonds_re_pr_t(bonds, atoms, all_prot, all_coming):
    out=defaultdict(dict)
    for bon in bonds:
        if bon[0] in atoms and bon[1] in atoms:
            if bon[0] in all_prot and bon[1] in all_prot:
                out[str(bon[0]) + ',' + str(bon[1])] = {'re': (bon[0], bon[1], bon[2].order),'pr': (bon[0], bon[1], None)}
            elif bon[0] in all_coming and bon[1] in all_coming:
                out[str(bon[0]) + ',' + str(bon[1])] = {'re': (bon[0], bon[1], None),'pr': (bon[0], bon[1], bon[2].p_order)}
            else:
                out[str(bon[0]) + ',' + str(bon[1])] = {'re':(bon[0], bon[1], bon[2].order), 'pr' :(bon[0], bon[1],bon[2].p_order)}
            #out[str(bon[0]) + ',' + str(bon[1])] = {'re':(bon[0], bon[1], bon[2].order), 'pr' :(bon[0], bon[1],bon[2].p_order)}
    return out

def bonds_re_pr(bonds, atoms):
    out=defaultdict(dict)
    for bon in bonds:
        if bon[0] in atoms or bon[1] in atoms:
            out[str(bon[0]) + ',' + str(bon[1])] = {'re':(bon[0], bon[1], bon[2].order), 'pr' :(bon[0], bon[1],bon[2].p_order)}
           # out[str(bon[0])+','+str(bon[1])]={'re':(bon[0], bon[1],bon[2]._reactant),'pr' :(bon[0],bon[1],bon[2]._product)}
    return out

def add_at(source, mol_source, out_source):
    for x in source:
        out_source.add_atom(mol_source.atom(x).copy(), x, charge = mol_source.atom(x).charge)
    return out_source

def enumeration_cgr(reaction):

    cgrs = ~reaction

    all_prot = diff_atoms(reaction.reactants, reaction.products)
    all_coming = diff_atoms(reaction.products, reaction.reactants)

    prot_list = []
    coming_list = []

    if all_prot:
        prot_gr = cgrs.substructure(all_prot).split()
        prot_list = [set(x) for x in prot_gr]
    if all_coming:
        coming_gr = cgrs.substructure(all_coming).split()
        coming_list = [set(x) for x in coming_gr]

    all_prot_come = prot_list + coming_list
    other_list = list(cgrs.centers_list)

    for x in all_prot_come:
        unio = []

        for i, y in enumerate(other_list):
            if set(x).intersection(y):
                unio.append(i)

        if len(unio)>1:
            qq = set()
            for i in reversed(unio):
                qq.update(other_list[i])
                other_list.pop(i)
            other_list.append(list(qq))
        del(unio)

    del (prot_list, coming_list)

    cycles = []
    for x in reaction.reactants:
        cycles.extend(x.sssr)
    for x in reaction.products:
        for y in x.sssr:
            if y not in cycles:
                cycles.append(y)
# объединение реакционных центров при циклизации в общий для этих рц цикл
    for y in cycles:
        kk=cgrs.substructure(y)
        ept = []
        new_ind = []
        unite = []
        if all(z[2].order == 4 for z in list(kk.bonds())) and any(z[2].p_order != 4 for z in list(kk.bonds()))\
                or all(z[2].p_order == 4 for z in list(kk.bonds())) and any(z[2].order != 4 for z in list(kk.bonds())):
            unite.extend(y)
        else:
            for x in other_list:
                #if set(x).intersection(kk)  and len(set(x).intersection(kk))>1:
                if len(set(x).intersection(kk))>1:
                    ept.append(set(x).intersection(kk))

            if len(ept) >= 2:
                for x in ept:
                    for i, p in enumerate(x):
                        for i2, m in enumerate(x):
                            if i!=i2:
                                if any((m == mp[0] and p == mp[1]) or (p==mp[0] and m == mp[1]) for mp in kk.bonds()) and (kk.bond(m,p).order == None or kk.bond(m,p).p_order == None):
                                    unite.extend([m,p])

        if unite:
            for i3, zop in enumerate(other_list):
                if set(zop).intersection(unite):
                    new_ind.append(i3)
        if len(new_ind)>1:
            y = []
            new_ind.reverse()
            for x in new_ind:
                y.extend(other_list[x])
                other_list.pop(x)
            other_list.append(y)
        del(ept, new_ind, unite,kk)
#конец объединения
    del (cycles, all_prot_come) #  хреновая оптимизация

    if 1 < len(other_list):

        #variants_reaction = []
        reactants_to_work = reduce(or_, reaction.reactants)
        product_to_work = reduce(or_, reaction.products)
        list_rc = []

        const_at = atom_re_pr(set(reactants_to_work).difference(all_prot).difference(cgrs.center_atoms),reactants_to_work, product_to_work)
        const_bond = bonds_re_pr(cgrs.bonds(), const_at.keys())

        prot_atoms = {}
        prot_bonds = {}
        comming_atoms = {}
        comming_bonds = {}
        for big_rc in other_list:

            if all_prot:
                prot_atoms = atom_re_pr(prot_come(big_rc, all_prot, prot_gr), reactants_to_work, product_to_work)
                prot_bonds = bonds_re_pr(cgrs.bonds(), prot_atoms.keys())

            if all_coming:
                comming_atoms = atom_re_pr(prot_come(big_rc, all_coming, coming_gr), reactants_to_work, product_to_work)
                comming_bonds = bonds_re_pr(cgrs.bonds(), comming_atoms.keys())

            atom_t = atom_re_pr(set(big_rc), reactants_to_work, product_to_work)
            t_bond = bonds_re_pr_t(cgrs.bonds(), atom_t.keys(), all_prot, all_coming)

            list_rc.append([atom_t, prot_atoms, comming_atoms, t_bond, prot_bonds, comming_bonds])


        for e, t_rc in enumerate(list_rc):
            for state in list(product([1, 0], repeat=len(list_rc) - 1)):

                new_reactant = MoleculeContainer()
                new_product = MoleculeContainer()
                new_all_bonds_reactants = []
                new_all_bonds_products =[]

                add_at(const_at, reactants_to_work, new_reactant)
                add_at(const_at, product_to_work, new_product)

                for x in const_bond.values():
                    new_all_bonds_reactants.append(x['re'])
                    new_all_bonds_products.append(x['pr'])

                for x, y in t_rc[0].items():
                    if y['re']:
                        new_reactant.add_atom(y['re'].copy(), x, charge = reactants_to_work.atom(x).charge)
                    if y['pr']:
                        new_product.add_atom(y['pr'].copy(), x, charge = product_to_work.atom(x).charge)

                for x, y in t_rc[1].items():
                    new_reactant.add_atom(y['re'].copy(), x, charge = reactants_to_work.atom(x).charge )  # атомы уходящей группы

                for x, y in t_rc[2].items():
                    new_product.add_atom(y['pr'].copy(), x, charge = product_to_work.atom(x).charge) # приходящие атомы

                for x in t_rc[5].values():
                    new_all_bonds_products.append(x['pr'])

                for x in t_rc[3].values():
                    if x['re'][2]:
                        new_all_bonds_reactants.append(x['re'])
                    if x['pr'][2]:
                        new_all_bonds_products.append(x['pr'])

                for x in t_rc[4].values():
                    new_all_bonds_reactants.append(x['re'])

                for s, big_rc in zip(state, list_rc[:e] + list_rc[e+1:]):
                    if s == 0:

                        for x, y in big_rc[0].items():
                            if y['re']:
                                new_reactant.add_atom(y['re'].copy(), x, charge = reactants_to_work.atom(x).charge)
                                new_product.add_atom(y['re'].copy(), x, charge = reactants_to_work.atom(x).charge)

                        for x, y in big_rc[1].items():
                            new_reactant.add_atom(y['re'].copy(), x, charge = reactants_to_work.atom(x).charge)
                            new_product.add_atom(y['re'].copy(), x, charge = reactants_to_work.atom(x).charge)

                        for x in big_rc[3].values():
                            if x['re'][2]:
                                new_all_bonds_reactants.append(x['re'])
                                new_all_bonds_products.append(x['re'])

                        for x in big_rc[4].values():
                            new_all_bonds_reactants.append(x['re'])
                            new_all_bonds_products.append(x['re'])

                    if s == 1:

                        for x, y in big_rc[0].items():
                            if y['pr']:
                                new_reactant.add_atom(y['pr'].copy(), x, charge = product_to_work.atom(x).charge)
                                new_product.add_atom(y['pr'].copy(), x, charge = product_to_work.atom(x).charge)

                        for x, y in big_rc[2].items():
                            new_reactant.add_atom(y['pr'].copy(), x, charge = product_to_work.atom(x).charge)
                            new_product.add_atom(y['pr'].copy(), x, charge = product_to_work.atom(x).charge)

                        for x in big_rc[3].values():
                            if x['pr'][2]:
                                new_all_bonds_reactants.append(x['pr'])
                                new_all_bonds_products.append(x['pr'])

                        for x in big_rc[5].values():
                            new_all_bonds_reactants.append(x['pr'])
                            new_all_bonds_products.append(x['pr'])

                [new_reactant.add_bond(*x) for x in new_all_bonds_reactants]
                [new_product.add_bond(*x) for x in new_all_bonds_products ]
             #   [new_product.add_bond(*x) for x in new_all_bonds_products if not new_product.has_edge(x[0], x[1])]

                new_reaction = ReactionContainer(reactants = (new_reactant,), products =(new_product,))
                new_reaction.meta.update(reaction.meta)
               # variants_reaction.append(new_reaction)
                yield new_reaction
    else:
        yield reaction