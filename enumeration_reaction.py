from CGRtools.files import RDFread, RDFwrite, MRVread
from CGRtools.algorithms import sssr
from collections import defaultdict
import copy
from itertools import product
from CGRtools.containers import ReactionContainer, MoleculeContainer, CGRContainer
from CGRtools.preparer import CGRpreparer


def big_mol(mols):
    to_work= MoleculeContainer()
    for x in mols:
        to_work._node.update(x._node)
        to_work._adj.update(x._adj)
    return to_work

def prot_come(rc, prot_or_come, gr):
    out_set = set()
    if set(rc).intersection(prot_or_come):
        for g in gr:
            if set(rc).intersection(g):
                out_set.update(g)
    return out_set

def atom_re_pr(set_atom, mol_source_re, mol_source_pr):
    out = defaultdict(dict)
    for x in set_atom:
        if x in mol_source_re.atoms_numbers:
            re = mol_source_re.atom(x)
        else:
            re = None
        if x in mol_source_pr.atoms_numbers:
            pr = mol_source_pr.atom(x)
        else:
            pr = None
        out[x] ={'re':re, 'pr' : pr}
    return out

def bonds_re_pr(bonds, atoms):
    out=defaultdict(dict)
    for bon in bonds:
        if set([bon[0], bon[1]]).intersection(atoms):
            out[str(bon[0])+','+str(bon[1])]={'re':(bon[0], bon[1], bon[2]._reactant),'pr' :(bon[0],bon[1],bon[2]._product)}
    return out

def add_at(source, mol_source, out_source):
    for x in source:
        if x not in out_source:
            out_source.add_atom(mol_source.atom(x), x)
    return out_source

def enumeration_cgr(reaction, all_prot, all_coming ):

    cgrs = ~reaction

    pp = cgrs.center_atoms

    # all_prot = diff_atoms(reaction.reactants, reaction.products)
    # all_coming = diff_atoms(reaction.products, reaction.reactants)
    # all_cycl = []
    # for x in reaction.reagents:
    #     all_cycl.extend(x.sssr)
    # for x in reaction.products:
    #     all_cycl.extend(x.sssr)
    prot_gr = cgrs.substructure(all_prot).split()
    coming_gr = cgrs.substructure(all_coming).split()

    prot_list = [set(x) for x in prot_gr]
    coming_list = [set(x) for x in coming_gr]

    all = []
    all.extend(prot_list)
    all.extend(coming_list)
    other_list = cgrs.centers_list

    for x in all:
        index_list = []
        for i, y in enumerate(other_list):
            if x.intersection(y):
                index_list.append(i)
                x.update(y)
        if index_list:
            index_list.reverse()
            uu = []
            for i in index_list:
                uu.extend(other_list.pop(i))
            other_list.append(uu)

    cycles = []
    for x in reaction.reactants:
        cycles.extend(x.sssr)
    for x in reaction.products:
        for y in x.sssr:
            if y not in cycles:
                cycles.extend(y)
    #for y in cgrs.split():
       # for kk in y.sssr:
    for y in cycles:

        for kk in cgrs.substructure(y):
            ept = []
            new_ind = []
            un_in_cy = []
            for x in other_list:
                if set(x).intersection(kk)  and len(set(x).intersection(kk))>1:
                    ept.append(set(x).intersection(kk))

            if len(ept) >= 2:
                hei = y.substructure(kk)
                for x in ept:
                    for i, p in enumerate(x):
                        for i2, m in enumerate(x):
                            if i!=i2:
                                if hei.has_edge(m,p) and (hei.bond(m,p).order == None or hei.bond(m,p).p_order == None):
                                    un_in_cy.extend([m,p])

            if un_in_cy:
                for i3, zop in enumerate(other_list):
                    if set(zop).intersection(un_in_cy):
                        new_ind.append(i3)
            if len(new_ind)>1:
                y = []
                new_ind.reverse()
                for x in new_ind:
                    y.extend(other_list[x])
                    other_list.pop(x)
                other_list.append(y)

    if 1 < len(other_list):

        variants_reaction = []
        reactants_to_work = big_mol(reaction.reactants)
        product_to_work = big_mol(reaction.products)

        list_rc = []

        const_at = atom_re_pr(set(reactants_to_work).difference(all_prot).difference(cgrs.center_atoms),
                              reactants_to_work, product_to_work)
        const_bond = bonds_re_pr(cgrs.bonds(), const_at.keys())

        for big_rc in other_list:

            atom_t = atom_re_pr(set(big_rc).difference(all_prot).difference(all_coming), reactants_to_work, product_to_work)
            prot_atoms = atom_re_pr(prot_come(big_rc, all_prot, prot_gr), reactants_to_work, product_to_work)
            comming_atoms = atom_re_pr(prot_come(big_rc, all_coming, coming_gr), reactants_to_work, product_to_work)

            prot_bonds = bonds_re_pr(cgrs.bonds(), prot_atoms.keys())
            comming_bonds = bonds_re_pr(cgrs.bonds(), comming_atoms.keys())
            t_bond = bonds_re_pr(cgrs.bonds(), atom_t.keys())

            list_rc.append([atom_t, prot_atoms, comming_atoms, t_bond, prot_bonds, comming_bonds])

        for e , t_rc in enumerate(list_rc):
            for state in list(product([1, 0], repeat=len(list_rc) - 1)):

                new_reactant = MoleculeContainer()
                new_product = MoleculeContainer()
                new_all_bonds_reactants = []
                new_all_bonds_products =[]

                new_reactant = add_at(const_at,reactants_to_work, new_reactant)
                new_product = add_at(const_at,product_to_work, new_product)

                for x, y in t_rc[0].items():
                    new_reactant.add_atom(y['re'], x)
                    new_product.add_atom(y['pr'], x)

                for x, y in t_rc[1].items():
                    new_reactant.add_atom(y['re'], x)  # атомы уходящей группы

                for x, y in t_rc[2].items():
                    new_product.add_atom(y['pr'], x) # приходящие атомы

                for x in const_bond.values():
                    new_all_bonds_reactants.append(x['re'])
                    new_all_bonds_products.append(x['pr'])

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
                            new_reactant.add_atom(y['re'], x)
                            new_product.add_atom(y['re'], x)

                        for x, y in big_rc[1].items():
                            new_reactant.add_atom(y['re'], x)
                            new_product.add_atom(y['re'], x)

                        for x in big_rc[3].values():
                            if x['re'][2]:
                                new_all_bonds_reactants.append(x['re'])
                                new_all_bonds_products.append(x['re'])

                        for x in big_rc[4].values():
                            new_all_bonds_reactants.append(x['re'])
                            new_all_bonds_products.append(x['re'])

                    if s == 1:

                        for x, y in big_rc[0].items():
                            new_reactant.add_atom(y['pr'], x)
                            new_product.add_atom(y['pr'], x)

                        for x, y in big_rc[2].items():
                            new_reactant.add_atom(y['pr'], x)
                            new_product.add_atom(y['pr'], x)

                        for x in big_rc[3].values():
                            if x['pr'][2]:
                                new_all_bonds_reactants.append(x['pr'])
                                new_all_bonds_products.append(x['pr'])

                        for x in big_rc[5].values():
                            new_all_bonds_reactants.append(x['pr'])
                            new_all_bonds_products.append(x['pr'])

                [new_reactant.add_bond(*x) for x in new_all_bonds_reactants if not new_reactant.has_edge(x[0], x[1])]
                [new_product.add_bond(*x) for x in new_all_bonds_products if not new_product.has_edge(x[0], x[1])]
                variants_reaction.append(ReactionContainer(reactants = (new_reactant,), products =(new_product,)))
        return variants_reaction
    return [reaction]