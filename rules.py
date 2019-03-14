from CGRtools.files import RDFread, RDFwrite, MRVread, SDFwrite
from collections import defaultdict
import copy
from itertools import product
from CGRtools.containers import ReactionContainer, MoleculeContainer, CGRContainer
from enumeration_reaction import enumeration_cgr
from cycl import  cycl
from constructor import constructor
import pickle
#fg_fg = defaultdict()
# def enumeration_cgr():
#
#     variant_reaction = ReactionContainer()
#
#     mol_chenge_r = MoleculeContainer()
#     mol_chenge_p = MoleculeContainer()
#
#     for rc in split_centers[e].edges():
#
#         for atom_t in rc:
#
#             for ml_re, ml_pr in zip(reaction.reactants, reaction.products):
#
#                 if atom_t in ml_re:
#
#                     if atom_t not in  mol_chenge_r:
#                         mol_chenge_r.add_atom(ml_re.atom(atom_t), atom_t)
#
#                     for m in ml_re._adj[atom_t].items():
#
#                         if m[0] not in mol_chenge_r:
#                             mol_chenge_r.add_atom(ml_re.atom(m[0]), m[0])
#
#                         mol_chenge_r.add_bond(atom_t, m[0], m[1])
#
#                 if atom_t in ml_pr:
#
#                     if atom_t not in mol_chenge_p:
#                         mol_chenge_p.add_atom(ml_pr.atom(atom_t), atom_t)
#
#                     for m in ml_pr._adj[atom_t].items():
#
#                         if m[0] not in mol_chenge_p:
#                             mol_chenge_p.add_atom(ml_pr.atom(m[0]), m[0])
#
#                         mol_chenge_p.add_bond(atom_t, m[0], m[1])
#
#     for atom_c in constant:
#         for ml_re, ml_pr in zip(reaction.reactants, reaction.products):
#
#             if atom_c in ml_re:
#
#                 if atom_c not in mol_chenge_r:
#                     mol_chenge_r.add_atom(ml_re.atom(atom_c), atom_c)
#
#                 for m in ml_re._adj[atom_c].items():
#
#                     if m[0] not in mol_chenge_r:
#                         mol_chenge_r.add_atom(ml_re.atom(m[0]), m[0])
#                     mol_chenge_r.add_bond(atom_c, m[0], m[1])
#
#             if atom_c in ml_pr:
#
#                 if atom_c not in mol_chenge_p:
#                     mol_chenge_p.add_atom(ml_pr.atom(atom_c), atom_c)
#
#                 for m in ml_pr._adj[atom_c].items():
#
#                     if m[0] not in mol_chenge_p:
#                         mol_chenge_p.add_atom(ml_pr.atom(m[0]), m[0])
#                     mol_chenge_p.add_bond(atom_c, m[0], m[1])
#
#     for s, big_rc in zip(state, split_centers[:e] + split_centers[e + 1:]):
#         #mol_const = MoleculeContainer()
#
#         for rc in big_rc.edges():
#
#             if s == 1:
#
#                 for atom_p in rc:
#
#                     for ml_pr in reaction.products:
#
#                         if atom_p in ml_pr:
#
#                             if atom_p not in mol_chenge_r:
#                                 mol_chenge_r.add_atom(ml_pr.atom(atom_p), atom_p)
#                             if atom_p not in mol_chenge_p:
#                                 mol_chenge_p.add_atom(ml_pr.atom(atom_p), atom_p)
#
#                             for m in  ml_pr._adj[atom_p].items():
#
#                                 if m[0] in rc:
#
#                                     if m[0] not in mol_chenge_r:
#                                         mol_chenge_r.add_atom(ml_pr.atom(m[0]), m[0])
#                                     mol_chenge_r.add_bond(atom_p,m[0],m[1])
#
#                                     if m[0] not in mol_chenge_p:
#                                         mol_chenge_p.add_atom(ml_pr.atom(m[0]), m[0])
#                                     mol_chenge_p.add_bond(atom_p, m[0], m[1])
#
#
#             if s == 0:
#
#                 for atom_r in rc:
#
#                     for ml_re in reaction.reactants:
#
#                         if atom_r in ml_re:
#                             if atom_r not in mol_chenge_r:
#                                 mol_chenge_r.add_atom(ml_re.atom(atom_r), atom_r)
#                             if atom_r not in mol_chenge_p:
#                                 mol_chenge_p.add_atom(ml_re.atom(atom_r), atom_r)
#
#                             for m in ml_re._adj[atom_r].items():
#
#                                 if m[0] in rc:
#
#                                     if m[0] not in mol_chenge_r:
#                                         mol_chenge_r.add_atom(ml_re.atom(m[0]), m[0])
#                                     mol_chenge_r.add_bond(atom_r, m[0], m[1])
#
#                                     if m[0] not in mol_chenge_p:
#                                         mol_chenge_p.add_atom(ml_re.atom(m[0]), m[0])
#                                     mol_chenge_p.add_bond(atom_r, m[0], m[1])
#
#     variant_reaction.reactants.append(mol_chenge_r)
#     variant_reaction.products.append(mol_chenge_p)
        # for atom_c in constant:
        #     for ml_re in reaction.reactants:
        #
        #         if atom_c in ml_re:
        #             mol_const.add_atom(ml_re.atom(atom_c), atom_c)
        #
        #             for m in ml_re._adj[atom_c].items():
        #
        #                 if m[0] not in mol_const:
        #                     mol_const.add_atom(ml_re.atom(m[0]), m[0])
        #                     mol_const.add_bond(atom_c, m[0], m[1])
        #
        # for rc in big_rc.edges():
        #
        #     if s == 1:
        #
        #         for atom_p in rc:
        #
        #             for ml_pr in reaction.products:
        #
        #                 if atom_p in ml_pr:
        #                     mol_const.add_atom(ml_pr.atom(atom_p),
        #                                        atom_p)
        #
        #                     for m in ml_pr._adj[atom_p].items():
        #
        #                         if m[0] in rc:
        #
        #                             if m[0] not in mol_const:
        #                                 mol_const.add_atom(
        #                                     ml_pr.atom(m[0]), m[0])
        #                                 mol_const.add_bond(atom_p, m[0],
        #                                                    m[1])
        #
        #
        #                                 # if new_cgr.edges[rc].get('p_bond'):
        #                                 #     new_cgr.edges[rc]['s_bond'] = new_cgr.edges[rc]['p_bond']# order= bond
        #                                 # else:
        #                                 #     new_cgr.remove_edge(*rc)
        #                                 #     for x in prot_gr:
        #                                 #         if set(x).intersection(rc):
        #                                 #             dl.extend(x)
        #
        #     if s == 0:
        #
        #         for atom_r in rc:
        #             variant_reaction
        #             for ml_re in reaction.reactants:
        #
        #                 if atom_r in ml_re:
        #                     mol_const.add_atom(ml_re.atom(atom_r),
        #                                        atom_r)
        #
        #                     for m in ml_re._adj[atom_r].items():
        #
        #                         if m[0] in rc:
        #
        #                             if m[0] not in mol_const:
        #                                 mol_const.add_atom(
        #                                     ml_re.atom(m[0]), m[0])
        #                                 mol_const.add_bond(atom_r, m[0],
        #                                                    m[1])

        # for big_rc in split_centers[e]:
        #     for rc in big_rc.edges():
        #
        #         for atom_t in rc:
        #             for ml_re in reaction.reactants:
        #
        #                 if atom_t in ml_re:
        #                     mol_chenge_r.add_atom(ml_re.atom(atom_t), atom_t)
        #
        #                     for m in ml_re._adj[atom_r].items():
        #
        #                         if m[0] in rc:
        #
        #                             if m[0] not in mol_chenge_r:
        #                                 mol_chenge_r.add_atom(ml_re.atom(m[0]), m[0])
        #                                 mol_chenge_r.add_bond(atom_r, m[0], m[1])
        #
        #             for ml_pr in reaction.products:
        #
        #                 if atom_t in ml_pr:
        #                     mol_chenge_p.add_atom(ml_pr.atom(atom_t), atom_t)
        #
        #                     for m in ml_pr._adj[atom_t].items():
        #
        #                         if m[0] in rc:
        #
        #                             if m[0] not in mol_chenge_p:
        #                                 mol_chenge_r.add_atom(ml_pr.atom(m[0]), m[0])
        #                                 mol_chenge_r.add_bond(atom_r, m[0], m[1])
                # if new_cgr.edges[rc].get('s_bond'):
                #     new_cgr.edges[rc]['p_bond'] = new_cgr.edges[rc]['s_bond']
                # else:
                #     new_cgr.remove_edge(*rc)
                #     for x in coming_gr:
                #         if set(x).intersection(rc):
                #             dl.extend(x)
    #
    # [new_cgr.delete_atom(x) for x in set(dl)]
# def enumeration_cgr():
#
#     cgrs = ~reaction
#
#     #cgrs.reset_query_marks()
#     # constant = set(x for x in cgrs).difference(prot).difference(coming).difference(set(center_atoms))
#     #split_centers = [cgrs.substructure(x) for x in cgrs.get_centers_list()]
#
#     if len(cgrs.get_centers_list()) > 1:
#         variant_reaction = copy.deepcopy(reaction)
#
#         reagent_mols = defaultdict(list)
#         #reagent_bonds = defaultdict(list)
#         product_mols = defaultdict(list)
#         #product_bonds = defaultdict(list)
#         #center_atoms = cgrs.get_centers_list()
#
#         prot = set(x for x in variant_reaction.reactants for x in x).difference(x for x in variant_reaction.products for x in x)
#         coming = set(x for x in variant_reaction.products for x in x).difference(x for x in variant_reaction.reactants for x in x)
#
#         prot_gr = cgrs.substructure(prot).split()
#         coming_gr = cgrs.substructure(coming).split()
#         conformity = {}
#
#         prot_atoms = defaultdict(list)
#         prot_bonds = defaultdict(list)
#         coming_atoms = defaultdict(list)
#         coming_bonds = defaultdict(list)
#
#         # for molecule in variant_reaction.reactants:
#         #     for x in molecule:
#         #         if x in prot:
#         #             prot_atoms[molecule]
#         for molecule in variant_reaction.reactants:
#             for x in molecule:
#                 for molecule2 in  variant_reaction.products:
#                     if x in molecule2:
#                         conformity[x]=((molecule, molecule.node[x], molecule.adj[x]), (molecule2, molecule2.node[x], molecule2.adj[x]))
#                         # словарь атомов
#         for molecule in variant_reaction.reactants:
#
#             for i in cgrs.get_centers_list():
#                 if set(molecule).intersection(i):
#                     atoms_prot = []
#                     for y in prot_gr:
#                         if set(molecule).intersection(i).intersection(y):
#                             atoms_prot.extend(y)
#                     #reagent_mols[molecule].append([conformity[x][1][0] for x in set(molecule).intersection(i).difference(atoms_prot)])
#                     # поведение атомов кора реагентов в продуктах
#                     reagent_mols[molecule].append(set(molecule).intersection(i).difference(atoms_prot))
#                     #reagent_bonds[molecule].append([conformity[x][1][1] for x in set(molecule).intersection(i).difference(atoms_prot)])
#                     # связи атомов кора реагентов в продуктах
#                     for x in set(molecule).intersection(i).difference(atoms_prot):
#                         prot_atoms[x].append([molecule.node[x] for x in atoms_prot])
#                         # уходящая группа для атомов кора в продуктах
#                         prot_bonds[x].append([molecule.adj[x] for x in atoms_prot])
#                         # связи уходящей группы для атомов кора в продуктах
#                    # kek = set(molecule).intersection(i).difference(atoms_prot)
#                     # for molecule2 in variant_reaction.products:
#                     #     if set(molecule).intersection(i).difference(atoms_prot).intersection(molecule2):
#                     #reagent_mols[molecule].append(([molecule2[x] for x in set(molecule).intersection(i).difference(atoms_prot)], atoms_prot))
#                     #reagent_mols[molecule].append((set(molecule).intersection(i).difference(atoms_prot), atoms_prot))
#                     #reagent_mols[molecule].append(([molecule[x] for x in set(molecule).intersection(i)], [molecule[x] for x in atoms_prot]))
#                     # reagent_mols[molecule].append(([conformity[x][1] for x in set(molecule).intersection(i).difference(atoms_prot)], atoms_prot))
#                     # prot_atoms[molecule].append(())
#
#         for molecule in variant_reaction.products:
#
#             for i in cgrs.get_centers_list():
#                 if set(molecule).intersection(i):
#                     atoms_come = []
#                     for y in coming_gr:
#                         if set(molecule).intersection(i).intersection(y):
#                             atoms_come.extend(y)
#
#                     #product_mols[molecule].append([conformity[x][0][0] for x in set(molecule).intersection(i).difference(atoms_come)])
#                     # поведение атомов кора продуктов в реагентах
#                     product_mols[molecule].append(set(molecule).intersection(i).difference(atoms_come))
#                    # product_bonds[molecule].append([conformity[x][0][1] for x in set(molecule).intersection(i).difference(atoms_come)])
#                     # связи атомов кора продуктов в реагента
#                     for x in set(molecule).intersection(i).difference(atoms_come):
#                         coming_atoms[x].append([molecule.node[x] for x in atoms_come])
#                         # приходящая группа для атомов кора в реагентах
#                         coming_bonds[x].append([molecule.adj[x] for x in atoms_come])
#                         # связи приходящей группы для атомов кора в реагентах
#
#             # for i in cgrs.get_centers_list():
#             #     if set(molecule).intersection(i):
#             #         atoms_come = []
#             #         for y in coming_gr:
#             #             if set(molecule).intersection(i).intersection(y):
#             #                 atoms_come.extend(y)
#             #         product_mols[molecule].append((set(molecule).intersection(i), atoms_come))
#                     #product_mols[molecule].append(([molecule[x] for x in set(molecule).intersection(i)], [molecule[x] for x in atoms_come]))
#
#         for e, x in enumerate(cgrs.get_centers_list()):
#             for state in list(product([1, 0], repeat=len(cgrs.get_centers_list()) - 1)):
#                 cc = cgrs.get_centers_list()
#                 c1c = cgrs.get_centers_list()[:e] + cgrs.get_centers_list()[e+1:]
#                 for s, big_rc in zip(state, cgrs.get_centers_list()[:e] + cgrs.get_centers_list()[e:]):
#                     if s == 0:
#                         for x in variant_reaction.products:
#                             for y in product_mols[x]:
#                                 if y.intersection(big_rc):
#                                     for z in y:
#                                         adj = conformity[z][0][1]
#                                         de = [(u, v) for u, nbrs in adj.items() for v in nbrs]
#                                         x.update(edges=de, nodes=adj)
#                                         if z in prot_atoms:
#                                             for pr_at , pr_bon in zip (prot_atoms[z] ,prot_bonds[z]):
#                                                 for adj in pr_bon:
#                                                     de = [(u, v) for u, nbrs in adj.items() for v in nbrs]
#                                                     x.update(edges=de, nodes=adj)
#                             # for x in reagent_mols[molecule]:
#                             #     if set(big_rc).intersection(x[0]):
#                             #         for p in variant_reaction.products:
#                             #             if set(p).intersection(x[0]):
#                             #                  mlpp = [p[z] for z in set(p).intersection(x[0])]
#                     if s == 1:
#                         for x in variant_reaction.reactants:
#                             for y in reagent_mols[x]:
#                                 if y.intersection(big_rc):
#                                     for z in y:
#                                         adj = conformity[z][0][1]
#                                         de = [(u, v) for u, nbrs in conformity[z][0][1].items() for v in nbrs]
#                                         x.update(edges=de, nodes=adj)
#                                         if z in coming_atoms:
#                                             for com_at, com_bon in zip(coming_atoms[z], coming_bonds[z]):
#                                                 for adj in com_bon:
#                                                     de = [(u, v) for u, nbrs in adj.items() for v in nbrs]
#                                                     x.update(edges=de, nodes=adj)
        # for cgr in cgrs.split():
        #     if set(center_atoms).intersection(cgr):
        #         for e, x in enumerate(cgrs.get_centers_list()):
        #             for state in list(product([1, 0], repeat=len(split_centers) - 1)):
        #                 agraculture = 1






                    # for s, big_rc in zip(state, split_centers[:e] + split_centers[e + 1:]):
                    #     for rc in big_rc.edges():
                    #         if s == 1:
                    #             for x in prot_gr:
                    #                 if set(rc).intersection(x):
                    #                     for y in variant_reaction.reactants:
                    #                         if set(y).intersection(x):

                    #kk = set(rc).intersection(prot)


            #     for x in vr:
            #         if kk <= vr:
            #             [vr.pop(x) for x in kk]
            #     #variant_reaction =
            # if s == 0:
            #     1
    #             if s == 1:
    # fg_file.write(variant_reaction)
    # return variant_reaction
# _______________________________________ready_to_use____________________________________________________________________
# def enumeration_cgr():
#
#     cgrs = ~reaction
#
#     if 1 < len(cgrs.centers_list):
#
#         variants_reaction = []
#         reagent_to_work = MoleculeContainer()
#         product_to_work = MoleculeContainer()
#         const_bond_reagent = []
#         const_bond_product = []
#
#         all_prot = set(x for x in reaction.reactants for x in x).difference(x for x in reaction.products for x in x)
#         prot_gr = cgrs.substructure(all_prot).split()
#         all_coming = set(x for x in reaction.products for x in x).difference(x for x in reaction.reactants for x in x)
#         coming_gr = cgrs.substructure(all_coming).split()
#
#         for molecule in reaction.reactants:
#             reagent_to_work._node.update(molecule._node)
#             reagent_to_work._adj.update(molecule._adj)
#         for molecule in reaction.products:
#             product_to_work._node.update(molecule._node)
#             product_to_work._adj.update(molecule._adj)
#
#         list_rc = []
#         const_at = set(reagent_to_work).difference(all_prot).difference(all_coming).difference(cgrs.center_atoms)
#         #const_at = set(reagent_to_work).difference(all_prot).difference(all_coming)
#         for bon in cgrs.bonds():
#             #if set([bon[0], bon[1]]).issubset(const_at) and not set([bon[0], bon[1]]).issubset((cgrs.center_atoms)): # ? and not set([bon[0], bon[1]]).issubset((cgrs.center_atoms))
#             #if set([bon[0], bon[1]]).issubset(const_at) and not set([bon[0], bon[1]]).issubset((cgrs.center_atoms)): # Issubset or intersection
#             if set([bon[0], bon[1]]).issubset(const_at):
#                 const_bond_reagent.append((bon[0],bon[1],bon[2]._reactant))
#                 const_bond_product.append((bon[0],bon[1],bon[2]._product))
#         for big_rc in cgrs.centers_list:
#
#             atom_t = {}
#             prot_atoms = {}
#             comming_atoms = {}
#             t_bond = []
#             prot_bonds = []
#             comming_bonds = []
#
#             for at in set(big_rc).difference(all_prot).difference(all_coming):
#                 atom_t[at]=(reagent_to_work.atom(at),product_to_work.atom(at))
#
#             if set(big_rc).intersection(all_prot):
#                 for prot in prot_gr:
#                     if set(big_rc).intersection(prot):
#                         for at in prot:
#                             prot_atoms[at] = reagent_to_work.atom(at)
#
#             if set(big_rc).intersection(all_coming):
#                 for come in coming_gr:
#                     if set(big_rc).intersection(come):
#                         for at in come:
#                             comming_atoms[at] = product_to_work.atom(at)
#
#             for bon in cgrs.bonds():
#                 if set(prot_atoms.keys()).intersection([bon[0], bon[1]]):
#                     prot_bonds.append((bon[0],bon[1],bon[2]._reactant))
#                 elif set(comming_atoms.keys()).intersection([bon[0], bon[1]]):
#                     comming_bonds.append((bon[0], bon[1], bon[2]._product))
#                 elif set([bon[0],bon[1]]).intersection((set(big_rc).difference(all_prot).difference(all_coming))): #INTERSECTION OR ISSUBSET
#                     t_bond.append(((bon[0],bon[1], bon[2]._reactant), (bon[0],bon[1], bon[2]._product)))
#             list_rc.append([atom_t, prot_atoms, comming_atoms, t_bond, prot_bonds, comming_bonds])
#
#             #list_rc.append([
#             # [0] atom_t[at]=([0] reagent_to_work.atom(at),[1] product_to_work.atom(at)),
#             # [1] prot_atoms[at] = reagent_to_work.atom(at),
#             # [2] comming_atoms[at] = product_to_work.atom(at),
#             # [3]t_bond = [((bon[0],bon[1], bon[2]._reagent), (bon[0],bon[1], bon[2]._product))],
#             # [4] prot_bonds = [(bon[0],bon[1],bon[2]._reagent)],
#             # [5] comming_bonds = [(bon[0], bon[1], bon[2]._product)]
#             # ])
#
#         for e , t_rc in enumerate(list_rc):
#             for state in list(product([1, 0], repeat=len(list_rc) - 1)):
#                 new_reagent = MoleculeContainer()
#                 new_product = MoleculeContainer()
#                 new_all_bonds_reactants = []
#                 new_all_bonds_products =[]
#
#                 [new_reagent.add_atom(reagent_to_work.atom(x), x) for x in const_at  if x not in new_reagent] # неизменные атомы
#                 [new_reagent.add_atom(list_rc[e][0][x][0], x) for x in list_rc[e][0] if x not in new_reagent] # атомы трансформации реагентные
#                 [new_reagent.add_atom(list_rc[e][1][x], x) for x in list_rc[e][1] if x not in new_reagent]  # атомы уходящей группы
#
#                 [new_all_bonds_reactants.append(x) for x in const_bond_reagent if x[2]]
#                 [new_all_bonds_reactants.append(x[0]) for x in list_rc[e][3] if x[0][2]]
#                 [new_all_bonds_reactants.append(x) for x in list_rc[e][4]]
#
#                 # [new_reagent.add_bond(*x) for x in const_bond_reagent if x[2]]  # связи неизменных
#                 # [new_reagent.add_bond(*x[0]) for x in list_rc[e][3] if x[0][2] and not new_reagent.has_edge(x[0][0],x[0][1])]
#                 # [new_reagent.add_bond(*x) for x in list_rc[e][4] if not new_reagent.has_edge(x[0],x[1])]
#
#                 [new_product.add_atom(product_to_work.atom(x), x) for x in const_at if x not in new_product] # неизменные атомы
#                 [new_product.add_atom(list_rc[e][0][x][1], x) for x in list_rc[e][0] if x not in new_product]  # атомы трансформации продуктов
#                 [new_product.add_atom(list_rc[e][2][x], x) for x in list_rc[e][2] if x not in new_product]  # приходящие атомы
#
#                 [new_all_bonds_products.append(x) for x in const_bond_product if x[2]]
#                 [new_all_bonds_products.append(x[1]) for x in list_rc[e][3] if x[1][2]]
#                 [new_all_bonds_products.append(x) for x in list_rc[e][5] if x[2]]
#
#                 # [new_product.add_bond(*x) for x in const_bond_product if x[2] and not new_product.has_edge(x[0],x[1])]
#                 # [new_product.add_bond(*x[1]) for x in list_rc[e][3] if x[1][2] and not new_product.has_edge(x[1][0],x[1][1])]
#                 # [new_product.add_bond(*x) for x in list_rc[e][5] if x[2] and not new_product.has_edge(x[0],x[1])]
#
#                 for s, big_rc in zip(state, list_rc[:e] + list_rc[e+1:]):
#                     if s == 0:
#
#                         [new_reagent.add_atom(big_rc[0][x][0], x) for x in big_rc[0] if x not in new_reagent]
#                         [new_product.add_atom(big_rc[0][x][0], x) for x in big_rc[0] if x not in new_product]
#                         [new_reagent.add_atom(big_rc[1][x], x) for x in big_rc[1] if x not in new_reagent]
#                         [new_product.add_atom(big_rc[1][x], x) for x in big_rc[1] if x not in new_product]
#
#                         [new_all_bonds_reactants.append(x[0]) for x in big_rc[3] if x[0][2]]
#                         [new_all_bonds_products.append(x[0]) for x in big_rc[3] if x[0][2]]
#                         [new_all_bonds_reactants.append(x) for x in big_rc[4] if x[2]]
#                         [new_all_bonds_products.append(x) for x in big_rc[4] if x[2]]
#
#                         # [new_reagent.add_bond(*x[0]) for x in big_rc[3] if x[0][2] and not new_reagent.has_edge(x[0][0],x[0][1])]
#                         # [new_product.add_bond(*x[0]) for x in big_rc[3] if x[0][2] and not new_product.has_edge(x[0][0],x[0][1])]
#                         # [new_reagent.add_bond(*x) for x in big_rc[4] if x[2] and not new_reagent.has_edge(x[0],x[1])]
#                         # [new_product.add_bond(*x) for x in big_rc[4] if x[2] and not new_product.has_edge(x[0],x[1])]
#
#                         # [new_reagent.add_atom(reagent_to_work.atom(x), x) for x in const_at if
#                         #  x not in new_reagent]  # неизменные атомы
#                         # [new_reagent.add_bond(*x) for x in const_bond_reagent if x[2] and not new_reagent.has_edge(x[0],x[1])]  # связи неизменных
#                         #
#                         # [new_product.add_atom(product_to_work.atom(x), x) for x in const_at if
#                         #  x not in new_product]  # неизменные атомы
#                         # [new_product.add_bond(*x) for x in const_bond_reagent if
#                         #  x[2] and not new_product.has_edge(x[0], x[1])]
#
#                     if s == 1:
#
#                         [new_reagent.add_atom(big_rc[0][x][1], x) for x in big_rc[0] if x not in new_reagent]
#                         [new_product.add_atom(big_rc[0][x][1], x) for x in big_rc[0] if x not in new_product]
#                         [new_reagent.add_atom(big_rc[2][x], x) for x in big_rc[2] if x not in new_reagent]
#                         [new_product.add_atom(big_rc[2][x], x) for x in big_rc[2] if x not in new_product]
#
#                         [new_all_bonds_reactants.append(x[1]) for x in big_rc[3] if x[1][2]]
#                         [new_all_bonds_products.append(x[1]) for x in big_rc[3] if x[1][2]]
#                         [new_all_bonds_reactants.append(x) for x in big_rc[5] if x[2]]
#                         [new_all_bonds_products.append(x) for x in big_rc[5] if x[2]]
#
#                         # [new_reagent.add_bond(*x[1]) for x in big_rc[3] if x[1][2] and not new_reagent.has_edge(x[1][0],x[1][1])]
#                         # [new_product.add_bond(*x[1]) for x in big_rc[3] if x[1][2] and not new_product.has_edge(x[1][0],x[1][1])]
#                         # [new_reagent.add_bond(*x) for x in big_rc[5] if x[2] and not new_reagent.has_edge(x[0],x[1])]
#                         # [new_product.add_bond(*x) for x in big_rc[5] if x[2] and not new_product.has_edge(x[0],x[1])]
#
#                 [new_reagent.add_bond(*x) for x in new_all_bonds_reactants if not new_reagent.has_edge(x[0], x[1])]
#                 [new_product.add_bond(*x) for x in new_all_bonds_products if not new_product.has_edge(x[0], x[1])]
#                 variants_reaction.append(ReactionContainer(reactants = (new_reagent,), products =(new_product,)))
#
#         return variants_reaction
#     else:
#         return [reaction]
#________________________________________________________________________________________________________________________
# def diff_atoms(mols_1, mols_2):
#     set_1 = set()
#     set_2 = set()
#     for mols in mols_1:
#         for at in mols:
#             set_1.add(at)
#     for mols in mols_2:
#         for at in mols:
#             set_2.add(at)
#     return set_1.difference(set_2)
#
# def big_mol(mols):
#     to_work= MoleculeContainer()
#     for x in mols:
#         to_work._node.update(x._node)
#         to_work._adj.update(x._adj)
#     return to_work
#
# def prot_come(rc, prot_or_come, gr):
#     out_set = set()
#     if set(rc).intersection(prot_or_come):
#         for g in gr:
#             if set(rc).intersection(g):
#                 out_set.update(g)
#     return out_set
#
# def atom_re_pr(set_atom, mol_source_re, mol_source_pr):
#     out = defaultdict(dict)
#     for x in set_atom:
#         if x in mol_source_re.atoms_numbers:
#             re = mol_source_re.atom(x)
#         else:
#             re = None
#         if x in mol_source_pr.atoms_numbers:
#             pr = mol_source_pr.atom(x)
#         else:
#             pr = None
#         out[x] ={'re':re, 'pr' : pr}
#     return out
#
# def bonds_re_pr(bonds, atoms):
#     out=defaultdict(dict)
#     for bon in bonds:
#         if set([bon[0], bon[1]]).intersection(atoms):
#             out[str(bon[0])+','+str(bon[1])]={'re':(bon[0], bon[1], bon[2]._reactant),'pr' :(bon[0],bon[1],bon[2]._product)}
#             #out.append(((bon[0], bon[1], bon[2]._reactant), (bon[0], bon[1], bon[2]._product)))
#     return out
#
# def add_at(source, mol_source, out_source):
#     for x in source:
#         if x not in out_source:
#             out_source.add_atom(mol_source.atom(x), x)
#     return out_source
#
# def enumeration_cgr(reaction):
#
#     cgrs = ~reaction
#
#     all_prot = diff_atoms(reaction.reactants, reaction.products)
#     all_coming = diff_atoms(reaction.products, reaction.reactants)
#
#     prot_gr = cgrs.substructure(all_prot).split()
#     coming_gr = cgrs.substructure(all_coming).split()
#
#     united_prot = defaultdict(list)
#     united_come = defaultdict(list)
#     united = []
#     other_list = []
#
#     for x in cgrs.centers_list:
#         flafg = 1
#         # for y in prot_gr:
#         #     if set(x).intersection(y):
#         #
#         #         flafg = 0
#         #         [united_prot[y].append(x) for x in x]
#         #
#         # for y in coming_gr:
#         #
#         #     if set(x).intersection(y):
#         #         flafg = 0
#         #         [united_come[y].append(x) for x in x]
#         for y, z in zip(prot_gr, coming_gr):
#             if set(x).intersection(y):
#                 flafg = 0
#                 [united_prot[y].append(x) for x in x]
#             if set(x).intersection(z):
#                 flafg = 0
#                 [united_come[z].append(x) for x in x]
#         if flafg:
#             other_list.append(x)
#
#
#
#     for y in united_prot.values():
#         flafg = 1
#         for x in united_come.values():
#             if set(x).intersection(y):
#                 flafg=0
#                 x.extend(y)
#                 united.append(set(x))
#         if flafg:
#             other_list.append(y)
#     other_list.extend(united)
#
#     if 1 < len(other_list):
#
#         variants_reaction = []
#         reactants_to_work = big_mol(reaction.reactants)
#         product_to_work = big_mol(reaction.products)
#
#         list_rc = []
#
#         const_at = atom_re_pr(set(reactants_to_work).difference(all_prot).difference(cgrs.center_atoms),
#                               reactants_to_work, product_to_work)
#         const_bond = bonds_re_pr(cgrs.bonds(), const_at.keys())
#
#         for big_rc in other_list:
#
#             atom_t = atom_re_pr(set(big_rc).difference(all_prot).difference(all_coming), reactants_to_work, product_to_work)
#             prot_atoms = atom_re_pr(prot_come(big_rc, all_prot, prot_gr), reactants_to_work, product_to_work)
#             comming_atoms = atom_re_pr(prot_come(big_rc, all_coming, coming_gr), reactants_to_work, product_to_work)
#
#             prot_bonds = bonds_re_pr(cgrs.bonds(), prot_atoms.keys())
#             comming_bonds = bonds_re_pr(cgrs.bonds(), comming_atoms.keys())
#             t_bond = bonds_re_pr(cgrs.bonds(), atom_t.keys())
#
#             list_rc.append([atom_t, prot_atoms, comming_atoms, t_bond, prot_bonds, comming_bonds])
#
#         for e , t_rc in enumerate(list_rc):
#             for state in list(product([1, 0], repeat=len(list_rc) - 1)):
#
#                 new_reactant = MoleculeContainer()
#                 new_product = MoleculeContainer()
#                 new_all_bonds_reactants = []
#                 new_all_bonds_products =[]
#
#                 new_reactant = add_at(const_at,reactants_to_work, new_reactant)
#                 new_product = add_at(const_at,product_to_work, new_product)
#
#                 for x, y in t_rc[0].items():
#                     new_reactant.add_atom(y['re'], x)
#                     new_product.add_atom(y['pr'], x)
#
#                 for x, y in t_rc[1].items():
#                     new_reactant.add_atom(y['re'], x)  # атомы уходящей группы
#
#                 for x, y in t_rc[2].items():
#                     new_product.add_atom(y['pr'], x) # приходящие атомы
#
#                 for x in const_bond.values():
#                     new_all_bonds_reactants.append(x['re'])
#                     new_all_bonds_products.append(x['pr'])
#
#                 for x in t_rc[5].values():
#                     new_all_bonds_products.append(x['pr'])
#
#                 for x in t_rc[3].values():
#                     if x['re'][2]:
#                         new_all_bonds_reactants.append(x['re'])
#                     if x['pr'][2]:
#                         new_all_bonds_products.append(x['pr'])
#
#                 for x in t_rc[4].values():
#                     new_all_bonds_reactants.append(x['re'])
#
#                 # [new_reagent.add_bond(*x) for x in const_bond_reagent if x[2]]  # связи неизменных
#                 # [new_reagent.add_bond(*x[0]) for x in list_rc[e][3] if x[0][2] and not new_reagent.has_edge(x[0][0],x[0][1])]
#                 # [new_reagent.add_bond(*x) for x in list_rc[e][4] if not new_reagent.has_edge(x[0],x[1])]
#
#
#
#                 # [new_product.add_bond(*x) for x in const_bond_product if x[2] and not new_product.has_edge(x[0],x[1])]
#                 # [new_product.add_bond(*x[1]) for x in list_rc[e][3] if x[1][2] and not new_product.has_edge(x[1][0],x[1][1])]
#                 # [new_product.add_bond(*x) for x in list_rc[e][5] if x[2] and not new_product.has_edge(x[0],x[1])]
#
#                 for s, big_rc in zip(state, list_rc[:e] + list_rc[e+1:]):
#                     if s == 0:
#                         for x, y in big_rc[0].items():
#                             new_reactant.add_atom(y['re'], x)
#                             new_product.add_atom(y['re'], x)
#
#                         for x, y in big_rc[1].items():
#                             new_reactant.add_atom(y['re'], x)
#                             new_product.add_atom(y['re'], x)
#
#
#                         # [new_reagent.add_atom(big_rc[0][x][0], x) for x, y in big_rc[0].items() if x not in new_reagent]
#                         # [new_product.add_atom(big_rc[0][x][0], x) for x in big_rc[0] if x not in new_product]
#                         # [new_reagent.add_atom(big_rc[1][x], x) for x in big_rc[1] if x not in new_reagent]
#                         #[new_product.add_atom(big_rc[1][x], x) for x in big_rc[1] if x not in new_product]
#
#                         for x in big_rc[3].values():
#                             if x['re'][2]:
#                                 new_all_bonds_reactants.append(x['re'])
#                                 new_all_bonds_products.append(x['re'])
#
#                         for x in big_rc[4].values():
#                             new_all_bonds_reactants.append(x['re'])
#                             new_all_bonds_products.append(x['re'])
#
#                         # [new_all_bonds_reactants.append(x[0]) for x in big_rc[3] if x[0][2]]
#                         # [new_all_bonds_products.append(x[0]) for x in big_rc[3] if x[0][2]]
#                         # [new_all_bonds_reactants.append(x) for x in big_rc[4] if x[2]]
#                         # [new_all_bonds_products.append(x) for x in big_rc[4] if x[2]]
#
#                         # [new_reagent.add_bond(*x[0]) for x in big_rc[3] if x[0][2] and not new_reagent.has_edge(x[0][0],x[0][1])]
#                         # [new_product.add_bond(*x[0]) for x in big_rc[3] if x[0][2] and not new_product.has_edge(x[0][0],x[0][1])]
#                         # [new_reagent.add_bond(*x) for x in big_rc[4] if x[2] and not new_reagent.has_edge(x[0],x[1])]
#                         # [new_product.add_bond(*x) for x in big_rc[4] if x[2] and not new_product.has_edge(x[0],x[1])]
#
#                         # [new_reagent.add_atom(reagent_to_work.atom(x), x) for x in const_at if
#                         #  x not in new_reagent]  # неизменные атомы
#                         # [new_reagent.add_bond(*x) for x in const_bond_reagent if x[2] and not new_reagent.has_edge(x[0],x[1])]  # связи неизменных
#                         #
#                         # [new_product.add_atom(product_to_work.atom(x), x) for x in const_at if
#                         #  x not in new_product]  # неизменные атомы
#                         # [new_product.add_bond(*x) for x in const_bond_reagent if
#                         #  x[2] and not new_product.has_edge(x[0], x[1])]
#
#                     if s == 1:
#
#                         for x, y in big_rc[0].items():
#                             new_reactant.add_atom(y['pr'], x)
#                             new_product.add_atom(y['pr'], x)
#
#                         for x, y in big_rc[2].items():
#                             new_reactant.add_atom(y['pr'], x)
#                             new_product.add_atom(y['pr'], x)
#                         # [new_reagent.add_atom(big_rc[0][x][1], x) for x in big_rc[0] if x not in new_reagent]
#                         # [new_product.add_atom(big_rc[0][x][1], x) for x in big_rc[0] if x not in new_product]
#                         # [new_reagent.add_atom(big_rc[2][x], x) for x in big_rc[2] if x not in new_reagent]
#                         # [new_product.add_atom(big_rc[2][x], x) for x in big_rc[2] if x not in new_product]
#                         for x in big_rc[3].values():
#                             if x['pr'][2]:
#                                 new_all_bonds_reactants.append(x['pr'])
#                                 new_all_bonds_products.append(x['pr'])
#
#                         for x in big_rc[5].values():
#                             new_all_bonds_reactants.append(x['pr'])
#                             new_all_bonds_products.append(x['pr'])
#                         # [new_all_bonds_reactants.append(x[1]) for x in big_rc[3] if x[1][2]]
#                         # [new_all_bonds_products.append(x[1]) for x in big_rc[3] if x[1][2]]
#                         # [new_all_bonds_reactants.append(x) for x in big_rc[5] if x[2]]
#                         # [new_all_bonds_products.append(x) for x in big_rc[5] if x[2]]
#
#                         # [new_reagent.add_bond(*x[1]) for x in big_rc[3] if x[1][2] and not new_reagent.has_edge(x[1][0],x[1][1])]
#                         # [new_product.add_bond(*x[1]) for x in big_rc[3] if x[1][2] and not new_product.has_edge(x[1][0],x[1][1])]
#                         # [new_reagent.add_bond(*x) for x in big_rc[5] if x[2] and not new_reagent.has_edge(x[0],x[1])]
#                         # [new_product.add_bond(*x) for x in big_rc[5] if x[2] and not new_product.has_edge(x[0],x[1])]
#
#                 [new_reactant.add_bond(*x) for x in new_all_bonds_reactants if not new_reactant.has_edge(x[0], x[1])]
#                 [new_product.add_bond(*x) for x in new_all_bonds_products if not new_product.has_edge(x[0], x[1])]
#                 variants_reaction.append(ReactionContainer(reactants = (new_reactant,), products =(new_product,)))
#         return variants_reaction
#     return [reaction]

# def cycl(new_reaction):
#
#     new_reaction.reset_query_marks()
#     cycles = []
#     # prot = set(new_reaction.reagents).difference(new_reaction.products)
#     # coming = set(new_reaction.products).difference(new_reaction.reagents)
#
#     new_cgr = ~new_reaction
#
#     new_cgr.reset_query_marks()
#     # prot_gr = new_cgr.substructure(prot).split()
#     # coming_gr = new_cgr.substructure(coming).split()
#     try:
#         cycles = new_cgr.sssr
#     except:
#         'No path between 1 and 29'
#     # if not cycles:                                             #  для проверки возникшей ошибки с отсутствующим циклом
#     #     with SDFwrite('/home/ravil/Desktop/Toster_peka.sdf',extralabels = True) as peka:
#     #         peka.write(new_cgr)
#     multiple_b_at = []
#     center_new_cgr = set(new_cgr.center_atoms)
#     hetero_atoms=[]
#     for x in center_new_cgr:
#         for y in new_cgr[x]:
#             if new_cgr.node[y].element != 'C' and new_cgr.node[y].element != 'H' and y not in center_new_cgr:
#                 hetero_atoms.append(y)
#     if center_new_cgr:
#
#         if cycles:
#             usefull_cycles = []
#             dinamic_cyc = []
#             non_dinamic_cyc = []
#
#             for set_cyc in cycles:
#                 if len(set_cyc) < 5 and center_new_cgr.intersection(set_cyc):  # малые циклы 3 , 4 ноды
#                     usefull_cycles.extend(set_cyc)
#
#                 if center_new_cgr.intersection(set_cyc) :
#
#                     if all(new_cgr.nodes[x].hybridization == 4 for x in set_cyc) or all(new_cgr.nodes[x].p_hybridization == 4 for x in set_cyc):    # hybridiztion
#                         usefull_cycles.extend(set_cyc)
#
#                 # for x in new_cgr.subgraph(set_cyc).edges():
#                 #     if new_cgr.subgraph(set_cyc).edges[x].order == None:
#                 #         pp = 1
#                 if any(new_cgr.subgraph(set_cyc).edges[x].order == None for x in new_cgr.subgraph(set_cyc).edges()):  # если замыкание цикла, если нет то цикл не интересен
#
#                     if len(set(set_cyc).intersection(center_new_cgr)) > 2:  # нахождение динамических циклов
#                         dinamic_cyc.extend(set_cyc)
#
#                     elif len(set(set_cyc).intersection(center_new_cgr)) == 2:  # если нет динамических циклов нахлждение не динамических циклов , с выбором самого малого цикла
#                         if not non_dinamic_cyc:
#                             non_dinamic_cyc = set_cyc
#                         else:
#                             if len(non_dinamic_cyc) > len(set_cyc):
#                                 non_dinamic_cyc = set_cyc
#             if dinamic_cyc:
#                 usefull_cycles.extend(dinamic_cyc)
#
#             elif non_dinamic_cyc:
#                 usefull_cycles.extend(non_dinamic_cyc)
#
#             center_new_cgr.update(usefull_cycles)
#
#         for atom in set(center_new_cgr):
#             for x in new_cgr[atom]:
#                 if new_cgr[atom][x].order  and new_cgr[atom][x].order  > 1 or new_cgr[atom][x].p_order and new_cgr[atom][x].p_order > 1:
#                     multiple_b_at.append(x)
#         center_new_cgr.update(multiple_b_at)
#
#         # for x in prot_gr:
#         #     if center_new_cgr.intersection(x):
#         #         center_new_cgr.update(x)
#         # for x in coming_gr:
#         #     if center_new_cgr.intersection(x):
#         #         center_new_cgr.update(x)
#
#         center_new_cgr.update(hetero_atoms)
#
#     return center_new_cgr, new_reaction

# def constructor(center_atoms, new_reaction,fg_fg):
#     new_cgr = ~new_reaction
#     c_a = new_cgr.center_atoms
#     prot = set(new_reaction.reagents).difference(new_reaction.products)
#     coming = set(new_reaction.products).difference(new_reaction.reagents)
#     prot_gr = new_cgr.substructure(prot).split()
#     coming_gr = new_cgr.substructure(coming).split()
#
#     all_re = MoleculeContainer()
#     all_pr =MoleculeContainer()
#     for x in prot_gr:
#         if center_atoms.intersection(x):
#             center_atoms.update(x)
#     for x in coming_gr:
#         if center_atoms.intersection(x):
#             center_atoms.update(x)
#
#     for mol in new_reaction.reactants:
#         all_re._node.update(mol._node)
#         all_re._adj.update(mol._adj)
#     for mol in new_reaction.products:
#         all_pr._node.update(mol._node)
#         all_pr._adj.update(mol._adj)
#     #fg = new_cgr.substructure(center_atoms)
#    # fg_center = fg.substructure(set(fg).difference(prot).difference(coming))
#     for_dell_re = set(all_re).difference(center_atoms)
#     for_dell_pr =set(all_pr).difference(center_atoms)
#     [all_re.delete_atom(x) for x in for_dell_re]
#     [all_pr.delete_atom(x) for x in  for_dell_pr]
#     # if center_atoms:
#     #
#     #     for x in prot:
#     #         if x in re:
#     #             fg.nodes[x].pop('s_hyb', 0)
#     #             fg.nodes[x].pop('p_hyb', 0)
#     #     for x in coming:
#     #         if x in fg:
#     #             fg.nodes[x].pop('s_hyb', 0)
#     #             fg.nodes[x].pop('p_hyb', 0)
#
#     flag = 0
#     for x in all_re.split():
#         if set(c_a).intersection(x):
#             flag += 1
#
#     fg = ReactionContainer(reactants=(all_re,) , products=(all_pr,))
#
#     if flag == 2:
#         fg.meta['composition'] = '2'
#     elif flag == 1:
#         fg.meta['composition'] = '1'
#     elif flag > 2:
#         fg.meta['composition'] = 'multi'
#     str_graph = str(fg)
#
#     if str_graph not in fg_fg:
#         fg_fg[str_graph] = fg
#         fg_fg[str_graph].meta['id'] = [n]
#     # else:
#     #     fg_fg[str_graph].meta['id'].append(n)
#     return fg_fg
def diff_atoms(mols_1, mols_2):
    set_1 = set()
    set_2 = set()
    for mols in mols_1:
        for at in mols:
            set_1.add(at)
    for mols in mols_2:
        for at in mols:
            set_2.add(at)
    return set_1.difference(set_2)

with  open('/home/ravil/reaxys.0001.pickle' , 'rb') as data, \
    RDFwrite('/home/ravil/Desktop/gg.rdf', 'w') as ffffffff,\
    RDFwrite('/home/ravil/Desktop/Error.rdf') as err:
#MRVread ('/home/ravil/Desktop/query_graphs_and_other/base/db_task/readliner_out/dart_standart_0.mrv', remap = False) as reaction_file,\
#db_task/readliner_out/diol.rdf
# db.rdf       RDFread('/home/ravil/Desktop/345.rdf') as reaction_file,\
    reaction_file = pickle.load(data)
    fg_fg = {}
    for n, reaction in enumerate(reaction_file, start = 1):
        if n < 48:
            continue
        if n > 48:
            break

        print(n)
        all_prot = diff_atoms(reaction.reactants, reaction.products)
        all_coming = diff_atoms(reaction.products, reaction.reactants)

        reaction.meta['id'] = n
        try:
            cgrs = ~reaction
        except ValueError:
            continue

        ffffffff.write(reaction)
        reaction.meta.clear()
        if cgrs.center_atoms:
            pomidorka = enumeration_cgr(reaction, all_prot, all_coming)

            for new_reaction in pomidorka:
                new_reaction.standardize()
                new_reaction.meta['id'] = n
                ffffffff.write(new_reaction)
                constructor(*cycl(new_reaction),fg_fg, all_prot, all_coming , n)

with RDFwrite('/home/ravil/Desktop/grg.rdf', 'w') as ffffffff:
    for x in fg_fg.values():
        ffffffff.write(x)




