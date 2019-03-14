from CGRtools.files import RDFread, RDFwrite, MRVread
from CGRtools.algorithms import sssr
from collections import defaultdict
import copy
from itertools import product
from CGRtools.containers import ReactionContainer, MoleculeContainer, CGRContainer
from CGRtools.preparer import CGRpreparer

def cycl(new_reaction):

    new_reaction.reset_query_marks()


    new_cgr = ~new_reaction
    new_cgr.reset_query_marks()
    # prot_gr = new_cgr.substructure(prot).split()
    # coming_gr = new_cgr.substructure(coming).split()

    cycles = []
    for x in new_reaction.reactants:
        cycles.extend(x.sssr)
    for x in new_reaction.products:
        for y in x.sssr:
            if y not in cycles:
                cycles.extend(y)

    multiple_b_at = []
    center_new_cgr = set(new_cgr.center_atoms)
    hetero_atoms=[]
    for x in center_new_cgr:
        for y in new_cgr[x]:
            if new_cgr.node[y].element != 'C' and new_cgr.node[y].element != 'H' and y not in center_new_cgr:
                hetero_atoms.append(y)
    if center_new_cgr:

        if cycles:
            usefull_cycles = []
            dinamic_cyc = []
            non_dinamic_cyc = []

            for set_cyc in cycles:
                if len(set_cyc) < 5 and center_new_cgr.intersection(set_cyc):  # малые циклы 3 , 4 ноды
                    usefull_cycles.extend(set_cyc)

                if center_new_cgr.intersection(set_cyc) :

                    if all(new_cgr.nodes[x].hybridization == 4 for x in set_cyc) or all(new_cgr.nodes[x].p_hybridization == 4 for x in set_cyc):    # hybridiztion
                        usefull_cycles.extend(set_cyc)

                if any(new_cgr.subgraph(set_cyc).edges[x].order == None for x in new_cgr.subgraph(set_cyc).edges()):  # если замыкание цикла, если нет то цикл не интересен

                    if len(set(set_cyc).intersection(center_new_cgr)) > 2:  # нахождение динамических циклов
                        dinamic_cyc.extend(set_cyc)

                    elif len(set(set_cyc).intersection(center_new_cgr)) == 2:  # если нет динамических циклов нахлждение не динамических циклов , с выбором самого малого цикла
                        if not non_dinamic_cyc:
                            non_dinamic_cyc = set_cyc
                        else:
                            if len(non_dinamic_cyc) > len(set_cyc):
                                non_dinamic_cyc = set_cyc
            if dinamic_cyc:
                usefull_cycles.extend(dinamic_cyc)

            elif non_dinamic_cyc:
                usefull_cycles.extend(non_dinamic_cyc)

            center_new_cgr.update(usefull_cycles)

        for atom in set(center_new_cgr):
            for x in new_cgr[atom]:
                if new_cgr[atom][x].order  and new_cgr[atom][x].order  > 1 or new_cgr[atom][x].p_order and new_cgr[atom][x].p_order > 1:
                    multiple_b_at.append(x)
        center_new_cgr.update(multiple_b_at)

        center_new_cgr.update(hetero_atoms)

    return center_new_cgr, new_reaction