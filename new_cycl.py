from CGRtools.containers import  CGRContainer
def cycl(new_reaction):

    new_cgr = ~new_reaction

    cycles = []
    for x in new_reaction.reactants:
        cycles.extend(x.sssr)
    for x in new_reaction.products:
        for y in x.sssr:
            if y not in cycles:
                cycles.append(y)

    multiple_b_at = []
    hetero_atoms=[]

    center_new_cgr = set(new_cgr.center_atoms)

    if center_new_cgr:
        for x in center_new_cgr:
            for y in list(new_cgr.bonds()):
                if x == y[0]:
                    if new_cgr.atom(y[1]).atomic_symbol != 'C' and new_cgr.atom(y[1]).atomic_symbol != 'H' and y[1] not in center_new_cgr:
                        hetero_atoms.append(y[1])

        if cycles:
            usefull_cycles = []
            dinamic_cyc = []
            non_dinamic_cyc = []

            for set_cyc in cycles:
                if len(set_cyc) < 5 and center_new_cgr.intersection(set_cyc):  # малые циклы 3 , 4 ноды
                    usefull_cycles.extend(set_cyc)

                if center_new_cgr.intersection(set_cyc) :

                    if all(new_cgr.atom(x).hybridization == 4 for x in set_cyc) or all(new_cgr.atom(x).p_hybridization == 4 for x in set_cyc):    # hybridiztion
                        usefull_cycles.extend(set_cyc)

                if any(x[2].order == None for x in list(new_cgr.substructure(set_cyc).bonds())):  # если замыкание цикла, если нет то цикл не интересен

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

            center_new_cgr.update(set(new_cgr.augmented_substructure(center_new_cgr, deep=1)).difference(center_new_cgr))
            center_new_cgr.update(usefull_cycles)

        for atom in set(center_new_cgr):
            for x in list(new_cgr.bonds()):
                if x[0] == atom and x[1] not in center_new_cgr and (x[2].order and x[2].order > 1 or x[2].p_order and x[2].p_order > 1):
                    multiple_b_at.append(x[1])

        center_new_cgr.update(multiple_b_at)

        center_new_cgr.update(hetero_atoms)

    return center_new_cgr, new_reaction