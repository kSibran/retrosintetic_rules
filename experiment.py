from MAVR import Mavr
from CGRtools.files import RDFread
theclass = Mavr()

with RDFread('/home/ravil/Desktop/Projects/retro/example.rdf') as reaction_file:
    for reaction in reaction_file:
        cgr = ~reaction
        print(cgr.centers_list)

        zakzak = theclass.enumerate_reaction(reaction)
        for x in zakzak:
            cc =1


def diff_atoms(mols_1, mols_2):
    set_1 = set()
    set_2 = set()
    [set_1.update(set(mol)) for mol in mols_1]
    [set_2.update(set(mol)) for mol in mols_2]
    return set_1.difference(set_2)

cgrs = ~reaction
#объединение реакционных центров по уходящим и приходящим группам
all_prot = diff_atoms(reaction.reactants, reaction.products)
all_coming = diff_atoms(reaction.products, reaction.reactants)

all_prot_come = all_prot + all_coming
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
#конец объединения

# объединение реакционных центров при циклизации в общий для этих рц цикл
cycles = []
for x in reaction.reactants:
    cycles.extend(x.sssr)
for x in reaction.products:
    for y in x.sssr:
        if y not in cycles:
            cycles.append(y)

for y in cycles:
    kk = cgrs.substructure(y)
    ept = []
    new_ind = []
    unite = []
    if all(z[2].order == 4 for z in list(kk.bonds())) and any(z[2].p_order != 4 for z in list(kk.bonds())) \
            or all(z[2].p_order == 4 for z in list(kk.bonds())) and any(z[2].order != 4 for z in list(kk.bonds())):
        unite.extend(y)        # при ароматизации и деароматизации весь цикл это реакционный центр
    else:
        for x in other_list:
            if len(set(x).intersection(kk)) > 1:
                ept.append(set(x).intersection(kk)) # пересечение рц в цикле

        if len(ept) >= 2:
            for x in ept:
                for i, p in enumerate(x):
                    for i2, m in enumerate(x):
                        if i != i2:
                            if any((m == mp[0] and p == mp[1]) or (p == mp[0] and m == mp[1]) for mp in kk.bonds()) and (kk.bond(m, p).order == None or kk.bond(m, p).p_order == None):
                                unite.extend([m, p]) # объединение информмации о всех атомах учавствующих в образовании или разрыве цикла

    if unite:
        for i3, zop in enumerate(other_list):
            if set(zop).intersection(unite):
                new_ind.append(i3)
    if len(new_ind) > 1:
        y = []
        new_ind.reverse()
        for x in new_ind:
            y.extend(other_list[x])
            other_list.pop(x)
        other_list.append(y) #объединение всех атомов реакционных центров в один реакционный центр участвующих в образовании, разрыве или ароматизации цикла
    del (ept, new_ind, unite, kk)
# конец объединения