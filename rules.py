from CGRtools.files import RDFread, RDFwrite, MRVread, SDFwrite
from collections import defaultdict
import copy
from itertools import product
from CGRtools.containers import ReactionContainer, MoleculeContainer, CGRContainer
from enumeration_reaction import enumeration_cgr
from cycl import  cycl
from constructor import constructor
import pickle

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

with  open('/home/ravil/reaxys.0000.pickle' , 'rb') as data,\
    RDFwrite('/home/ravil/Desktop/gg.rdf', 'w') as ffffffff,\
    RDFwrite('/home/ravil/Desktop/Error.rdf') as err: #'/home/ravil/Desktop/query_graphs_and_other/base/db_task/readliner_out/diol.rdf'

    reaction_file = pickle.load(data)
    fg_fg = {}
    for n, reaction in enumerate(reaction_file, start = 1):
        print(n)
        reac = (reaction.reactants,reaction.products)
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




