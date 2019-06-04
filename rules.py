from CGRtools.files import RDFread, RDFwrite
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

with  RDFread('/home/ravil/2_5_reactions.rdf') as reaction_file,\
    RDFwrite('/home/ravil/Desktop/Error.rdf') as err: #'/home/ravil/Desktop/query_graphs_and_other/base/db_task/readliner_out/diol.rdf'

    fg_fg = {}
    for n, reaction in enumerate(reaction_file, start = 1):
        print(n)
        all_prot = diff_atoms(reaction.reactants, reaction.products)
        all_coming = diff_atoms(reaction.products, reaction.reactants)

        reaction.meta['id'] = n
        try:
            cgrs = ~reaction
        except ValueError:
            continue

        if cgrs.center_atoms:

            print('len_fg_fg = ' + str(len(fg_fg)))
            perebor = enumeration_cgr(reaction, all_prot, all_coming)
            for new_reaction in perebor:
                new_reaction.standardize()
                new_reaction.reset_query_marks()
                new_reaction.meta['id'] = n
                constructor(*cycl(new_reaction),fg_fg, all_prot, all_coming , n)


with open('/home/ravil/Desktop/retro_rules_with_labels', 'wb') as fgfgfg:
    pickle.dump(fg_fg, fgfgfg)




