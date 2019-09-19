from CGRtools.files import RDFread, RDFwrite
from enumeration_reaction import enumeration_cgr
from new_cycl import  cycl
from constructor import constructor
import pickle


det = False
with RDFread():
    fg_fg = {}
    for n, reaction in enumerate(reaction_file, start = 1):
        print(n)

        # if n != 58925:
        #     continue
        # err.write(reaction)
        if n != 180:
            continue
        reaction.meta['id'] = n
        try:
            cgrs = ~reaction
        except ValueError:
            continue

        if cgrs.center_atoms:
            v= cgrs.center_atoms
            if any(x.is_radical or x.p_is_radical for _, x in cgrs.atoms()):
                continue
            if cgrs.center_atoms:
                print('len_fg_fg = ' + str(len(fg_fg)))
                perebor = enumeration_cgr(reaction)

                for new_reaction in perebor:
                    new_reaction.standardize()
                   # new_reaction.reset_query_marks()
                    new_reaction.meta['id'] = n

                    if not constructor(*cycl(new_reaction),fg_fg,n):
                        print('COMPOSITION IS None '+ str(n))
                        det = True
                        break
            if det :
                break
with RDFwrite('/home/ravil/Desktop/Projects/retro/rules_09_19.rdf') as fg:
    for x in fg_fg.values():
        fg.write(x)
