from CGRtools.files import RDFread, RDFwrite
from enumeration_reaction import enumeration_cgr
from new_cycl import  cycl
from constructor import constructor

class Mavr:
    def enumerate_reaction(self,data):
        return enumeration_cgr(data)
    #def c
__all__ = ['Mavr']