import sys
import numpy as np
from pymses.utils import constants as C
from pymses.filters import CellsToPoints
from functools import partial as partial

GRAV_COMPAT = True

"""
* To select column i in a dataset created by pymses, use (..., i) vs (:,i)
(:,i) is correct for slicing, but not for raytracing
whereas (..., i) is more general
* 
"""
class GenericQuantity(object):
    name = None
    amrFields = None
    _func = None
    _factorFunc = None
    unit = None
    texName = None
    texUnitName = None

    @classmethod
    def process(cls, ro, dset=None, unit=None, verbose=False):
        if dset == None:
            source = ro.amr_source(cls.amrFields,
                                   grav_compat=GRAV_COMPAT,
                                   verbose=verbose)
            dset = CellsToPoints(source).flatten()
        return cls.get_calc(ro, unit)(dset)

    @classmethod
    def get_calc(cls, ro, unit=None):
        """
        Returns a function f evaluating a Pointdataset.
        This function f returns : _func(dset)*factor
        """
        if unit is None:
            unit = cls.unit
        # varplot = partial(lambda dset,func,factor: func(dset)*factor,
        #                      func=cls._func,
        #                      factor=cls.factorFunc(ro).express(unit))
        varplot = lambda dset: cls._func(dset) * cls.factorFunc(ro).express(unit)
        return varplot

    @classmethod
    def label(cls):
        if cls.unit is C.none:
            return r"$%s$" % cls.texName
        else:
            return r"$%s \ \left( \ %s \ \right)$" %(cls.texName, cls.texUnitName)

##########################################################################
# Scalar decorator
##########################################################################
def absDecorator(quantity):
    class abs(quantity):
        name = quantity.name+"_abs"
        texName = r"\vert "+quantity.texName+r" \vert"
        _func = lambda dset: np.abs(quantity._func(dset))
        _func = staticmethod(_func)
    return abs

def ScalarDecorator(quantity):
    class ScalarQuantity(quantity):
        pass
    ScalarQuantity.abs = absDecorator(quantity)
    return ScalarQuantity

##########################################################################
# Vector Decorator
##########################################################################
class xyzDecorator:
    dict = {'x':0, 'y':1, 'z':2}
    def __init__(self, item):
        self.item = item
    def __call__(self, quantity):
        @ScalarDecorator
        class xyz(quantity):
            name = quantity.name+"_%s" % self.item
            texName = r"\rm "+quantity.texName+"_{%s}" % self.item
            _func = lambda dset: quantity._func(dset)[..., self.dict[self.item]]
            _func = staticmethod(_func)
        return xyz

class dotDecorator:
    def __init__(self, u):
        self.u = u
    def __call__(self, quantity):
        @ScalarDecorator
        class proj(quantity):
            name = quantity.name+"_u"
            texName = r"\bf "+quantity.texName+r" \cdot u"
            _func = lambda dset: np.dot(quantity._func(dset), self.u)
            _func = staticmethod(_func)
        return proj

def normDecorator(quantity):
    @ScalarDecorator
    class norm(quantity):
        name = quantity.name+"_norm"
        texName = r"\bf \Vert "+quantity.texName+r" \Vert"
        _func = lambda dset: np.linalg.norm(quantity._func(dset), axis=-1)
        _func = staticmethod(_func)
    return norm

def VectorDecorator(quantity):
    class VectorQuantity(quantity):
        pass
    VectorQuantity._original = quantity
    VectorQuantity.dot = classmethod(lambda cls, u: dotDecorator(u)(cls._original))
    VectorQuantity.x = xyzDecorator('x')(quantity)
    VectorQuantity.y = xyzDecorator('y')(quantity)
    VectorQuantity.z = xyzDecorator('z')(quantity)
    VectorQuantity.norm = normDecorator(quantity)
    VectorQuantity.texName = r"\bf "+quantity.texName
    return VectorQuantity

##########################################################################
# MultivaluedDecorator
##########################################################################
def totDecorator(quantity):
    @ScalarDecorator
    class tot(quantity):
        name = quantity.name+"_tot"
        texName = r"\rm "+quantity.texName+"^{tot}"
        _func = lambda dset: np.sum(quantity._func(dset), axis=-1)
        _func = staticmethod(_func)
    return tot

class grpDecorator:
    def __init__(self, i):
        self.i = i
    def __call__(self, quantity):
        @ScalarDecorator
        class igrp(quantity):
            name = quantity.name+"_%d" % self.i
            texName = r"\rm "+quantity.texName+r"^{%d}" % self.i
            _func = lambda dset: quantity._func(dset)[..., self.i]
            _func = staticmethod(_func)
        return igrp

def MultiValuedDecorator(quantity):
    class MultiValuedQuantity(quantity):
        pass
    MultiValuedQuantity._original = quantity
    MultiValuedQuantity.grp = classmethod(lambda cls, i: grpDecorator(cls._original))
    MultiValuedQuantity.tot = totDecorator(quantity)
    return MultiValuedQuantity

##########################################################################
#
##########################################################################
def get_amrfields_to_read(*fields):
    fields = [field for field in fields if field is not None]
    amr_fields_to_read = []
    for field in fields:
        amr_fields_to_read.extend(field.amrFields)
    amr_fields_to_read = list(set(amr_fields_to_read))
    return amr_fields_to_read

@ScalarDecorator
class rho(GenericQuantity):
    name = "rho"
    amrFields = ["rho"]
    _func = staticmethod(lambda dset: dset["rho"])
    factorFunc = staticmethod(lambda ro: ro.info["unit_density"])
    unit = C.g/C.cm**3
    texName = r"\rm \rho"
    texUnitName = r"\rm g.cm^{-3}"

@ScalarDecorator
class T(GenericQuantity):
    name = "T"
    amrFields = ["rho","P"]
    _func = staticmethod(lambda dset: dset["P"]/dset["rho"])
    @staticmethod
    def factorFunc(ro):
        return ro.info["unit_temperature"] * ro.info["mu_gas"]
    unit = C.K
    texName = r"\rm T"
    texUnitName = r"\rm K"

@ScalarDecorator
class mass(GenericQuantity):
    name = "mass"
    amrFields = ["rho"]
    @staticmethod
    def _func(dset):
        if 'size' in dset.fields.keys():
            mass = dset["rho"] * dset["size"]**3
        elif 'level' in dset.fields.keys():
            mass = dset["rho"] * 1./2**(3*dset["level"])
        else:
            sys.exit(1)
        return mass
    @staticmethod
    def factorFunc(ro):
        return ro.info["unit_density"] * ro.info["unit_length"]**3
    unit = C.Msun
    texName = r"\rm mass"
    texUnitName = r"\rm M_{\odot}"

@ScalarDecorator
class Mjeans(GenericQuantity):
    name = "Mjeans"
    amrFields = ["P","rho"]
    _func = staticmethod(lambda dset: dset["P"]**1.5 / dset["rho"]**2)
    @staticmethod
    def factorFunc(ro):
        coeff = np.pi**2.5/6
        num = (ro.info["unit_pressure"]/C.G)**1.5
        den = ro.info["unit_density"]**2
        return coeff * num / den
    unit = C.Msun
    texName = r"\rm M_{jeans}"
    texUnitName = r"\rm M_{\odot}"

@ScalarDecorator
class FracMjeans(GenericQuantity):
    name = "FracMjeans"
    amrFields = ["P","rho"]
    @staticmethod
    def _func(dset):
        if 'size' in dset.fields.keys():
            mass = dset["rho"] * dset["size"]**3
        elif 'level' in dset.fields.keys():
            mass = dset["rho"] * 1./2**(3*dset["level"])
        else:
            sys.exit(1)
        Mjeans = dset["P"]**1.5 / dset["rho"]**2
        return mass/Mjeans
    @staticmethod
    def factorFunc(ro):
        coeff = np.pi**2.5/6
        num = (ro.info["unit_pressure"]/C.G)**1.5
        den = ro.info["unit_density"]**2
        Mjeans = coeff * num / den
        mass = ro.info["unit_density"] * ro.info["unit_length"]**3
        return mass / Mjeans
    unit = C.none
    texName = r"\rm \frac{m}{M_{jeans}}"
    texUnitName = r""

@ScalarDecorator
class phi(GenericQuantity):
    name = "phi"
    amrFields = ["phi"]
    _func = staticmethod(lambda dset: dset['phi'])
    @staticmethod
    def factorFunc(ro):
        return ro.info["unit_velocity"] / ro.info["unit_time"]
    unit = C.m/C.s**2
    texName = r"\phi"
    texUnitName = r"\rm m.s^{-2}"

@ScalarDecorator
class Eint(GenericQuantity):
    name = "eint"
    amrFields = ["eint"]
    _func = staticmethod(lambda dset: dset['eint'])
    @staticmethod
    def factorFunc(ro):
        return ro.info["unit_density"] * ro.info["unit_velocity"]**2
    unit = C.erg/C.cm**3
    texName = r"e_{int}"
    texUnitName = r"\rm J"

@VectorDecorator
class vel(GenericQuantity):
    name = "vel"
    amrFields = ["vel"]
    _func = staticmethod(lambda dset: dset['vel'])
    factorFunc = staticmethod(lambda ro: ro.info["unit_velocity"])
    unit = C.km/C.s
    texName = r"v"
    texUnitName = r"\rm km.s^{-1}"

@VectorDecorator
class B(GenericQuantity):
    name = "B"
    amrFields = ["Bl","Br"]
    _func = staticmethod(lambda dset: (dset['Bl']+dset['Br'])/2.)
    @staticmethod
    def factorFunc(ro):
        coeff = 4*np.pi
        factor = coeff * ro.info['unit_density']*ro.info['unit_velocity']**2
        return factor.express(C.erg/C.cm**3)**0.5 * C.Gauss
    unit = C.Gauss
    texName = r"B"
    texUnitName = r"\rm G"

@VectorDecorator
class g(GenericQuantity):
    name = "g"
    amrFields = ["g"]
    @staticmethod
    def _func(dset):
        return dset['g']
    @staticmethod
    def factorFunc(ro):
        return ro.info["unit_velocity"] / ro.info["unit_time"]
    unit = C.m/C.s**2
    texName = r"g"
    texUnitName = r"\rm m.s^{-2}"

@MultiValuedDecorator
class Er(GenericQuantity):
    name = "Er"
    amrFields = ['Er']
    @staticmethod
    def _func(dset):
        return dset['Er']
    @staticmethod
    def factorFunc(ro):
        return ro.info["unit_density"]*ro.info["unit_velocity"]**2
    unit = C.erg/C.cm**3
    texName = r"E_{r}"
    texUnitName = r"\rm erg.cm^{-3}"

def PositionDecorator(quantity):
    class Position(quantity):
        @staticmethod
        def t(x):
            @VectorDecorator
            class T(quantity):
                _func = staticmethod(lambda dset: quantity._func(dset) - np.array(x))
            return T
    return Position

@VectorDecorator
@PositionDecorator
class pos(GenericQuantity):
    name = "pos"
    amrFields = []
    _func = staticmethod(lambda dset: dset.points)
    factorFunc = staticmethod(lambda ro: ro.info["unit_length"])
    unit = C.au
    texName = r"r"
    texUnitName = r"\rm AU"
