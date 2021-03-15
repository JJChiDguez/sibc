from struct import pack, unpack
from sidh.montgomery.curve import MontgomeryCurve
from sidh.montgomery.isogeny import Isogeny as MontgomeryIsogeny
from sidh.constants import parameters
from sidh.common import attrdict