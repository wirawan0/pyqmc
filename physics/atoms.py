# $Id: atoms.py,v 1.2 2010-09-07 15:10:56 wirawan Exp $
#
# pyqmc.physics.atoms module
# Created: 20100903
# Wirawan Purwanto
#
# This module is part of PyQMC project.
#
# Information about atoms
#

# Rigged with the help of Wikipedia,
# http://en.wikipedia.org/wiki/List_of_elements_by_symbol
# - taken from the source code (see pyqmc/_scratch subdir)
# - rigged with _rig_atoms.py

ATOM_LIST = """\
  1 H   hydrogen
  2 He  helium
  3 Li  lithium
  4 Be  beryllium
  5 B   boron
  6 C   carbon
  7 N   nitrogen
  8 O   oxygen
  9 F   fluorine
 10 Ne  neon
 11 Na  sodium natrium
 12 Mg  magnesium
 13 Al  aluminium aluminum
 14 Si  silicon
 15 P   phosphorus phosphorous
 16 S   sulfur sulphur
 17 Cl  chlorine
 18 Ar  argon
 19 K   potassium kalium
 20 Ca  calcium
 21 Sc  scandium
 22 Ti  titanium
 23 V   vanadium
 24 Cr  chromium
 25 Mn  manganese
 26 Fe  iron
 27 Co  cobalt
 28 Ni  nickel
 29 Cu  copper
 30 Zn  zinc
 31 Ga  gallium
 32 Ge  germanium
 33 As  arsenic
 34 Se  selenium
 35 Br  bromine
 36 Kr  krypton
 37 Rb  rubidium
 38 Sr  strontium
 39 Y   yttrium
 40 Zr  zirconium
 41 Nb  niobium
 42 Mo  molybdenum
 43 Tc  technetium
 44 Ru  ruthenium
 45 Rh  rhodium
 46 Pd  palladium
 47 Ag  silver
 48 Cd  cadmium
 49 In  indium
 50 Sn  tin
 51 Sb  antimony
 52 Te  tellurium
 53 I   iodine
 54 Xe  xenon
 55 Cs  caesium cesium
 56 Ba  barium
 57 La  lanthanum
 58 Ce  cerium
 59 Pr  praseodymium
 60 Nd  neodymium
 61 Pm  promethium
 62 Sm  samarium
 63 Eu  europium
 64 Gd  gadolinium
 65 Tb  terbium
 66 Dy  dysprosium
 67 Ho  holmium
 68 Er  erbium
 69 Tm  thulium
 70 Yb  ytterbium
 71 Lu  lutetium
 72 Hf  hafnium
 73 Ta  tantalum
 74 W   tungsten
 75 Re  rhenium
 76 Os  osmium
 77 Ir  iridium
 78 Pt  platinum
 79 Au  gold
 80 Hg  mercury
 81 Tl  thallium
 82 Pb  lead
 83 Bi  bismuth
 84 Po  polonium
 85 At  astatine
 86 Rn  radon
 87 Fr  francium
 88 Ra  radium
 89 Ac  actinium
 90 Th  thorium
 91 Pa  protactinium
 92 U   uranium
 93 Np  neptunium
 94 Pu  plutonium
 95 Am  americium
 96 Cm  curium
 97 Bk  berkelium
 98 Cf  californium
 99 Es  einsteinium
100 Fm  fermium
101 Md  mendelevium
102 No  nobelium
103 Lr  lawrencium
104 Rf  rutherfordium
105 Db  dubnium
106 Sg  seaborgium
107 Bh  bohrium
108 Hs  hassium
109 Mt  meitnerium
110 Ds  darmstadtium
111 Rg  roentgenium
112 Cn  copernicium
113 Uut ununtrium
114 Uuq ununquadium
115 Uup ununpentium
116 Uuh ununhexium
117 Uus ununseptium
118 Uuo ununoctium
"""

class atom:
  """Representation of an atomic species.

  Useful fields are:
  . no
  . symb
  . name()
  . names
  """
  def __init__(self, no, symb, names):
    self.no = no
    self.symb = symb
    self.names = names
  def name(self):
    return self.names[0]
  def __repr__(self):
    return "atom(%d,'%s',%s)" % (self.no, self.symb, self.names)

def parse_atom_list(ATOM_LIST=ATOM_LIST):
  """Internal routine to parse the atom list for the first time."""
  atom_list_by_no = {}
  atom_list_by_symbol = {}
  atom_list_by_name = {}
  atom_list_by_whatever = {}

  for L in ATOM_LIST.split("\n"):
    #print L
    fld = L.split()
    if len(fld) < 3: continue
    no = int(fld[0])
    symb = fld[1]
    names = fld[2:]

    at = atom(no, symb, names)
    atom_list_by_no[no] = at
    atom_list_by_symbol[symb] = at
    atom_list_by_whatever[no] = at
    atom_list_by_whatever[symb] = at
    for n in names:
      atom_list_by_name[n] = at
      atom_list_by_whatever[n] = at

  return (atom_list_by_no, atom_list_by_symbol, atom_list_by_name, atom_list_by_whatever)


(atom_list_by_no,
 atom_list_by_symbol,
 atom_list_by_name,
 atom_list_by_whatever) = parse_atom_list()

def get(symb=None):
  try:
    return atom_list_by_whatever[symb]
  except:
    pass
  # Try three more things
  try:
    no = int(symb)
    return atom_list_by_no[no]
  except:
    pass
  try:
    return atom_list_by_name[symb.lower()]
  except:
    pass
  try:
    if len(symb) > 1:
      return atom_list_by_symbol[symb[0].upper()+symb[1:].lower()]
    else:
      return atom_list_by_symbol[symb.upper()]
  except:
    raise KeyError, "Unknown atomic symbol: %s" % symb

