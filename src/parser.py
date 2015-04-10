#!/usr/bin/env python

import sys
from pyparsing import *

# --- Define Grammars for the parser (pyparsing)
integer       = Word(nums)            ### '0123456789'
StrangeName   = Word(printables)      ### '0123456789abc...wxyzABC...WXYZ!"#$%&\\' ()*+,-./:;<=>?@[\\]^_`{|}~'
integer       = Word(nums)            ### '0123456789'
number        = Word(alphanums)       ### 'abc...wyzABC...WXYZ0123456789'
decimalNumber = Combine((Optional(Literal("-"))+Optional(integer)+Optional(Literal("."))+integer))
name          = Word(alphas)          ### 'abc...wxyzABC...WXYZ'


endLine       = Literal("\n")
end           = Literal("\n").suppress()  # go to the end of the line, and suppress it, same as EOL?
EOL           = LineEnd().suppress()      # go to the end of the line, and suppress it, same as end?
all           = SkipTo(end)               # go to the end of the line, match the next line
ligneParLigne = ZeroOrMore(SkipTo('\n').setResultsName("ligne")).setResultsName("lignes")

element = oneOf( """H He Li Be B C N O F Ne Na Mg Al Si P S Cl 
Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge
As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag
Cd In Sn Sb Te I Xe Cs Ba Lu Hf Ta W Re Os
Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Lr Rf
Db Sg Bh Hs Mt Ds Rg Uub Uut Uuq Uup Uuh Uus
Uuo La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm
Yb Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No""" )

