#! /usr/bin/env python
"""
Creates an alpha-helix given an input sequence 
"""

from BuildStructure import placePeptide

seq = "AAA"
placePeptide(seq, [(-57, -47)] * len(seq), rotlib="Richardson.mode")


