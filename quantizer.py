#!/usr/bin/env python

import sys
import os
import gzip
import argparse

# Quantizers used for testing of this program
def F_(qsc):
    return ord('_') - 33

def FR(qsc):
    return ord('R') - 33

# Quantizers
def F10(qsc):
    return 10

def Q2(qsc):
    if qsc <= 7:
        return 5
    else:
        return 15

def Q4(qsc):
    if qsc <= 7:
        return 5
    elif qsc <= 13:
        return 12
    elif qsc <= 19:
        return 18
    else:
        return 24

def Q8(qsc):
    if qsc <= 6:
        return 5
    elif qsc <= 11:
        return 10
    elif qsc <= 16:
        return 15
    elif qsc <= 21:
        return 20
    elif qsc <= 26:
        return 25
    elif qsc <= 31:
        return 30
    elif qsc <= 36:
        return 35
    else:
        return 40


def Q0(qsc):
    return qsc

def QuantFromName(qname):
    if qname == "Q2":
        return Q2
    elif qname == "Q4":
        return Q4
    elif qname == "Q8":
        return Q8
    elif qname == "Q0":
        return Q0
    elif qname == "F10":
        return F10
    elif qname == "F_":
        return F_
    elif qname == "FR":
        return FR
    else:
        raise SystemError("Error: Undefined quantizer " + qname + "\n")


# Split a read into four components
def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

# Test if file is in zip format
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


##################################
# Main program
##################################
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Input fastq file, either in plain text formtat or gzip compressed.")
parser.add_argument("Qprimary", help="Quantizer used in non-repetitive regions. One of Q2, Q4, Q8, or F10.")
parser.add_argument("Qrun", help="Quantizer used near repetitive regions. One of Q2, Q4, Q8, or F10. Specify the same as Qprimary to use a fixed quantizer.")
parser.add_argument("output_file", help="Output fastq file in plain text formtat.")
args = parser.parse_args()

# Define specific quantizers from command line parameters
Qprimary = QuantFromName(args.Qprimary)
Qrun = QuantFromName(args.Qrun)

# All qscores within distanceToRun from a run of the same base
# or a run of the same dinucleotide
# are quantized with Qrun.
# The other qscores are quantized with Qprimary
distanceToRun = 5

# The minimum number of repetitions of the same base that is considered a run
runTheshold = 5

# The minimum number of repetitions of the same pair of bases
# that is considered a dinucleotide run
diRunTheshold = 4

if not os.path.exists(args.input_file):
    raise SystemError("Error: File does not exist\n")

o_file=open(args.output_file,mode="w")

if is_gz_file(args.input_file):
    fh = gzip.open(args.input_file, 'rt')
else:
    fh = open(args.input_file, 'r')

lines = []
for line in fh:
    lines.append(line.rstrip())
    if len(lines) == 4:
        read = process(lines)
        qual = read['quality']
        bases = read['sequence']
        
        newqual=""
        lastChar = ""
        secondLastChar = ""
        runlength = 1
        diRunlength = 0
        qualLength = len(qual)
        
        # indexes i goes thorough the base sequence positions,
        # with a lag of distanceToRun + runTheshold - 1 with respect to j.
        # When a run is detected in position j, Qrun is applied to position i
        # until i reaches position top, which is set to distanceToRun positions
        # ahead from the detected run.
        i = -distanceToRun - runTheshold + 1
        top = -1
        for j in range(qualLength + distanceToRun + runTheshold - 1):
            # If j lies within read limits
            if j < qualLength:
                if bases[j] == lastChar:
                    runlength += 1
                else:
                    runlength = 1
                if bases[j] == secondLastChar:
                    diRunlength += 1
                else:
                    diRunlength = 0
                secondLastChar = lastChar
                lastChar = bases[j]
                if runlength >= runTheshold or diRunlength + 2 >= 2*diRunTheshold:
                    top = j + distanceToRun
            # If i lies within read limits
            if i >= 0:
                qsc=ord(qual[i])-33
                if i <= top:
                    nqsc = Qrun(qsc)
                else:
                    nqsc = Qprimary(qsc)
                newqual = newqual + chr(nqsc+33)
            i += 1
        read['quality'] = newqual
        for key in ['name', 'sequence', 'optional', 'quality']:
            o_file.write(read[key] + "\n")
        lines = []
o_file.close()
