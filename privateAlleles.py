#!/usr/bin/python

from comline import ComLine
from popmap import Popmap
from vcf import VCF

import sys



def main():
	input = ComLine(sys.argv[1:])
	pops = Popmap(input.args.popmap)
	vcf = VCF(input.args.vcf, pops)

main()

raise SystemExit
