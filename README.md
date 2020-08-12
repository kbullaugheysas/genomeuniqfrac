# genomeuniqfrac

## Overview

Utility to build an index of kmers that are unique in the genome. This is useful for programs that only trust kmers that exist in only one copy anywhere in the genome. It's rather inefficient memory/speed wise.

## Usage

usage: genomeuniqfrac [options]
  -input string
    	Input file (required)
  -k int
    	kmer size
  -output string
    	Output file (optional)
