package main

/* This script reads in a file and records the coordinates of unique kmers. */

import (
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"strings"
	"unicode/utf8"

	"github.com/pierrec/lz4"
)

type Args struct {
	K      int
	Input  string
	Output string
}

var args = Args{}
var linesPrinted int
var unique int
var kmerToIndex map[string]int
var kmerCounts []int

func init() {
	log.SetFlags(0)
	flag.IntVar(&args.K, "k", 0, "kmer size")
	flag.StringVar(&args.Input, "input", "", "Input file (required)")
	flag.StringVar(&args.Output, "output", "", "Output file (optional)")

	flag.Usage = func() {
		log.Println("usage: genomeuniqfrac [options]")
		flag.PrintDefaults()
	}
}

func complement(r rune) rune {
	if r == 'A' {
		return 'T'
	}
	if r == 'T' {
		return 'A'
	}
	if r == 'G' {
		return 'C'
	}
	if r == 'C' {
		return 'G'
	}
	if r == 'N' {
		return 'N'
	}
	if r == 'a' {
		return 't'
	}
	if r == 't' {
		return 'a'
	}
	if r == 'g' {
		return 'c'
	}
	if r == 'c' {
		return 'g'
	}
	log.Fatalf("Invalid base '%s'", r)
	return r
}

func revcomp(s string) string {
	n := utf8.RuneCountInString(s)
	o := make([]rune, n)
	i := n - 1
	for _, c := range s {
		o[i] = complement(c)
		i--
	}
	return string(o)
}

func main() {
	flag.Parse()

	if args.K == 0 || args.Input == "" {
		flag.Usage()
		os.Exit(1)
	}
	if args.Output == "" {
		log.Println("no output file specified, will print stats only")
	} else {
		if !strings.HasSuffix(args.Output, ".lz4") {
			log.Fatal("-output file must end on '.lz4'")
		}
	}

	log.Println("reading input")
	data, err := ioutil.ReadFile(args.Input)
	if err != nil {
		log.Fatalf("failed read read file: %v", err)
	}
	log.Println("bytes:", len(data))

	log.Println("converting input")
	sequence := strings.ToUpper(string(data))

	log.Println("reverse complementing")
	sequenceRC := revcomp(sequence)

	log.Println("processing input")
	kmerToIndex = make(map[string]int)
	kmerCounts = make([]int, 0)
	for i := 0; i+args.K <= len(sequence); i++ {
		kmer := sequence[i:(i + args.K)]
		idx, found := kmerToIndex[kmer]
		if !found {
			idx = len(kmerToIndex)
			kmerToIndex[kmer] = idx
			kmerCounts = append(kmerCounts, 0)
		}
		kmerCounts[idx] += 1
		// Reverse complement
		kmerRC := sequenceRC[i:(i + args.K)]
		idx, found = kmerToIndex[kmerRC]
		if !found {
			idx = len(kmerToIndex)
			kmerToIndex[kmerRC] = idx
			kmerCounts = append(kmerCounts, 0)
		}
		// Sequences that are their own reverse complement would get a count of
		// two, even if that sequence only happens one place in the genome. Of
		// course the orientation is ambiguous, but at least the location is
		// distinct. In this case we don't double count.
		if kmer != kmerRC {
			kmerCounts[idx] += 1
		}
	}

	log.Println("finding unique kmers")
	var fp *os.File
	var zw *lz4.Writer
	if args.Output != "" {
		var err error
		fp, err = os.Create(args.Output)
		if err != nil {
			log.Fatalf("failed to open output file %s: %v", args.Output, err)
		}
		defer fp.Close()
		zw = lz4.NewWriter(fp)
		defer zw.Close()
	}
	for i := 0; i+args.K <= len(sequence); i++ {
		outputKmer(i, sequence[i:(i+args.K)], zw)
		// TODO: I think it's unnecessary to output the index of the reverse
		// complement sequence at this offset because that could be found by
		// looking up the index of the reverse complement of a kmer when using
		// the index.
		outputKmer(-i, sequenceRC[i:(i+args.K)], zw)
	}
	if zw != nil {
		err = zw.Flush()
		if err != nil {
			log.Fatalf("Failed to flush output: %v", err)
		}
	}

	log.Println("kmers:", len(kmerCounts))
	log.Println("unique:", unique)
	log.Println("lines:", linesPrinted)
}

func outputKmer(coord int, kmer string, w io.Writer) {
	idx := kmerToIndex[kmer]
	if kmerCounts[idx] == 1 {
		unique += 1
		if w != nil {
			// Only output kmers that are ATGC
			okay := true
			for j := 0; j < args.K; j++ {
				b := kmer[j]
				if b != 'A' && b != 'T' && b != 'C' && b != 'G' {
					okay = false
					break
				}
			}
			if okay {
				fmt.Fprintf(w, "%d\t%s\n", coord, kmer)
				linesPrinted += 1
			}
		}
	}
}

// END
