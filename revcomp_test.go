package main

import "testing"

func BenchmarkReverseComplement(b *testing.B) {
	for i := 0; i < b.N; i++ {
		revcomp("CTACGAGAGGACCTATTGTCTGCAGACGGGTAACGCTAAGCATTCCCTCATGATATTGAG")
	}
}
