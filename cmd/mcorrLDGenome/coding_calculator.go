// Copyright 2022 Asher Preska Steinberg
//
// Initially modified from coding_calculator.go of mcorr
// https://github.com/kussell-lab/mcorr
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package main

import (
	"github.com/kussell-lab/biogo/seq"
	"github.com/kussell-lab/ncbiftp/taxonomy"
	"sync"
)

// Calculator define a interface for calculating correlations.
type Calculator interface {
	CalcP2(a Alignment, others ...Alignment) (corrResults CorrResults)
}

// CodingCalculator for calculating coding sequences.
type CodingCalculator struct {
	CodingTable   *taxonomy.GeneticCode
	MaxCodonLen   int
	CodonOffset   int
	CodonPosition int
	Synonymous    bool
}

// NewCodingCalculator return a CodingCalculator
func NewCodingCalculator(codingTable *taxonomy.GeneticCode, maxCodonLen, codonOffset int, codonPosition int, synonymous bool) *CodingCalculator {
	return &CodingCalculator{
		CodingTable:   codingTable,
		MaxCodonLen:   maxCodonLen,
		CodonOffset:   codonOffset,
		CodonPosition: codonPosition,
		Synonymous:    synonymous,
	}
}

// CalcP2 calculate P2
func (cc *CodingCalculator) CalcP2(a Alignment, others ...Alignment) CorrResults {
	results := calcP2Coding(a, cc.CodonOffset, cc.CodonPosition, cc.MaxCodonLen, cc.CodingTable, cc.Synonymous)
	return CorrResults{ID: a.ID, Results: results}
}

func calcP2Coding(aln Alignment, codonOffset, codonPosition, maxCodonLen int, codingTable *taxonomy.GeneticCode, synonymous bool) (results []CorrResult) {
	codonSequences := [][]Codon{}
	for _, s := range aln.Sequences {
		codons := extractCodons(s, codonOffset)
		codonSequences = append(codonSequences, codons)
	}

	var wg sync.WaitGroup
	//tell wg how many threads are running concurrently
	wg.Add(len(codonSequences[0]))
	//mutex lock for safe access to results
	var mutex = &sync.Mutex{}
	//was previously len(maxCodonLen)
	//so now the maximum separation between the two sites can be the length of gene
	for l := 0; l < len(codonSequences[0]); l++ {
		//spawn threads for loop iterations
		go func(l int) {
			defer wg.Done()
			for i := 0; i+l < len(codonSequences[0]); i++ {
				totalP2 := 0.0
				totaln := 0
				codonPairs := []CodonPair{}
				j := i + l
				for _, cc := range codonSequences {
					if i+l < len(cc) {
						codonPairs = append(codonPairs, CodonPair{A: cc[i], B: cc[j]})
					}
				}

				multiCodonPairs := [][]CodonPair{}
				if synonymous {
					multiCodonPairs = synonymousSplit(codonPairs, codingTable)
				} else {
					multiCodonPairs = append(multiCodonPairs, codonPairs)
				}
				for _, codonPairs := range multiCodonPairs {
					if len(codonPairs) >= 2 {
						nc := doubleCodons(codonPairs, codonPosition)
						xy, n := nc.P11(0)
						totalP2 += xy
						totaln += n

					}
				}
				if totaln > 0 {
					res1 := CorrResult{
						x_pos:    i * 3,
						Lag:      l * 3,
						P11:      totalP2 / float64(totaln),
						totalP11: totalP2,
						N:        totaln,
						Type:     "P2",
					}
					mutex.Lock()
					results = append(results, res1)
					mutex.Unlock()
				}
			}
		}(l)
	}

	wg.Wait()

	return
}

// CorrResult stores a correlation result.
type CorrResult struct {
	x_pos    int
	Lag      int
	P11      float64
	totalP11 float64 //will remove later
	N        int
	P1a      float64 //probability of difference at site a
	Na       int     //number of site As
	P1b      float64 //probability of difference at site
	Nb       int     // number of site Bs
	Type     string
}

// CorrResults stores a list of CorrResult with an gene ID.
type CorrResults struct {
	ID      string
	Results []CorrResult
}

//doubleCodons collects all codon pairs (where a codon pair is a codon at position i and i+l)
//from the set of sequences and adds them into a covariance matrix for the set of sequences
func doubleCodons(codonPairs []CodonPair, codonPosition int) *NuclCov {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	c := NewNuclCov(alphabet)
	for _, codonPair := range codonPairs {
		a := codonPair.A[codonPosition]
		b := codonPair.B[codonPosition]
		//add a value at site a and b (i and i+l)
		c.Add(a, b)
	}
	return c
}

// Codon is a byte list of length 3
type Codon []byte

// CodonSequence is a sequence of codons.
type CodonSequence []Codon

// CodonPair is a pair of Codons.
type CodonPair struct {
	A, B Codon
}

// extractCodons return a list of codons from a DNA sequence.
func extractCodons(s seq.Sequence, offset int) (codons []Codon) {
	for i := offset; i+3 <= len(s.Seq); i += 3 {
		c := s.Seq[i:(i + 3)]
		codons = append(codons, c)
	}
	return
}

// synonymousSplit split a list of codon pairs into multiple
// synonymous pairs. You check which AAs each codon in the two site pair produces,
// then add it to the multiCodonPair list at an index corresponding to the AAs
//the two sites produce
func synonymousSplit(codonPairs []CodonPair, codingTable *taxonomy.GeneticCode) (multiCodonPairs [][]CodonPair) {
	aaList := []string{}
	for _, codonPair := range codonPairs {
		// check gap.
		containsGap := false
		for _, codon := range []Codon{codonPair.A, codonPair.B} {
			for i := 0; i < 3; i++ {
				if codon[i] == '-' || codon[i] == 'N' {
					containsGap = true
					break
				}
			}
		}
		if containsGap {
			continue
		}

		codonA := string(codonPair.A)
		codonB := string(codonPair.B)
		a := codingTable.Table[codonA]
		b := codingTable.Table[codonB]
		ab := string([]byte{a, b})
		index := -1
		for i := 0; i < len(aaList); i++ {
			if aaList[i] == ab {
				index = i
			}
		}
		if index == -1 {
			index = len(aaList)
			aaList = append(aaList, ab)
			multiCodonPairs = append(multiCodonPairs, []CodonPair{})
		}

		multiCodonPairs[index] = append(multiCodonPairs[index], codonPair)
	}

	return
}
