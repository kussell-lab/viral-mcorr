// Copyright 2022 Asher Preska Steinberg
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
	"fmt"
	"github.com/kussell-lab/biogo/seq"
	"github.com/kussell-lab/mcorr"
	"github.com/kussell-lab/ncbiftp/taxonomy"
	"gopkg.in/cheggaaa/pb.v2"
	"os"
	"strings"
	"sync"
)

// calcQsAll calculates Qs at all positions for a given alignment and writes it to the outputcsv

func calcKsAll(seqMap map[string][]Codon, seqpairs [][]string, codonOffset, codonPosition int,
	codingTable *taxonomy.GeneticCode, synonymous bool, outFile string, numDigesters int, bar *pb.ProgressBar) {
	//numDigesters := 20
	codonSequences := [][]Codon{}

	for _, s := range seqMap {
		codonSequences = append(codonSequences, s)
	}
	done := make(chan struct{})

	pairChan := makeSeqPairChan(done, seqMap, seqpairs)
	//start a fixed number of go routines
	c := make(chan map[string]mcorr.CorrResult)
	var wg sync.WaitGroup
	for i := 0; i < numDigesters; i++ {
		wg.Add(1)
		go calcKsPair(done, pairChan, c, synonymous, codingTable, codonPosition, i, bar, &wg)
	}

	go func() {
		wg.Wait()
		close(c)
	}()
	//end of pipeline; write files ...
	for res := range c {
		writeCsvOut(outFile, res)
	}

}

//writeCsvOut writes results to the output csv
func writeCsvOut(outFile string, results map[string]mcorr.CorrResult) {
	f, err := os.OpenFile(outFile, os.O_APPEND|os.O_WRONLY, 0600)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	for pairID, res := range results {
		f.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n",
			res.Lag, res.Mean, res.Variance, res.N, res.Type, pairID))
	}
}

//calcQs calculates Qs for a given lag across all initial positions
func calcKsPair(done <-chan struct{}, pairChan <-chan SeqPair, resChan chan<- map[string]mcorr.CorrResult,
	synonymous bool, codingTable *taxonomy.GeneticCode, codonPosition int, id int, bar *pb.ProgressBar, wg *sync.WaitGroup) {
	defer wg.Done()
	//fmt.Printf("Worker %d starting\n", id)
	//var results []CorrResult
	for seqPair := range pairChan {
		//fmt.Printf("lag %d starting \n", l)
		KsResMap := mapKsRes(seqPair, synonymous, codingTable, codonPosition)
		select {
		case resChan <- KsResMap:
			//lag := 3 * l
			//fmt.Printf("\rlag %d done", lag)
			if bar != nil {
				bar.Add(1)
			}
		case <-done:
			return
		}
	}
	//fmt.Printf("Worker %d done\n", id)
}

//mapKsRes calculates pairwise synonymous distance for a sequence pair
func mapKsRes(seqPair SeqPair, synonymous bool,
	codingTable *taxonomy.GeneticCode, codonPos int) map[string]mcorr.CorrResult {
	ksMap := make(map[string]mcorr.CorrResult)
	pairID := seqPair.genomeName1 + "_vs_" + seqPair.genomeName2
	seq1 := seqPair.genome1
	seq2 := seqPair.genome2
	d := 0.0
	t := 0
	for k := 0; k < len(seq1) && k < len(seq2); k++ {
		c1 := seq1[k]
		c2 := seq2[k]
		a1, found1 := codingTable.Table[string(c1)]
		a2, found2 := codingTable.Table[string(c2)]
		if found1 && found2 && a1 == a2 {
			b1 := seq1[k]
			b2 := seq2[k]

			good := true
			if synonymous {
				d1, found1 := codingTable.Table[string(c1)]
				d2, found2 := codingTable.Table[string(c2)]
				if found1 && found2 && d1 == d2 {
					good = true
				} else {
					good = false
				}
			}
			if good {
				var codonPositions []int
				if codonPos < 0 || codonPos > 2 {
					codonPositions = []int{0, 1, 2}
				} else {
					codonPositions = append(codonPositions, codonPos)
				}
				for _, codonP := range codonPositions {
					if c1[codonP] != c2[codonP] {
						if b1[codonP] != b2[codonP] {
							d++
						}
					}
					t++
				}
			}
		}
	}
	//collect the goods only if we've actually computed something  ...
	if t > 0 {
		ks := mcorr.CorrResult{}
		ks.Lag = 0
		ks.Mean = d / float64(t)
		ks.N = t
		ks.Type = "Ks"
		ksMap[pairID] = ks
	}

	return ksMap
}

// Codon is a byte list of length 3
type Codon []byte

// extractCodons return a list of codons from a DNA sequence.
func extractCodons(s seq.Sequence, offset int) (codons []Codon) {
	for i := offset; i+3 <= len(s.Seq); i += 3 {
		c := s.Seq[i:(i + 3)]
		codons = append(codons, c)
	}
	return
}

//makeSeqPairChan returns a channel of sequence pairs
func makeSeqPairChan(done <-chan struct{}, seqMap map[string][]Codon, seqpairs [][]string) <-chan SeqPair {
	SeqPairChan := make(chan SeqPair)
	go func() {
		defer close(SeqPairChan)
		for _, seqpair := range seqpairs {
			seqName1 := seqpair[0]
			seqName2 := seqpair[1]
			seq1 := seqMap[seqName1]
			seq2 := seqMap[seqName2]
			pairSeqs := SeqPair{seqName1, seq1,
				seqName2, seq2}
			select {
			case SeqPairChan <- pairSeqs:
			case <-done:
			}

		}
	}()
	return SeqPairChan
}

//SeqPair pair of sequences to be analyzed
type SeqPair struct {
	genomeName1 string
	genome1     []Codon
	genomeName2 string
	genome2     []Codon
}

func getNames(s string) (geneName, genomeName string) {
	terms := strings.Split(s, " ")
	//this is for the helicobacter test files
	//geneName = terms[0]
	//genomeName = terms[1]
	//this is the genomeName for the MSA files assembled from ReferenceAlignmentGenerator
	geneName = terms[0]
	genomeName = terms[2]
	return
}
