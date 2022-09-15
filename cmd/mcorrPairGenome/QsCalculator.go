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

//calcQsAll calculates Qs at all positions for a given alignment and writes it to the output csv

func calcQsAll(seqMap map[string][]Codon, seqpairs [][]string, codonOffset, codonPosition int, maxCodonLen int,
	codingTable *taxonomy.GeneticCode, synonymous bool, outFile string, numDigesters int, bar *pb.ProgressBar) {
	//numDigesters := 20
	codonSequences := [][]Codon{}

	for _, s := range seqMap {
		codonSequences = append(codonSequences, s)
	}
	done := make(chan struct{})

	pairChan := makeSeqPairChan(done, seqMap, seqpairs)
	//start a fixed number of go routines
	c := make(chan map[string]mcorr.CorrResults)
	var wg sync.WaitGroup
	for i := 0; i < numDigesters; i++ {
		wg.Add(1)
		go calcQsPair(done, pairChan, c, synonymous, codingTable, codonPosition,
			maxCodonLen, i, bar, &wg)
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
func writeCsvOut(outFile string, corrResMap map[string]mcorr.CorrResults) {
	f, err := os.OpenFile(outFile, os.O_APPEND|os.O_WRONLY, 0600)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	for pairID, corrRes := range corrResMap {
		results := corrRes.Results
		//save d_sample ...
		var ds float64
	Loop:
		for _, res := range results {
			switch lag := res.Lag; {
			case lag == 0 && res.Mean == 0:
				//stop writing this result if ds = 0
				res.Type = "Ks"
				f.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, pairID))
				break Loop
			case lag == 0:
				res.Type = "Ks"
				f.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, pairID))
				ds = res.Mean
			default:
				f.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean/ds, res.Variance, res.N, res.Type, pairID))
			}
		}
		//save d_sample ...
		//var ds float64
		//for _, res := range results {
		//	if res.Lag == 0 {
		//		res.Type = "Ks"
		//		f.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, pairID))
		//		ds = res.Mean
		//	} else {
		//		f.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean/ds, res.Variance, res.N, res.Type, pairID))
		//	}
		//}
	}
}

//calcQsPair calculates Qs for a given pair across all positions
func calcQsPair(done <-chan struct{}, pairChan <-chan SeqPair, resChan chan<- map[string]mcorr.CorrResults,
	synonymous bool, codingTable *taxonomy.GeneticCode, codonPosition int,
	maxCodonLen int, id int, bar *pb.ProgressBar, wg *sync.WaitGroup) {
	defer wg.Done()
	//fmt.Printf("Worker %d starting\n", id)
	//var results []CorrResult
	for seqPair := range pairChan {
		//fmt.Printf("lag %d starting \n", l)
		QsResMap := mapQsRes(seqPair, synonymous, codingTable, codonPosition, maxCodonLen)
		select {
		case resChan <- QsResMap:
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

//mapQsRes calculates correlation profiles for a sequence pair
func mapQsRes(seqPair SeqPair, synonymous bool, codingTable *taxonomy.GeneticCode,
	codonPos int, maxCodonLen int) map[string]mcorr.CorrResults {
	QsMap := make(map[string]mcorr.CorrResults)

	//name pairs
	var pairID string
	if seqPair.genomeName1 > seqPair.genomeName2 {
		pairID = seqPair.genomeName2 + "_vs_" + seqPair.genomeName1
	} else {
		pairID = seqPair.genomeName1 + "_vs_" + seqPair.genomeName2
	}
	//pairID := seqPair.genomeName1 + "_vs_" + seqPair.genomeName2
	seq1 := seqPair.genome1
	seq2 := seqPair.genome2
	//initiate result collection
	crRes := mcorr.CorrResults{ID: pairID}
	for l := 0; l < maxCodonLen; l++ {
		d := 0.0
		t := 0
		for k := 0; k < len(seq1)-l; k++ {
			c1 := seq1[k]
			c2 := seq2[k]
			a1, found1 := codingTable.Table[string(c1)]
			a2, found2 := codingTable.Table[string(c2)]
			if found1 && found2 && a1 == a2 {
				b1 := seq1[k+l]
				b2 := seq2[k+l]

				good := true
				if synonymous {
					d1, found1 := codingTable.Table[string(b1)]
					d2, found2 := codingTable.Table[string(b2)]
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
		cr := mcorr.CorrResult{}
		cr.Lag = l * 3
		cr.Mean = d / float64(t)
		cr.N = t
		if l == 0 {
			cr.Type = "Ks"
		} else {
			cr.Type = "P2"
		}

		crRes.Results = append(crRes.Results, cr)
	}
	QsMap[pairID] = crRes

	return QsMap
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
