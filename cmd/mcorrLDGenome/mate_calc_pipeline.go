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
	"github.com/kussell-lab/ncbiftp/taxonomy"
	"math"
	"sync"
)

// calcQsAll calculates Qs at all positions for a given alignment and writes it to the outputcsv

func calcQsMatesAll(seqMap1, seqMap2 map[string][]Codon, codonOffset, codonPosition, minCodonLen int, maxCodonLen int,
	codingTable *taxonomy.GeneticCode, synonymous bool, outFile string, numDigesters int) {
	//get our two lists of codon sequences
	var cs1 []CodonSequence
	var cs2 []CodonSequence
	//get our codon sequences
	for _, s := range seqMap1 {
		cs1 = append(cs1, s)
	}

	for _, s := range seqMap2 {
		cs2 = append(cs2, s)
	}
	done := make(chan struct{})

	lagChan := startLagChan(done, minCodonLen, maxCodonLen, cs1)
	//start a fixed number of go routines
	c := make(chan map[pos_key]CorrResult)
	var wg sync.WaitGroup
	fmt.Printf("starting probability calculations ...\n")
	for i := 0; i < numDigesters; i++ {
		wg.Add(1)
		go calcQsMates(done, lagChan, c, cs1, cs2, synonymous, codingTable, codonPosition, i, &wg)
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

//calcQsMates calculates Qs for a given lag across all initial positions
func calcQsMates(done <-chan struct{}, lagChan <-chan int, resChan chan<- map[pos_key]CorrResult, cs1, cs2 []CodonSequence, synonymous bool,
	codingTable *taxonomy.GeneticCode, codonPosition int, id int, wg *sync.WaitGroup) {
	defer wg.Done()
	//fmt.Printf("Worker %d starting\n", id)
	//var results []CorrResult
	for l := range lagChan {
		//fmt.Printf("lag %d starting \n", l)
		corrResMap := mapCorrResMates(cs1, cs2, synonymous, codingTable, codonPosition, l)
		select {
		case resChan <- corrResMap:
			lag := 3 * l
			fmt.Printf("\rlag %d done", lag)
		case <-done:
			return
		}
	}
	//fmt.Printf("Worker %d done\n", id)
}

func mapCorrResMates(cs1, cs2 []CodonSequence, synonymous bool,
	codingTable *taxonomy.GeneticCode, codonPosition int, l int) map[pos_key]CorrResult {
	corrResMap := make(map[pos_key]CorrResult)
	//loop through initial positions for a given lag
	for i := 0; i+l < len(cs1[0]); i++ {
		//collect P2
		totalP2 := 0.0
		totaln := 0
		//collect probability of difference at site a (Pa) and b (Pb)
		totalPa := 0.0
		totalPb := 0.0
		j := i + l
		if _, ok := corrResMap[pos_key{i, j}]; !ok {
			//get codonPairs, which are two codons on the same sequence
			//separated by a distance i+l
			cpList1 := extractCodonPairs(cs1, i, j, codingTable, synonymous)
			cpList2 := extractCodonPairs(cs2, i, j, codingTable, synonymous)
			for _, cp1 := range cpList1 {
				nc1, nc1A, nc1B := doubleCodonsAll(cp1, codonPosition)
				for _, cp2 := range cpList2 {
					nc2, nc2A, nc2B := doubleCodonsAll(cp2, codonPosition)
					if synonymous {
						aa1 := translateCodonPair(cp1[0], codingTable)
						aa2 := translateCodonPair(cp2[0], codingTable)
						if aa1 == aa2 {
							xy, n := nc1.MateP11(nc2, 0)
							xx, _ := nc1A.MateP11(nc2A, 0)
							yy, _ := nc1B.MateP11(nc2B, 0)
							totalP2 += xy
							totalPa += xx
							totalPb += yy
							totaln += n
						}
					} else {

					}
				}
			}
			if totaln > 0 {
				res1 := CorrResult{
					x_pos:    i * 3,
					Lag:      l * 3,
					P11:      totalP2 / float64(totaln),
					totalP11: totalP2,
					P1a:      totalPa / float64(totaln),
					P1b:      totalPb / float64(totaln),
					N:        totaln,
					Type:     "P2",
				}
				//results = append(results, res1)
				corrResMap[pos_key{i, j}] = res1
			} else {
				//fill if there were no sequence pairs for the position
				//totaln = 0
				res1 := CorrResult{
					x_pos:    i * 3,
					Lag:      l * 3,
					P11:      math.NaN(),
					totalP11: totalP2,
					P1a:      math.NaN(),
					P1b:      math.NaN(),
					N:        totaln,
					Type:     "P2",
				}
				corrResMap[pos_key{i, j}] = res1
			}
		}
	}
	return corrResMap
}

func extractCodonPairs(codonSequences []CodonSequence, i, j int,
	codingTable *taxonomy.GeneticCode, synonymous bool) [][]CodonPair {
	codonPairs := []CodonPair{}
	for _, cc := range codonSequences {
		if i < len(cc) && j < len(cc) {
			pair := CodonPair{A: cc[i], B: cc[j]}
			codonPairs = append(codonPairs, pair)
		}
	}

	if synonymous {
		return synonymousSplit(codonPairs, codingTable)
	}

	return [][]CodonPair{codonPairs}
}

func translateCodonPair(cp CodonPair, codingTable *taxonomy.GeneticCode) string {
	a := codingTable.Table[string(cp.A)]
	b := codingTable.Table[string(cp.B)]
	return string([]byte{a, b})
}

//startLagChan returns a channel of lags
func startLagChan(done <-chan struct{}, minCodonLen int, maxCodonLen int, codonSequences []CodonSequence) <-chan int {
	lagChan := make(chan int)
	var maxlag int
	var minlag int
	if maxCodonLen == 0 {
		maxlag = len(codonSequences[0])
		minlag = 0
	} else {
		maxlag = maxCodonLen
		minlag = minCodonLen
	}
	go func() {
		defer close(lagChan)
		for l := minlag; l < maxlag; l++ {
			select {
			case lagChan <- l:
			case <-done:
			}
		}
	}()
	return lagChan
}
