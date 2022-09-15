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
	"github.com/kussell-lab/mcorr"
	"github.com/kussell-lab/ncbiftp/taxonomy"
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
	c := make(chan mcorr.CorrResult)
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
	//end of pipeline; make a map then write to file

	resMap := make(map[int]mcorr.CorrResult)
	for res := range c {
		resMap[res.Lag] = res
	}

	WriteResults(resMap, maxCodonLen, outFile)

}

//calcQsMates calculates Qs for a given lag across all initial positions
func calcQsMates(done <-chan struct{}, lagChan <-chan int, resChan chan<- mcorr.CorrResult, cs1, cs2 []CodonSequence, synonymous bool,
	codingTable *taxonomy.GeneticCode, codonPosition int, id int, wg *sync.WaitGroup) {
	defer wg.Done()
	//fmt.Printf("Worker %d starting\n", id)
	//var results []CorrResult
	for l := range lagChan {
		//fmt.Printf("lag %d starting \n", l)
		corrRes := calcCorrResMates(cs1, cs2, synonymous, codingTable, codonPosition, l)
		select {
		case resChan <- corrRes:
			lag := 3 * l
			fmt.Printf("\rlag %d done", lag)
		case <-done:
			return
		}
	}
	//fmt.Printf("Worker %d done\n", id)
}

func calcCorrResMates(cs1, cs2 []CodonSequence, synonymous bool,
	codingTable *taxonomy.GeneticCode, codonPosition int, l int) (corrRes mcorr.CorrResult) {
	//corrResMap := make(map[pos_key]CorrResult)
	//loop through initial positions for a given lag
	//collect P2
	totalP2 := 0.0
	totaln := 0
	for i := 0; i+l < len(cs1[0]); i++ {

		j := i + l
		//get codonPairs, which are two codons on the same sequence
		//separated by a distance i+l
		cpList1 := extractCodonPairs(cs1, i, j, codingTable, synonymous)
		cpList2 := extractCodonPairs(cs2, i, j, codingTable, synonymous)
		for _, cp1 := range cpList1 {
			nc1 := doubleCodons(cp1, codonPosition)
			for _, cp2 := range cpList2 {
				nc2 := doubleCodons(cp2, codonPosition)
				if synonymous {
					aa1 := translateCodonPair(cp1[0], codingTable)
					aa2 := translateCodonPair(cp2[0], codingTable)
					if aa1 == aa2 {
						xy, n := nc1.MateP11(nc2, 0)
						totalP2 += xy
						totaln += n
					}
				} else {

				}
			}
		}
		if totaln > 0 {
			corrRes = mcorr.CorrResult{
				Lag:  l * 3,
				Mean: totalP2 / float64(totaln),
				N:    totaln,
				Type: "P2",
			}
		}
	}
	return corrRes
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
