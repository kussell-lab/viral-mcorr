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
	"math"
	"os"
	"sync"
)

// calcQsAll calculates Qs at all positions for a given alignment and writes it to the outputcsv

func calcQsAll(seqMap map[string][]Codon, codonOffset, codonPosition, minCodonLen int, maxCodonLen int,
	codingTable *taxonomy.GeneticCode, synonymous bool, outFile string, numDigesters int) {
	//numDigesters := 20
	codonSequences := [][]Codon{}
	//for _, s := range aln.Sequences {
	//	codons := extractCodons(s, codonOffset)
	//	codonSequences = append(codonSequences, codons)
	//}

	for _, s := range seqMap {
		codonSequences = append(codonSequences, s)
	}
	done := make(chan struct{})

	lagChan := makeLagChan(done, minCodonLen, maxCodonLen, codonSequences)
	//start a fixed number of go routines
	c := make(chan map[pos_key]CorrResult)
	var wg sync.WaitGroup
	fmt.Printf("starting probability calculations ...\n")
	for i := 0; i < numDigesters; i++ {
		wg.Add(1)
		go calcQs(done, lagChan, c, codonSequences, synonymous, codingTable, codonPosition, i, &wg)
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

//makeLagChan returns a channel of lags
func makeLagChan(done <-chan struct{}, minCodonLen int, maxCodonLen int, codonSequences [][]Codon) <-chan int {
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

//calcQs calculates Qs for a given lag across all initial positions
func calcQs(done <-chan struct{}, lagChan <-chan int, resChan chan<- map[pos_key]CorrResult, codonSequences [][]Codon, synonymous bool,
	codingTable *taxonomy.GeneticCode, codonPosition int, id int, wg *sync.WaitGroup) {
	defer wg.Done()
	//fmt.Printf("Worker %d starting\n", id)
	//var results []CorrResult
	for l := range lagChan {
		//fmt.Printf("lag %d starting \n", l)
		corrResMap := mapCorrRes(codonSequences, synonymous, codingTable, codonPosition, l)
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

func mapCorrRes(codonSequences [][]Codon, synonymous bool,
	codingTable *taxonomy.GeneticCode, codonPosition int, l int) map[pos_key]CorrResult {
	corrResMap := make(map[pos_key]CorrResult)
	//loop through initial positions for a given lag
	for i := 0; i+l < len(codonSequences[0]); i++ {
		//collect P2
		totalP2 := 0.0
		totaln := 0
		//collect probability of difference at site a (Pa) and b (Pb)
		totalPa := 0.0
		totalPb := 0.0
		//totalna := 0
		//totalnb := 0
		codonPairs := []CodonPair{}
		j := i + l
		if _, ok := corrResMap[pos_key{i, j}]; !ok {
			//get codonPairs, which are two codons on the same sequence
			//separated by a distance i+l
			for _, cc := range codonSequences {
				if i+l < len(cc) {
					codonPairs = append(codonPairs, CodonPair{A: cc[i], B: cc[j]})
				}
			}
			//now split the codonPairs into different sets of codon pairs
			//corresponding to the amino acids they produce at site i and i+l
			//this is the multiCodonPair list
			multiCodonPairs := [][]CodonPair{}
			if synonymous {
				multiCodonPairs = synonymousSplit(codonPairs, codingTable)
			} else {
				multiCodonPairs = append(multiCodonPairs, codonPairs)
			}
			for _, codonPairs := range multiCodonPairs {
				if len(codonPairs) >= 2 {
					//
					nc, ncA, ncB := doubleCodonsAll(codonPairs, codonPosition)
					xy, n := nc.P11(0)
					xx, _ := ncA.P11(0)
					yy, _ := ncB.P11(0)
					totalP2 += xy
					totalPa += xx
					totalPb += yy
					totaln += n
					//totalPa += xx
					//totalna += nA
					//totalPb += yy
					//totalnb += nB

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

//doubleCodonsAll collects all codon pairs (where a codon pair is a codon at position i and i+l)
//from the set of sequences and adds them into a covariance matrix for the set of sequences
//and for positions i and i+l (referred to as a and b)
func doubleCodonsAll(codonPairs []CodonPair, codonPosition int) (c, ca, cb *mcorr.NuclCov) {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	c = mcorr.NewNuclCov(alphabet)
	ca = mcorr.NewNuclCov(alphabet)
	cb = mcorr.NewNuclCov(alphabet)
	for _, codonPair := range codonPairs {
		a := codonPair.A[codonPosition]
		b := codonPair.B[codonPosition]
		//add a value at site a and b (i and i+l)
		c.Add(a, b)
		ca.Add(a, a)
		cb.Add(b, b)
	}
	return c, ca, cb
}

//pos_key for corrResMap
type pos_key struct {
	pos_x int
	lag   int
}

//initCsvOut initializes the output csv
func initCsvOut(outFile string) {
	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	w.WriteString("# x: the initial position of the probability\n")
	w.WriteString("# l: the distance between two genomic positions\n")
	w.WriteString("# P11: joint probability of difference\n")
	w.WriteString("# P1a: probability of difference at site x\n")
	w.WriteString("# P1b: probability of difference at site x+l\n")
	w.WriteString("# n: the total number of seq pairs used for calculation\n")
	w.WriteString("# t: the type of result: ds is for d_sample, and Qs is for joint probability\n")
	w.WriteString("# g: the gene name.\n")
	w.WriteString("# pos: position of gene on the genome.\n")

	w.WriteString("x,l,P11,P1a,P1b,n,t,g,pos\n")
	w.Close()
}

//writeCsvOut writes results to the output csv
func writeCsvOut(outFile string, results map[pos_key]CorrResult) {
	f, err := os.OpenFile(outFile, os.O_APPEND|os.O_WRONLY, 0600)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	for _, res := range results {
		if res.Lag == 0 {
			res.Type = "ds"
		} else {
			res.Type = "Qs"
		}
		f.WriteString(fmt.Sprintf("%d,%d,%g,%g,%g,%d,%s,%s,%s\n",
			res.x_pos, res.Lag, res.P11, res.P1a, res.P1b, res.N,
			res.Type, "all CDS", "n/a"))
	}
}
