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
	"github.com/boltdb/bolt"
	"github.com/kussell-lab/mcorr"
	"github.com/kussell-lab/ncbiftp/taxonomy"
	"gopkg.in/cheggaaa/pb.v2"
	"os"
	"strconv"
	"sync"
)

//calcQsAll calculates all lags in a multithreaded fashion, good for large datasets ...
func calcQsAll(db *bolt.DB, codonOffset, codonPosition, maxCodonLen, numCodons int,
	codingTable *taxonomy.GeneticCode, synonymous bool, outFile string, numDigesters int, bar *pb.ProgressBar) {
	done := make(chan struct{})
	lagChan := makeLagChan(done, maxCodonLen)

	c := make(chan mcorr.CorrResult)
	var wg sync.WaitGroup
	for i := 0; i < numDigesters; i++ {
		wg.Add(1)
		go calcQs(done, lagChan, c, db, synonymous, codingTable,
			codonPosition, numCodons, i, bar, &wg)
	}

	go func() {
		wg.Wait()
		close(c)
	}()

	//end of pipeline; dump into a map temporarily, then write to a csv file....
	resMap := make(map[int]mcorr.CorrResult)
	for res := range c {
		resMap[res.Lag] = res
	}
	//now write to a csv file ....
	writeCsvOut(outFile, resMap, maxCodonLen)

}

//makeLagChan returns a channel of lags
func makeLagChan(done <-chan struct{}, maxCodonLen int) <-chan int {
	lagChan := make(chan int)
	go func() {
		defer close(lagChan)
		for l := 0; l < maxCodonLen; l++ {
			select {
			case lagChan <- l:
			case <-done:
			}
		}
	}()
	return lagChan
}

//calcQs calculates Qs for a given lag across all positions
func calcQs(done <-chan struct{}, lagChan <-chan int, resChan chan<- mcorr.CorrResult,
	db *bolt.DB, synonymous bool, codingTable *taxonomy.GeneticCode, codonPosition, numCodons int,
	id int, bar *pb.ProgressBar, wg *sync.WaitGroup) {
	defer wg.Done()
	for l := range lagChan {
		corrRes := calcCorrRes(db, codonPosition, l, numCodons, codingTable, synonymous)
		select {
		case resChan <- corrRes:
			if bar != nil {
				bar.Add(1)
			}
		case <-done:
			return
		}
	}
}

func calcCorrRes(db *bolt.DB, codonPosition, lag int, numCodons int,
	codingTable *taxonomy.GeneticCode, synonymous bool) (corrRes mcorr.CorrResult) {
	//l is lag
	l := lag
	totalP2 := 0.0
	totaln := 0
	for i := 0; i+l < numCodons; i++ {
		codonPairs := []CodonPair{}
		j := i + l
		var codonsA []Codon
		var codonsB []Codon
		codonsA = getCodons(db, i)
		codonsB = getCodons(db, j)

		for p, codonA := range codonsA {
			codonPairs = append(codonPairs, CodonPair{A: codonA, B: codonsB[p]})
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
				nc := doubleCodons(codonPairs, codonPosition)
				xy, n := nc.P11(0)
				totalP2 += xy
				totaln += n

			}
		}
	}

	switch totaln > 0 {
	case l != 0:
		corrRes = mcorr.CorrResult{
			Lag:  l * 3,
			Mean: totalP2 / float64(totaln),
			N:    totaln,
			Type: "P2",
		}
	case l == 0:
		corrRes = mcorr.CorrResult{
			Lag:  l * 3,
			Mean: totalP2 / float64(totaln),
			N:    totaln,
			Type: "Ks",
		}
	}
	//if totaln > 0 {
	//	corrRes = mcorr.CorrResult{
	//		Lag:  l * 3,
	//		Mean: totalP2 / float64(totaln),
	//		N:    totaln,
	//		Type: "P2",
	//	}
	//}

	return corrRes
}

func getCodons(db *bolt.DB, pos int) (codons []Codon) {
	fn := func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte("codons"))
		codonPos := strconv.Itoa(pos)
		v := b.Get([]byte(codonPos))
		codons = BytestoCodons(v, 0)

		return nil
	}

	db.View(fn)

	return
}

// BytestoCodons convert bytes back to codons
func BytestoCodons(s []byte, offset int) (codons []Codon) {
	for i := offset; i+3 <= len(s); i += 3 {
		c := s[i:(i + 3)]
		codons = append(codons, c)
	}
	return
}

//initCsvOut initializes the output csv
func initCsvOut(outFile string) {
	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	w.WriteString("# l: the distance between two genomic positions\n")
	w.WriteString("# m: the mean value of the correlation profile\n")
	w.WriteString("# v: the variance of the correlation profile\n")
	w.WriteString("# n: the total number of seq pairs used for calculation\n")
	w.WriteString("# t: the type of result: Ks is for d_sample and P2 is for correlation profile\n")
	w.WriteString("# b: the bootstrap number (all means used all alignments).\n")

	w.WriteString("l,m,v,n,t,b\n")
	w.Close()
}

//writeCsvOut writes results to the output csv
func writeCsvOut(outFile string, results map[int]mcorr.CorrResult, maxCodonLen int) {
	f, err := os.OpenFile(outFile, os.O_APPEND|os.O_WRONLY, 0600)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	//get ds so we can divide through to get conditional probability
	dsRes := results[0]
	ds := dsRes.Mean
	var Pab float64
	i := 0
	for i < maxCodonLen*3 {
		res := results[i]
		if i != 0 {
			Pab = res.Mean / ds
		} else {
			Pab = res.Mean
		}

		f.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n",
			res.Lag, Pab, res.Variance, res.N, res.Type, "all"))
		i = i + 3
	}
}
