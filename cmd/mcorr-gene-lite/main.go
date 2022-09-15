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

// script written by Asher Preska Steinberg (apsteinberg@nyu.edu)
import (
	"fmt"
	"github.com/boltdb/bolt"
	"github.com/kussell-lab/biogo/seq"
	"github.com/kussell-lab/mcorr"
	"github.com/kussell-lab/ncbiftp/taxonomy"
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v2"
	"io"
	"log"
	"math/rand"
	"os"
	"runtime"
	"strings"
	"time"
)

// global variables.
func main() {
	app := kingpin.New("mcorr-gene-lite", "Calculate mutation correlation from bacterial single gene alignments stored as boltdb files.")
	app.Version("v20211015")

	dbFile := app.Arg("in", "Alignment stored as boltdb file.").Required().String()
	outPrefix := app.Arg("out", "Output prefix.").Required().String()

	maxl := app.Flag("max-corr-length", "Maximum distance of correlation (base pairs)").Default("300").Int()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	//numBoot := app.Flag("num-boot", "Number of bootstrapping on alleles").Default("0").Int()
	numDigesters := app.Flag("num-threads", "number of threads").Default("8").Int()
	showProgress := app.Flag("show-progress", "Show progress").Bool()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	//start timer

	start := time.Now()

	if *ncpu <= 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)

	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3

	// show progress bar?
	var bar *pb.ProgressBar
	if *showProgress {
		//max := getNumberOfAlignments(*alnFile)
		bar = pb.StartNew(maxCodonLen)
		defer bar.Finish()
	}

	synonymous := true
	codonPos := 3
	codonOffset := 0

	//initialize output csv
	outFile := *outPrefix + ".csv"
	initCsvOut(outFile)

	//open the database for usage, and determine the number of codons which are stored in there
	db, err := bolt.Open(*dbFile, 0600, &bolt.Options{Timeout: 1 * time.Second})
	if err != nil {
		log.Fatal(err)
	}
	defer db.Close()
	//count of codons
	numCodons := 0
	db.View(func(tx *bolt.Tx) error {
		// Assume bucket exists and has keys
		b := tx.Bucket([]byte("codons"))

		c := b.Cursor()
		for k, _ := c.First(); k != nil; k, _ = c.Next() {
			//fmt.Printf("key=%s, value=%s\n", k, v)
			numCodons++
		}

		return nil
	})
	fmt.Printf("Total number of codons: %d\n", numCodons)
	calcQsAll(db, codonOffset, codonPos-1, maxCodonLen, numCodons,
		codingTable, synonymous, outFile, *numDigesters, bar)

	//total time to complete ...
	duration := time.Since(start)
	fmt.Println("Time to calculate correlation profiles:", duration)
}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string
	Sequences []seq.Sequence
}

// calcSingleClade calculate correlation functions in a single cluster of sequence.
func calcSingleClade(alnChan chan Alignment, calculator Calculator) (corrResChan chan mcorr.CorrResults) {
	corrResChan = make(chan mcorr.CorrResults)
	done := make(chan bool)

	ncpu := runtime.GOMAXPROCS(0)
	for i := 0; i < ncpu; i++ {
		go func() {
			for aln := range alnChan {
				if len(aln.Sequences) > 1 {
					results := calculator.CalcP2(aln)
					corrResChan <- results
				}
			}
			done <- true
		}()
	}

	go func() {
		defer close(corrResChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()
	return
}

// setNumThreads sets number of CPU cores for using.
// if ncpu == 0, we will used all core avaible.
func setNumThreads(ncpu int) {
	if ncpu == 0 {
		ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(ncpu)
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(file string) (alnChan chan Alignment) {
	alnChan = make(chan Alignment)
	go func() {
		defer close(alnChan)

		c := readXMFA(file)
		for alignment := range c {
			alnID := strings.Split(alignment[0].Id, " ")[0]
			alnChan <- Alignment{ID: alnID, Sequences: alignment}
		}
	}()

	return
}

// bootstrapAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func bootstrapAlignments(file string, numBoot int) (alnChan chan Alignment) {
	alnChan = make(chan Alignment)
	go func() {
		defer close(alnChan)

		c := readXMFA(file)
		var alignment []seq.Sequence
		alignment = <-c
		numSeqs := len(alignment)
		//send the first alignment along the channel ...
		alnChan <- Alignment{ID: "all", Sequences: alignment}
		for i := 0; i < numBoot; i++ {
			bootstrap := bootstrapSeqs(alignment, numSeqs)
			id := fmt.Sprintf("boot_%d", i)
			alnChan <- Alignment{ID: id, Sequences: bootstrap}
		}
	}()

	return
}

func bootstrapSeqs(alignment []seq.Sequence, numSeqs int) (bootstrap []seq.Sequence) {
	for i := 0; i < numSeqs; i++ {
		//pick a random sequence from the original alignment
		j := rand.Intn(numSeqs)
		seq := alignment[j]
		//add it to our bootstrap
		bootstrap = append(bootstrap, seq)
	}
	return bootstrap
}

// getNumberOfAlignments return total number of alignments in a xmfa file.
func getNumberOfAlignments(file string) (count int) {
	c := readXMFA(file)
	for a := range c {
		if len(a) >= 2 {
			count++
		}
	}
	return
}

// readXMFA reads a xmfa format file and returns a channel of []seq.Sequence.
func readXMFA(file string) chan []seq.Sequence {
	c := make(chan []seq.Sequence)
	go func() {
		defer close(c)

		f := mustOpen(file)
		defer f.Close()

		rd := seq.NewXMFAReader(f)
		for {
			a, err := rd.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			if len(a) >= 2 {
				c <- a
			}
		}
	}()
	return c
}

// mustOpen is a helper function to open a file.
// and panic if error occurs.
func mustOpen(file string) (f *os.File) {
	var err error
	f, err = os.Open(file)
	if err != nil {
		panic(err)
	}
	return
}

// WriteResults writes correlation results of the original sample and the bootstraps to a .csv file
func WriteResults(corrResChan chan mcorr.CorrResults, outFile string) {

	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString("# l: the distance between two genomic positions\n")
	w.WriteString("# m: the mean value of correlatio profile\n")
	w.WriteString("# v: the variance of correlation profile\n")
	w.WriteString("# n: the total number of alignments used for calculation\n")
	w.WriteString("# t: the type of result: Ks is for d_sample, and P2 is for correlation profile\n")
	w.WriteString("# b: the bootstrap number (all means used all alignments).\n")

	w.WriteString("l,m,v,n,t,b\n")
	for corrRes := range corrResChan {
		//for _, bs := range bootstraps {
		//	results := bs.Results()
		results := corrRes.Results
		bootID := corrRes.ID
		//save d_sample ...
		var ds float64
		for _, res := range results {
			if res.Lag == 0 {
				res.Type = "Ks"
				w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean, res.Variance, res.N, res.Type, bootID))
				ds = res.Mean
			} else {
				w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n", res.Lag, res.Mean/ds, res.Variance, res.N, res.Type, bootID))
			}

		}
	}
}
