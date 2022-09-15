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
	"github.com/kussell-lab/biogo/seq"
	"github.com/kussell-lab/ncbiftp/taxonomy"
	"gopkg.in/alecthomas/kingpin.v2"
	"io"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

// global variables.
func main() {
	app := kingpin.New("mcorrLDGenome", "Calculate mutation correlation across CDS regions from an XMFA file.")
	app.Version("v20210513")

	alnFile := app.Arg("aln", "Alignment file in XMFA format.").Required().String()
	outPrefix := app.Arg("out", "Output prefix.").Required().String()
	minl := app.Flag("min-corr-length", "min distance of correlation (base pairs)").Default("0").Int()
	maxl := app.Flag("max-corr-length", "Max distance of correlation (base pairs)").Default("0").Int()
	mateAln := app.Flag("mate-aln", "Second alignment for calculating between clades").Default("").String()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	mates := app.Flag("between-clades", "just calculate correlations between pairs from different xmfa files").Default("false").Bool()
	numDigesters := app.Flag("num-threads", "number of threads").Default("50").Int()
	//showProgress := app.Flag("show-progress", "Show progress").Bool()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	if *ncpu <= 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)

	//timer
	start := time.Now()

	// prepare calculator.
	//var calculator Calculator
	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3
	minCodonLen := *minl / 3

	synonymous := true
	codonPos := 3
	codonOffset := 0

	//make one giant alignment of all CDS regions ...
	// make a list of starting positions so we know how to order things and
	// dump the alignments into a map with starting positions as keys
	var startSlice []int
	c := readAlignments(*alnFile)
	for a := range c {
		startPos := getStartPos(a)
		startSlice = append(startSlice, startPos)
	}
	//sort the slice numerically
	sort.Ints(startSlice)

	// now go through and make a map of codon sequences
	// where we add on codons to the end of each strain sequence
	fmt.Print("fetching CDS regions\n")
	//the keys are the strain names and the values are the codon sequences for each
	//seqMap := make(map[string][]Codon)
	var seqMap map[string][]Codon
	var seqMap1 map[string][]Codon
	if *mateAln != "" {
		if *mates {
			seqMap = makeSeqMap(startSlice, *alnFile, codonOffset)
			seqMap1 = makeSeqMap(startSlice, *mateAln, codonOffset)
		} else {
			seqMap = combinedSeqMap(startSlice, *alnFile, *mateAln, codonOffset)
		}
	} else {
		seqMap = makeSeqMap(startSlice, *alnFile, codonOffset)
	}

	numSeqs := len(seqMap)
	//get total number of codons
	var numCodons int
	for _, seq := range seqMap {
		numCodons = len(seq)
		break
	}
	fmt.Print("done fetching CDS regions\n")
	if *mates {
		numSeqs = numSeqs + len(seqMap1)
		fmt.Printf("total number of strains: %d\n", numSeqs)
	} else {
		fmt.Printf("total number of strains: %d\n", numSeqs)
	}
	fmt.Printf("total number of codons: %d\n", numCodons)

	//initialize output csv
	outFile := *outPrefix + ".csv"
	initCsvOut(outFile)
	if *mates {
		calcQsMatesAll(seqMap, seqMap1, codonOffset, codonPos-1, minCodonLen, maxCodonLen, codingTable, synonymous, outFile, *numDigesters)
	} else {
		calcQsAll(seqMap, codonOffset, codonPos-1, minCodonLen, maxCodonLen, codingTable, synonymous, outFile, *numDigesters)
	}

	duration := time.Since(start)
	fmt.Println("Time to calculate LD:", duration)
	//corrResChan := calcSingleClade(alnChan, calculator)
	//what's in the json is actually Qs NOT P2!
	//resChan := mcorr.PipeOutCorrResults(corrResChan, *outPrefix+".json")
	//division by d_sample or P2 is not until here!!!
	//CollectWrite(corrResChan, *outPrefix+".csv")
}

func makeSeqMap(startSlice []int, alnFile string, codonOffset int) (seqMap map[string][]Codon) {
	seqMap = make(map[string][]Codon)
	for _, i := range startSlice {
		a := getGene(alnFile, i)
		addCodons(a, seqMap, codonOffset)
	}
	return seqMap
}

func combinedSeqMap(startSlice []int, alnFile, mateAln string, codonOffset int) (seqMap map[string][]Codon) {
	seqMap = make(map[string][]Codon)
	for _, i := range startSlice {
		// get the gene alignment from the first file
		aln1 := getGene(alnFile, i)
		//get the gene alignment from the second file
		aln2 := getGene(mateAln, i)
		// add on the alignment 2 sequences onto alignment 1
		aln1.Sequences = append(aln1.Sequences, aln2.Sequences...)
		addCodons(aln1, seqMap, codonOffset)
	}
	return seqMap
}

func getGene(alnFile string, startCodon int) (gene Alignment) {
	c := readAlignments(alnFile)
	for a := range c {
		startPos := getStartPos(a)
		if startCodon == startPos {
			gene = a
			break
		}
	}
	return gene
}

func getStartPos(aln Alignment) int {
	genomePos := aln.genePos
	terms := strings.Split(genomePos, "+")
	startPos, _ := strconv.Atoi(terms[0])
	return startPos
}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string //gene ID
	genePos   string // position of gene on the genome
	Sequences []seq.Sequence
}

// calcSingleClade calculate correlation functions in a single cluster of sequence.
func calcSingleClade(alnChan chan Alignment, calculator Calculator) (corrResChan chan CorrResults) {
	corrResChan = make(chan CorrResults)
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
			header := strings.Split(alignment[0].Id, " ")
			alnID := header[0]
			genomePos := header[1]
			alnChan <- Alignment{ID: alnID, genePos: genomePos, Sequences: alignment}
		}
	}()

	return
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
			// can have 1 sequence in both xmfas
			if len(a) >= 1 {
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

func getNames(s string) (geneName, genomeName string) {
	terms := strings.Split(s, " ")
	//this is the genomeName for the MSA files assembled from ReferenceAlignmentGenerator
	geneName = terms[0]
	genomeName = terms[2]
	return
}

//addCodons adds codons to each strain sequence in the sequence map
func addCodons(a Alignment, seqMap map[string][]Codon, codonOffset int) {
	//add codons to the map for each sequence
	var wg sync.WaitGroup
	//tell wg how many threads will be running concurrently
	wg.Add(len(a.Sequences))
	//mutex lock for safe access to seqMap
	var mutex = &sync.Mutex{}
	alnSeqs := a.Sequences
	//give duplicate strain names a suffix to differentiate
	duplicates := make(map[string]int)
	for j := 0; j < len(alnSeqs); j++ {
		go func(j int) {
			defer wg.Done()
			s := alnSeqs[j]
			codons := extractCodons(s, codonOffset)
			_, seqName := getNames(s.Id)
			//check if the strain name has been done ...
			mutex.Lock()
			i, found := duplicates[seqName]
			//update the count for the strain
			duplicates[seqName]++
			mutex.Unlock()
			if found {
				id := strconv.Itoa(i)
				seqName = seqName + "_" + id
			}
			mutex.Lock()
			seqMap[seqName] = append(seqMap[seqName], codons...)
			mutex.Unlock()
		}(j)
	}
	wg.Wait()
}
