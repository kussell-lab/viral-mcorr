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

//script by Asher Preska Steinberg (apsteinberg@nyu.edu)
import (
	"fmt"
	"github.com/kussell-lab/biogo/seq"
	"github.com/kussell-lab/ncbiftp/taxonomy"
	"gonum.org/v1/gonum/stat/combin"
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v2"
	"io"
	//"math/bits"
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
	app := kingpin.New("mcorrPairGenome", "Calculate position-averaged correlation profiles across whole genome for strain pairs from sequence alignments in XMFA format.")
	app.Version("v20211013")

	alnFile := app.Arg("in", "Alignment file in XMFA format.").Required().String()
	outPrefix := app.Arg("out", "Output prefix.").Required().String()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	//numBoot := app.Flag("num-boot", "Number of bootstrapping on genes").Default("1000").Int()
	maxl := app.Flag("max-corr-length", "Maximum distance of correlation (base pairs)").Default("300").Int()
	mateAln := app.Flag("mate-aln", "Second alignment").Default("").String()
	numDigesters := app.Flag("num-threads", "number of threads").Default("8").Int()
	showProgress := app.Flag("show-progress", "Show progress").Bool()

	kingpin.MustParse(app.Parse(os.Args[1:]))

	if *ncpu <= 0 {
		*ncpu = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*ncpu)

	if *numDigesters > *ncpu {
		*numDigesters = *ncpu
	}

	//timer
	start := time.Now()

	// prepare calculator.
	//var calculator Calculator
	codingTable := taxonomy.GeneticCodes()["11"]
	maxCodonLen := *maxl / 3

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
	//initialize output csv
	outFile := *outPrefix + ".csv"
	initCsvOut(outFile)

	// now go through and make a map of codon sequences
	// where we add on codons to the end of each strain sequence
	fmt.Print("fetching CDS regions\n")
	//the keys are the strain names and the values are the codon sequences for each
	//var seqMap map[string][]Codon
	//add additional sequences if necessary
	var seqMap map[string][]Codon
	if *mateAln == "" {
		seqMap = makeSeqMap(startSlice, *alnFile, codonOffset)
	} else {
		seqMap = combinedSeqMap(startSlice, *alnFile, *mateAln, codonOffset)
	}

	//make a map of all sequence names and get all pairs
	seqNameMap := make(map[int]string)
	i := 0
	for seqName, _ := range seqMap {
		seqNameMap[i] = seqName
		i = i + 1
	}

	numSeqs := i
	pairList := combin.Combinations(numSeqs, 2)

	var seqpairs [][]string
	for _, pair := range pairList {
		seq1 := seqNameMap[pair[0]]
		seq2 := seqNameMap[pair[1]]
		seqpair := []string{seq1, seq2}
		seqpairs = append(seqpairs, seqpair)
	}
	//
	//test := []string{"ABC", "B","C"}
	//testpairs := Combinations(test, 2)
	//fmt.Println(testpairs)
	//seqpairs := Combinations(seqNames, 2)

	numpairs := len(pairList)
	fmt.Println(numpairs, "of pairwise corr profiles to compute")
	// show progress bar
	var bar *pb.ProgressBar
	if *showProgress {
		//max := maxCodonLen
		bar = pb.StartNew(numpairs)
		defer bar.Finish()
	}
	calcQsAll(seqMap, seqpairs, codonOffset, codonPos-1, maxCodonLen,
		codingTable, synonymous, outFile, *numDigesters, bar)

	//time it
	duration := time.Since(start)
	fmt.Println("Time to calculate pairwise Ks:", duration)
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

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string //gene ID
	genePos   string // position of gene on the genome
	Sequences []seq.Sequence
}

//initCsvOut initializes the output csv
func initCsvOut(outFile string) {
	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	w.WriteString("l,m,v,n,t,b\n")
	w.Close()
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

func getStartPos(aln Alignment) int {
	genomePos := aln.genePos
	terms := strings.Split(genomePos, "+")
	startPos, _ := strconv.Atoi(terms[0])
	return startPos
}
