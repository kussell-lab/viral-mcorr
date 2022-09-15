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

//script written by Asher Preska Steinberg (apsteinberg@nyu.edu)

import (
	"bufio"
	"fmt"
	"github.com/boltdb/bolt"
	"github.com/kussell-lab/biogo/seq"
	"gopkg.in/alecthomas/kingpin.v2"
	"io"
	"log"
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
	app := kingpin.New("makeGeneDB", "Convert XMFA file into boltdb files for each gene.")
	app.Version("v20211003")

	alnFile := app.Arg("aln", "Alignment file in XMFA format.").Required().String()
	//outPrefix := app.Arg("out", "Output prefix.").Required().String()
	//minl := app.Flag("min-corr-length", "min distance of correlation (base pairs)").Default("0").Int()
	//maxl := app.Flag("max-corr-length", "Max distance of correlation (base pairs)").Default("0").Int()
	//mateAln := app.Flag("mate-aln", "Second alignment for calculating between clades").Default("").String()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	//mates := app.Flag("between-clades", "just calculate correlations between pairs from different xmfa files").Default("false").Bool()
	//numDigesters := app.Flag("num-threads", "number of threads").Default("50").Int()
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
	//codingTable := taxonomy.GeneticCodes()["11"]
	//maxCodonLen := *maxl / 3
	//minCodonLen := *minl / 3
	//
	//synonymous := true
	//codonPos := 3
	codonOffset := 0

	// make a list of starting positions so we know how to order things and
	// dump the alignments into a map with starting positions as keys
	var startSlice []int
	c := readAlignments(*alnFile)

	for a := range c {
		startPos, _ := getStartStop(a)
		startSlice = append(startSlice, startPos)
	}
	//make the codon databases ...

	//sort the slice numerically
	sort.Ints(startSlice)
	print(startSlice)
	//get the number of codons ....

	// now go through and make a map of codon sequences
	// where we add on codons to the end of each strain sequence
	fmt.Print("fetching CDS regions\n")
	//make codon databases
	numSeqs, numCodons, startCodons, _ := makeCodonDB(startSlice, *alnFile, codonOffset)

	fmt.Print("done fetching CDS regions\n")

	writeCodonStartList(startCodons)

	fmt.Printf("total number of strains: %d\n", numSeqs)

	fmt.Printf("total number of codons: %d\n", numCodons)

	duration := time.Since(start)
	fmt.Println("Time to make gene db files:", duration)
}

func makeSeqMap(startSlice []int, alnFile string, codonOffset int) (seqMap map[string][]Codon) {
	seqMap = make(map[string][]Codon)
	for _, i := range startSlice {
		a, _ := getGene(alnFile, i)
		addCodons(a, seqMap, codonOffset)
	}
	return seqMap
}

func combinedSeqMap(startSlice []int, alnFile, mateAln string, codonOffset int) (seqMap map[string][]Codon) {
	seqMap = make(map[string][]Codon)
	for _, i := range startSlice {
		// get the gene alignment from the first file
		aln1, _ := getGene(alnFile, i)
		//get the gene alignment from the second file
		aln2, _ := getGene(mateAln, i)
		// add on the alignment 2 sequences onto alignment 1
		aln1.Sequences = append(aln1.Sequences, aln2.Sequences...)
		addCodons(aln1, seqMap, codonOffset)
	}
	return seqMap
}

func getGene(alnFile string, startCodon int) (gene Alignment, stopPos int) {
	fmt.Printf("retrieving codons starting at " + strconv.Itoa(startCodon) + "\n")
	c := readAlignments(alnFile)
	for a := range c {
		var startPos int
		startPos, stopPos = getStartStop(a)
		if startCodon == startPos {
			gene = a
			break
		}
	}
	return gene, stopPos
}

// getStartStop get start/stop positions for each gene
func getStartStop(aln Alignment) (int, int) {
	genomePos := aln.genePos
	terms := strings.Split(genomePos, "+")
	startPos, _ := strconv.Atoi(terms[0])
	stopPos, _ := strconv.Atoi(terms[1])
	return startPos, stopPos
}

// Alignment is an array of mutliple sequences with same length.
type Alignment struct {
	ID        string //gene ID
	genePos   string // position of gene on the genome
	Sequences []seq.Sequence
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

//makeStrainMap maps codons to strains
func makeStrainMap(a Alignment, codonOffset int) (strainMap map[string][]Codon) {
	strainMap = make(map[string][]Codon)
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
			strainMap[seqName] = append(strainMap[seqName], codons...)
			mutex.Unlock()
		}(j)
	}
	wg.Wait()
	return
}

// makeCodonMap assembles the codons of each gene into a map where keys are codon positions
// and the values are all of the codons for each sequence
// in the order defined by strain list
func makeCodonMap(strainMap map[string][]Codon, startCodon int, strainList []string) (codonMap map[int][]Codon) {
	//add codons for each sequences to a map where the keys are codon positions and the values
	// are the codon sequence for each position
	codonMap = make(map[int][]Codon)
	for _, strain := range strainList {
		codonSeq := strainMap[strain]
		for k, cc := range codonSeq {
			codonPos := startCodon + k
			codonMap[codonPos] = append(codonMap[codonPos], cc)
		}
	}
	return
}

// makeCodonDB makes a keymap database for each codon position
//and return a list of databases ...
func makeCodonDB(startSlice []int, alnFile string, codonOffset int) (numSeqs int, numCodons int, startCodons []int, dbMap map[int]*bolt.DB) {
	startCodon := 0
	var codonMap map[int][]Codon
	dbMap = make(map[int]*bolt.DB)
	var strainList []string
	for _, startPos := range startSlice {
		a, stopPos := getGene(alnFile, startPos)
		strainMap := makeStrainMap(a, codonOffset)
		if startCodon == 0 {
			for strain, _ := range strainMap {
				strainList = append(strainList, strain)
			}
		}
		codonMap = makeCodonMap(strainMap, startCodon, strainList)
		//make the name of the database
		dbName := strconv.Itoa(startCodon)
		db := createDB(dbName + ".db")
		createBucket(db, "codons")
		loadCodons(db, "codons", codonMap)
		dbMap[startCodon] = db
		db.Close()
		//add to start codon slice
		startCodons = append(startCodons, startCodon)
		//calculate the first codon position in the next db
		//get the length of the cds in codons ...
		cdslen := (stopPos - startPos + 1) / 3
		startCodon = startCodon + cdslen
		//fmt.Print(strconv.Itoa(startCodon)+"\n")

	}
	//startCodons = append(startCodons, startCodon)
	for _, v := range codonMap {
		numSeqs = len(v)
	}
	numCodons = startCodon
	return
}

type codonDB struct {
	db         *bolt.DB
	startCodon int
}

// createDB creates a bolt db.
func createDB(dbFile string) *bolt.DB {
	db, err := bolt.Open(dbFile, 0600, nil)
	if err != nil {
		log.Fatal(err)
	}
	return db
}

// createBucket creates a bucket.
func createBucket(db *bolt.DB, bucketName string) {
	fn := func(tx *bolt.Tx) error {
		_, err := tx.CreateBucketIfNotExists([]byte(bucketName))
		return err
	}

	err := db.Update(fn)
	if err != nil {
		log.Fatal(err)
	}
}

func loadCodons(db *bolt.DB, bucketName string, codonMap map[int][]Codon) {
	fn := func(tx *bolt.Tx) error {
		b := tx.Bucket([]byte(bucketName))
		for pos, codon := range codonMap {
			//codon position
			codonPos := strconv.Itoa(pos)
			//if pos == 5000 {
			//	print("mammajamma")
			//}
			//all of the alleles of the codon for that position
			codonBytes := CodonstoBytes(codon, 0)
			err := b.Put([]byte(codonPos), codonBytes)
			if err != nil {
				return err
			}
		}

		return nil
	}

	err := db.Update(fn)
	if err != nil {
		log.Fatal(err)
	}
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

// CodonstoBytes convert codons to bytes
func CodonstoBytes(codons []Codon, offset int) (s []byte) {
	for _, cc := range codons {
		s = append(s, cc...)
	}
	return
}

//writeCodonStartList writes codon starts to a list so we know the names of the boltdb files
func writeCodonStartList(sampledata []int) {
	file, err := os.Create("gene_boltdb_list.txt")

	if err != nil {
		log.Fatalf("failed creating file: %s", err)
	}

	datawriter := bufio.NewWriter(file)

	for _, data := range sampledata {
		_, _ = datawriter.WriteString(strconv.Itoa(data) + "\n")
	}

	datawriter.Flush()
	file.Close()
}
