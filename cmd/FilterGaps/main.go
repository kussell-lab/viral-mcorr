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
	"gopkg.in/alecthomas/kingpin.v2"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"
)

func main() {

	app := kingpin.New("FilterGaps", "Filters gene alignments with >10% gaps (or user defined threshold)")
	app.Version("v20210305")
	alnFile := app.Arg("master_MSA", "multi-sequence alignment file for all genes").Required().String()
	outdir := app.Arg("outdir", "output directory for the filtered MSA").Required().String()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	numFilters := app.Flag("threads", "Number of alignments to process at a time (default: 8)").Default("8").Int()
	cutoff := app.Flag("cutoff", "cutoff percentage; default is 10%").Default("10").Int()
	fillgaps := app.Flag("fill-gaps", "fill gene sequences that didn't make the cutoff with dashes as a placeholder").Default("False").Bool()
	kingpin.MustParse(app.Parse(os.Args[1:]))
	if *ncpu == 0 {
		*ncpu = runtime.NumCPU()
	}

	runtime.GOMAXPROCS(*ncpu)

	//alnFile := "/Volumes/aps_timemachine/recombo/APS168_gapfiltered/test_3"
	//outdir := "/Volumes/aps_timemachine/recombo/APS168_gapfiltered/gapfiltered"
	//numFilters := 4
	//cutoff := 99
	//timer
	start := time.Now()
	makeFilteredMSA(*outdir, *alnFile, *cutoff)
	done := make(chan struct{})
	//read in alignments
	alignments, errc := readAlignments(done, *alnFile)
	//start a fixed number of goroutines to read alignments and split into core/flex
	c := make(chan Alignment)
	//convert to fraction
	gapCutoff := float64(*cutoff) / 100
	var wg sync.WaitGroup
	for i := 0; i < *numFilters; i++ {
		wg.Add(1)
		go FilterGappedAlns(done, alignments, c, gapCutoff, *fillgaps, i, &wg)
	}

	go func() {
		wg.Wait()
		close(c)
	}()
	//end of pipeline; write files
	for filteredAln := range c {
		if len(filteredAln.Sequences) > 0 {
			writeMSA(filteredAln, *outdir, *alnFile, *cutoff)
		}
	}
	if err := <-errc; err != nil { // HLerrc
		panic(err)
	}
	//add the number of core and flex to the bottom of the spreadsheet

	duration := time.Since(start)
	fmt.Println("Time to filter gapped alignments:", duration)
}

// readAlignments reads sequence alignment from a extended Multi-FASTA file,
// and return a channel of alignment, which is a list of seq.Sequence
func readAlignments(done <-chan struct{}, file string) (<-chan Alignment, <-chan error) {
	alignments := make(chan Alignment)
	errc := make(chan error, 1)
	go func() {
		defer close(alignments)

		f, err := os.Open(file)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		xmfaReader := seq.NewXMFAReader(f)
		numAln := 0
		for {
			alignment, err := xmfaReader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			if len(alignment) > 0 {
				numAln++
				alnID := strings.Split(alignment[0].Id, " ")[0]
				select {
				case alignments <- Alignment{alnID, numAln, alignment}:
					fmt.Printf("\rRead %d alignments.", numAln)
					fmt.Printf("\r alignment ID: %s", alnID)
				case <-done:
					fmt.Printf(" Total alignments %d\n", numAln)
				}
			}
		}
		errc <- err
	}()
	return alignments, errc
}

// Alignment is an array of multiple sequences with same length.
type Alignment struct {
	ID        string
	num       int
	Sequences []seq.Sequence
}

// FilterGappedAlns reads gene alignments from the master MSA, filters out sequences with >10% gaps,
// or whatever cutoff the user sets
// then sends these processed results on alnChan until either the master MSA or done channel is closed.
func FilterGappedAlns(done <-chan struct{}, alignments <-chan Alignment, filteredAlns chan<- Alignment,
	cutoff float64, fillgaps bool, id int, wg *sync.WaitGroup) {
	//threshold := float64(cutoff)/100
	defer wg.Done()
	//fmt.Printf("Worker %d starting\n", id)
	for aln := range alignments { // HLpaths
		//collect those sequences with =< 10% gaps
		var filteredSeqs []seq.Sequence
		for _, s := range aln.Sequences {
			//count gaps in the gene alignment
			gaps := countGaps(s)
			//gene alignment length
			seqLength := float64(len(s.Seq))
			percentGaps := gaps / seqLength
			//add sequences that were below the cutoff
			if percentGaps <= cutoff {
				filteredSeqs = append(filteredSeqs, s)
			} else if fillgaps == true {
				//if fillgaps is true add in dashes for sequences that didn't make the cut
				var blankseq seq.Sequence
				blankseq = s
				//fill the seq with dashes
				gapFiller(blankseq)
				//append to filteredSeqs
				filteredSeqs = append(filteredSeqs, blankseq)
			}

		}
		//just include those gene alignments with <=10% gaps
		var filteredAln Alignment
		filteredAln = Alignment{aln.ID, aln.num, filteredSeqs}

		select {
		case filteredAlns <- filteredAln:
		case <-done:
			return
		}
	}
	//fmt.Printf("Worker %d done\n", id)

}

// countGaps counts the number of gaps in a gene sequence
func countGaps(s seq.Sequence) (NumGaps float64) {
	for i := 0; i < len(s.Seq); i++ {
		b := s.Seq[i]
		if b == '-' || b == 'N' {
			NumGaps++
		}
	}
	return
}

// gapFiller fills a sequence which didn't make the cutoff with dashes
func gapFiller(s seq.Sequence) {
	for i := 0; i < len(s.Seq); i++ {
		s.Seq[i] = '-'
	}
	return
}

//makeFilteredMSA makes the outdir and initializes the MSA files for core and flexible genomes
func makeFilteredMSA(outdir string, alnFile string, cutoff int) {
	if _, err := os.Stat(outdir); os.IsNotExist(err) {
		os.Mkdir(outdir, 0755)
	}
	terms := strings.Split(alnFile, "/")
	alnFileName := terms[len(terms)-1]
	//take care of file extensions ...
	terms = strings.Split(alnFileName, ".")
	var MSAname string
	if len(terms) > 1 {
		alnFileName = terms[0]
		extension := terms[1]
		cut := strconv.Itoa(cutoff)
		MSAname = alnFileName + "_leq" + cut + "_gaps." + extension
	} else {
		cut := strconv.Itoa(cutoff)
		MSAname = alnFileName + "_leq" + cut + "_gaps"
	}

	MSA := filepath.Join(outdir, MSAname)
	f, err := os.Create(MSA)
	check(err)
	f.Close()
	f, err = os.Create(MSA)
	check(err)
	f.Close()
}

//check for errors
func check(e error) {
	if e != nil {
		panic(e)
	}
}

//writeMSA write the gene to the correct MSA (core or flex)
func writeMSA(c Alignment, outdir string, alnFile string, cutoff int) {
	terms := strings.Split(alnFile, "/")
	alnFileName := terms[len(terms)-1]
	//take care of file extensions ...
	terms = strings.Split(alnFileName, ".")
	var MSAname string
	if len(terms) > 1 {
		alnFileName = terms[0]
		extension := terms[1]
		cut := strconv.Itoa(cutoff)
		MSAname = alnFileName + "_leq" + cut + "_gaps." + extension
	} else {
		cut := strconv.Itoa(cutoff)
		MSAname = alnFileName + "_leq" + cut + "_gaps"
	}
	MSA := filepath.Join(outdir, MSAname)
	//f, err := os.OpenFile(MSA, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0600)
	//f, err := os.Create(MSA)
	f, err := os.OpenFile(MSA, os.O_APPEND|os.O_WRONLY, 0600)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	aln := c
	for _, s := range aln.Sequences {
		f.WriteString(">" + s.Id + "\n")
		f.Write(s.Seq)
		f.WriteString("\n")
	}
	f.WriteString("=\n")
}
