// Copyright 2022 Asher Preska Steinberg
//
// Initially modified from nucl_cov.go of mcorr
// https://github.com/kussell-lab/mcorr
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
	"bytes"
	"fmt"
)

// NuclCov contains covariance of nucleotide acid in a DNA sequence.
type NuclCov struct {
	//the doublet matrix contains a position for each of the 16 combos
	//that site A and B (i and i+l) could be for a given sequence
	// (e.g., AA, AT, AC, AG is the first row)
	//a 1 is filled in at one of the 16 indices of the matrix for each
	//sequence which has the combination
	Doublets []int
	Alphabet []byte
}

// NewNuclCov return a NuclCov given the alphabet.
func NewNuclCov(alphabet []byte) *NuclCov {
	sizeOfAlphabet := len(alphabet)
	nc := NuclCov{Alphabet: alphabet}
	nc.Doublets = make([]int, sizeOfAlphabet*sizeOfAlphabet)
	return &nc
}

// Add insert a pair of nucliotide acids.
// It returns error when the nucliotide acid is not in the alphabet.
//the doublet matrix contains a position for each of the 16 combinations
//site A and B (i and i+l) could be for a given sequence
func (nc *NuclCov) Add(a, b byte) error {
	//the alphabet is {A,T,C,G} so indexA or B will be
	//indexA = 0 for A, 1 for T, 2 for C, 3 for G, and -1 if there's nothing there
	//indexA is for site i, indexB is for site i+l
	indexA := bytes.IndexByte(nc.Alphabet, a)
	indexB := bytes.IndexByte(nc.Alphabet, b)
	sizeOfAlphabet := len(nc.Alphabet)
	if indexA >= 0 && indexB >= 0 {
		//so the doublet is essential a 4x4 matrix, with a 1 filled in
		//at the spot of the matrix where this combination would be (AA, AT, AC, AG, TC, TG, etc)
		nc.Doublets[indexA*sizeOfAlphabet+indexB]++
		return nil
	}

	var err error
	if indexA < 0 && indexB < 0 {
		err = fmt.Errorf("%c and %c are not in Alphabet: %s", a, b, string(nc.Alphabet))
	} else if indexA < 0 {
		err = fmt.Errorf("%c is not in Alphabet: %s", a, string(nc.Alphabet))
	} else {
		err = fmt.Errorf("%c is not in Alphabet: %s", b, string(nc.Alphabet))
	}

	return err
}

// Count returns the total number of pairs.
func (nc *NuclCov) Count() int {
	n := 0
	for _, a := range nc.Doublets {
		n += a
	}
	return n
}

// P00 returns the probability of 00.
func (nc *NuclCov) P00(minAlleleNum int) (xy float64, n int) {
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > minAlleleNum {
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > minAlleleNum {
					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
			xy += float64(nc.Doublets[i] * (nc.Doublets[i] - 1) / 2)
		}
	}
	return
}

// P11 returns the probability of 11 (difference at both sites)
// this goes through the doublet matrix, and calculates each possible doublet pair permutation
// (e.g., AA*AT
// doublet pairs which emit no substitution at one or both sites are counted
// in the sequence pair count n, but not the covariance (xy)
// (e.g., AA*AT, emits no substitution at site A, so goes into n but not xy;
//
func (nc *NuclCov) P11(minAlleleNum int) (xy float64, n int) {
	sizeOfAlphabet := len(nc.Alphabet)
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > minAlleleNum {
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > minAlleleNum {
					c := float64(nc.Doublets[i] * nc.Doublets[j])
					if i%sizeOfAlphabet != j%sizeOfAlphabet && i/sizeOfAlphabet != j/sizeOfAlphabet {
						xy += c
					}
					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
		}
	}
	return
}

// MateP11 calculate covariance between two clusters.
func (nc *NuclCov) MateP11(nc2 *NuclCov, minAlleleNum int) (xy float64, n int) {
	sizeOfAlphabet := len(nc.Alphabet)
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > minAlleleNum {
			for j := 0; j < len(nc2.Doublets); j++ {
				if i != j && nc2.Doublets[j] > minAlleleNum {
					c := float64(nc.Doublets[i] * nc2.Doublets[j])
					if i%sizeOfAlphabet != j%sizeOfAlphabet && i/sizeOfAlphabet != j/sizeOfAlphabet {
						xy += c
					}
				}
			}
		}
	}
	n1 := 0
	n2 := 0
	for i := 0; i < len(nc.Doublets); i++ {
		n1 += nc.Doublets[i]
		n2 += nc2.Doublets[i]
	}
	n = n1 * n2
	return
}

// MateP00 calculate covariance between two clusters.
func (nc *NuclCov) MateP00(nc2 *NuclCov, minAlleleNum int) (xy float64, n int) {
	n1, n2 := 0, 0
	for i := 0; i < len(nc.Doublets); i++ {
		xy += float64(nc.Doublets[i] * nc2.Doublets[i])
		n1 += nc.Doublets[i]
		n2 += nc2.Doublets[i]
	}
	n = n1 * n2
	return
}

// Append another NuclCov.
func (nc *NuclCov) Append(nc2 *NuclCov) error {
	// Check alphabet
	diffAlphabetError := fmt.Errorf("Different alphabet %s, %s", string(nc.Alphabet), string(nc2.Alphabet))
	if len(nc.Alphabet) != len(nc2.Alphabet) {
		return diffAlphabetError
	}
	for i, a := range nc.Alphabet {
		b := nc2.Alphabet[i]
		if a != b {
			return diffAlphabetError
		}
	}

	for i := 0; i < len(nc.Doublets); i++ {
		nc.Doublets[i] += nc2.Doublets[i]
	}

	return nil
}

// covXY returns the joint probability of events at xy, xx, and yy
//this won't work ... remember it iterates over sequence pairs so ii and jj are just the same sequence pair
//what you want to do is feed in multiple sets of doublets, one for xx, one for yy and one for xy
//you loop through the doublets, which are 4x4 matrices of ATCG combinations for a sequence pair
func (nc *NuclCov) covXY(minAlleleNum int) (xy float64, xx float64, yy float64, n int) {
	//this is 4 (ATCG is our alphabet)
	sizeOfAlphabet := len(nc.Alphabet)
	//loop through sequence pairs and calculate the covariance matrix
	for i := 0; i < len(nc.Doublets); i++ {
		if nc.Doublets[i] > minAlleleNum {
			//start at i + 1 because this then calculates the upper triangle matrix
			//and skips diagonal (which is self vs self)
			for j := i + 1; j < len(nc.Doublets); j++ {
				if nc.Doublets[j] > minAlleleNum {
					//joint event count for site x and y
					c := float64(nc.Doublets[i] * nc.Doublets[j])
					//joint event count for site x
					cXX := float64(nc.Doublets[i] * nc.Doublets[i])
					//joint event count for site y
					cYY := float64(nc.Doublets[j] * nc.Doublets[j])
					if i%sizeOfAlphabet != j%sizeOfAlphabet && i/sizeOfAlphabet != j/sizeOfAlphabet {
						xy += c
						xx += cXX
						yy += cYY
					}

					n += nc.Doublets[i] * nc.Doublets[j]
				}
			}
			n += nc.Doublets[i] * (nc.Doublets[i] - 1) / 2
		}
	}
	return
}
