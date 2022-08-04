package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/bio/sketches"
	"github.com/twotwotwo/sorts/sortutil"
)

func main() {
	if len(os.Args) != 6 {
		fmt.Printf("usage: %s <fasta file with ONE sequence> <k> <m> <scale> <mutations: [12]>\n", filepath.Base(os.Args[0]))
		os.Exit(0)
	}

	// go tool pprof -http=:8080 cpu.pprof
	// defer profile.Start(profile.CPUProfile, profile.ProfilePath(".")).Stop()

	// go tool trace -http=:8080 trace.out
	// defer profile.Start(profile.TraceProfile, profile.ProfilePath(".")).Stop()

	// go tool pprof -http=:8080 mem.pprof
	// defer profile.Start(profile.MemProfile, profile.MemProfileRate(1), profile.ProfilePath(".")).Stop()
	// defer profile.Start(profile.MemProfile, profile.ProfilePath(".")).Stop()

	file := os.Args[1]
	k, err := strconv.Atoi(os.Args[2])
	if err != nil || k <= 0 {
		checkError(fmt.Errorf("<k> should be positive interger: %s", os.Args[2]))
	}
	m, err := strconv.Atoi(os.Args[3])
	if err != nil || m <= 0 {
		checkError(fmt.Errorf("<m> should be positive interger: %s", os.Args[3]))
	}
	scale, err := strconv.Atoi(os.Args[4])
	if err != nil || m <= 0 {
		checkError(fmt.Errorf("<scale> should be positive interger: %s", os.Args[4]))
	}
	nMutations, err := strconv.Atoi(os.Args[5])
	if err != nil || m <= 0 {
		checkError(fmt.Errorf("<mutations> should be positive interger: %s", os.Args[5]))
	}
	if nMutations > 2 {
		checkError(fmt.Errorf("<mutations> should be in range of [1, 2]: %s", os.Args[5]))
	}

	outfh := bufio.NewWriter(os.Stdout)
	defer outfh.Flush()

	fmt.Fprintf(outfh, "# file: %s\n", file)
	fmt.Fprintf(outfh, "# k-mer: %d\n", k)
	fmt.Fprintf(outfh, "# m-mer: %d\n", m)
	fmt.Fprintf(outfh, "# scale: %d\n", scale)
	fmt.Fprintf(outfh, "# mutations: %d\n", nMutations)

	fmt.Fprintf(outfh, "pos\tstatus\tkmer\thash\tmutant\tmpos\tref\talt\tmpos2\tref2\talt2\n")

	outfh.Flush()

	fastxReader, err := fastx.NewDefaultReader(file)
	if err != nil {
		checkError(fmt.Errorf("%s: %s", file, err))
	}

	var i, end int
	var record *fastx.Record

	for {
		record, err = fastxReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(err)
			break
		}

		if len(record.Seq.Seq) < k {
			continue
		}

		record.Seq.Seq = bytes.ToUpper(record.Seq.Seq)

		end = len(record.Seq.Seq) - k

		var wg sync.WaitGroup
		ch := make(chan *KmerResult, runtime.NumCPU())
		done := make(chan int)

		go func() {
			var id int
			var ok bool
			buf := make(map[int]*KmerResult, 1024)

			var mutant Mutant
			for r := range ch {
				if id == r.Pos {
					if len(r.Mutants) == 0 {
						fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
							r.Pos+1, 0, r.Kmer, r.Hash, "", "", "", "", "", "", "")
					} else {
						for _, mutant = range r.Mutants {
							switch len(mutant.Changes) {
							case 1:
								fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%02d\t%c\t%c\t\t\t\n",
									r.Pos+1, 1, r.Kmer, r.Hash, mutant.Kmer,
									mutant.Changes[0].Pos+1, r.Kmer[mutant.Changes[0].Pos], mutant.Changes[0].Base)
							case 2:
								fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%02d\t%c\t%c\t%02d\t%c\t%c\n",
									r.Pos+1, 1, r.Kmer, r.Hash, mutant.Kmer,
									mutant.Changes[0].Pos+1, r.Kmer[mutant.Changes[0].Pos], mutant.Changes[0].Base,
									mutant.Changes[1].Pos+1, r.Kmer[mutant.Changes[1].Pos], mutant.Changes[1].Base)
							}
						}
					}

					id++
					continue
				}

				buf[r.Pos] = r

				if r, ok = buf[id]; ok {
					if len(r.Mutants) == 0 {
						fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
							r.Pos+1, 0, r.Kmer, r.Hash, "", "", "", "", "", "", "")
					} else {
						for _, mutant = range r.Mutants {
							switch len(mutant.Changes) {
							case 1:
								fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%02d\t%c\t%c\t\t\t\n",
									r.Pos+1, 1, r.Kmer, r.Hash, mutant.Kmer,
									mutant.Changes[0].Pos+1, r.Kmer[mutant.Changes[0].Pos], mutant.Changes[0].Base)
							case 2:
								fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%02d\t%c\t%c\t%02d\t%c\t%c\n",
									r.Pos+1, 1, r.Kmer, r.Hash, mutant.Kmer,
									mutant.Changes[0].Pos+1, r.Kmer[mutant.Changes[0].Pos], mutant.Changes[0].Base,
									mutant.Changes[1].Pos+1, r.Kmer[mutant.Changes[1].Pos], mutant.Changes[1].Base)
							}
						}
					}

					delete(buf, id)
					id++
				}
			}
			if len(buf) > 0 {
				ids := make([]int, len(buf))
				i := 0
				for id = range buf {
					ids[i] = id
					i++
				}
				sortutil.Ints(ids)
				for _, id = range ids {
					r := buf[id]

					if len(r.Mutants) == 0 {
						fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
							r.Pos+1, 0, r.Kmer, r.Hash, "", "", "", "", "", "", "")
					} else {
						for _, mutant = range r.Mutants {
							switch len(mutant.Changes) {
							case 1:
								fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%02d\t%c\t%c\t\t\t\n",
									r.Pos+1, 1, r.Kmer, r.Hash, mutant.Kmer,
									mutant.Changes[0].Pos+1, r.Kmer[mutant.Changes[0].Pos], mutant.Changes[0].Base)
							case 2:
								fmt.Fprintf(outfh, "%d\t%d\t%s\t%d\t%s\t%02d\t%c\t%c\t%02d\t%c\t%c\n",
									r.Pos+1, 1, r.Kmer, r.Hash, mutant.Kmer,
									mutant.Changes[0].Pos+1, r.Kmer[mutant.Changes[0].Pos], mutant.Changes[0].Base,
									mutant.Changes[1].Pos+1, r.Kmer[mutant.Changes[1].Pos], mutant.Changes[1].Base)
							}
						}
					}
				}
			}

			done <- 1
		}()

		for i = 0; i <= end; i++ {
			wg.Add(1)

			go func(i int, _seq *seq.Seq) {
				switch nMutations {
				case 1:
					ch <- checkKmer(i, _seq, k, m, scale)
				case 2:
					ch <- checkKmer2(i, _seq, k, m, scale)
				}
				wg.Done()
			}(i, record.Seq.SubSeq(i+1, i+k))
		}
		wg.Wait()
		close(ch)
		<-done

		break
	}
}

var basesUpper = [4]byte{'A', 'C', 'G', 'T'}
var basesLower = [4]byte{'a', 'c', 'g', 't'}

type Change struct {
	Pos  int
	Base byte
}
type Mutant struct {
	Kmer    []byte
	Changes []Change
}
type KmerResult struct {
	Pos     int
	Kmer    []byte
	Hash    uint64
	Mutants []Mutant
}

func checkKmer2(pos int, kmer *seq.Seq, k int, m int, scale int) *KmerResult {
	var iter *sketches.Iterator
	var err error

	iter, err = sketches.NewSimHashIterator(kmer, k, m, scale, true, false)
	if err != nil {
		checkError(fmt.Errorf("%s: %s", kmer, err))
	}
	hash0, _ := iter.NextSimHash()

	var i, j, ii, jj int
	var b, B, bb, BB byte
	var hash uint64
	mutants := make([]Mutant, 0, 8)
	for j, B = range kmer.Seq {
		for jj = j + 1; jj < k; jj++ {
			BB = kmer.Seq[jj]
			for i, b = range basesUpper {
				if b == B {
					continue
				}
				for ii, bb = range basesUpper {
					if bb == BB {
						continue
					}

					mutant := kmer.Clone2()
					mutant.Seq[j] = basesLower[i]
					mutant.Seq[jj] = basesLower[ii]

					iter, err = sketches.NewSimHashIterator(mutant, k, m, scale, true, false)
					if err != nil {
						checkError(fmt.Errorf("%s: %s", mutant, err))
					}
					hash, _ = iter.NextSimHash()
					iter.NextSimHash() // evoke sync.Pool.Put

					if hash != hash0 {
						continue
					}

					// fmt.Printf("%s\t%s\t%02d\t%c\t%c\t%d\n", kmer.Seq, mutant.Seq, j, B, b, hash)
					mutants = append(mutants, Mutant{
						Kmer: mutant.Seq,
						Changes: []Change{
							{Pos: j, Base: b},
							{Pos: jj, Base: bb},
						},
					})
				}
			}
		}
	}

	return &KmerResult{
		Pos:     pos,
		Kmer:    kmer.Seq,
		Hash:    hash0,
		Mutants: mutants,
	}
}

func checkKmer(pos int, kmer *seq.Seq, k int, m int, scale int) *KmerResult {
	var iter *sketches.Iterator
	var err error

	iter, err = sketches.NewSimHashIterator(kmer, k, m, scale, true, false)
	if err != nil {
		checkError(fmt.Errorf("%s: %s", kmer, err))
	}
	hash0, _ := iter.NextSimHash()

	var i, j int
	var b, B byte
	var hash uint64
	mutants := make([]Mutant, 0, 8)
	for j, B = range kmer.Seq {
		for i, b = range basesUpper {
			if b == B {
				continue
			}

			mutant := kmer.Clone2()
			mutant.Seq[j] = basesLower[i]

			iter, err = sketches.NewSimHashIterator(mutant, k, m, scale, true, false)
			if err != nil {
				checkError(fmt.Errorf("%s: %s", mutant, err))
			}
			hash, _ = iter.NextSimHash()
			iter.NextSimHash() // evoke sync.Pool.Put

			if hash != hash0 {
				continue
			}

			// fmt.Printf("%s\t%s\t%02d\t%c\t%c\t%d\n", kmer.Seq, mutant.Seq, j, B, b, hash)
			mutants = append(mutants, Mutant{
				Kmer: mutant.Seq,
				Changes: []Change{
					{Pos: j, Base: b},
				},
			})
		}
	}

	return &KmerResult{
		Pos:     pos,
		Kmer:    kmer.Seq,
		Hash:    hash0,
		Mutants: mutants,
	}
}

func checkError(err error) {
	if err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err)
		os.Exit(-1)
	}
}
