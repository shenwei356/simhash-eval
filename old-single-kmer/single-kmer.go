package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"

	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/bio/sketches"
)

func main() {
	if len(os.Args) != 5 {
		fmt.Printf("usage: %s <file> <k> <m> <scale>\n", filepath.Base(os.Args[0]))
		os.Exit(0)
	}

	file := os.Args[1]
	k, err := strconv.Atoi(os.Args[2])
	if err != nil || k <= 0 {
		checkError(fmt.Errorf("<k> should be positive interger: %s", os.Args[2]))
	}
	m, err := strconv.Atoi(os.Args[3])
	if err != nil || m <= 0 {
		checkError(fmt.Errorf("<m> should be positive interger: %s", os.Args[3]))
	}
	s, err := strconv.Atoi(os.Args[4])
	if err != nil || m <= 0 {
		checkError(fmt.Errorf("<s> should be positive interger: %s", os.Args[4]))
	}

	fmt.Printf("file: %s\nk-mer: %d\nm-mer: %d\nscale: %d\n\n", file, k, m, s)

	fastxReader, err := fastx.NewDefaultReader(file)
	if err != nil {
		checkError(fmt.Errorf("%s: %s", file, err))
	}

	var code uint64
	var ok bool
	var iter *sketches.Iterator
	var record *fastx.Record
	var idx int

	data := make(map[uint64]map[string]interface{})

	for {
		record, err = fastxReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(err)
			break
		}

		iter, err = sketches.NewSimHashIterator(record.Seq, k, m, s, true, false)
		if err != nil {
			checkError(fmt.Errorf("%s: %s", record.Name, err))
		}

		for {
			code, ok = iter.NextSimHash()
			if !ok {
				break
			}

			idx = iter.Index()
			if _, ok = data[code]; !ok {
				data[code] = map[string]interface{}{string(record.Seq.Seq[idx : idx+k]): struct{}{}}
			} else {
				data[code][string(record.Seq.Seq[idx:idx+k])] = struct{}{}
			}
		}
	}

	outfh := bufio.NewWriter(os.Stdout)
	for hash, seqs := range data {
		if len(seqs) == 1 {
			continue
		}

		fmt.Fprintf(outfh, ">%d %d\n", hash, len(seqs))
		for s := range seqs {
			fmt.Fprintf(outfh, "%s\n", s)
		}
		fmt.Fprintf(outfh, "\n\n")
	}
	outfh.Flush()
}

func checkError(err error) {
	if err != nil {
		fmt.Printf("%s\n", err)
		os.Exit(-1)
	}
}
