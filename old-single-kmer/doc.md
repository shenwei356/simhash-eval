    for i in $(seq 1 21); do \
        for b in a c g t; do \
            seqkit mutate --quiet -p $i:$b t.fna; \
        done; \
    done \
    | seqkit rmdup -s -i \
    > t.fa

    go run single-kmer.go t.fa  21 5 6 > t.txt
