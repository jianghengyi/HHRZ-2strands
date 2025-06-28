## Build readsend
- install Golang
- go into Scan_read_ends
- execute: go build

## Finding read bump
```sh
sh batch_fastp_process.sh target_file post_fix
sh batch_align_sort_human.sh prefix_name_file postfix read_count
sh batch_extract_reads.sh target_dir target_file input_bed
sh batch_readbump_test.sh target_dir target_file observing_window_file lib_type
```

## Finding top 3 expression tissues
```sh
ruby top3_expression.rb target_gens.txt profile.tsv sample.tsv
```