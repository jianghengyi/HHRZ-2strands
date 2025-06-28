package main

import (
	"bufio"
	"io"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

// for observing window record
type Observing_window struct {
	Rname  string // reference genome name, eg. chr1
	Start  int64
	Stop   int64
	Site   int64 // for example cut site, Ligation site
	Strand byte
}

func (w *Observing_window) Equal(x *Observing_window) bool {
	if w.Rname != x.Rname || w.Start != x.Start || w.Stop != x.Stop {
		return false
	}
	if w.Site != x.Site || w.Strand != x.Strand {
		return false
	}
	return true
}

func Digest_observing_window_line(line string, fix bool) *Observing_window {
	fields := strings.Split(line, "\t")
	if len(fields) < 5 {
		return nil
	}
	var result Observing_window
	//

	// todo: some test of the fields

	result.Rname = fields[0]
	result.Start, _ = strconv.ParseInt(fields[1], 0, 64)
	result.Stop, _ = strconv.ParseInt(fields[2], 0, 64)
	result.Site, _ = strconv.ParseInt(fields[3], 0, 64)
	result.Strand = fields[4][0] // bed file with no score
	if result.Strand != '+' && result.Strand != '-' {
		result.Strand = fields[5][0] // normal bed file
	}

	if result.Strand == '-' && fix {
		// fix lgt cut meaning
		result.Site += 1
	}
	return &result
}

func Read_Observing_windows(fh io.Reader, fix bool) []*Observing_window {
	fileScanner := bufio.NewScanner(fh)
	fileScanner.Split(bufio.ScanLines)
	var records []*Observing_window

	for fileScanner.Scan() {
		line := fileScanner.Text()
		record := Digest_observing_window_line(line, fix)
		in_array := false
		for _, r := range records {
			if r.Equal(record) {
				in_array = true
				break
			}
		}

		if record != nil && !in_array {
			records = append(records, record)
		}
	}
	sort.Sort(Window_by_order(records))
	return records
}

type Window_by_order []*Observing_window

func (w Window_by_order) Len() int {
	return len(w)
}

func Extract_number(s string) int {
	num_pattern := regexp.MustCompile(`(\d+)`)
	m := num_pattern.FindStringSubmatch(s)
	if len(m) > 0 {
		digit, err := strconv.Atoi(m[1])
		if err != nil {
			return 0
		}
		return digit
	}
	if strings.Contains(s, "X") || strings.Contains(s, "x") {
		return 999998
	}
	if strings.Contains(s, "Y") || strings.Contains(s, "y") {
		return 999999
	}
	return 0

}

func (w Window_by_order) Less(i, j int) bool {
	if w[i].Rname == w[j].Rname {
		return w[i].Start < w[j].Start
	}
	a := Extract_number(w[i].Rname)
	b := Extract_number(w[j].Rname)
	if a == 0 && b != 0 || a != 0 && b == 0 {
		return w[i].Rname < w[j].Rname
	}
	return a < b
}
func (w Window_by_order) Swap(i, j int) {
	w[i], w[j] = w[j], w[i]
}

func Range_include(a, b, c, d int64) bool {
	//a<b, c<d
	// [a,b] include c , d
	if a <= c && b >= c || a <= d && b >= d {
		return true
	}
	return false
}
