package main

import (
	"bufio"
	"errors"
	"io"
	"os"
	"regexp"
	"strconv"
	"strings"
)

/*
	sample record
SRR13367164.3129510	99	NC_000002.12	196	40	4M2I70M	=	226	106
GAAGAAAGAGAGCTTGGATTTTCAGAGTGCAATGGGGAAGAGTACAGAAGCAAGTGGAAACTGAGCCTCATGGAAA
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEE
AS:i:-21	XN:i:0	XM:i:2	XO:i:1	XG:i:2	NM:i:4	MD:Z:3A3A66	YS:i:-5	YT:Z:CP
*/

type Read struct {
	Qname         string
	Flag          int
	Rname         string
	Pos           int64
	Mapq          int
	Cigar         string
	Rnext         string
	Pnext         int64
	Tlen          int
	Seq           string
	Qual          string
	Strand        byte
	Signiture     string
	Signiture_rna string
	Cartoon       string
}

func (r *Read) Rname_() string {
	return r.Rname
}
func (r *Read) Pos_() int64 {
	return r.Pos
}
func (r *Read) Signiture_() string {
	return r.Signiture
}
func (r *Read) Signiture_set(sg string) {
	r.Signiture = sg
}
func (r *Read) Signiture_rna_set(sg string) {
	r.Signiture_rna = sg
}

func (r *Read) Signiture_rna_() string {
	return r.Signiture_rna
}
func (r *Read) Cartoon_() string {
	return r.Cartoon
}

func (r *Read) Cartoon_set(ct string) {
	r.Cartoon = ct
}

func (r *Read) Name_() string {
	return r.Qname
}

func (r *Read) Seq_() string {
	return r.Seq
}

func (r *Read) CIGAR_() string {
	return r.Cigar
}
func (r *Read) Start_() int64 {
	return r.Pos
}

func (r *Read) Stop_() int64 {
	// here should be changed later
	cigar_length, err := Cigar_to_genomic_length(r.Cigar)
	if err != nil {
		return r.Pos + int64(r.Tlen) - 1
	}
	return r.Pos + int64(cigar_length) - 1
}

func (r *Read) Strand_mapped() byte {
	if r == nil {
		return '?'
	}
	if r.Flag&0x10 != 0 {
		return '-'
	}
	return '+'
}

func (r *Read) Source() int {

	if r.Flag&0x40 != 0 {
		return 1
	} else if r.Flag&0x80 != 0 {
		return 2
	}
	return 0
}

func (r *Read) Valid_end_pos() int64 {
	if r == nil {
		return 0
	}
	strand := r.Strand_mapped()
	switch strand {
	case '+':
		return r.Pos
	case '-':
		return r.Stop_()
	default:
		return r.Pos
	}
}

func (r *Read) Represent(lib_type int) byte {
	return Represent(r, lib_type)
}

func (r *Read) Infer_break_site(lib_type int) *Break_site {
	return Infer_break_site(r, lib_type)
}

// Calculate the genomic length from a CIGAR string
func Cigar_to_genomic_length(cigar string) (int, error) {
	// Define regex to match CIGAR operations
	re := regexp.MustCompile(`(\d+)([MIDNSHP=X])`)

	// Find all matches
	matches := re.FindAllStringSubmatch(cigar, -1)

	// Initialize length
	genomicLength := 0

	// Iterate over matches and calculate length
	for _, match := range matches {
		// Extract the length and operation
		length, err := strconv.Atoi(match[1])
		if err != nil {
			return 0, err
		}
		op := match[2]
		// Only add length for operations that consume reference sequence
		switch op {
		case "M", "D", "N":
			genomicLength += length
		}
	}
	return genomicLength, nil
}

func Digest_sam_line(line string, store_seq bool) *Read {
	var record = new(Read)
	fields := strings.Split(line, "\t")
	if len(fields) < 10 {
		return nil
	}

	record.Qname = fields[0]
	record.Flag, _ = strconv.Atoi(fields[1])
	record.Rname = fields[2]
	record.Pos, _ = strconv.ParseInt(fields[3], 0, 64)
	record.Mapq, _ = strconv.Atoi(fields[4])
	record.Cigar = fields[5]
	record.Rnext = fields[6]
	record.Pnext, _ = strconv.ParseInt(fields[7], 0, 64)
	record.Tlen, _ = strconv.Atoi(fields[8])
	if store_seq {
		record.Seq = fields[9]
		record.Qual = fields[10]
	}
	return record
}

func Read_sam_as_reads(fh_sam io.Reader) []Abstract_read {
	fileScanner := bufio.NewScanner(fh_sam)
	fileScanner.Split(bufio.ScanLines)
	var records []Abstract_read
	for fileScanner.Scan() {
		line := fileScanner.Text()
		if len(line) > 0 {
			if line[0] == '@' {
				continue
			}
		}

		record := Digest_sam_line(line, true)
		if record != nil {
			records = append(records, Abstract_read(record))
		}
	}
	return records
}

func Read_sam_file(file_name string) ([]Abstract_read, error) {
	fh_sam, err := os.Open(file_name)
	if err != nil {
		return nil, err
	}
	defer fh_sam.Close()
	collection := Read_sam_as_reads(fh_sam)
	if len(collection) == 0 {
		return nil, errors.New("SAM file:" + file_name + " contains not valid line")
	}
	return collection, nil
}

func Pack_reads_into_pairs(collection []*Read) []*Read_pair {
	pairs := make(map[string]*Read_pair)
	name_order := []int{}
	for i, read := range collection {
		pair, ok := pairs[read.Qname]
		n := read.Source()
		if !ok {
			var new_pair Read_pair

			if n == 1 {
				new_pair.Name = read.Qname
				new_pair.Rname = read.Rname
				new_pair.Read1 = read
				new_pair.Read2 = nil
			} else if n == 2 {
				new_pair.Read1 = nil
				new_pair.Read2 = read
			}
			pairs[read.Qname] = &new_pair
			name_order = append(name_order, i)
		} else {
			if n == 1 {
				pair.Read1 = read
			} else if n == 2 {
				pair.Read2 = read
			}
		}
	}
	var result []*Read_pair
	for _, i := range name_order {
		pair := pairs[collection[i].Qname]
		result = append(result, pair)
	}
	return result
}

// read pair
type Read_pair struct {
	Name  string
	Read1 *Read
	Read2 *Read
	Rname string // reference genome name
}

func (pair *Read_pair) Start() int64 {
	var t int64
	if pair.Read1 != nil {
		t = pair.Read1.Valid_end_pos()
	}
	if pair.Read2 != nil {
		tt := pair.Read2.Valid_end_pos()
		if tt < t {
			return tt
		} else {
			return t
		}
	}
	return -1
}

func (pair *Read_pair) Stop() int64 {
	var t int64
	if pair.Read1 != nil {
		t = pair.Read1.Valid_end_pos()
	}
	if pair.Read2 != nil {
		tt := pair.Read2.Valid_end_pos()
		if tt > t {
			return tt
		} else {
			return t
		}
	}
	return -1
}

func (pair *Read_pair) Strand_mapped(lib_type int) byte {
	// mean RNA is on which Strand
	if lib_type == 1 {
		if pair.Read2 != nil {
			return pair.Read2.Strand_mapped()
		} else if pair.Read1 != nil {
			tmp_strand := pair.Read1.Strand_mapped()
			if tmp_strand == '+' {
				return '-'
			} else if tmp_strand == '-' {
				return '+'
			} else {
				return '?'
			}
		}
	} else if lib_type == 2 {
		if pair.Read1 != nil {
			return pair.Read1.Strand_mapped()
		} else if pair.Read2 != nil {
			tmp_strand := pair.Read2.Strand_mapped()
			if tmp_strand == '+' {
				return '-'
			} else if tmp_strand == '-' {
				return '+'
			} else {
				return '?'
			}
		}

	}
	return '?'
}

func (pair *Read_pair) To_line() string {
	var builder strings.Builder
	builder.WriteString(pair.Name + "\t")
	if pair.Read1 != nil {
		builder.WriteString(pair.Read1.Seq + "\t")
		builder.WriteString(string(pair.Read1.Strand_mapped()) + "\t")
		builder.WriteString(strconv.Itoa(pair.Read1.Source()) + "\t")
		builder.WriteString(strconv.FormatInt(pair.Read1.Pos, 10) + "\t")
	} else {
		builder.WriteString("\t\t\t\t")
	}

	if pair.Read2 != nil {
		builder.WriteString(pair.Read2.Seq + "\t")
		builder.WriteString(string(pair.Read2.Strand_mapped()) + "\t")
		builder.WriteString(strconv.Itoa(pair.Read2.Source()) + "\t")
		builder.WriteString(strconv.FormatInt(pair.Read2.Pos, 10) + "\t")
	} else {
		builder.WriteString("\t\t\t\t")
	}
	builder.WriteString(strconv.FormatInt(pair.Start(), 10) + "\t")
	builder.WriteString(strconv.FormatInt(pair.Stop(), 10))

	return builder.String()
}

func (pair *Read_pair) Infer_break_sites(lib_type int) []*Break_site {
	// only for type 1 and type 2 library constructed
	// when we talk about RNA break site, it means the RNA breaks "after" that site
	// and we know RNA direction is 5' -> 3'
	// when A read is mapped to the plus(+) Strand: the underlying possible break site can be
	//    1) one base before the 5' end( read.Start-1)
	//    2) the base of the '3 end(read.Stop)
	// when A read is mapped to the minus(-) Strand: the underlying break possible break site can be
	//    1) one base before the 5' end( read.Stop+1)
	//    2) the base of the 3'end( read.Start)
	var result []*Break_site
	if pair.Read1 != nil && pair.Read2 != nil {
		bs1 := new(Break_site)
		bs2 := new(Break_site)
		bs1.Rname = pair.Rname
		bs2.Rname = pair.Rname
		if pair.Strand_mapped(lib_type) == '+' {
			//             v         v
			//      ------->--------->---
			// -------------++++++++++--------
			// -------------------------------

			bs1.Strand = '+'
			bs1.Pos = pair.Start() - 1
			bs2.Strand = '+'
			bs2.Pos = pair.Stop()

		} else if pair.Strand_mapped(lib_type) == '-' {
			// -------------++++++++++--------
			// -------------------------------
			//            ------<----<----
			//                  ^     ^
			bs1.Strand = '-'
			bs1.Pos = pair.Start()
			bs2.Strand = '-'
			bs2.Pos = pair.Stop() + 1
		}
		result = append(result, bs1)
		result = append(result, bs2)
	} else {
		if pair.Read1 != nil {
			bs1 := pair.Read1.Infer_break_site(lib_type)
			if bs1 != nil {
				result = append(result, bs1)
			}
		} else if pair.Read2 != nil {
			bs1 := pair.Read2.Infer_break_site(lib_type)
			if bs1 != nil {
				result = append(result, bs1)
			}
		}
	}
	return result
}

func Select_paired_Reads(paired_reads []*Read_pair, window *Observing_window) []*Read_pair {
	var selected_reads []*Read_pair
	for _, pair := range paired_reads {
		if pair.Rname != window.Rname {
			continue
		}
		if Range_include(window.Start, window.Stop, pair.Start(), pair.Stop()) {
			selected_reads = append(selected_reads, pair)
		}
	}
	return selected_reads
}

func Select_reads(reads []*Read, window *Observing_window) []*Read {
	var selected_reads []*Read
	for _, read := range reads {
		if Range_include(window.Start, window.Stop, read.Valid_end_pos(), read.Valid_end_pos()) {
			selected_reads = append(selected_reads, read)
		}
	}
	return selected_reads
}

func Check_paired_Reads(paired_reads []*Read_pair, window *Observing_window, lib_type int) ([]*Read_pair, int, int) { // passed, pumps
	// lib_type ==1 ==> RNA=read2
	// lib_type ==2 ==> RNA=read1
	// lib_type ==0 ==> unstranded
	var selected_reads []*Read_pair
	passed := 0
	bumped := 0
	// var strand_mapped byte
	for _, pair := range paired_reads {
		if pair.Rname != window.Rname {
			continue
		}

		// if lib_type == 1 {
		// 	strand_mapped = pair.Read2.Strand_mapped()
		// } else if lib_type == 2 {
		// 	strand_mapped = pair.Read1.Strand_mapped()
		// } else {
		// 	strand_mapped = window.Strand
		// }
		// if strand_mapped == window.Strand {
		// 	// the read of the RNA is one the same strand as observing windows
		// 	if pair.Start()-1 == window.Site || pair.Stop() == window.Site {
		// 		bumped += 1
		// 	} else if pair.Start() < window.Site && pair.Stop() > window.Site {
		// 		passed += 1
		// 		selected_reads = append(selected_reads, pair)
		// 	}
		// }

		if pair.Start() < window.Site && pair.Stop() > window.Site {
			passed += 1
			continue
		}

		candidates := pair.Infer_break_sites(lib_type)
		for _, candi := range candidates {
			if candi.Strand != window.Strand {
				continue
			}
			if candi.Pos == window.Site {
				bumped += 1
				selected_reads = append(selected_reads, pair)
			}
		}

	}

	return selected_reads, passed, bumped
}

func Scan_paired_reads_for_break(paired_reads []*Read_pair, lib_type int, threshold int) string {
	// when we talk about RNA break site, it means the RNA breaks "after" that site
	// and we know RNA direction is 5' -> 3'
	// when A read is mapped to the plus(+) Strand: the underlying possible break site can be
	//    1) one base before the 5' end( read.Start-1)
	//    2) the base of the '3 end(read.Stop)
	// when A read is mapped to the minus(-) Strand: the underlying break possible break site can be
	//    1) one base before the 5' end( read.Stop+1)
	//    2) the base of the 3'end( read.Start)

	var counter = make(map[string]int)
	// var candidates []*Break_site
	for _, pair := range paired_reads {
		candi := pair.Infer_break_sites(lib_type)
		for _, cc := range candi {
			k := cc.Rname + "\t" + strconv.FormatInt(cc.Pos, 10) + "\t" + string(cc.Strand)
			_, ok := counter[k]
			if !ok {
				counter[k] = 1
			} else {
				counter[k] += 1
			}
		}
	}
	var builder strings.Builder
	for k, v := range counter {
		if v >= threshold {
			builder.WriteString(k + "\t" + strconv.Itoa(v))
			builder.WriteString("\n")
		}
	}
	return builder.String()
}
