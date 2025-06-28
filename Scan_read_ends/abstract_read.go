package main

import (
	"fmt"
	"regexp"
	"strings"
)

type Abstract_read interface {
	Source() int
	Strand_mapped() byte
	Valid_end_pos() int64
	Pos_() int64
	Rname_() string
	Name_() string
	Start_() int64
	Stop_() int64
	Seq_() string
	CIGAR_() string
	Signiture_() string
	Signiture_rna_() string
	Signiture_set(string)
	Signiture_rna_set(string)
	Cartoon_() string
	Cartoon_set(string)
}

func Represent(r Abstract_read, lib_type int) byte {
	read_source := r.Source()
	if lib_type == 1 { // RNA on the read2, cDNA on read1
		if read_source == 1 {
			return 'c'
		} else if read_source == 2 {
			return 'r'
		}
	} else if lib_type == 2 { // RNA on the read1, cDNA on read2
		if read_source == 1 {
			return 'r'
		} else if read_source == 2 {
			return 'c'
		}
	}
	return '?'
}

func Infer_break_site(r Abstract_read, lib_type int) *Break_site {

	read_represent := Represent(r, lib_type)
	if read_represent == '?' {
		return nil
	}
	bs1 := new(Break_site)
	bs1.Rname = r.Rname_()
	if r.Strand_mapped() == '+' { // read mapping to the + strand
		if read_represent == 'r' { // read representing RNA strand
			bs1.Strand = '+'
			// test CIGAR missing
			bs1.Pos = r.Pos_() - 1 // Pos-1, one base before 5'
			bs1.Cartoon = "|>>>>"  // the RNA is on the plus trand, this read is included because RNA 5 end
			bs1.Signiture = "^" + r.Seq_()[0:6]
			bs1.Signiture_rna = "|" + r.Seq_()[0:6]
		} else if read_represent == 'c' {

			// test CIGAR missing, and can not be tested,
			// because no reference to it pair
			bs1.Strand = '-'
			bs1.Pos = r.Pos_() // pos
			bs1.Cartoon = "|<<<<"
			read_seq := r.Seq_()[0:6]
			bs1.Signiture = "^" + read_seq
			bs1.Signiture_rna = Reverse_complement(&read_seq) + "|"
		}

	} else if r.Strand_mapped() == '-' {
		// this is quite over-simpified
		if read_represent == 'r' {
			bs1.Strand = '-'
			bs1.Pos = r.Stop_() + 1
			bs1.Cartoon = "<<<<|"
			read_seq := r.Seq_()[len(r.Seq_())-6 : len(r.Seq_())]
			bs1.Signiture = read_seq + "$"
			bs1.Signiture_rna = "|" + Reverse_complement(&read_seq)
		} else if read_represent == 'c' {
			bs1.Strand = '+'
			bs1.Pos = r.Stop_()
			bs1.Cartoon = ">>>>|"
			read_seq := r.Seq_()[len(r.Seq_())-6 : len(r.Seq_())]
			bs1.Signiture = read_seq + "$"
			bs1.Signiture_rna = read_seq + "|"
		}
	}
	return bs1
}

func Infer_break_site_unstranded(r Abstract_read, gene_strand byte) []*Break_site {

	var result []*Break_site
	if r.Strand_mapped() == '+' {
		switch gene_strand {
		case '+':
			bs_temp := new(Break_site)
			bs_temp.Rname = r.Rname_()
			bs_temp.Strand = '+'
			bs_temp.Pos = r.Start_() - 1
			read_seq := r.Seq_()[0:6]
			bs_temp.Signiture = "^" + read_seq
			bs_temp.Signiture_rna = "|" + read_seq
			bs_temp.Cartoon = "|>>>>"
			result = append(result, bs_temp)
		case '-':
			bs_temp := new(Break_site)
			bs_temp.Rname = r.Rname_()
			bs_temp.Strand = '-'
			bs_temp.Pos = r.Start_()
			read_seq := r.Seq_()[0:6]
			bs_temp.Signiture = "^" + read_seq
			bs_temp.Cartoon = "|<<<<"
			bs_temp.Signiture_rna = Reverse_complement(&read_seq) + "|"
			result = append(result, bs_temp)
		default:
			bs1 := new(Break_site)
			bs2 := new(Break_site)
			bs1.Rname = r.Rname_()
			bs2.Rname = r.Rname_()
			bs1.Pos = r.Start_() - 1
			bs1.Strand = '+'
			read_seq := r.Seq_()[0:6]
			bs1.Signiture = "^" + read_seq
			bs2.Pos = r.Start_()
			bs2.Strand = '-'
			bs2.Signiture = "^" + read_seq
			result = []*Break_site{bs1, bs2}
		}

	} else if r.Strand_mapped() == '-' {
		switch gene_strand {
		case '+':
			bs_temp := new(Break_site)
			bs_temp.Rname = r.Rname_()
			bs_temp.Strand = '+'
			bs_temp.Pos = r.Stop_()
			read_seq := r.Seq_()[len(r.Seq_())-6 : len(r.Seq_())]
			bs_temp.Signiture = read_seq + "$"
			bs_temp.Signiture_rna = Reverse_complement(&read_seq) + "|"
			bs_temp.Cartoon = ">>>>|"
			result = append(result, bs_temp)
		case '-':
			bs_temp := new(Break_site)
			bs_temp.Rname = r.Rname_()
			bs_temp.Strand = '-'
			bs_temp.Pos = r.Stop_() + 1
			read_seq := r.Seq_()[len(r.Seq_())-6 : len(r.Seq_())]
			bs_temp.Signiture = read_seq + "$"
			bs_temp.Signiture_rna = "|" + Reverse_complement(&read_seq)
			bs_temp.Cartoon = "<<<<|"
			result = append(result, bs_temp)

		default:
			bs1 := new(Break_site)
			bs2 := new(Break_site)
			bs1.Pos = r.Stop_() + 1
			bs1.Strand = '-'
			read_seq := r.Seq_()[len(r.Seq_())-6 : len(r.Seq_())] // seq in on the same in the genome
			bs1.Signiture = read_seq + "$"
			bs2.Pos = r.Stop_()
			bs2.Strand = '+'
			bs2.Signiture = read_seq + "$"
			result = []*Break_site{bs1, bs2}
		}
	}
	return result
}

func Check_break_site(reads []Abstract_read, window *Observing_window, lib_type int, collect_reads bool) ([]Abstract_read, []Abstract_read, int, int, int, int, int) {
	var selected_reads []Abstract_read
	var minus_bumped_reads []Abstract_read
	all := 0
	passed := 0
	covered := 0
	bumped := 0
	bumped_minus := 0
	// var strand_mapped byte

	for _, read := range reads {
		if read.Rname_() != window.Rname {
			continue
		}

		if read.Stop_() < window.Start || read.Start_() > window.Stop {
			continue
		}
		all += 1

		if read.Start_() < window.Site && read.Stop_() > window.Site {
			covered += 1
		}

		if read.Start_() >= window.Start && read.Start_() <= window.Stop || read.Stop_() >= window.Start && read.Stop_() <= window.Stop {
			passed += 1
		}

		candi := Infer_break_site(read, lib_type)
		if candi != nil {
			if candi.Pos == window.Site {
				if candi.Strand != window.Strand {
					bumped_minus += 1
					if collect_reads{
						minus_bumped_reads = append(minus_bumped_reads, read)
					}					
					continue
				}
				bumped += 1
				read.Signiture_set(candi.Signiture)
				read.Signiture_rna_set(candi.Signiture_rna)
				read.Cartoon_set(candi.Cartoon)
				if collect_reads{
					selected_reads = append(selected_reads, read)
				}
				
			}
		}
	}
	return selected_reads, minus_bumped_reads, all, covered, passed, bumped, bumped_minus
}

func Scan_break_site(reads []Abstract_read, window *Observing_window, lib_type int) (map[string]*Break_site_stats, int) {
	all := 0
	var result = make(map[string]*Break_site_stats)
	for _, read := range reads {
		if read.Rname_() != window.Rname {
			continue
		}
		// -------------------------------------------
		//           |                 |
		//    ------                     ____

		if read.Stop_() < window.Start || read.Start_() > window.Stop {
			// fmt.Println(window.Start, window.Stop, read.Start_(), read.Stop_())
			continue
		}
		all += 1

		candi := Infer_break_site(read, lib_type)
		if candi != nil {
			stat, ok := result[candi.Str()]
			if ok {
				stat.Add_read(read)
			} else {
				stat_new := new(Break_site_stats)
				stat_new.Set_Break_site(candi)
				stat_new.Add_read(read)
				result[candi.Str()] = stat_new
			}
		}
	}
	return result, all
}

func Count_reads(reads []Abstract_read, window *Observing_window) (int, int, int) {
	all := 0
	passed := 0
	covered := 0
	for _, read := range reads {
		if read.Rname_() != window.Rname {
			continue
		}

		if read.Stop_() < window.Start || read.Start_() > window.Stop {
			continue
		}
		all += 1

		if read.Start_() < window.Site && read.Stop_() > window.Site {
			covered += 1
		}

		if read.Start_() >= window.Start && read.Start_() <= window.Stop || read.Stop_() >= window.Start && read.Stop_() <= window.Stop {
			passed += 1
		}
	}
	return all, covered, passed
}

func Reads_to_lines(reads []Abstract_read) []string {
	var result []string
	front_pattern := regexp.MustCompile(`^\d+M`)
	tail_pattern := regexp.MustCompile(`\d+M$`)
	for _, read := range reads {
		var cigar_pass string
		if strings.HasPrefix(read.Signiture_(), "^") {
			// regexp.Match("H.* ", "Hello World!")
			if front_pattern.MatchString(read.CIGAR_()) {
				cigar_pass = "PASS"
			} else {
				cigar_pass = "FAIL"
			}
		} else if strings.HasSuffix(read.Signiture_(), "$") {
			if tail_pattern.MatchString(read.CIGAR_()) {
				cigar_pass = "PASS"
			} else {
				cigar_pass = "FAIL"
			}
		}
		r := fmt.Sprintf("\t%s/%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s", read.Name_(), read.Source(),
			read.Seq_(), read.CIGAR_(),
			read.Start_(), read.Stop_(),
			string(read.Strand_mapped()),
			read.Signiture_(),
			read.Signiture_rna_(),
			read.Cartoon_(),
			cigar_pass,
		)
		result = append(result, r)
	}
	return result
}

// for cut site
type Break_site struct {
	Rname         string
	Strand        byte
	Pos           int64
	Cartoon       string
	Signiture     string
	Signiture_rna string
}

func (bs *Break_site) Str() string {
	if bs == nil {
		return "NA"
	}
	s := fmt.Sprintf("%s\t%d\t%s", bs.Rname, bs.Pos, string(bs.Strand))
	return s
}

type Break_site_stats struct {
	Site  *Break_site
	Reads []Abstract_read
}

func (collection *Break_site_stats) Set_Break_site(bs *Break_site) {
	collection.Site = bs
}

func (collection *Break_site_stats) Add_read(read Abstract_read) {
	collection.Reads = append(collection.Reads, read)
}

func (collection *Break_site_stats) Read_count() int {
	return len(collection.Reads)
}

// for unstranded library construction
type Frag_end struct {
	Pos  int64
	Type int // 3 or 5
}

func Frag_end_equal(a *Frag_end, b *Frag_end) bool {
	if a == nil || b == nil {
		return false
	}
	if a.Pos == b.Pos && a.Type == b.Type {
		return true
	}
	return false
}

func Read_to_frag_end(r Abstract_read) *Frag_end {
	if r == nil {
		return nil
	}
	if r.Strand_mapped() == '+' {
		return &Frag_end{Pos: r.Start_(), Type: 5}
	} else if r.Strand_mapped() == '-' {
		return &Frag_end{Pos: r.Stop_(), Type: 3}
	}
	return nil
}

func Cut_site_to_frag_ends(window *Observing_window) []*Frag_end {
	if window.Strand == '+' {
		frag_end1 := Frag_end{Pos: window.Site, Type: 3}
		frag_end2 := Frag_end{Pos: window.Site + 1, Type: 5}
		return []*Frag_end{&frag_end1, &frag_end2}
	} else if window.Strand == '-' {
		frag_end1 := Frag_end{Pos: window.Site, Type: 3}
		frag_end2 := Frag_end{Pos: window.Site - 1, Type: 5}
		return []*Frag_end{&frag_end1, &frag_end2}
	}
	return nil
}

func Unstranded_read_check_break_site(reads []Abstract_read, window *Observing_window, ob_window_is_on_gene bool, collect_reads bool) ([]Abstract_read, []Abstract_read, int, int, int, int, int) {
	// end_sites := Cut_site_to_frag_ends(window)
	var selected_reads []Abstract_read
	var minus_bumped_reads []Abstract_read
	var all int
	var passed int
	var covered int
	var bumped int
	var minus_bumps int

	for _, read := range reads {
		if read.Rname_() != window.Rname {
			continue
		}

		if read.Stop_() < window.Start || read.Start_() > window.Stop {
			// fmt.Println(window.Start, window.Stop, read.Start_(), read.Stop_())
			continue
		}
		all += 1

		if read.Start_() < window.Site && read.Stop_() > window.Site {
			covered += 1

		}

		if read.Start_() >= window.Start && read.Start_() <= window.Stop || read.Stop_() >= window.Start && read.Stop_() <= window.Stop {
			passed += 1
		}

		var test_strand byte
		test_strand = '?'
		if ob_window_is_on_gene {
			test_strand = window.Strand
		}

		bss := Infer_break_site_unstranded(read, test_strand)

		for _, bs := range bss {
			if bs.Pos == window.Site {
				if bs.Strand == window.Strand {
					read.Signiture_set(bs.Signiture)
					if collect_reads{
						selected_reads = append(selected_reads, read)
					}
					bumped += 1
				} else {
					minus_bumps += 1
					if collect_reads{
						minus_bumped_reads = append(minus_bumped_reads, read)
					}
				}
			}
		}

	}
	return selected_reads, minus_bumped_reads, all, covered, passed, bumped, minus_bumps
}

func complement(s byte) byte {
	switch {
	case s == 'C', s == 'c':
		return 'G'
	case s == 'T', s == 't':
		return 'A'
	case s == 'A', s == 'a':
		return 'T'
	case s == 'G', s == 'g':
		return 'C'
	case s == 'N', s == 'n':
		return 'N'
	default:
		return 'X'
	}
}

func Reverse_complement(s *string) string {
	r := []byte(*s)
	for i, j := 0, len(r)-1; i < len(r)/2; i, j = i+1, j-1 {

		r[i], r[j] = complement(r[j]), complement(r[i])
	}
	if len(r)%2 == 1 {
		t := r[len(r)/2]
		r[len(r)/2] = complement(t)
	}
	return string(r)
}
