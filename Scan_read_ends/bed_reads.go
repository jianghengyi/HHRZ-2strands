package main

import (
	"bufio"
	"errors"
	"io"
	"os"
	"strconv"
	"strings"
)

// example record
// 1	64886274	64966365	SRR15503362.42773667/1	255	+
// 1	64886276	64966367	SRR15503362.11705087/2	255	-

type Bed_read struct {
	Rname         string
	Start         int64
	Stop          int64
	Name          string
	Score         int
	Strand        byte
	Signiture     string
	Signiture_rna string
	Cartoon       string
}

func (r *Bed_read) Rname_() string {
	if r == nil {
		return "NA"
	}
	return r.Rname
}
func (r *Bed_read) Seq_() string {
	return "."
}
func (r *Bed_read) CIGAR_() string {
	return "."
}
func (r *Bed_read) Signiture_() string {
	return r.Signiture
}
func (r *Bed_read) Signiture_rna_() string {
	return r.Signiture_rna
}
func (r *Bed_read) Signiture_set(sg string) {
	r.Signiture = sg
}
func (r *Bed_read) Signiture_rna_set(sg string) {
	r.Signiture_rna = sg
}

func (r *Bed_read) Cartoon_() string{
	return r.Cartoon
}

func (r *Bed_read) Cartoon_set(ct string) {
	r.Cartoon=ct
}


func (r *Bed_read) Name_() string {
	if r == nil {
		return "NA"
	}
	return r.Name
}

func (r *Bed_read) Pos_() int64 {
	if r == nil {
		return -1
	}
	return r.Start
}
func (r *Bed_read) Start_() int64 {
	if r == nil {
		return -1
	}
	return r.Start
}

func (r *Bed_read) Stop_() int64 {
	if r == nil {
		return -1
	}
	return r.Stop
}

func (r *Bed_read) Strand_mapped() byte {
	if r == nil {
		return '?'
	}
	return r.Strand
}

func Digest_to_bed_read(line string) *Bed_read {
	fields := strings.Split(line, "\t")
	if len(fields) < 5 {
		return nil
	}
	record := new(Bed_read)
	record.Rname = fields[0]
	record.Start, _ = strconv.ParseInt(fields[1], 10, 64)
	record.Start += 1 // in bed file coordinate start from 0
	record.Stop, _ = strconv.ParseInt(fields[2], 10, 64)
	// in bed file stop is not included in range, offseting the 0 based coordinate
	record.Name = fields[3]
	record.Score, _ = strconv.Atoi(fields[4])
	record.Strand = fields[5][0]
	return record
}

func (r *Bed_read) Source() int {
	if r == nil {
		return 0
	}
	tt := strings.Split(r.Name, "/")
	if len(tt) < 2 {
		return 0
	}
	v, _ := strconv.Atoi(tt[1])
	return v
}

func (r *Bed_read) Valid_end_pos() int64 {
	// strand ='+'
	// +-----

	if r == nil {
		return -1
	}
	if r.Strand == '+' {
		return r.Start
	} else if r.Strand == '-' {
		return r.Stop
	}
	return -1
}

func (r *Bed_read) Seq_name() string {
	if r == nil {
		return "NA"
	}
	tt := strings.Split(r.Name, "/")
	if len(tt) < 2 {
		return "NA"
	}
	return tt[0]
}

func Read_bed_as_reads(fh_bed io.Reader) []Abstract_read {
	fileScanner := bufio.NewScanner(fh_bed)
	fileScanner.Split(bufio.ScanLines)
	var records []Abstract_read
	for fileScanner.Scan() {
		line := fileScanner.Text()
		if len(line) > 0 {
			if line[0] == '@' {
				continue
			}
			if line[0] == '#' {
				continue
			}
		}

		record := Digest_to_bed_read(line)
		if record != nil {
			records = append(records, Abstract_read(record))
		}
	}
	return records
}

func Read_bed_file(file_name string) ([]Abstract_read, error) {
	fh_bed, err := os.Open(file_name)
	if err != nil {
		return nil, err
	}
	defer fh_bed.Close()
	collection := Read_bed_as_reads(fh_bed)
	if len(collection) == 0 {
		return nil, errors.New("Bed file:" + file_name + " contains not valid line")
	}
	return collection, nil
}