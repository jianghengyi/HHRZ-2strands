package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"sort"
	"strings"
	"sync"
	"time"
)

var target_file = flag.String("f", "", "the target sam file, which is exported from bam file")
var target_type = flag.String("t", "sam", "the target file type, which can be sam or bed")
var window_file = flag.String("w", "", "the bed file that describe the observing window")
var fix_ob_window = flag.Bool("b", false, "fix observing window notion")
var lib_type = flag.Int("l", -1, "the library type, 1 for RNA=read2, 2 for RNA=read1, 0 for unstranded")
var report_threshold = flag.Int("s", 0, "threshold to report for calc command")
var repress_reads_printing = flag.Bool("r", false, "repress hit reads printing, default false")
var minus_bumps_printing = flag.Bool("n", false, "to print minus bump reads, default false")
var cmd = flag.String("c", "test", "command: can be test and scan and print")
var thread_n = flag.Int("j", 2, "threads to open when running test")
var cpu_cores = flag.Int("u", 2, "CPU core to use to runing test")
var force_ob_window_on_gene_strand = flag.Bool("g", false, "when the library type is 0, assume observing windows strand is the gene strand")

type Test_job struct {
	Lib_type   int
	Collection []Abstract_read
	Ob_window  *Observing_window
}

type line_window_order []string

func (w line_window_order) Len() int {
	return len(w)
}
func (w line_window_order) Less(i, j int) bool {
	cells_i := strings.Split(w[i], "\t")
	cells_j := strings.Split(w[j], "\t")
	if len(cells_i) < 4 || len(cells_j) < 4 {
		return w[i] < w[j]
	}
	a := Extract_number(cells_i[0])
	b := Extract_number(cells_j[0])

	if a == 0 && b != 0 || a != 0 && b == 0 {
		// when Rname =X or Rname = Y
		// compare literal string value
		return cells_i[0] < cells_j[0]
	}

	if a == b {
		// no X or Y chromo compare chromo
		// when Rname is eaqual compare Start
		c := Extract_number(cells_i[1])
		d := Extract_number(cells_j[1])
		return c < d
	}
	return a < b
}
func (w line_window_order) Swap(i, j int) {
	w[i], w[j] = w[j], w[i]
}

func test_worker(jobs <-chan *Test_job, results chan<- string, wg *sync.WaitGroup, report_threshold int, collect_reads bool) {
	wg.Add(1)
	defer wg.Done()
	window_printed := false
	var output_builder strings.Builder
	var selected_reads, minus_bumps []Abstract_read
	var all, covered, passed, bumped, bumped_minus int
	for j := range jobs {
		output_builder.Reset()
		if j.Lib_type == 0 {
			selected_reads, minus_bumps, all, covered, passed, bumped, bumped_minus = Unstranded_read_check_break_site(j.Collection, j.Ob_window, *force_ob_window_on_gene_strand, collect_reads)
		} else {
			selected_reads, minus_bumps, all, covered, passed, bumped, bumped_minus = Check_break_site(j.Collection, j.Ob_window, j.Lib_type, collect_reads)
		}
		observing_window := j.Ob_window

		if len(selected_reads) > report_threshold {
			window_printed = true
			fmt.Fprintf(&output_builder, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", observing_window.Rname, observing_window.Start,
				observing_window.Stop, observing_window.Site, string(observing_window.Strand), all, covered, passed, bumped, bumped_minus)
			if *repress_reads_printing {
				results <- output_builder.String()
				continue
			}
			if len(selected_reads) > 0 {
				read_lines := Reads_to_lines(selected_reads)
				for _, line := range read_lines {
					if strings.Contains(line, "PASS") {
						fmt.Fprintf(&output_builder, "%s\t*hit\n", line)
					} else {
						fmt.Fprintf(&output_builder, "%s\t?hit\n", line)
					}
				}
			}
		}
		if !*minus_bumps_printing {
			results <- output_builder.String()
			continue
		}
		if len(minus_bumps) > report_threshold {
			if !window_printed {
				fmt.Fprintf(&output_builder, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", observing_window.Rname, observing_window.Start,
					observing_window.Stop, observing_window.Site, string(observing_window.Strand), all, covered, passed, bumped, bumped_minus)
			}
			if len(minus_bumps) > 0 {
				read_lines := Reads_to_lines(minus_bumps)
				for _, line := range read_lines {
					fmt.Fprintf(&output_builder, "%s\t*minus_bump\n", line)
				}
			}
		}

		results <- output_builder.String()
	}
}

func count_worker(jobs <-chan *Test_job, results chan<- string, wg *sync.WaitGroup) {
	wg.Add(1)
	defer wg.Done()
	all := 0
	covered := 0
	passed := 0
	var ob_window *Observing_window
	var output_builder strings.Builder
	for j := range jobs {
		output_builder.Reset()
		ob_window = j.Ob_window
		all, covered, passed = Count_reads(j.Collection, ob_window)
		fmt.Fprintf(&output_builder, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n", ob_window.Rname, ob_window.Start,
			ob_window.Stop, ob_window.Site, string(ob_window.Strand), all, covered, passed)
		results <- output_builder.String()
	}
}

func test_reads_in_ob_window(thread_n int, ob_windows []*Observing_window, collection []Abstract_read, lib_type int, is_testing bool) {
	jobs := make(chan *Test_job, thread_n)
	results := make(chan string, thread_n)
	finished := 0
	var wg sync.WaitGroup
	var out_lines []string

	// open a routine to collect result
	wg.Add(1)
	go func() {
		defer wg.Done()
		for {
			r := <-results
			out_lines = append(out_lines, r)
			finished++
			if finished >= len(ob_windows) {
				close(jobs)
				close(results)
				break
			}
		}
		time.Sleep(1 * time.Second)
		sort.Sort(line_window_order(out_lines))
		for _, l := range out_lines {
			fmt.Print(l)
		}

	}()

	// open workers
	if is_testing {
		//testing
		for i := 0; i < thread_n; i++ {
			go test_worker(jobs, results, &wg, *report_threshold, !*repress_reads_printing)
		}
	} else {
		//counting
		for i := 0; i < thread_n; i++ {
			go count_worker(jobs, results, &wg)
		}
	}

	for _, observing_window := range ob_windows {
		jb := Test_job{
			Lib_type:   lib_type,
			Ob_window:  observing_window,
			Collection: collection,
		}
		jobs <- &jb
	}
	wg.Wait()
}

func main() {
	flag.Parse()
	usage := `usage: readsends -w observing_window.tsv -t sam -f condition1.sam -c test -l 1

options:
	-f  the target sam file or bed file
	-t  the target file type, which can be sam or bed
	-w  the tsv file that describe the observing window
	-c  command, can be: test, scan, count, default is test
	-l  libary type , can be: 1,2,0
	-s  the threshold to report when using the print list command
	-r  repress hit reads printing, default false
	-n  to print minus bump reads, default false
	-j  threads to run simultaneously
	-u  cpu cores to use
	-b  to fix cut coordinate on the - strand if the coordinate is one base before the real cut_site
	-g  when the library type is 0, assume observing windows strand is the gene strand

Descriptions(version 0.1.2):
1. break site notion: 
	in this program, break site means after this site the RNA breaks
2. command:
	scan: to include possible breaks, regardless of input cut site infomation
	test: to test if the input cut site can bump into the found read end, default
3. liberary types:
	1:  stranded liberary, read1 sequence = cDNA, read2 sequence = RNA
	2:  stranded liberary, read1 sequence = RNA, read2 sequence = cDNA
	0:  unstranded liberary
	use infer_experiment.py to test which library the dataset being tested belongs to
	if the result is like blew, the lib_type is 1
	    Fraction of reads explained by "1++,1--,2+-,2-+": 0.0049	
    	raction of reads explained by "1+-,1-+,2++,2--": 0.9000
	if the result is like blew, the lib_type is 2
		Fraction of reads explained by "1++,1--,2+-,2-+": 0.9000	
    	raction of reads explained by "1+-,1-+,2++,2--": 0.0046
4. Observing Window file:
	using normal bed file is ok when the name field in the field is just the break site position
	the cordinate is 1-based, each line is in the order of the following, deliminated by tab
	chromosome  start   stop  strand  cut_site 100
	`

	if *thread_n > runtime.NumCPU() {
		*thread_n = runtime.NumCPU()
	}
	if *cpu_cores > runtime.NumCPU() {
		*cpu_cores = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(*cpu_cores)

	fh_window, err := os.Open(*window_file)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Failed to open file:%s", *window_file)
		fmt.Println(usage)
		os.Exit(1)
	}

	defer fh_window.Close()

	var ob_windows = Read_Observing_windows(fh_window, *fix_ob_window)
	if len(ob_windows) == 0 {
		fmt.Fprintf(os.Stderr, "Observing windows file:%s contains not valid line", *window_file)
		fmt.Println(usage)
		os.Exit(1)
	}
	if *lib_type == -1 {
		fmt.Fprintln(os.Stderr, "library type is not set, please set it to 0, 1, or 2")
		fmt.Println(usage)
		os.Exit(1)
	}
	// var collection_sam []*Read
	var collection []Abstract_read
	if *target_type == "sam" {
		collection, err = Read_sam_file(*target_file)
		if err != nil {
			fmt.Println("reading sam file error:" + err.Error())
			os.Exit(1)
		}
		// collection = Sam_slice_to_Ab(collection_sam)
	} else if *target_type == "bed" {
		collection, err = Read_bed_file(*target_file)
		if err != nil {
			fmt.Println("reading bed file error:" + err.Error())
			os.Exit(1)
		}
	}

	if len(collection) == 0 {
		fmt.Println("target file contains no valid record")
		os.Exit(0)
	}

	switch *cmd {
	case "test":
		fmt.Println("Rname\tStart\tStop\tSite\tStrand\tAll\tCovered\tPasses\tBumps\tBump_minus")
		// limit the max channels to prevent system fail

		// ch_buffer_size := *thread_n
		test_reads_in_ob_window(*thread_n, ob_windows, collection, *lib_type, true)

	case "scan":
		fmt.Println("Rname\tStart\tStop\tSite\tStrand\tAll\tCovered\tPasses\tBumps\tBump_minus")
		w := new(Observing_window)
		for _, observing_window := range ob_windows {
			w.Rname = observing_window.Rname
			w.Start = observing_window.Start
			w.Stop = observing_window.Stop
			w.Strand = observing_window.Strand
			for i := observing_window.Start; i < observing_window.Stop+1; i++ {
				w.Site = i
				test_reads_in_ob_window(*thread_n, []*Observing_window{w}, collection, *lib_type, true)
			}
		}
	case "count":
		fmt.Println("Rname\tStart\tStop\tSite\tStrand\tAll\tCovered\tPasses")
		test_reads_in_ob_window(*thread_n, ob_windows, collection, 0, false)
	default:
		fmt.Printf("command \"%s\" unknown", *cmd)
	}

}
