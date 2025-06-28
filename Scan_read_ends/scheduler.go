package main

import (
	"sync"
	"sync/atomic"
	"time"
)

type Job interface {
	Run() Job_result
	Id() int
}

type Job_result interface {
	Process()
}

type Scheduler struct {
	Jobs        chan Job
	Results     chan Job_result
	Accepted    int64
	Finished    int64
	Stop_signal int32
	Wg          *sync.WaitGroup
}

func (sch *Scheduler) Stop() {
	for {
		// wait
		if atomic.LoadInt64(&sch.Accepted) == atomic.LoadInt64(&sch.Finished) {
			atomic.AddInt32(&sch.Stop_signal, 1)
			close(sch.Jobs)
			close(sch.Results)
			break
		}
		time.Sleep(100 * time.Millisecond)
	}
	sch.Wg.Wait()
}

func (sch *Scheduler) Report() (int64, int64) {
	return sch.Accepted, sch.Finished
}

func (sch *Scheduler) Start(n int) {
	sch.Wg = new(sync.WaitGroup)
	sch.Jobs = make(chan Job, n)
	sch.Accepted = 0
	sch.Results = make(chan Job_result, n)
	for i := 0; i < n; i++ {
		sch.Wg.Add(1)
		go func() {
			// this is worker
			defer sch.Wg.Done()
			for j := range sch.Jobs {
				atomic.AddInt64(&sch.Accepted, 1)
				sch.Results <- j.Run()
			}
		}()
	}

	sch.Wg.Add(1)
	go func() {
		// this the for collecting result
		defer sch.Wg.Done()
		for {
			if sch.Stop_signal > 0 {
				break
			}
			r := <-sch.Results
			if r != nil {
				r.Process()
				atomic.AddInt64(&sch.Finished, 1)
			}

		}
	}()
}

func (sch *Scheduler) Feed(jb Job) {
	sch.Jobs <- jb
}

// // define job
type Read_ends_job struct {
	Ab_reads  []Abstract_read
	Ob_window *Observing_window
	Lib_type  int
}

type Read_ends_result struct {
	Bumps       []Abstract_read
	Minus_bumps []Abstract_read
	All         int
	Covered     int
	Passed      int
	Bump_count  int
	Minus_count int
}

func (*Read_ends_result) Process() {

}

// func (j *Read_ends_job) Run() Job_result {
// 	// ([]Abstract_read, []Abstract_read, int, int, int, int, int)
// 	selected_reads, minus_bumped_reads, all, covered, passed, bumped, bumped_minus := Check_break_site(j.Ab_reads, j.Ob_window, j.Lib_type)
// 	return &Read_ends_result{selected_reads, minus_bumped_reads, all, covered, passed, bumped, bumped_minus}

// }
