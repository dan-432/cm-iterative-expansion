// Command cmiterativeexpansion iteratively expands a covariance model:
// search a fixed sequence pool, align whatever hits, rebuild the model,
// and repeat until the included-hit set stops changing or a round limit
// is reached.
package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"os"
	"time"

	"cmiterativeexpansion/src/infernal"
	"cmiterativeexpansion/src/pipeline"
)

const banner = `
   ╭─────────────────────────────────────────────╮
   │             cmIterativeExpansion             │
   │     search -> align -> build -> calibrate    │
   ╰─────────────────────────────────────────────╯
              5'-N
                   \
                    G≡C
                    C≡G
                    A=U
                    U-A
                   /
              3'-N
`

func main() {
	cmPath := flag.String("cm", "", "path to the starting covariance model (required)")
	fastaPath := flag.String("fasta", "", "path to the FASTA sequence pool to search each round; expected to already contain the sequences the seed model was built from (required)")
	outDir := flag.String("out", "cie_out", "directory for per-round output (alignments, models)")
	maxRounds := flag.Int("max-rounds", 10, "maximum number of expansion rounds")
	minNewHits := flag.Int("min-new-hits", 1, "a round changing the included-hit set by fewer than this many sequences is treated as converged")
	evalue := flag.Float64("evalue", 0.01, "cmsearch inclusion E-value threshold (-E)")
	disableFilters := flag.Bool("disable-filters", true, "disable cmsearch's HMM acceleration filters (--max) so weak true homologs aren't filtered out before CM scoring")
	cpu := flag.Int("cpu", 4, "CPU threads for cmsearch/cmalign/cmcalibrate")
	seed := flag.Int("seed", 0, "cmcalibrate random seed (0 leaves Infernal's own default)")
	structured := flag.Bool("structured", true, "cmbuild --hand (true, preserve consensus structure) vs --noss (false, infer from scratch)")
	quiet := flag.Bool("quiet", false, "suppress the banner and per-round progress output")
	flag.Parse()

	if *cmPath == "" || *fastaPath == "" {
		fmt.Fprintln(os.Stderr, "usage: cmiterativeexpansion -cm model.cm -fasta seqs.fasta [flags]")
		flag.PrintDefaults()
		os.Exit(2)
	}

	if !*quiet {
		fmt.Print(banner)
		printInputSummary(*cmPath, *fastaPath, *outDir, *maxRounds, *minNewHits, *evalue, *cpu, *structured, *disableFilters, *seed)
	}

	cfg := pipeline.Config{
		CMPath:     *cmPath,
		FastaPath:  *fastaPath,
		OutDir:     *outDir,
		MaxRounds:  *maxRounds,
		MinNewHits: *minNewHits,
		Structured: *structured,
		Searcher: infernal.Searcher{
			CPU:            *cpu,
			EValueCutoff:   *evalue,
			DisableFilters: *disableFilters,
		},
		Aligner:    infernal.Aligner{CPU: *cpu},
		Builder:    infernal.Builder{},
		Calibrator: infernal.Calibrator{CPU: *cpu, Seed: *seed},
	}

	if !*quiet {
		start := time.Now()
		cfg.Progress = func(msg string) {
			fmt.Printf("  [%6s] %s\n", time.Since(start).Round(time.Second), msg)
		}
	}

	result, err := pipeline.Run(context.Background(), cfg)
	if err != nil {
		log.Fatalf("cmiterativeexpansion: %v", err)
	}

	printSummary(result)
}

func printInputSummary(cmPath, fastaPath, outDir string, maxRounds, minNewHits int, evalue float64, cpu int, structured, disableFilters bool, seed int) {
	nSeqs, err := countFastaRecords(fastaPath)
	seqCount := fmt.Sprintf("%d", nSeqs)
	if err != nil {
		seqCount = "unknown (couldn't pre-read FASTA: " + err.Error() + ")"
	}

	fmt.Println("inputs")
	fmt.Printf("  seed model     : %s\n", cmPath)
	fmt.Printf("  sequence pool  : %s (%s sequences)\n", fastaPath, seqCount)
	fmt.Printf("  output dir     : %s\n", outDir)
	fmt.Println("settings")
	fmt.Printf("  max rounds     : %d\n", maxRounds)
	fmt.Printf("  convergence    : stop when round-to-round change < %d sequences\n", minNewHits)
	fmt.Printf("  E-value cutoff : %g\n", evalue)
	fmt.Printf("  cpu threads    : %d\n", cpu)
	fmt.Printf("  cmbuild mode   : %s\n", map[bool]string{true: "--hand (preserve structure)", false: "--noss (infer from scratch)"}[structured])
	fmt.Printf("  cmsearch filter: %s\n", map[bool]string{true: "disabled (--max, full CM scoring)", false: "default HMM filters"}[disableFilters])
	if seed != 0 {
		fmt.Printf("  calibrate seed : %d\n", seed)
	}
	fmt.Println()
	fmt.Println("progress")
}

func countFastaRecords(path string) (int, error) {
	f, err := os.Open(path)
	if err != nil {
		return 0, err
	}
	defer f.Close()
	records, err := infernal.ParseFASTA(f)
	if err != nil {
		return 0, err
	}
	return len(records), nil
}

func printSummary(result pipeline.Result) {
	fmt.Println()
	fmt.Println("summary")
	for _, r := range result.Rounds {
		fmt.Printf("  round %2d: %4d included  (+%d / -%d)  -> %s\n",
			r.Round, len(r.IncludedIDs), len(r.NewIDs), len(r.LostIDs), r.CMPath)
	}
	if result.Converged {
		fmt.Printf("\nconverged after %d rounds\n", len(result.Rounds))
	} else {
		fmt.Printf("\nstopped at max-rounds without converging\n")
	}
	fmt.Printf("final model: %s\n", result.FinalCM)
}