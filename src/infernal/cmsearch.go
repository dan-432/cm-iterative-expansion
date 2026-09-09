package infernal

import (
	"bufio"
	"context"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// Hit is one reported row from cmsearch's --tblout output.
type Hit struct {
	TargetName string
	Query      string
	SeqFrom    int
	SeqTo      int
	Strand     string
	Score      float64
	Evalue     float64
	Inc        string // "!" (included) or "?" (marginal) per --tblout
}

// Searcher wraps Infernal's cmsearch for homology discovery: given a
// model and a sequence pool, which sequences (and which regions of them)
// does the model consider hits.
type Searcher struct {
	// BinPath is the path to the cmsearch binary. If empty, "cmsearch" is
	// resolved from PATH at call time.
	BinPath string

	// CPU is passed as cmsearch's --cpu. Defaults to 1 if unset.
	CPU int

	// EValueCutoff is passed as cmsearch's -E, the reporting E-value
	// threshold. If zero, Infernal's own default (10) applies.
	EValueCutoff float64

	// DisableFilters, if true, passes --max to cmsearch, disabling the
	// HMM-based acceleration filter pipeline so every candidate region is
	// scored by the full CM. Slower, but avoids the filter silently
	// discarding a true but weakly-scoring homolog -- worth it for
	// motif discovery on modest-sized sequence pools.
	DisableFilters bool

	// KeepTempFiles, if true, leaves the temporary --tblout file on disk
	// instead of deleting it (useful for debugging a parsing mismatch by
	// hand against real cmsearch output).
	KeepTempFiles bool
}

func (s Searcher) bin() string { return resolveBin(s.BinPath, "cmsearch") }

func (s Searcher) cpu() int {
	if s.CPU <= 0 {
		return 1
	}
	return s.CPU
}

// Search runs cmsearch with cmPath against every sequence in fastaPath and
// returns every reported hit row (one per hit region -- a single sequence
// can appear multiple times if the model matches it in more than one
// place). Use BestHitPerTarget to collapse to one representative hit per
// sequence.
func (s Searcher) Search(ctx context.Context, cmPath, fastaPath string) ([]Hit, error) {
	tblFile, err := os.CreateTemp("", "cie-tblout-*.tbl")
	if err != nil {
		return nil, fmt.Errorf("cmsearch: creating temp tblout file: %w", err)
	}
	tblPath := tblFile.Name()
	tblFile.Close()
	defer func() {
		if !s.KeepTempFiles {
			os.Remove(tblPath)
		}
	}()

	args := []string{"--cpu", strconv.Itoa(s.cpu()), "--noali", "--tblout", tblPath}
	if s.EValueCutoff > 0 {
		args = append(args, "-E", strconv.FormatFloat(s.EValueCutoff, 'g', -1, 64))
	}
	if s.DisableFilters {
		args = append(args, "--max", "--smxsize", "1024")
	}
	args = append(args, cmPath, fastaPath)

	if out, err := runCapture(ctx, s.bin(), args...); err != nil {
		return nil, fmt.Errorf("cmsearch against %s: %w\ncmsearch output:\n%s", cmPath, err, string(out))
	}

	hits, err := parseTbloutHits(tblPath)
	if err != nil {
		return nil, fmt.Errorf("cmsearch against %s: parsing --tblout: %w", cmPath, err)
	}
	return hits, nil
}

// BestHitPerTarget collapses possibly-multiple hit rows per sequence down
// to the single highest-scoring region per target name.
func BestHitPerTarget(hits []Hit) map[string]Hit {
	best := make(map[string]Hit, len(hits))
	for _, h := range hits {
		existing, ok := best[h.TargetName]
		if !ok || h.Score > existing.Score {
			best[h.TargetName] = h
		}
	}
	return best
}

// --tblout is a fixed, documented 18-column layout (see the Infernal
// user guide's "tabular output" section) that doesn't change release to
// release, unlike cmscore's plain-text table.
const (
	tbloutTargetNameCol = 0
	tbloutQueryNameCol  = 2
	tbloutSeqFromCol    = 7
	tbloutSeqToCol      = 8
	tbloutStrandCol     = 9
	tbloutScoreCol      = 14
	tbloutEvalueCol     = 15
	tbloutIncCol        = 16
)

func parseTbloutHits(path string) ([]Hit, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var hits []Hit
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) <= tbloutIncCol {
			continue // malformed/short line, not a real data row
		}

		seqFrom, err := strconv.Atoi(fields[tbloutSeqFromCol])
		if err != nil {
			continue
		}
		seqTo, err := strconv.Atoi(fields[tbloutSeqToCol])
		if err != nil {
			continue
		}
		score, err := strconv.ParseFloat(fields[tbloutScoreCol], 64)
		if err != nil {
			continue
		}
		evalue, err := strconv.ParseFloat(fields[tbloutEvalueCol], 64)
		if err != nil {
			continue
		}

		hits = append(hits, Hit{
			TargetName: fields[tbloutTargetNameCol],
			Query:      fields[tbloutQueryNameCol],
			SeqFrom:    seqFrom,
			SeqTo:      seqTo,
			Strand:     fields[tbloutStrandCol],
			Score:      score,
			Evalue:     evalue,
			Inc:        fields[tbloutIncCol],
		})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return hits, nil
}
