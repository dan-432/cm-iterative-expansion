// Package pipeline implements the iterative expansion loop: search a
// fixed sequence pool with a covariance model, align whatever hits, build
// a new model from the alignment, calibrate it, and repeat until the set
// of included hits stops changing (or a round limit is hit).
package pipeline

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"sort"

	"cmiterativeexpansion/src/infernal"
)

// Config controls a single iterative expansion run.
type Config struct {
	// CMPath is the starting covariance model.
	CMPath string

	// FastaPath is the sequence pool searched every round. It's expected
	// to already contain the sequences used to build the seed model, so
	// no separate alignment-merge step against the original alignment is
	// needed: each round's alignment is built fresh from every currently
	// included hit (which will always include the original members, since
	// they trivially still hit their own model).
	FastaPath string

	// OutDir is where per-round CMs, alignments, and tblout files are
	// written, one subdirectory per round.
	OutDir string

	// MaxRounds caps how many search/align/build/calibrate iterations run.
	MaxRounds int

	// MinNewHits: a round whose included-hit set changes (new hits gained
	// plus hits lost, relative to the previous round) by fewer than this
	// many sequences is treated as converged and stops the loop. Defaults
	// to 1 (any change at all keeps iterating) if left at zero.
	MinNewHits int

	// Structured is passed through to cmbuild: true means --hand
	// (preserve consensus columns/structure from the cmalign-produced
	// alignment), false means --noss (infer from scratch).
	Structured bool

	Searcher   infernal.Searcher
	Aligner    infernal.Aligner
	Builder    infernal.Builder
	Calibrator infernal.Calibrator

	// Progress, if non-nil, is called with short human-readable status
	// updates as the run proceeds (round start, each stage completing,
	// round summary). Optional -- nil means no progress reporting, which
	// is what tests use to keep output quiet.
	Progress func(msg string)
}

func (cfg Config) progress(format string, args ...interface{}) {
	if cfg.Progress != nil {
		cfg.Progress(fmt.Sprintf(format, args...))
	}
}

// RoundResult records what happened in a single iteration.
type RoundResult struct {
	Round         int
	CMPath        string
	AlignmentPath string
	IncludedIDs   []string // best-hit sequence IDs included this round, sorted
	NewIDs        []string // included this round but not the previous round
	LostIDs       []string // included the previous round but not this one
}

// Result is the full record of an expansion run.
type Result struct {
	Rounds    []RoundResult
	Converged bool
	FinalCM   string
}

// Run executes the iterative expansion loop described by cfg.
func Run(ctx context.Context, cfg Config) (Result, error) {
	if cfg.MaxRounds <= 0 {
		return Result{}, fmt.Errorf("pipeline: MaxRounds must be > 0")
	}
	minNewHits := cfg.MinNewHits
	if minNewHits <= 0 {
		minNewHits = 1
	}

	if err := os.MkdirAll(cfg.OutDir, 0o755); err != nil {
		return Result{}, fmt.Errorf("pipeline: creating out dir: %w", err)
	}

	byID, err := readFastaByID(cfg.FastaPath)
	if err != nil {
		return Result{}, fmt.Errorf("pipeline: reading input FASTA: %w", err)
	}

	var result Result
	currentCM := cfg.CMPath
	var previousIDs map[string]bool

	for round := 1; round <= cfg.MaxRounds; round++ {
		roundDir := filepath.Join(cfg.OutDir, fmt.Sprintf("round_%02d", round))
		if err := os.MkdirAll(roundDir, 0o755); err != nil {
			return result, fmt.Errorf("pipeline: round %d: creating round dir: %w", round, err)
		}

		cfg.progress("round %d: searching %s against %s", round, cfg.FastaPath, currentCM)
		hits, err := cfg.Searcher.Search(ctx, currentCM, cfg.FastaPath)
		if err != nil {
			return result, fmt.Errorf("pipeline: round %d: cmsearch: %w", round, err)
		}
		best := infernal.BestHitPerTarget(hits)
		cfg.progress("round %d: %d sequences hit above threshold", round, len(best))

		includedIDs := make([]string, 0, len(best))
		for id := range best {
			includedIDs = append(includedIDs, id)
		}
		sort.Strings(includedIDs)
		if len(includedIDs) == 0 {
			return result, fmt.Errorf("pipeline: round %d: cmsearch returned no hits against %s", round, currentCM)
		}

		newIDs, lostIDs := diffIDs(previousIDs, includedIDs)

		extracted := make([]infernal.FastaRecord, 0, len(includedIDs))
		for _, id := range includedIDs {
			rec, ok := byID[id]
			if !ok {
				return result, fmt.Errorf("pipeline: round %d: hit %q not found in input FASTA %s", round, id, cfg.FastaPath)
			}
			h := best[id]
			region, err := infernal.ExtractRegion(rec, h.SeqFrom, h.SeqTo, h.Strand)
			if err != nil {
				return result, fmt.Errorf("pipeline: round %d: extracting hit region for %q: %w", round, id, err)
			}
			extracted = append(extracted, region)
		}

		cfg.progress("round %d: aligning %d hit regions against %s", round, len(extracted), currentCM)
		alignmentPath := filepath.Join(roundDir, "alignment.stk")
		if err := cfg.Aligner.Align(ctx, currentCM, extracted, alignmentPath); err != nil {
			return result, fmt.Errorf("pipeline: round %d: cmalign: %w", round, err)
		}

		cfg.progress("round %d: rebuilding model from alignment", round)
		nextCM := filepath.Join(roundDir, "model.cm")
		if err := cfg.Builder.Build(ctx, alignmentPath, nextCM, cfg.Structured); err != nil {
			return result, fmt.Errorf("pipeline: round %d: cmbuild: %w", round, err)
		}

		cfg.progress("round %d: calibrating rebuilt model", round)
		if err := cfg.Calibrator.Calibrate(ctx, nextCM); err != nil {
			return result, fmt.Errorf("pipeline: round %d: cmcalibrate: %w", round, err)
		}

		result.Rounds = append(result.Rounds, RoundResult{
			Round:         round,
			CMPath:        nextCM,
			AlignmentPath: alignmentPath,
			IncludedIDs:   includedIDs,
			NewIDs:        newIDs,
			LostIDs:       lostIDs,
		})

		currentCM = nextCM
		previousIDs = toSet(includedIDs)

		cfg.progress("round %d: complete -- %d included (%d new, %d lost)",
			round, len(includedIDs), len(newIDs), len(lostIDs))

		if round > 1 && len(newIDs) < minNewHits && len(lostIDs) < minNewHits {
			cfg.progress("converged after round %d", round)
			result.Converged = true
			result.FinalCM = currentCM
			return result, nil
		}
	}

	cfg.progress("reached max rounds (%d) without converging", cfg.MaxRounds)
	result.FinalCM = currentCM
	return result, nil
}

func readFastaByID(path string) (map[string]infernal.FastaRecord, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	records, err := infernal.ParseFASTA(f)
	if err != nil {
		return nil, err
	}

	byID := make(map[string]infernal.FastaRecord, len(records))
	for _, r := range records {
		byID[r.ID()] = r
	}
	return byID, nil
}

// diffIDs compares a sorted current ID list against the previous round's
// included set (nil on round 1, in which case every current ID counts as
// new). Both returned slices are sorted.
func diffIDs(previous map[string]bool, current []string) (newIDs, lostIDs []string) {
	currentSet := toSet(current)
	for _, id := range current {
		if !previous[id] {
			newIDs = append(newIDs, id)
		}
	}
	for id := range previous {
		if !currentSet[id] {
			lostIDs = append(lostIDs, id)
		}
	}
	sort.Strings(newIDs)
	sort.Strings(lostIDs)
	return newIDs, lostIDs
}

func toSet(ids []string) map[string]bool {
	set := make(map[string]bool, len(ids))
	for _, id := range ids {
		set[id] = true
	}
	return set
}