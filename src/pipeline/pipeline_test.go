package pipeline

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"cmiterativeexpansion/src/infernal"
)

const fixtureDir = "../../tests/fixtures"

func writeTestFasta(t *testing.T, dir string) string {
	t.Helper()
	path := filepath.Join(dir, "seqs.fasta")
	content := ">seqA\nACGUACGUACGU\n>seqB\nACGUACGUACGUACGU\n>seqC\nACGU\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("writing test FASTA: %v", err)
	}
	return path
}

func TestRunConvergesWithFakeInfernal(t *testing.T) {
	dir := t.TempDir()
	fastaPath := writeTestFasta(t, dir)
	seedCM := filepath.Join(dir, "seed.cm")
	if err := os.WriteFile(seedCM, []byte("INFERNAL1/a\nNAME seed\n//\n"), 0o644); err != nil {
		t.Fatalf("writing seed CM: %v", err)
	}

	cfg := Config{
		CMPath:     seedCM,
		FastaPath:  fastaPath,
		OutDir:     filepath.Join(dir, "out"),
		MaxRounds:  5,
		MinNewHits: 1,
		Structured: true,
		Searcher:   infernal.Searcher{BinPath: filepath.Join(fixtureDir, "fake_cmsearch")},
		Aligner:    infernal.Aligner{BinPath: filepath.Join(fixtureDir, "fake_cmalign")},
		Builder:    infernal.Builder{BinPath: filepath.Join(fixtureDir, "fake_cmbuild")},
		Calibrator: infernal.Calibrator{BinPath: filepath.Join(fixtureDir, "fake_cmcalibrate")},
	}

	result, err := Run(context.Background(), cfg)
	if err != nil {
		t.Fatalf("Run error: %v", err)
	}
	if !result.Converged {
		t.Fatalf("expected convergence, got %d rounds without converging: %+v", len(result.Rounds), result.Rounds)
	}
	// fake_cmsearch's hit set is round-independent, so this should
	// converge on round 2 (round 1 establishes the baseline, which
	// counts every ID as "new" by definition).
	if len(result.Rounds) != 2 {
		t.Fatalf("len(result.Rounds) = %d, want 2", len(result.Rounds))
	}

	first, second := result.Rounds[0], result.Rounds[1]
	if len(first.IncludedIDs) != 3 {
		t.Fatalf("round 1 IncludedIDs = %v, want 3 sequences", first.IncludedIDs)
	}
	if len(first.NewIDs) != 3 {
		t.Fatalf("round 1 NewIDs = %v, want all 3 (baseline round)", first.NewIDs)
	}
	if len(second.NewIDs) != 0 || len(second.LostIDs) != 0 {
		t.Fatalf("round 2 should show no change, got NewIDs=%v LostIDs=%v", second.NewIDs, second.LostIDs)
	}
	if result.FinalCM != second.CMPath {
		t.Fatalf("FinalCM = %q, want %q", result.FinalCM, second.CMPath)
	}
	if !strings.HasSuffix(result.FinalCM, "model.cm") {
		t.Fatalf("FinalCM = %q, want it to end in model.cm", result.FinalCM)
	}
}

func TestRunStopsAtMaxRoundsWithoutConverging(t *testing.T) {
	// A MinNewHits of zero sequences (i.e. requiring the hit set to be
	// perfectly identical AND non-empty on the very first comparison)
	// combined with only a single allowed round means Run should stop
	// at MaxRounds rather than ever reaching round > 1 to check
	// convergence at all.
	dir := t.TempDir()
	fastaPath := writeTestFasta(t, dir)
	seedCM := filepath.Join(dir, "seed.cm")
	os.WriteFile(seedCM, []byte("INFERNAL1/a\nNAME seed\n//\n"), 0o644)

	cfg := Config{
		CMPath:     seedCM,
		FastaPath:  fastaPath,
		OutDir:     filepath.Join(dir, "out"),
		MaxRounds:  1,
		MinNewHits: 1,
		Structured: true,
		Searcher:   infernal.Searcher{BinPath: filepath.Join(fixtureDir, "fake_cmsearch")},
		Aligner:    infernal.Aligner{BinPath: filepath.Join(fixtureDir, "fake_cmalign")},
		Builder:    infernal.Builder{BinPath: filepath.Join(fixtureDir, "fake_cmbuild")},
		Calibrator: infernal.Calibrator{BinPath: filepath.Join(fixtureDir, "fake_cmcalibrate")},
	}

	result, err := Run(context.Background(), cfg)
	if err != nil {
		t.Fatalf("Run error: %v", err)
	}
	if result.Converged {
		t.Fatalf("did not expect convergence with MaxRounds=1")
	}
	if len(result.Rounds) != 1 {
		t.Fatalf("len(result.Rounds) = %d, want 1", len(result.Rounds))
	}
}

func TestDiffIDs(t *testing.T) {
	previous := toSet([]string{"a", "b", "c"})
	newIDs, lostIDs := diffIDs(previous, []string{"b", "c", "d"})
	if len(newIDs) != 1 || newIDs[0] != "d" {
		t.Fatalf("newIDs = %v, want [d]", newIDs)
	}
	if len(lostIDs) != 1 || lostIDs[0] != "a" {
		t.Fatalf("lostIDs = %v, want [a]", lostIDs)
	}
}

func TestDiffIDsFirstRound(t *testing.T) {
	newIDs, lostIDs := diffIDs(nil, []string{"a", "b"})
	if len(newIDs) != 2 {
		t.Fatalf("newIDs = %v, want both a and b treated as new on round 1", newIDs)
	}
	if len(lostIDs) != 0 {
		t.Fatalf("lostIDs = %v, want none", lostIDs)
	}
}
