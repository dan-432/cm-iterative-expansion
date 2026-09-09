package infernal

import (
	"context"
	"fmt"
	"os"
	"strconv"
)

// Aligner wraps Infernal's cmalign binary.
type Aligner struct {
	// BinPath is the path to the cmalign binary. If empty, "cmalign" is
	// resolved from PATH at call time.
	BinPath string

	// CPU is passed as cmalign's --cpu. Defaults to 1 if unset.
	CPU int
}

func (a Aligner) bin() string { return resolveBin(a.BinPath, "cmalign") }

func (a Aligner) cpu() int {
	if a.CPU <= 0 {
		return 1
	}
	return a.CPU
}

// Align aligns records to cmPath and writes the result as a Stockholm
// alignment at alignmentPath.
func (a Aligner) Align(ctx context.Context, cmPath string, records []FastaRecord, alignmentPath string) error {
	if len(records) == 0 {
		return fmt.Errorf("infernal: Align called with empty sequence set")
	}

	tmpFasta, err := os.CreateTemp("", "cie-align-*.fasta")
	if err != nil {
		return fmt.Errorf("cmalign: creating temp FASTA: %w", err)
	}
	tmpPath := tmpFasta.Name()
	defer os.Remove(tmpPath)

	if err := WriteFASTA(tmpFasta, records); err != nil {
		tmpFasta.Close()
		return fmt.Errorf("cmalign: writing temp FASTA: %w", err)
	}
	if err := tmpFasta.Close(); err != nil {
		return fmt.Errorf("cmalign: closing temp FASTA: %w", err)
	}

	args := []string{
		"-o", alignmentPath,
		"--outformat", "Stockholm",
		"--cpu", strconv.Itoa(a.cpu()),
		cmPath,
		tmpPath,
	}
	out, err := runCapture(ctx, a.bin(), args...)
	if err != nil {
		return fmt.Errorf("cmalign against %s: %w\ncmalign output:\n%s", cmPath, err, string(out))
	}
	return nil
}
