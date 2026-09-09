package infernal

import (
	"context"
	"fmt"
)

// Builder wraps Infernal's cmbuild binary.
type Builder struct {
	// BinPath is the path to the cmbuild binary. If empty, "cmbuild" is
	// resolved from PATH at call time.
	BinPath string
}

func (b Builder) bin() string { return resolveBin(b.BinPath, "cmbuild") }

// Build constructs a covariance model at cmPath from the Stockholm
// alignment at alignmentPath. If structured is true, --hand is passed,
// preserving the alignment's existing consensus columns and SS_cons
// (the case for a cmalign-produced alignment against an existing model).
// If false, --noss is passed, telling cmbuild to infer consensus columns
// and structure from scratch.
func (b Builder) Build(ctx context.Context, alignmentPath, cmPath string, structured bool) error {
	args := []string{}
	if structured {
		args = append(args, "--hand")
	} else {
		args = append(args, "--noss")
	}
	args = append(args, "--informat", "stockholm", cmPath, alignmentPath)

	out, err := runCapture(ctx, b.bin(), args...)
	if err != nil {
		return fmt.Errorf("cmbuild from %s: %w\ncmbuild output:\n%s", alignmentPath, err, string(out))
	}
	return nil
}
