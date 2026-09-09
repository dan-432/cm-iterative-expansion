package infernal

import (
	"context"
	"fmt"
	"strconv"
)

// Calibrator wraps Infernal's cmcalibrate binary.
type Calibrator struct {
	// BinPath is the path to the cmcalibrate binary. If empty,
	// "cmcalibrate" is resolved from PATH at call time.
	BinPath string

	CPU int

	// Seed is passed as cmcalibrate's --seed. Zero (Go's default) omits
	// the flag entirely and leaves Infernal to use its own default seed,
	// rather than assuming a seed is always supplied.
	Seed int
}

func (c Calibrator) bin() string { return resolveBin(c.BinPath, "cmcalibrate") }

func (c Calibrator) cpu() int {
	if c.CPU <= 0 {
		return 1
	}
	return c.CPU
}

func (c Calibrator) Calibrate(ctx context.Context, cmPath string) error {
	args := []string{"--cpu", strconv.Itoa(c.cpu())}
	if c.Seed != 0 {
		args = append(args, "--seed", strconv.Itoa(c.Seed))
	}
	args = append(args, cmPath)

	out, err := runCapture(ctx, c.bin(), args...)
	if err != nil {
		return fmt.Errorf("cmcalibrate for %s: %w\ncmcalibrate output:\n%s", cmPath, err, string(out))
	}
	return nil
}
