# cm-iterative-expansion
Tool developed as part of my honours project. Used infernal tools to itertively expand a covariance model by aligning new sequences discovered by cmsearch on a provided sequence dataset

## Relationship to cmsimile's `infernal` package

`src/infernal/` here is a deliberate fork of the `infernal` package in
[cmsimile](../cmsimile), not a shared dependency. The two tools use
Infernal for different purposes and that difference shows up directly in
what each `infernal` package needs to do:

- **cmsimile** scores a fixed pool of sequences (including
  dinucleotide-shuffled nulls) against a CM for statistical comparison.
  Its `cmsearch` wrapper (`Scorer`) is built to guarantee a real score for
  *every* sequence in the pool, however weak, and depends on cmsimile's
  own `seqpool` package for that.
- **cmiterativeexpansion** searches a sequence database to *discover*
  homologs above a threshold, extracts the hit regions, aligns them, and
  rebuilds the model. Its `cmsearch` wrapper (`Searcher`) reports hit
  coordinates and E-values rather than forcing a score for every input
  sequence, and has no need for `seqpool` at all — it works directly with
  a lightweight `FastaRecord` type.

Both packages wrap the same four Infernal binaries and share a family
resemblance (same `runCapture`/`resolveBin` subprocess pattern, same
tblout-parsing approach), but `cmbuild`/`cmalign`/`cmcalibrate` here take
different parameters and `cmsearch` here returns a different shape of
result entirely. Given the package is a thin wrapper around Infernal's
CLI — which isn't shifting under either tool — a real fork with a small
amount of duplicated boilerplate was judged less overhead than a shared
module dependency or a `go.work` link between two actively-developed
repos. If `cmiterativeexpansion`'s needs from `cmsearch`/`cmalign` ever
converge back with cmsimile's, that's a reasonable point to reconsider.

One functional difference worth flagging explicitly: cmsimile's
`Calibrator.Calibrate` reads `ctx.Value("seed")` and type-asserts it to
`int` unconditionally, which panics if the context doesn't carry that
key. This package's `Calibrator` instead takes an explicit `Seed int`
field (zero omits `--seed` entirely and lets Infernal use its own
default) — worth porting back to cmsimile if that code path is still live
there.

## Approach

Given a seed covariance model and a FASTA pool of DNA/RNA sequences
(expected to already contain whatever sequences the seed model was built
from), each round does:

1. **Search** -- `cmsearch` the full pool against the current model
   (`--max` to disable HMM acceleration filters, so a weak-but-real
   homolog isn't discarded before ever reaching full CM scoring; `-E` to
   set the inclusion threshold).
2. **Extract** -- for every sequence with a hit, pull out just the hit
   region using cmsearch's reported coordinates, reverse-complementing
   minus-strand hits so everything going into the alignment is in a
   consistent orientation.
3. **Align** -- `cmalign` those extracted regions against the current
   model.
4. **Rebuild** -- `cmbuild` a new model from that alignment (`--hand` by
   default, preserving the consensus structure the model already
   defines).
5. **Calibrate** -- `cmcalibrate` the new model so the next round's
   `cmsearch` E-values are meaningful.

No merge against the original seed alignment is needed: because the input
pool already contains the seed sequences, they trivially re-hit their own
descendant model every round, so each round's fresh alignment naturally
keeps including them.

The loop stops when either:
- the set of included hits changes by fewer than `-min-new-hits`
  sequences from the previous round (**converged**), or
- `-max-rounds` is reached without converging.

Every round's alignment and model are kept on disk under
`-out/round_NN/`, so a run can be inspected or resumed from any
intermediate point by hand.

## envs/

`envs/cie_env.yaml` describes the conda environment (Infernal, providing `cmemit`/`cmsearch`/`cmstat`/`cmfetch`; Go toolchain). Build with:

```
conda env create -f envs/cie_env.yaml
conda activate cie
```


## Building

Requires Go 1.22+.

```bash
# from the repo root
go build -o bin/cmiterativeexpansion ./src/
```

This produces a single static binary at `bin/cmiterativeexpansion` --
Infernal itself (`cmsearch`, `cmalign`, `cmbuild`, `cmcalibrate`) is
resolved from `$PATH` at run time and isn't bundled, so it needs to be
installed/loaded (e.g. via the `cie_env.yaml` conda environment) wherever
the binary runs.

Run the test suite (uses the fake Infernal binaries under
`tests/fixtures/` for the pipeline-level tests, so it doesn't need a real
Infernal install):

```bash
go vet ./...
go test ./...
```

Basic usage once built:

```bash
bin/cmiterativeexpansion \
  -cm seed.cm \
  -fasta sequences.fasta \
  -out cie_out \
  -max-rounds 10 \
  -evalue 0.01 \
  -cpu 8
```

Run `bin/cmiterativeexpansion -h` for the full flag list, or pass
`-quiet` to suppress the banner/progress output and just get the final
summary (useful for scripted/batch runs).