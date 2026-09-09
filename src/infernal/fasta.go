package infernal

import (
	"bufio"
	"fmt"
	"io"
	"strings"
)

// FastaRecord is a single parsed FASTA entry: the header line (without the
// leading '>') and the concatenated, line-unwrapped sequence.
type FastaRecord struct {
	Header string
	Seq    string
}

// ID returns the first whitespace-delimited token of the header -- the
// sequence identifier cmsearch's --tblout reports in its "target name"
// column, as distinct from any free-text description that may follow it
// on the header line.
func (r FastaRecord) ID() string {
	if i := strings.IndexAny(r.Header, " \t"); i >= 0 {
		return r.Header[:i]
	}
	return r.Header
}

// ParseFASTA parses simple FASTA text into a slice of records. It
// tolerates wrapped sequence lines and blank lines, and does not
// validate alphabet.
func ParseFASTA(r io.Reader) ([]FastaRecord, error) {
	scanner := bufio.NewScanner(r)
	// Genomic flanking regions can exceed the default bufio.Scanner
	// token size on a single wrapped line in unusual cases; give it
	// plenty of headroom.
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 16*1024*1024)

	var records []FastaRecord
	var cur *FastaRecord
	var seqBuilder strings.Builder

	flush := func() {
		if cur != nil {
			cur.Seq = seqBuilder.String()
			records = append(records, *cur)
		}
	}

	for scanner.Scan() {
		line := strings.TrimRight(scanner.Text(), "\r\n")
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, ">") {
			flush()
			header := strings.TrimSpace(strings.TrimPrefix(line, ">"))
			cur = &FastaRecord{Header: header}
			seqBuilder.Reset()
			continue
		}
		if cur == nil {
			return nil, fmt.Errorf("infernal: FASTA data before any header line: %q", line)
		}
		seqBuilder.WriteString(strings.TrimSpace(line))
	}
	flush()

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("infernal: reading FASTA: %w", err)
	}
	return records, nil
}

// WriteFASTA writes records back out in plain (unwrapped) FASTA. Infernal
// doesn't care about line wrapping, so this favours simplicity over
// producing conventionally-wrapped 60/70-column output.
func WriteFASTA(w io.Writer, records []FastaRecord) error {
	bw := bufio.NewWriter(w)
	for _, rec := range records {
		if _, err := fmt.Fprintf(bw, ">%s\n%s\n", rec.Header, rec.Seq); err != nil {
			return fmt.Errorf("infernal: writing FASTA record %q: %w", rec.Header, err)
		}
	}
	return bw.Flush()
}

var complement = map[byte]byte{
	'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A',
	'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'u': 'a',
	'N': 'N', 'n': 'n',
}

// ReverseComplement returns the reverse complement of a DNA/RNA sequence.
// Bases outside the standard IUPAC A/C/G/T/U (case-insensitive) are passed
// through unchanged rather than erroring, since ambiguity codes show up
// legitimately in real genomic sequence.
func ReverseComplement(seq string) string {
	out := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		b := seq[len(seq)-1-i]
		if c, ok := complement[b]; ok {
			out[i] = c
		} else {
			out[i] = b
		}
	}
	return string(out)
}

// ExtractRegion returns the subsequence of rec spanning a cmsearch hit,
// given the hit's 1-based seqFrom/seqTo coordinates and strand exactly as
// reported in --tblout. For minus-strand hits, cmsearch reports seqFrom >
// seqTo (the hit read 3'->5' along the plus strand), so the region
// [seqTo, seqFrom] is extracted from the plus strand and reverse
// complemented.
func ExtractRegion(rec FastaRecord, seqFrom, seqTo int, strand string) (FastaRecord, error) {
	start, end := seqFrom, seqTo
	if strand == "-" {
		start, end = seqTo, seqFrom
	}
	if start < 1 || end < start || end > len(rec.Seq) {
		return FastaRecord{}, fmt.Errorf(
			"infernal: region %d-%d (strand %q) out of bounds for %q (length %d)",
			seqFrom, seqTo, strand, rec.ID(), len(rec.Seq))
	}

	sub := rec.Seq[start-1 : end]
	if strand == "-" {
		sub = ReverseComplement(sub)
	}

	return FastaRecord{
		Header: fmt.Sprintf("%s/%d-%d", rec.ID(), seqFrom, seqTo),
		Seq:    sub,
	}, nil
}
