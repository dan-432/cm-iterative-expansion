package infernal

import (
	"os"
	"strings"
	"testing"
)

func TestParseFASTA(t *testing.T) {
	text := ">h1\nAC\nTG\n\n>h2 some description\nUUUU\n"
	recs, err := ParseFASTA(strings.NewReader(text))
	if err != nil {
		t.Fatalf("ParseFASTA error: %v", err)
	}
	if len(recs) != 2 {
		t.Fatalf("records len = %d, want 2", len(recs))
	}
	if recs[0].Header != "h1" || recs[0].Seq != "ACTG" {
		t.Fatalf("record 0 = %+v, want Header=h1 Seq=ACTG", recs[0])
	}
	if recs[1].ID() != "h2" || recs[1].Seq != "UUUU" {
		t.Fatalf("record 1 ID/Seq = %q/%q, want h2/UUUU", recs[1].ID(), recs[1].Seq)
	}
}

func TestParseFASTARejectsDataBeforeHeader(t *testing.T) {
	if _, err := ParseFASTA(strings.NewReader("ACGT\n>h\nTT\n")); err == nil {
		t.Fatalf("expected error for sequence data before header")
	}
}

func TestFastaRecordID(t *testing.T) {
	r := FastaRecord{Header: "seq123 some free-text description"}
	if got := r.ID(); got != "seq123" {
		t.Fatalf("ID() = %q, want seq123", got)
	}
}

func TestWriteFASTARoundTrip(t *testing.T) {
	var buf strings.Builder
	in := []FastaRecord{{Header: "a", Seq: "ACGT"}, {Header: "b desc", Seq: "TTTT"}}
	if err := WriteFASTA(&buf, in); err != nil {
		t.Fatalf("WriteFASTA error: %v", err)
	}
	out, err := ParseFASTA(strings.NewReader(buf.String()))
	if err != nil {
		t.Fatalf("re-parsing written FASTA: %v", err)
	}
	if len(out) != 2 || out[0].Seq != "ACGT" || out[1].ID() != "b" {
		t.Fatalf("round trip mismatch: %+v", out)
	}
}

func TestReverseComplement(t *testing.T) {
	got := ReverseComplement("ACGTN")
	want := "NACGT"
	if got != want {
		t.Fatalf("ReverseComplement(ACGTN) = %q, want %q", got, want)
	}
}

func TestExtractRegionPlusStrand(t *testing.T) {
	rec := FastaRecord{Header: "seq1", Seq: "AAACCCGGGTTT"}
	got, err := ExtractRegion(rec, 4, 6, "+")
	if err != nil {
		t.Fatalf("ExtractRegion error: %v", err)
	}
	if got.Seq != "CCC" || got.Header != "seq1/4-6" {
		t.Fatalf("got %+v, want Seq=CCC Header=seq1/4-6", got)
	}
}

func TestExtractRegionMinusStrand(t *testing.T) {
	// cmsearch reports minus-strand hits with seqFrom > seqTo.
	rec := FastaRecord{Header: "seq1", Seq: "AAACCCGGGTTT"}
	got, err := ExtractRegion(rec, 9, 7, "-")
	if err != nil {
		t.Fatalf("ExtractRegion error: %v", err)
	}
	// Plus-strand region 7-9 is "GGG"; reverse complement is "CCC".
	if got.Seq != "CCC" {
		t.Fatalf("got Seq=%q, want CCC", got.Seq)
	}
}

func TestExtractRegionOutOfBounds(t *testing.T) {
	rec := FastaRecord{Header: "seq1", Seq: "ACGT"}
	if _, err := ExtractRegion(rec, 1, 100, "+"); err == nil {
		t.Fatalf("expected out-of-bounds error")
	}
}

func TestParseTbloutHits(t *testing.T) {
	tmp, err := os.CreateTemp("", "tblout-*.txt")
	if err != nil {
		t.Fatalf("CreateTemp: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "" +
		"# target name         accession query name           accession mdl mdl from   mdl to       seq from         seq to strand trunc pass   gc  bias  score   E-value inc description of target\n" +
		"# ------------------- --------- -------------------- --------- --- -------- -------- -------------- -------------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------\n" +
		"seqA     -         motif1               RF00001    cm        1      119              1            121      +    no    1 0.61   0.0  108.2   1.5e-27 ! -\n" +
		"seqA      -         motif1               RF00001    cm        1       73        2451092        2451021      -    no    1 0.55   0.0   64.3   2.1e-09 ! -\n" +
		"seqB         -         motif1               RF00001    cm       10      150            45            189      +   5'    2 0.42   0.2   45.1     0.004 ? -\n" +
		"# [ok]"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("write temp tblout: %v", err)
	}
	if err := tmp.Close(); err != nil {
		t.Fatalf("close temp tblout: %v", err)
	}

	hits, err := parseTbloutHits(tmp.Name())
	if err != nil {
		t.Fatalf("parseTbloutHits error: %v", err)
	}
	if len(hits) != 3 {
		t.Fatalf("len(hits) = %d, want 3", len(hits))
	}

	best := BestHitPerTarget(hits)
	if len(best) != 2 {
		t.Fatalf("len(best) = %d, want 2", len(best))
	}
	if best["seqA"].Score != 108.2 {
		t.Fatalf("seqA best score = %v, want 108.2", best["seqA"].Score)
	}
	if best["seqA"].SeqFrom != 1 || best["seqA"].SeqTo != 121 || best["seqA"].Strand != "+" {
		t.Fatalf("seqA best hit coords = %d-%d (%s), want 1-121 (+)",
			best["seqA"].SeqFrom, best["seqA"].SeqTo, best["seqA"].Strand)
	}
	if best["seqB"].Inc != "?" {
		t.Fatalf("seqB Inc = %q, want ?", best["seqB"].Inc)
	}
}
