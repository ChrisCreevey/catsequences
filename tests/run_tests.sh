#!/usr/bin/env bash
# tests/run_tests.sh - regression test suite for catsequences
#
# Run from the repository root:
#   cc catsequences.c -o catsequences -lm
#   bash tests/run_tests.sh
#
# Exits 0 on success, non-zero if any test fails.

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${REPO}/catsequences"
EXPECTED="${REPO}/tests/expected"
TMPDIR_RUN="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_RUN"' EXIT

PASS=0
FAIL=0

# ── helper: compare two files, report result ─────────────────────
check() {
    local label="$1" got="$2" want="$3"
    if diff -q "$got" "$want" > /dev/null 2>&1; then
        echo "PASS: $label"
        PASS=$(( PASS + 1 ))
    else
        echo "FAIL: $label"
        diff "$want" "$got" | head -20
        FAIL=$(( FAIL + 1 ))
    fi
}

# ── test 1: basic concatenation via list file ─────────────────────
OUT="${TMPDIR_RUN}/basic"
"$BIN" -q -o "$OUT" "${REPO}/list"
check "basic FASTA output"      "${OUT}.fas"             "${EXPECTED}/basic_aln.fas"
check "basic partitions"        "${OUT}.partitions.txt"  "${EXPECTED}/basic_partitions.txt"
check "basic occupancy summary" "${OUT}.occupancy.txt"   "${EXPECTED}/basic_occupancy.txt"

# ── test 2: direct FASTA files as positional arguments ───────────
# Run from REPO so relative filenames match the expected partition output
OUT="${TMPDIR_RUN}/direct"
( cd "$REPO" && "$BIN" -q -o "$OUT" test1.txt test2.txt test3.txt )
check "direct-args FASTA"      "${OUT}.fas"            "${EXPECTED}/direct_aln.fas"
check "direct-args partitions" "${OUT}.partitions.txt" "${EXPECTED}/direct_partitions.txt"

# ── test 3: partition formats ─────────────────────────────────────
OUT="${TMPDIR_RUN}/raxml"
"$BIN" -q -f raxml  -o "$OUT" "${REPO}/list"
check "RAxML partition format"  "${OUT}.partitions.txt" "${EXPECTED}/raxml_partitions.txt"

OUT="${TMPDIR_RUN}/iqtree"
"$BIN" -q -f iqtree -o "$OUT" "${REPO}/list"
check "IQ-TREE partition format" "${OUT}.partitions.txt" "${EXPECTED}/iqtree_partitions.txt"

OUT="${TMPDIR_RUN}/nexpart"
"$BIN" -q -f nexus  -o "$OUT" "${REPO}/list"
check "NEXUS partition format"  "${OUT}.partitions.txt" "${EXPECTED}/nexus_part_partitions.txt"

# ── test 4: alignment output formats ─────────────────────────────
OUT="${TMPDIR_RUN}/nexus"
"$BIN" -q -F nexus  -o "$OUT" "${REPO}/list"
check "NEXUS alignment format"  "${OUT}.nex" "${EXPECTED}/nexus_aln.nex"

OUT="${TMPDIR_RUN}/phylip"
"$BIN" -q -F phylip -o "$OUT" "${REPO}/list"
check "PHYLIP alignment format" "${OUT}.phy" "${EXPECTED}/phylip_aln.phy"

# ── test 5: occupancy CSV ─────────────────────────────────────────
OUT="${TMPDIR_RUN}/occ"
"$BIN" -q -O -o "$OUT" "${REPO}/list"
check "occupancy matrix CSV" "${OUT}.occupancy.csv" "${EXPECTED}/occ_matrix.csv"

# ── test 6: missing-data character ────────────────────────────────
OUT="${TMPDIR_RUN}/missn"
"$BIN" -q -m N -o "$OUT" "${REPO}/list"
check "missing char N" "${OUT}.fas" "${EXPECTED}/missing_n.fas"

# ── test 7: species include filter ───────────────────────────────
OUT="${TMPDIR_RUN}/inc"
"$BIN" -q --include "${EXPECTED}/include.txt" -o "$OUT" "${REPO}/list"
check "species include filter" "${OUT}.fas" "${EXPECTED}/include_aln.fas"

# ── test 8: species exclude filter ───────────────────────────────
OUT="${TMPDIR_RUN}/exc"
"$BIN" -q --exclude "${EXPECTED}/exclude.txt" -o "$OUT" "${REPO}/list"
check "species exclude filter" "${OUT}.fas" "${EXPECTED}/exclude_aln.fas"

# ── test 9: pre-flight catches all missing files ──────────────────
printf 'test1.txt\nmissing_a.txt\ntest3.txt\nmissing_b.txt\n' \
    > "${TMPDIR_RUN}/badlist.txt"
ERROUT=$( cd "$REPO" && "$BIN" -o "${TMPDIR_RUN}/junk" \
              "${TMPDIR_RUN}/badlist.txt" 2>&1 || true )
NMISSING=$(echo "$ERROUT" | grep -c "Cannot open" || true)
if [ "$NMISSING" -eq 2 ]; then
    echo "PASS: pre-flight reports all missing files"
    PASS=$(( PASS + 1 ))
else
    echo "FAIL: pre-flight should report 2 missing files, got: $NMISSING"
    echo "$ERROUT"
    FAIL=$(( FAIL + 1 ))
fi

# ── test 10: duplicate handling – warn (default) ──────────────────
printf '>sp1\nAAAA\n>sp1\nTTTT\n>sp2\nCCCC\n' > "${TMPDIR_RUN}/dup.fasta"
OUT="${TMPDIR_RUN}/dup_warn"
STDOUT=$("$BIN" -q -o "$OUT" "${TMPDIR_RUN}/dup.fasta" 2>&1)
if echo "$STDOUT" | grep -q "Duplicate"; then
    echo "PASS: duplicate warn message printed"
    PASS=$(( PASS + 1 ))
else
    echo "FAIL: expected 'Duplicate' warning, got: $STDOUT"
    FAIL=$(( FAIL + 1 ))
fi

# ── test 11: duplicate handling – error ───────────────────────────
if ! "$BIN" -q --on-duplicate error \
       -o "${TMPDIR_RUN}/dup_err" "${TMPDIR_RUN}/dup.fasta" > /dev/null 2>&1; then
    echo "PASS: --on-duplicate error exits non-zero"
    PASS=$(( PASS + 1 ))
else
    echo "FAIL: --on-duplicate error should exit with failure"
    FAIL=$(( FAIL + 1 ))
fi

# ── test 12: duplicate handling – keep-last ───────────────────────
OUT="${TMPDIR_RUN}/dup_last"
"$BIN" -q --on-duplicate keep-last -o "$OUT" "${TMPDIR_RUN}/dup.fasta"
SP1SEQ=$(grep -A1 '^>sp1$' "${OUT}.fas" | tail -1)
if [ "$SP1SEQ" = "TTTT" ]; then
    echo "PASS: --on-duplicate keep-last keeps last sequence"
    PASS=$(( PASS + 1 ))
else
    echo "FAIL: keep-last expected 'TTTT', got '$SP1SEQ'"
    FAIL=$(( FAIL + 1 ))
fi

# ── test 13: glob expansion in list file ─────────────────────────
# Create a list file using a glob pattern and verify results match basic run
printf '%s/test*.txt\n' "$REPO" > "${TMPDIR_RUN}/globlist.txt"
OUT="${TMPDIR_RUN}/glob"
"$BIN" -q -o "$OUT" "${TMPDIR_RUN}/globlist.txt"
# The glob expands test1.txt test2.txt test3.txt in sorted order – same as basic
check "glob expansion in list file" "${OUT}.fas" "${EXPECTED}/basic_aln.fas"

# ── test 14: output prefix (-o) and missing char ─────────────────
OUT="${TMPDIR_RUN}/prefix_test/myprefix"
mkdir -p "${TMPDIR_RUN}/prefix_test"
"$BIN" -q -o "$OUT" "${REPO}/list"
if [ -f "${OUT}.fas" ] && [ -f "${OUT}.partitions.txt" ]; then
    echo "PASS: -o prefix creates correctly named files"
    PASS=$(( PASS + 1 ))
else
    echo "FAIL: -o prefix did not create expected files"
    FAIL=$(( FAIL + 1 ))
fi

# ── summary ───────────────────────────────────────────────────────
echo ""
echo "Results: ${PASS} passed, ${FAIL} failed"
[ "$FAIL" -eq 0 ]
