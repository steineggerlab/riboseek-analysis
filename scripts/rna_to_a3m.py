#!/usr/bin/env python3
"""
RNA to A3M converter using the same approach as AlphaFold3.

This script takes RNA sequences and generates A3M (MSA) files by searching
against user-specified databases using nhmmer, then realigning hits with
hmmalign for consistent MSA width.

Usage:
    python rna_to_a3m.py --input sequences.fasta --output_dir ./output \
        --db1 /path/to/database1.fasta --db2 /path/to/database2.fasta

Requirements:
    - nhmmer (from HMMER suite)
    - hmmbuild (from HMMER suite)
    - hmmalign (from HMMER suite)
"""

from __future__ import annotations

import argparse
import os
import re
import string
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterable, List, Optional, Sequence, Tuple

SHORT_SEQUENCE_CUTOFF = 50


@dataclass
class MsaToolResult:
    """Result from an MSA search tool."""
    target_sequence: str
    a3m: str
    e_value: float


def check_binary_exists(path: str, name: str) -> None:
    """Check if a binary exists at the given path."""
    if not os.path.exists(path):
        raise RuntimeError(f"{name} binary not found at {path}")


def create_query_fasta_file(sequence: str, path: str, linewidth: int = 80) -> None:
    """Create a FASTA file with the sequence."""
    with open(path, "w") as f:
        f.write(">query\n")
        i = 0
        while i < len(sequence):
            f.write(f"{sequence[i:(i + linewidth)]}\n")
            i += linewidth


def run_command(cmd: Sequence[str], cmd_name: str) -> subprocess.CompletedProcess:
    """Run a command and handle errors."""
    print(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        error_msg = f"{cmd_name} failed\nstdout:\n{e.stdout}\nstderr:\n{e.stderr}"
        raise RuntimeError(error_msg) from e
    return result


def lazy_parse_fasta_string(fasta_string: str) -> Iterable[tuple[str, str]]:
    """Parse a FASTA string and yield (sequence, description) tuples."""
    lines = fasta_string.strip().split("\n")
    current_desc = None
    current_seq = []

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_desc is not None:
                yield "".join(current_seq), current_desc
            current_desc = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line)

    if current_desc is not None:
        yield "".join(current_seq), current_desc


def parse_fasta(fasta_string: str) -> tuple[list[str], list[str]]:
    """Parse FASTA string and return sequences and descriptions."""
    sequences = []
    descriptions = []
    for seq, desc in lazy_parse_fasta_string(fasta_string):
        sequences.append(seq)
        descriptions.append(desc)
    return sequences, descriptions


def convert_stockholm_to_a3m(
    stockholm: IO[str],
    max_sequences: int | None = None,
    remove_first_row_gaps: bool = True,
) -> str:
    """Convert MSA in Stockholm format to A3M format."""
    descriptions = {}
    sequences = {}
    reached_max_sequences = False

    for line in stockholm:
        reached_max_sequences = max_sequences and len(sequences) >= max_sequences
        line = line.strip()
        if not line or line.startswith(("#", "//")):
            continue
        seqname, aligned_seq = line.split(maxsplit=1)
        if seqname not in sequences:
            if reached_max_sequences:
                continue
            sequences[seqname] = ""
        sequences[seqname] += aligned_seq

    if not sequences:
        return ""

    stockholm.seek(0)
    for line in stockholm:
        line = line.strip()
        if line[:4] == "#=GS":
            columns = line.split(maxsplit=3)
            seqname, feature = columns[1:3]
            value = columns[3] if len(columns) == 4 else ""
            if feature != "DE":
                continue
            if reached_max_sequences and seqname not in sequences:
                continue
            descriptions[seqname] = value
            if len(descriptions) == len(sequences):
                break

    a3m_sequences = {}
    query_sequence = next(iter(sequences.values()))
    for seqname, sto_sequence in sequences.items():
        if remove_first_row_gaps:
            a3m_sequences[seqname] = align_sequence_to_gapless_query(
                sto_sequence, query_sequence
            ).replace(".", "")
        else:
            a3m_sequences[seqname] = sto_sequence.replace(".", "")

    fasta_chunks = []
    for seqname, seq in a3m_sequences.items():
        fasta_chunks.append(f">{seqname} {descriptions.get(seqname, '')}")
        fasta_chunks.append(seq)

    return "\n".join(fasta_chunks) + "\n"


def align_sequence_to_gapless_query(sequence: str, query_sequence: str) -> str:
    """
    Align a sequence to a gapless query (A3M conversion).
    Converts insertions (relative to query) to lowercase.
    """
    result = []
    for query_char, seq_char in zip(query_sequence, sequence):
        if query_char == "-":
            if seq_char != "-":
                result.append(seq_char.lower())
        else:
            result.append(seq_char)
    return "".join(result)


class Hmmbuild:
    """Python wrapper for hmmbuild."""

    def __init__(self, binary_path: str, alphabet: str | None = None):
        self._binary_path = binary_path
        self._alphabet = alphabet
        check_binary_exists(self._binary_path, "hmmbuild")

    def build_profile_from_a3m(self, a3m: str) -> str:
        """Build HMM profile from A3M string."""
        lines = []
        for sequence, description in lazy_parse_fasta_string(a3m):
            sequence = re.sub("[a-z]+", "", sequence)
            lines.append(f">{description}\n{sequence}\n")
        msa = "".join(lines)

        with tempfile.TemporaryDirectory(dir="/mnt/scratch/sukhwan/") as tmp_dir:
            input_msa_path = os.path.join(tmp_dir, "query.msa")
            output_hmm_path = os.path.join(tmp_dir, "output.hmm")

            with open(input_msa_path, "w") as f:
                f.write(msa)

            cmd = [self._binary_path, "--informat", "afa"]
            if self._alphabet:
                cmd.append(f"--{self._alphabet}")
            cmd.extend([output_hmm_path, input_msa_path])

            run_command(cmd, "hmmbuild")

            with open(output_hmm_path) as f:
                hmm = f.read()

        return hmm


class Hmmalign:
    """Python wrapper for hmmalign."""

    def __init__(self, binary_path: str):
        self._binary_path = binary_path
        check_binary_exists(self._binary_path, "hmmalign")

    def align_sequences_to_profile(self, profile: str, sequences_a3m: str) -> str:
        """Align sequences to profile and return A3M string."""
        deletion_table = str.maketrans("", "", "-")
        sequences_no_gaps_a3m = []
        for seq, desc in lazy_parse_fasta_string(sequences_a3m):
            sequences_no_gaps_a3m.append(f">{desc}")
            sequences_no_gaps_a3m.append(seq.translate(deletion_table))
        sequences_no_gaps_a3m = "\n".join(sequences_no_gaps_a3m)

        with tempfile.TemporaryDirectory(dir="/mnt/scratch/sukhwan/") as tmp_dir:
            input_profile = os.path.join(tmp_dir, "profile.hmm")
            input_sequences = os.path.join(tmp_dir, "sequences.a3m")
            output_a3m_path = os.path.join(tmp_dir, "output.a3m")

            with open(input_profile, "w") as f:
                f.write(profile)

            with open(input_sequences, "w") as f:
                f.write(sequences_no_gaps_a3m)

            cmd = [
                self._binary_path,
                "-o", output_a3m_path,
                "--outformat", "A2M",
                input_profile,
                input_sequences,
            ]

            run_command(cmd, "hmmalign")

            with open(output_a3m_path, encoding="utf-8") as f:
                a3m = f.read()

        return a3m


class Nhmmer:
    """Python wrapper for nhmmer (RNA/DNA sequence search)."""

    def __init__(
        self,
        binary_path: str,
        hmmalign_binary_path: str,
        hmmbuild_binary_path: str,
        database_path: str,
        n_cpu: int = 8,
        e_value: float = 1e-3,
        max_sequences: int = 10000,
        alphabet: str = "rna",
    ):
        self._binary_path = binary_path
        self._hmmalign_binary_path = hmmalign_binary_path
        self._hmmbuild_binary_path = hmmbuild_binary_path
        self._database_path = database_path
        self._n_cpu = n_cpu
        self._e_value = e_value
        self._max_sequences = max_sequences
        self._alphabet = alphabet

        check_binary_exists(self._binary_path, "nhmmer")

    def query(self, target_sequence: str) -> MsaToolResult:
        """Query the database using nhmmer."""
        print(f"Querying database: {os.path.basename(self._database_path)}")
        print(f"Query length: {len(target_sequence)} nt")

        with tempfile.TemporaryDirectory(dir="/mnt/scratch/sukhwan/") as tmp_dir:
            input_a3m_path = os.path.join(tmp_dir, "query.a3m")
            output_sto_path = os.path.join(tmp_dir, "output.sto")
            Path(output_sto_path).touch()

            create_query_fasta_file(target_sequence, input_a3m_path)

            cmd_flags = [
                "-o", "/dev/null",
                "--noali",
                "--cpu", str(self._n_cpu),
                "-E", str(self._e_value),
            ]

            if self._alphabet:
                cmd_flags.append(f"--{self._alphabet}")

            cmd_flags.extend(["-A", output_sto_path])

            if self._alphabet == "rna" and len(target_sequence) < SHORT_SEQUENCE_CUTOFF:
                cmd_flags.extend(["--F3", "0.02"])

            cmd_flags.extend([input_a3m_path, self._database_path])

            cmd = [self._binary_path] + cmd_flags
            run_command(cmd, f"nhmmer ({os.path.basename(self._database_path)})")

            if os.path.getsize(output_sto_path) > 0:
                with open(output_sto_path) as f:
                    a3m_out = convert_stockholm_to_a3m(
                        f, max_sequences=self._max_sequences - 1
                    )

                if a3m_out.strip():
                    print(f"Aligning {len(a3m_out)} bytes of hits to query profile")

                    aligner = Hmmalign(self._hmmalign_binary_path)
                    target_sequence_fasta = f">query\n{target_sequence}\n"
                    profile_builder = Hmmbuild(
                        binary_path=self._hmmbuild_binary_path, alphabet=self._alphabet
                    )
                    profile = profile_builder.build_profile_from_a3m(target_sequence_fasta)
                    a3m_out = aligner.align_sequences_to_profile(
                        profile=profile, sequences_a3m=a3m_out
                    )
                    a3m_out = "".join([target_sequence_fasta, a3m_out])

                    a3m = "\n".join(
                        [f">{n}\n{s}" for s, n in lazy_parse_fasta_string(a3m_out)]
                    )
                else:
                    a3m = f">query\n{target_sequence}"
            else:
                a3m = f">query\n{target_sequence}"

        return MsaToolResult(
            target_sequence=target_sequence,
            e_value=self._e_value,
            a3m=a3m,
        )


class Msa:
    """Multiple Sequence Alignment container."""

    def __init__(
        self,
        query_sequence: str,
        sequences: Sequence[str],
        descriptions: Sequence[str],
        deduplicate: bool = True,
    ):
        if len(sequences) != len(descriptions):
            raise ValueError("Number of sequences and descriptions must match.")

        self.query_sequence = query_sequence

        if not deduplicate:
            self.sequences = list(sequences)
            self.descriptions = list(descriptions)
        else:
            self.sequences = []
            self.descriptions = []
            deletion_table = str.maketrans("", "", string.ascii_lowercase)
            unique_sequences = set()
            for seq, desc in zip(sequences, descriptions):
                sequence_no_deletions = seq.translate(deletion_table)
                if sequence_no_deletions not in unique_sequences:
                    unique_sequences.add(sequence_no_deletions)
                    self.sequences.append(seq)
                    self.descriptions.append(desc)

        self.sequences = self.sequences or [query_sequence]
        self.descriptions = self.descriptions or ["Original query"]

    @classmethod
    def from_a3m(
        cls,
        query_sequence: str,
        a3m: str,
        max_depth: int | None = None,
        deduplicate: bool = True,
    ) -> "Msa":
        """Parse A3M and build Msa object."""
        sequences, descriptions = parse_fasta(a3m)

        if max_depth is not None and 0 < max_depth < len(sequences):
            print(f"MSA cropped from depth {len(sequences)} to {max_depth}")
            sequences = sequences[:max_depth]
            descriptions = descriptions[:max_depth]

        return cls(
            query_sequence=query_sequence,
            sequences=sequences,
            descriptions=descriptions,
            deduplicate=deduplicate,
        )

    @classmethod
    def from_multiple_msas(cls, msas: Sequence["Msa"], deduplicate: bool = True) -> "Msa":
        """Merge multiple MSAs into one."""
        if not msas:
            raise ValueError("At least one MSA must be provided.")

        query_sequence = msas[0].query_sequence
        sequences = []
        descriptions = []

        for msa in msas:
            if msa.query_sequence != query_sequence:
                raise ValueError("Query sequences must match across all MSAs.")
            sequences.extend(msa.sequences)
            descriptions.extend(msa.descriptions)

        return cls(
            query_sequence=query_sequence,
            sequences=sequences,
            descriptions=descriptions,
            deduplicate=deduplicate,
        )

    @property
    def depth(self) -> int:
        return len(self.sequences)

    def to_a3m(self) -> str:
        """Return the MSA in A3M format."""
        a3m_lines = []
        for desc, seq in zip(self.descriptions, self.sequences):
            a3m_lines.append(f">{desc}")
            a3m_lines.append(seq)
        return "\n".join(a3m_lines) + "\n"


def get_rna_msa(
    target_sequence: str,
    db_paths: List[str],
    nhmmer_binary: str,
    hmmalign_binary: str,
    hmmbuild_binary: str,
    n_cpu: int = 8,
    e_value: float = 1e-3,
    max_sequences: int = 10000,
) -> Msa:
    """
    Get RNA MSA by searching multiple databases sequentially.

    Args:
        target_sequence: The RNA sequence to search for.
        db_paths: List of database paths to search.
        nhmmer_binary: Path to nhmmer binary.
        hmmalign_binary: Path to hmmalign binary.
        hmmbuild_binary: Path to hmmbuild binary.
        n_cpu: Number of CPUs for nhmmer.
        e_value: E-value threshold.
        max_sequences: Maximum sequences per database.

    Returns:
        Merged MSA from all databases.
    """
    msas = []
    for i, db_path in enumerate(db_paths):
        print(f"Searching database {i+1}/{len(db_paths)}: {os.path.basename(db_path)}")
        nhmmer = Nhmmer(
            binary_path=nhmmer_binary,
            hmmalign_binary_path=hmmalign_binary,
            hmmbuild_binary_path=hmmbuild_binary,
            database_path=db_path,
            n_cpu=n_cpu,
            e_value=e_value,
            max_sequences=max_sequences,
            alphabet="rna",
        )
        result = nhmmer.query(target_sequence)
        msa = Msa.from_a3m(
            query_sequence=target_sequence,
            a3m=result.a3m,
            deduplicate=False,
        )
        msas.append(msa)
        print(f"  -> Found {msa.depth} sequences")

    merged_msa = Msa.from_multiple_msas(msas, deduplicate=True)
    print(f"Merged MSA depth (deduplicated): {merged_msa.depth}")

    return merged_msa


def process_rna_sequences(
    input_fasta: str,
    output_dir: str,
    db_paths: List[str],
    nhmmer_binary: str,
    hmmalign_binary: str,
    hmmbuild_binary: str,
    n_cpu: int = 8,
    e_value: float = 1e-3,
    max_sequences: int = 10000,
) -> None:
    """
    Process multiple RNA sequences and generate A3M files.

    Args:
        input_fasta: Path to input FASTA file with RNA sequences.
        output_dir: Directory to write A3M output files.
        db_paths: List of database paths to search.
        nhmmer_binary: Path to nhmmer binary.
        hmmalign_binary: Path to hmmalign binary.
        hmmbuild_binary: Path to hmmbuild binary.
        n_cpu: Number of CPUs for nhmmer.
        e_value: E-value threshold.
        max_sequences: Maximum sequences per database.
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(input_fasta) as f:
        input_content = f.read()

    sequences_list, descriptions_list = parse_fasta(input_content)

    print(f"Found {len(sequences_list)} RNA sequence(s) to process")
    print(f"Searching against {len(db_paths)} database(s)")
    print("-" * 60)

    for i, (sequence, description) in enumerate(zip(sequences_list, descriptions_list)):
        seq_name = description.split()[0] if description else f"sequence_{i+1}"
        safe_name = re.sub(r"[^\w\-_.]", "_", seq_name)

        print(f"\n[{i+1}/{len(sequences_list)}] Processing: {seq_name}")
        print(f"Sequence length: {len(sequence)} nt")

        sequence = sequence.upper().replace("T", "U")

        msa = get_rna_msa(
            target_sequence=sequence,
            db_paths=db_paths,
            nhmmer_binary=nhmmer_binary,
            hmmalign_binary=hmmalign_binary,
            hmmbuild_binary=hmmbuild_binary,
            n_cpu=n_cpu,
            e_value=e_value,
            max_sequences=max_sequences,
        )

        output_path = os.path.join(output_dir, f"{safe_name}.a3m")
        with open(output_path, "w") as f:
            f.write(msa.to_a3m())

        print(f"Wrote: {output_path} ({msa.depth} sequences)")

    print("\n" + "=" * 60)
    print("Processing complete!")


def find_binary(name: str) -> str | None:
    """Try to find a binary in PATH."""
    import shutil
    return shutil.which(name)


def main():
    parser = argparse.ArgumentParser(
        description="Generate A3M files for RNA sequences using nhmmer (AlphaFold3-style).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with two databases
    python rna_to_a3m.py --input rna_sequences.fasta --output_dir ./a3m_output \\
        --db1 /path/to/rfam.fasta --db2 /path/to/rnacentral.fasta

    # With custom binary paths
    python rna_to_a3m.py --input rna_sequences.fasta --output_dir ./a3m_output \\
        --db1 /path/to/db1.fasta --db2 /path/to/db2.fasta \\
        --nhmmer /usr/local/bin/nhmmer \\
        --hmmalign /usr/local/bin/hmmalign \\
        --hmmbuild /usr/local/bin/hmmbuild

    # With custom parameters
    python rna_to_a3m.py --input rna_sequences.fasta --output_dir ./a3m_output \\
        --db1 /path/to/db1.fasta --db2 /path/to/db2.fasta \\
        --e_value 1e-5 --max_sequences 5000 --n_cpu 16
        """,
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file containing RNA sequences",
    )
    parser.add_argument(
        "--output_dir", "-o",
        required=True,
        help="Output directory for A3M files",
    )
    parser.add_argument(
        "--db1",
        required=True,
        help="Path to first database (FASTA format)",
    )
    parser.add_argument(
        "--db2",
        required=True,
        help="Path to second database (FASTA format)",
    )
    parser.add_argument(
        "--nhmmer",
        default=None,
        help="Path to nhmmer binary (default: search in PATH)",
    )
    parser.add_argument(
        "--hmmalign",
        default=None,
        help="Path to hmmalign binary (default: search in PATH)",
    )
    parser.add_argument(
        "--hmmbuild",
        default=None,
        help="Path to hmmbuild binary (default: search in PATH)",
    )
    parser.add_argument(
        "--n_cpu",
        type=int,
        default=8,
        help="Number of CPUs for nhmmer (default: 8)",
    )
    parser.add_argument(
        "--e_value",
        type=float,
        default=1e-3,
        help="E-value threshold for nhmmer (default: 1e-3)",
    )
    parser.add_argument(
        "--max_sequences",
        type=int,
        default=10000,
        help="Maximum sequences per database (default: 10000)",
    )

    args = parser.parse_args()

    nhmmer_binary = args.nhmmer or find_binary("nhmmer")
    hmmalign_binary = args.hmmalign or find_binary("hmmalign")
    hmmbuild_binary = args.hmmbuild or find_binary("hmmbuild")

    if not nhmmer_binary:
        print("ERROR: nhmmer not found. Specify with --nhmmer or add to PATH.")
        sys.exit(1)
    if not hmmalign_binary:
        print("ERROR: hmmalign not found. Specify with --hmmalign or add to PATH.")
        sys.exit(1)
    if not hmmbuild_binary:
        print("ERROR: hmmbuild not found. Specify with --hmmbuild or add to PATH.")
        sys.exit(1)

    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)
    if not os.path.exists(args.db1):
        print(f"ERROR: Database 1 not found: {args.db1}")
        sys.exit(1)
    if not os.path.exists(args.db2):
        print(f"ERROR: Database 2 not found: {args.db2}")
        sys.exit(1)

    print("=" * 60)
    print("RNA to A3M Converter (AlphaFold3-style)")
    print("=" * 60)
    print(f"Input:       {args.input}")
    print(f"Output:      {args.output_dir}")
    print(f"Database 1:  {args.db1}")
    print(f"Database 2:  {args.db2}")
    print(f"nhmmer:      {nhmmer_binary}")
    print(f"hmmalign:    {hmmalign_binary}")
    print(f"hmmbuild:    {hmmbuild_binary}")
    print(f"E-value:     {args.e_value}")
    print(f"Max seqs:    {args.max_sequences}")
    print(f"CPUs:        {args.n_cpu}")
    print("=" * 60)

    process_rna_sequences(
        input_fasta=args.input,
        output_dir=args.output_dir,
        db_paths=[args.db1, args.db2],
        nhmmer_binary=nhmmer_binary,
        hmmalign_binary=hmmalign_binary,
        hmmbuild_binary=hmmbuild_binary,
        n_cpu=args.n_cpu,
        e_value=args.e_value,
        max_sequences=args.max_sequences,
    )


if __name__ == "__main__":
    main()
