"""Report text files that mix CRLF and bare-LF line endings.

Most files in this repo are CRLF in the working tree, and an edit that writes
bare LF into one leaves it mixed. A file that is consistently CRLF or
consistently LF passes; binary files (any NUL byte) are skipped.

Usage:
    python tools/check_line_endings.py PATH [PATH ...]

With no paths, reads a Claude Code hook payload (JSON) from stdin and checks
its ``tool_input.file_path``. Exits 2 if any file is mixed, with the details
on stderr, so a PostToolUse hook shows them to Claude; exits 0 otherwise.
"""

import json
import sys
from pathlib import Path


def line_ending_counts(path):
    """Return ``(crlf, bare_lf)`` for a text file, or ``None`` to skip it."""
    try:
        data = Path(path).read_bytes()
    except OSError:
        return None
    if b"\0" in data:
        return None
    crlf = data.count(b"\r\n")
    return crlf, data.count(b"\n") - crlf


def paths_from_hook_payload(stream):
    try:
        payload = json.load(stream)
    except ValueError:
        return []
    file_path = payload.get("tool_input", {}).get("file_path")
    return [file_path] if file_path else []


def main(argv):
    paths = argv or paths_from_hook_payload(sys.stdin)
    mixed = []
    for path in paths:
        counts = line_ending_counts(path)
        if counts and counts[0] and counts[1]:
            mixed.append((path, *counts))
    for path, crlf, bare_lf in mixed:
        print(
            f"Mixed line endings in {path}: {crlf} CRLF, {bare_lf} bare LF. "
            "Restore the file's original ending on the edited lines.",
            file=sys.stderr,
        )
    return 2 if mixed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
