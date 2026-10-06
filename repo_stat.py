"""Contributor statistics for one or more git repositories.

Produces per-author commit and line counts with the provenance needed to
reproduce them, and reports the identity-folding problems that make naive
`git shortlog` numbers wrong.

Design notes:

  * The mailmap is passed via `git -c mailmap.file=...` rather than written
    into the repository, so no working tree is modified.

  * Git's mailmap has NO wildcard syntax. An entry of the form
    `Proper Name <proper@email> Commit Name <*>` matches a commit whose email
    is the literal string `*`, i.e. nothing, and is silently ignored. Folding
    is therefore driven off explicit email lists only.

  * Two reports surface identities the mailmap missed: an audit of unfolded
    authors matching the canonical person, and a list of any email address
    used under more than one display name. Populating the config from those
    two reports and re-running is the intended workflow.

  * Statistics are produced under two ref scopes: `--all` (every ref in the
    clone, including unmerged branches and any rebased duplicates) and the
    default branch alone (comparable to a forge's contributors page). The two
    disagree; both are written.

  * Line counts are reported raw and filtered. Vendored, generated, and bulk
    data files dominate raw counts and make them meaningless as a measure of
    authored work.

If you are running this for me: set REPO_PATHS below to the repositories you
want scanned, run `python repo_stats_generic.py`, and send back the CSV and
JSON files it writes plus whatever it printed to the terminal. Nothing else
needs configuring, and nothing needs re-running — the reports it prints are
notes for me, not tasks for you.

It scans every branch and tag present in the clone, then scans the default
branch on its own, and writes a separate file for each so the two can be
compared. It is read-only with respect to the repository: it checks nothing
out, creates and deletes no branches, and never writes into the repository or
its working tree. Its only writes are the output files, which go to
OUTPUT_DIR.

Note that it can only see refs that are actually in the clone. If the co-worker
running it has a shallow or single-branch clone, the all-branches pass will be
incomplete; a full `git clone` with `git fetch --all --tags` beforehand avoids
that.
"""

import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

REPO_PATHS = [
    ".",
]

# Directory the CSV and provenance files are written to.
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# The person the run is centred on. Their rank is echoed to the console after
# each scope, and the identity audit looks for unfolded variants of them.
CANONICAL_NAME = "Juliano Gianlupi"
CANONICAL_EMAIL = "julianoferrarigianlupi@gmail.com"

# Every email address the canonical person has committed under. Folding is by
# email; git cannot match on name alone, so a commit made under an address not
# listed here lands in its own row instead of being merged.
#
# todo: the identity audit at the end of each run prints any author whose name
#       or address looks like the canonical person but was not folded. Add
#       anything it surfaces here — internal repos may carry a corporate SSO
#       or forge-generated noreply address that is not yet on this list.
EMAIL_VARIANTS = [
    "julianoferrarigianlupi@gmail.com",
    "jferrarigian@unither.com",
    "35964356+JulianoGianlupi@users.noreply.github.com",
]

# Other contributors whose commits are split across several display names.
# Maps a chosen proper name to every address that should fold into it; the
# first address in each list is used as the proper address.
#
# Do not fold shared team, bot, or service accounts into an individual. Who
# used a shared account is not recoverable from commit metadata, and merging
# one into a person asserts something the data does not support. Leave such
# accounts as their own rows.
#
# todo: leave empty on a first pass. The shared-email report at the end of each
#       run lists every address appearing under more than one display name,
#       which is the raw material for filling this in later if the split
#       identities turn out to matter.
OTHER_IDENTITIES = {}

# Substrings used only to flag unfolded authors for manual review. Derived
# from CANONICAL_NAME when left empty; override for nicknames, transliterated
# spellings, or machine account names that would not match the tokens.
IDENTITY_HINTS = []

# Neutral ordering by default. A table intended to be shown to anyone else
# should not be sorted so that the canonical person occupies row one; the rank
# column and the console echo make their row easy to find without rearranging
# the table around them.
PIN_CANONICAL_FIRST = False

# Paths excluded from the filtered line counts. The fragments below are the
# generic offenders. Every repository has its own, so derive the rest from a
# per-file scan of the history:
#
#   git log --all --no-merges --numstat --format= \
#     | awk -F'\t' '$1 ~ /^[0-9]+$/ {s[$3]+=$1} END {for (f in s) print s[f], f}' \
#     | sort -rn | head -30
#
# todo: tune per repository if the filtered counts still look inflated. Run the
#       scan above, add whatever dominates the top of the list, repeat until
#       the top is authored source. Common culprits beyond the list below are
#       generated bindings and resource modules, rendered API documentation,
#       and large fixture or data files. Not required for a first pass — the
#       raw and filtered counts are both reported either way.
EXCLUDED_PATH_FRAGMENTS = [
    # Vendored third-party trees.
    "/third_party/",
    "/thirdparty/",
    "/3rdparty/",
    "/vendor/",
    "/external/",
    "/node_modules/",
    "/site-packages/",
    "/.venv/",
    # Build products and generated bindings.
    "/build/",
    "/dist/",
    "/target/",
    "/_generated/",
    "/generated/",
    "_pb2.py",
    "_pb2_grpc.py",
    "_wrap.cxx",
    "_wrap.cpp",
    ".min.js",
    ".min.css",
    # Dependency lockfiles: machine-written and often enormous.
    "package-lock.json",
    "yarn.lock",
    "pnpm-lock.yaml",
    "poetry.lock",
    "cargo.lock",
    "gemfile.lock",
    "composer.lock",
    # Shadow directories left by migrations from older version control.
    "/.svn/",
    "/cvs/",
]

# Extensions to ignore globally regardless of directory
EXCLUDED_EXTENSIONS = (
    # Data, Spreadsheets & Notebooks
    ".csv", ".tsv", ".ods", ".xlsx", ".xls", ".parquet", ".feather", ".json", ".ipynb",
    # Documents, Media & Design
    ".pdf", ".tif", ".tiff", ".png", ".jpg", ".jpeg", ".gif", ".svg", ".eps", ".psd", ".ai",
    # Generated / Minified Code & Maps
    ".min.js", ".min.css", ".map", ".bundle.js",
    # Databases & Dumps
    ".sqlite", ".sqlite3", ".db", ".sql", ".dump",
    # Archives, Binaries, Executables & Java Artifacts
    ".zip", ".tar", ".gz", ".7z", ".rar", ".exe", ".dll", ".so", ".dylib", ".bin", ".pif",
    ".class", ".jar", ".war", ".ear", ".jmod",
    # Logs & Lockfiles
    ".log", ".logs", ".lock",
)

# Extensions that are ordinary source elsewhere in a tree but are data when
# they appear under one of these directories.
DATA_DIR_FRAGMENTS = [
    "/fixtures/",
    "/testdata/",
    "/test_data/",
    "/snapshots/",
    "/__snapshots__/",
]
DATA_EXTENSIONS = [".txt", ".csv", ".tsv", ".json", ".xml", ".sql"]

# ---------------------------------------------------------------------------


def identity_hints():
    if IDENTITY_HINTS:
        return [h.lower() for h in IDENTITY_HINTS]
    return [t.lower() for t in CANONICAL_NAME.split() if len(t) > 2]


def run(cmd, cwd=None):
    """Run a command and return stripped stdout, or None on failure."""
    try:
        out = subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True, check=True
        )
        return out.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def write_mailmap():
    """Write a temporary mailmap outside the repository.

    Only the form git actually supports is emitted:
      Proper Name <proper@email> <commit@email>
    which folds any commit using `commit@email`, whatever name it carries.
    """
    fd, path = tempfile.mkstemp(suffix=".mailmap", text=True)
    with os.fdopen(fd, "w") as f:
        if CANONICAL_NAME and EMAIL_VARIANTS:
            primary = CANONICAL_EMAIL or EMAIL_VARIANTS[0]
            for email in EMAIL_VARIANTS:
                f.write(f"{CANONICAL_NAME} <{primary}> <{email}>\n")
        for proper_name, emails in OTHER_IDENTITIES.items():
            if not emails:
                continue
            primary = emails[0]
            for email in emails:
                f.write(f"{proper_name} <{primary}> <{email}>\n")
    return path


def is_excluded(path):
    lowered = "/" + path.lower()
    if lowered.endswith(EXCLUDED_EXTENSIONS):
        return True
    if any(frag in lowered for frag in EXCLUDED_PATH_FRAGMENTS):
        return True
    if any(frag in lowered for frag in DATA_DIR_FRAGMENTS) and any(
        lowered.endswith(ext) for ext in DATA_EXTENSIONS
    ):
        return True
    return False


def collect(repo_path, mailmap_path, ref_scope):
    """Walk the log once, accumulating per-author statistics.

    ref_scope is either "--all" or a branch name.
    """
    cmd = [
        "git",
        "-c",
        f"mailmap.file={mailmap_path}",
        "log",
        ref_scope,
        "--no-merges",
        "--ignore-submodules",
        "--use-mailmap",
        "--numstat",
        # %aN is mailmapped (so identities fold), %ae is the raw address (so
        # the emails column shows what was actually folded, rather than
        # echoing the proper address back). Unit separators delimit, so author
        # names containing pipes or tabs cannot corrupt the parse.
        "--format=C%x1f%aN%x1f%ae%x1f%at",
    ]

    stats = {}
    author = None
    commit_count = 0

    process = subprocess.Popen(
        cmd,
        cwd=repo_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    for line in process.stdout:
        line = line.rstrip("\n")
        if not line:
            continue

        if line.startswith("C\x1f"):
            _, name, email, ts = line.split("\x1f")
            author = name
            timestamp = int(ts)
            commit_count += 1
            if commit_count % 500 == 0:
                print(f"    {commit_count:,} commits...", end="\r")

            rec = stats.setdefault(
                author,
                {
                    "contributor name": author,
                    "emails": set(),
                    "number of commits": 0,
                    "number of lines added": 0,
                    "number of lines removed": 0,
                    "lines added (filtered)": 0,
                    "lines removed (filtered)": 0,
                    "binary files touched": 0,
                    "first_ts": timestamp,
                    "last_ts": timestamp,
                },
            )
            rec["emails"].add(email)
            rec["number of commits"] += 1
            rec["first_ts"] = min(rec["first_ts"], timestamp)
            rec["last_ts"] = max(rec["last_ts"], timestamp)
            continue

        if author is None:
            continue

        parts = line.split("\t")
        if len(parts) < 3:
            continue
        added, removed, path = parts[0], parts[1], parts[2]

        if added == "-" or removed == "-":
            # Binary file; numstat reports no line counts for these.
            stats[author]["binary files touched"] += 1
            continue
        if not (added.isdigit() and removed.isdigit()):
            continue

        stats[author]["number of lines added"] += int(added)
        stats[author]["number of lines removed"] += int(removed)
        if not is_excluded(path):
            stats[author]["lines added (filtered)"] += int(added)
            stats[author]["lines removed (filtered)"] += int(removed)

    process.wait()
    stderr = process.stderr.read().strip()
    if process.returncode != 0:
        raise RuntimeError(f"git log failed: {stderr}")

    print(f"    {commit_count:,} commits total.          ")
    return stats


def to_frame(stats):
    rows = []
    for data in stats.values():
        rows.append(
            {
                "contributor name": data["contributor name"],
                "emails": "; ".join(sorted(data["emails"])),
                "number of commits": data["number of commits"],
                "number of lines added": data["number of lines added"],
                "number of lines removed": data["number of lines removed"],
                "lines added (filtered)": data["lines added (filtered)"],
                "lines removed (filtered)": data["lines removed (filtered)"],
                "binary files touched": data["binary files touched"],
                "first commit date": datetime.fromtimestamp(
                    data["first_ts"], tz=timezone.utc
                ).strftime("%Y-%m-%d"),
                "last commit date": datetime.fromtimestamp(
                    data["last_ts"], tz=timezone.utc
                ).strftime("%Y-%m-%d"),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Neutral ordering: commits descending, name ascending to break ties.
    df = df.sort_values(
        by=["number of commits", "contributor name"],
        ascending=[False, True],
    ).reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))

    if PIN_CANONICAL_FIRST and CANONICAL_NAME:
        df["_pin"] = df["contributor name"] != CANONICAL_NAME
        df = df.sort_values(by="_pin", kind="stable").drop(columns=["_pin"])
        df = df.reset_index(drop=True)

    return df


def audit_identities(df):
    """Unfolded authors that look like the canonical person.

    Any row returned is a commit set the mailmap failed to capture, almost
    always because the address is missing from EMAIL_VARIANTS.
    """
    hints = identity_hints()
    if not hints:
        return []
    flagged = []
    for _, row in df.iterrows():
        if row["contributor name"] == CANONICAL_NAME:
            continue
        haystack = f"{row['contributor name']} {row['emails']}".lower()
        if any(hint in haystack for hint in hints):
            flagged.append(row.to_dict())
    return flagged


def shared_email_report(df):
    """Emails appearing under more than one display name.

    Each is an identity the mailmap has not folded. Populating
    OTHER_IDENTITIES from this list merges the remaining split contributors.
    """
    by_email = {}
    for _, row in df.iterrows():
        for email in row["emails"].split("; "):
            if email:
                by_email.setdefault(email, set()).add(row["contributor name"])
    return {e: sorted(n) for e, n in by_email.items() if len(n) > 1}


def default_branch(repo_path):
    head = run(
        ["git", "symbolic-ref", "--short", "refs/remotes/origin/HEAD"],
        cwd=repo_path,
    )
    if head:
        return head  # e.g. "origin/main"
    for candidate in ("origin/main", "origin/master", "main", "master"):
        if run(["git", "rev-parse", "--verify", candidate], cwd=repo_path):
            return candidate
    return "HEAD"


def provenance(repo_path, ref_scope, cmd_note):
    dirty = run(["git", "status", "--porcelain"], cwd=repo_path)
    return {
        "repository path": os.path.abspath(repo_path),
        "remote origin": run(
            ["git", "remote", "get-url", "origin"], cwd=repo_path
        ),
        "HEAD sha": run(["git", "rev-parse", "HEAD"], cwd=repo_path),
        "HEAD date": run(
            ["git", "log", "-1", "--format=%cI", "HEAD"], cwd=repo_path
        ),
        "ref scope": ref_scope,
        "working tree clean": (dirty == ""),
        "commits reachable from all refs": run(
            ["git", "rev-list", "--all", "--count"], cwd=repo_path
        ),
        "git version": run(["git", "--version"]),
        "run at (UTC)": datetime.now(timezone.utc).isoformat(),
        "command": cmd_note,
        "identities folded": {
            **({CANONICAL_NAME: EMAIL_VARIANTS} if CANONICAL_NAME else {}),
            **OTHER_IDENTITIES,
        },
        "note": (
            "Counts are by commit author, not committer, and exclude merge "
            "commits and submodules. Binary files are reported as a file "
            "count, not as lines, because numstat reports no line counts for "
            "them. The 'filtered' columns additionally exclude the vendored, "
            "generated, and bulk-data paths configured for this run. The "
            "'--all' scope counts commits reachable from any ref in the "
            "clone, which includes unmerged branches and any rebased "
            "duplicates of the same change; the default-branch scope does "
            "not. Author identities are folded by email address only; see "
            "the identity audit and shared-email report for any that were "
            "not folded."
        ),
    }


def report(df, label):
    if CANONICAL_NAME:
        me = df[df["contributor name"] == CANONICAL_NAME]
        if not me.empty:
            r = me.iloc[0]
            print(
                f"    {CANONICAL_NAME}: rank {r['rank']} of {len(df)}"
                f"  ·  {r['number of commits']} commits"
                f"  ·  {r['lines added (filtered)']} lines added (filtered),"
                f" {r['number of lines added']} raw"
                f"  ·  {r['first commit date']} to {r['last commit date']}"
            )
        else:
            print(
                f"    {CANONICAL_NAME}: no commits found under this ref scope."
            )

    if label != "all-refs":
        return

    flagged = audit_identities(df)
    if flagged:
        print(
            "\n  IDENTITY AUDIT — unfolded authors matching the canonical "
            "person:"
        )
        for row in flagged:
            print(
                f"    {row['contributor name']!r} <{row['emails']}>  "
                f"{row['number of commits']} commits  "
                f"{row['first commit date']} to {row['last commit date']}"
            )
        print(
            "  Note: these are unmerged identities. No action needed here; "
            "the raw rows are in the CSV.\n"
        )
    elif CANONICAL_NAME:
        print("  Identity audit: no unfolded variants found.")

    shared = shared_email_report(df)
    if shared:
        print(
            "\n  SHARED EMAILS — one address, several display names (each is "
            "an unfolded identity):"
        )
        for email, names in sorted(shared.items()):
            print(f"    {email}: {', '.join(names)}")
        print(
            "  Note: these contributors appear under more than one row. No "
            "action needed here; the raw rows are in the CSV.\n"
        )
    else:
        print(
            "  Shared-email report: no address used under more than one "
            "name.\n"
        )


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    mailmap_path = write_mailmap()
    try:
        for repo_path in REPO_PATHS:
            if not os.path.isdir(os.path.join(repo_path, ".git")):
                print(f"Skipping (not a git repository): {repo_path}")
                continue

            name = os.path.basename(os.path.abspath(repo_path)) or "repo"
            branch = default_branch(repo_path)
            print(f"\nRepository: {name}   default branch: {branch}")

            for label, scope in (
                ("all-refs", "--all"),
                ("default-branch", branch),
            ):
                print(f"  scope: {scope}")
                stats = collect(repo_path, mailmap_path, scope)
                df = to_frame(stats)
                if df.empty:
                    print("    No commits found under this ref scope.")
                    continue

                csv_path = os.path.join(OUTPUT_DIR, f"{name}_stats_{label}.csv")
                df.to_csv(csv_path, index=False)

                meta = provenance(
                    repo_path,
                    scope,
                    f"git log {scope} --no-merges --ignore-submodules "
                    f"--use-mailmap --numstat",
                )
                meta_path = os.path.join(
                    OUTPUT_DIR, f"{name}_provenance_{label}.json"
                )
                with open(meta_path, "w") as f:
                    json.dump(meta, f, indent=2)
                    f.write("\n")

                print(f"    wrote {csv_path}")
                print(f"    wrote {meta_path}")
                report(df, label)
    finally:
        os.remove(mailmap_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())

