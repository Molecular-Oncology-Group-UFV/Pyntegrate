"""
Static audit of the figure-generating scripts in this folder for Fig-7-class bugs.

The Fig 7 bug was a single unsuffixed `CHIP_COL` global consumed by all four panel
calls, so panels C/D plotted SSRP1 occupancy under SUPT16H titles. Nothing raised an
error. This script detects that signature and its neighbours by parsing the AST of
every script in the folder -- no pipeline code is imported and no figure is rendered.

Run:  python audit_figure_scripts.py [--table PATH]
Exit: 0 if no finding, 1 if any finding.
"""
import argparse, ast, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TABLE = os.path.join(os.path.dirname(HERE), "data", "final_table_filtered.tsv")
SUFFIX = re.compile(r"_(A|B|C|D|AB|CD|ABCD|1|2)$")
COLARG = re.compile(r"col$", re.I)


def subunit(s):
    """Which FACT subunit a string refers to. 'suptt16h' is the real (misspelled) header."""
    t = str(s).lower()
    has_ss = "ssrp1" in t
    has_su = "supt16h" in t or "suptt16h" in t
    if has_ss and has_su:
        return "BOTH"
    return "SSRP1" if has_ss else ("SUPT16H" if has_su else None)


def module_constants(tree, data_dir=None):
    """Module-level constants, with os.path.join(DATA_DIR, ...) resolved to real paths."""
    out = {}

    def ev(node):
        try:
            return ast.literal_eval(node)
        except Exception:
            pass
        # os.path.join(DATA_DIR, "sub", "file.bw")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "join":
            parts = []
            for a in node.args:
                if isinstance(a, ast.Constant) and isinstance(a.value, str):
                    parts.append(a.value)
                elif isinstance(a, ast.Name) and a.id == "DATA_DIR" and data_dir:
                    parts.append(data_dir)
                elif isinstance(a, ast.Name) and isinstance(out.get(a.id), str):
                    parts.append(out[a.id])
                else:
                    return "<expr>"
            return os.path.join(*parts)
        return "<expr>"

    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            out[node.targets[0].id] = ev(node.value)
    return out


def pipeline_calls(tree, consts):
    """Every *Pipeline(...) instantiation, with keyword args resolved to literals."""
    def resolve(node):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            return consts.get(node.id)
        try:
            return ast.literal_eval(node)
        except Exception:
            return None

    calls = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = node.func.attr if isinstance(node.func, ast.Attribute) else getattr(node.func, "id", "")
        if "Pipeline" not in name:
            continue
        args, raw = {}, {}
        for kw in node.keywords:
            raw[kw.arg] = kw.value
            val = resolve(kw.value)
            if val is not None:
                args[kw.arg] = val
        calls.append({"line": node.lineno, "name": name, "args": args, "raw": raw})
    return calls


def audit(path, headers, DATA_DIR=None):
    """Return a list of finding strings for one script."""
    src = open(path).read()
    tree = ast.parse(src)
    consts = module_constants(tree, data_dir=DATA_DIR)
    calls = pipeline_calls(tree, consts)
    findings = []   # (severity, message); severity is 'BUG' or 'WARN'

    # 1. Argument name implies one subunit, value is the other.
    for c in calls:
        for k, v in c["args"].items():
            if not isinstance(v, str):
                continue
            implied = "SSRP1" if re.search(r"ssrp1", k, re.I) else ("SUPT16H" if re.search(r"supt16h", k, re.I) else None)
            actual = subunit(v)
            if implied and actual and implied != actual:
                findings.append(("BUG", f"L{c['line']}: arg {k} implies {implied} but value is {actual}: {v!r}"))

    # 2. Panel suffixes mixed inside one call (CHIP_COL_A next to RNA_COL_B).
    #    A group suffix legitimately covers its members: CHIP_COL_AB alongside TITLE_A is
    #    correct, because panels A and B share that column and differ only in title/filter.
    #    Suffixes _1/_2 label the two operands of a correlation plot, not two panels.
    for c in calls:
        seen = {}
        for k, node in c["raw"].items():
            if isinstance(node, ast.Name):
                m = SUFFIX.search(node.id)
                if m and m.group(1) not in ("1", "2"):
                    seen.setdefault(m.group(1), []).append(node.id)
        letters = set(seen)
        incompatible = {a for a in letters
                        if not any(b != a and set(a) <= set(b) for b in letters)}
        if len(incompatible) > 1:
            findings.append(("BUG", f"L{c['line']}: mixes panel suffixes "
                            f"{ {k: v for k, v in seen.items() if k in incompatible} }"))

    # 3. THE FIG 7 SIGNATURE: one unsuffixed, subunit-specific column reused by >1 call.
    if len(calls) > 1:
        shared = {}
        for c in calls:
            for k, node in c["raw"].items():
                if not COLARG.search(k or ""):
                    continue
                key = node.id if isinstance(node, ast.Name) else (repr(node.value) if isinstance(node, ast.Constant) else None)
                if key and not SUFFIX.search(key):
                    shared.setdefault(key, []).append((c["line"], c["args"].get(k)))
        for key, uses in shared.items():
            if len(uses) < 2:
                continue
            vals = {subunit(v) for _, v in uses if v}
            if not (vals & {"SSRP1", "SUPT16H"}):
                continue
            col_sub = next(iter(vals & {"SSRP1", "SUPT16H"}))
            lines = [ln for ln, _ in uses]
            sharers = [c for c in calls if c["line"] in lines]
            # Sharing one subunit-specific column across panels is CORRECT when those
            # panels are all about that subunit and differ in some other way (Fig S2
            # shares 'SSRP1 FACT category' between its high and low SSRP1 panels).
            # It is the Fig 7 BUG when the sharing panels claim DIFFERENT subunits:
            # there, one SSRP1 occupancy column fed panels titled SSRP1 and SUPT16H.
            claims = {subunit(c["args"].get("title", "")) or subunit(c["args"].get("output_path", ""))
                      for c in sharers}
            claims.discard(None)
            if len(claims) > 1 or (claims and col_sub not in claims):
                findings.append(("BUG", f"{key}: {col_sub} column shared by {len(uses)} calls at lines {lines} "
                                f"whose titles claim {sorted(claims)} -- this is the Fig 7 signature"))

    # 4. Title/output claims one subunit, the data columns are the other.
    for c in calls:
        claimed = subunit(c["args"].get("title", "")) or subunit(c["args"].get("output_path", ""))
        actual = {subunit(v) for k, v in c["args"].items()
                  if isinstance(v, str) and v in headers and subunit(v)}
        stray = actual - {claimed} if claimed and claimed != "BOTH" else set()
        if stray:
            wrong = sorted({v for k, v in c["args"].items()
                            if isinstance(v, str) and v in headers and subunit(v) in stray})
            findings.append(("BUG", f"L{c['line']}: title/output claims {claimed} but these columns are "
                            f"{sorted(stray)}: {wrong}"))

    # 5. Column-shaped constants that are not real headers.
    if headers:
        for k, v in consts.items():
            if isinstance(v, str) and re.search(r"_COL$|_COL_", k) and v not in headers:
                findings.append(("BUG", f"{k} = {v!r} is not a column of the table"))

    # 6. Declared but never referenced (a renamed constant the calls no longer read).
    used = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    for k in consts:
        if re.search(r"_COL$|_COL_|_PATH$|_PATH_", k) and k not in used:
            findings.append(("WARN", f"{k} is declared but never used"))

    # 7. Output path outside the shared figures directory (stale-figure risk).
    for k, v in consts.items():
        if isinstance(v, str) and v.startswith("outputs/") and os.path.dirname(v) != "outputs/figs":
            findings.append(("BUG", f"{k} = {v!r} writes outside outputs/figs -- a rerun leaves the old file behind"))

    # 8. Declared input files that do not exist.
    for k, v in consts.items():
        for x in (v if isinstance(v, list) else [v]):
            if isinstance(x, str) and re.search(r"\.(bed|bed\.gz|bw|bigWig|tsv|csv|gtf|db)$", x) \
                    and not x.startswith("outputs/") and not os.path.exists(x):
                findings.append(("WARN", f"{k}: input file missing -- {x}"))

    # 9. Inline spec dicts (e.g. Fig 8A's VIOLIN_SPECS): every column inside one dict
    #    must refer to the same subunit as that entry's title.
    for name, val in consts.items():
        if not (isinstance(val, list) and val and all(isinstance(d, dict) for d in val)):
            continue
        for i, spec in enumerate(val):
            subs = {}
            for k, v in spec.items():
                if isinstance(v, str) and (not headers or v in headers):
                    s = subunit(v)
                    if s:
                        subs.setdefault(s, []).append(f"{k}={v!r}")
            claimed = subunit(spec.get("title", ""))
            if len(subs) > 1:
                findings.append(("BUG", f"{name}[{i}] ({spec.get('title','?')}): mixes subunits {dict(subs)}"))
            elif claimed and subs and claimed not in subs:
                findings.append(("BUG", f"{name}[{i}]: title claims {claimed} but columns are {list(subs)}"))

    return findings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=DEFAULT_TABLE)
    args = ap.parse_args()

    headers = set()
    if os.path.exists(args.table):
        with open(args.table) as fh:
            headers = set(fh.readline().rstrip("\n").split("\t"))
    else:
        print(f"note: table not found at {args.table}; header checks skipped\n")

    scripts = sorted(f for f in os.listdir(HERE) if f.endswith(".py") and f != os.path.basename(__file__))
    results = {s: audit(os.path.join(HERE, s), headers,
                        DATA_DIR=os.path.join(os.path.dirname(HERE), "data")) for s in scripts}

    nbug = nwarn = 0
    print("=" * 70)
    print("BUGS  (wrong data silently plotted, or a stale figure silently kept)")
    print("=" * 70)
    for s, found in results.items():
        bugs = [m for sev, m in found if sev == "BUG"]
        nbug += len(bugs)
        if bugs:
            print(f"  {s}")
            for m in bugs:
                print(f"      {m}")
    if not nbug:
        print("  none")

    print()
    print("=" * 70)
    print("WARNINGS  (reproducibility only -- absent inputs, unused constants)")
    print("=" * 70)
    for s, found in results.items():
        warns = [m for sev, m in found if sev == "WARN"]
        nwarn += len(warns)
        if warns:
            missing = [m for m in warns if "input file missing" in m]
            other = [m for m in warns if "input file missing" not in m]
            print(f"  {s}: {len(missing)} missing input file(s)" if missing else f"  {s}")
            for m in other:
                print(f"      {m}")
    if not nwarn:
        print("  none")

    print(f"\n{len(scripts)} scripts audited: {nbug} bug(s), {nwarn} warning(s).")
    return 1 if nbug else 0


if __name__ == "__main__":
    sys.exit(main())
