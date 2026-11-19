"""cofactor: compare two results files row by row (porting check).
Rows "<name> n= <n> max_rel=<x>  at <where>[ | verdict]". Identical = same row names, same n, same max_rel text.
Rows only in the new file are allowed if they are cofactor's own (EXCLUDED: row 0 failed..., REVERSE SKIPPED: ...).
usage: python3 compare_results.py <old> <new>   (exit 0 identical, 1 different)"""
import re
import sys

ROW = re.compile(r"^(.*?)\s+n=\s*(\d+)\s+max_rel=(\S+)")
OWN = ("EXCLUDED: row 0 failed", "REVERSE SKIPPED:")


def rows(path):
    out = {}
    for line in open(path):
        m = ROW.match(line.rstrip("\n"))
        if m and not line.startswith("MISMATCH"):
            out[m.group(1).strip()] = (m.group(2), m.group(3))
    return out


def main():
    old, new = rows(sys.argv[1]), rows(sys.argv[2])
    bad = 0
    for k in sorted(set(old) | set(new)):
        if k not in new:
            print(f"MISSING in new: {k}")
            bad += 1
        elif k not in old:
            if not k.startswith(OWN):
                print(f"EXTRA in new: {k}")
                bad += 1
        elif old[k] != new[k]:
            print(f"DIFF {k}: old n={old[k][0]} max_rel={old[k][1]}  new n={new[k][0]} max_rel={new[k][1]}")
            bad += 1
    print(f"compare: {len(old)} old rows, {len(new)} new rows, {bad} differences")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
