#!/usr/bin/env python3
"""Package Teeth3DS scans for ToothGroupNetwork inference (D81).

Selects the scans of one Teeth3DS part that are in the 3DTeethSeg'22 PRIVATE TEST list (the checkpoints were
trained on the challenge training lists, so only these are held out for the network), lays them out as
<case>/<case>_<jaw>.obj, and writes a zip under data/ (git-ignored; Teeth3DS is CC BY-NC-ND 4.0: keep it private).

    python3 tools/tgn/prepare_input.py data/data_part_6 data/tgn_input_part6.zip
"""
import os
import sys
import zipfile

SPLIT = "data/xctdy-osfstorage-archive/3DTeethSeg22_challenge_train_test_split/private-testing-set.txt"


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    part, out = sys.argv[1], sys.argv[2]
    with open(SPLIT) as f:
        held_out = {line.strip() for line in f if line.strip()}
    chosen = []
    for root, _, files in os.walk(part):
        for name in files:
            stem, ext = os.path.splitext(name)
            if ext == ".obj" and stem in held_out:
                chosen.append((stem, os.path.join(root, name)))
    chosen.sort()
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for stem, path in chosen:
            case = stem.split("_")[0]
            z.write(path, f"input/{case}/{stem}.obj")
    print(f"{len(chosen)} held-out scans -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
