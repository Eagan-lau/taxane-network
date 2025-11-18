#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import pandas as pd
from collections import OrderedDict, defaultdict

# ================= Configuration =================
ALL_CSV = "all.csv"                    # Column 2 = GeneID, Column 4 = Family; first row is header
GENE_CLUSTER_TXT = "gene_cluster.txt"  # Each line represents a group; elements are enclosed in quotes
OUT_CSV = "group_members_labeled.csv"  # Output table containing labeled group members
OUT_MISSING = "seeds_not_in_any_group.txt"  # Seed genes not found in any group
DEDUPLICATE = True   # If True, remove duplicated rows (same group + family + member)
# =================================================


def load_seed_family_map(path: str) -> pd.DataFrame:
    """
    Load all.csv and extract the seed gene ID (column 2) and family (column 4).
    Uses positional indexing to avoid dependency on column names.
    """
    df = pd.read_csv(path)

    # iloc[:,1] = GeneID; iloc[:,3] = Family
    sub = df.iloc[:, [1, 3]].copy()
    sub.columns = ["GeneID", "Family"]

    sub["GeneID"] = sub["GeneID"].astype(str).str.strip()
    sub["Family"] = sub["Family"].astype(str).str.strip()

    sub = sub.dropna(subset=["GeneID", "Family"])
    return sub


def parse_gene_clusters(path: str):
    """
    Parse gene_cluster.txt.  
    Each line corresponds to one group. Elements in the group are detected by:
        - Capturing any string enclosed in 'single' or "double" quotes
        - If no quoted elements are present, fallback to splitting the line by whitespace
    Returns:
        groups:      list of (group_id, [member1, member2, ...])
        elem2group:  mapping from element → group_id  
                     If an element appears in multiple groups, the *first* occurrence is used.
    """
    groups = []
    elem2group = {}

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f, start=1):
            # Capture strings enclosed by either ' or "
            elems = re.findall(r"""(['"])(.*?)\1""", line)

            # Extract the content part (elem[1])
            members = [m[1].strip() for m in elems if m[1].strip() != ""]

            # Fallback: whitespace split if no quoted strings found
            if not members:
                tokens = re.split(r"[\s,]+", line.strip())
                tokens = [t for t in tokens if t and not t.endswith(":")]
                members = tokens

            gid = f"group_{i}"

            # Remove duplicates while preserving order
            seen = OrderedDict()
            for m in members:
                if m not in seen:
                    seen[m] = True
            members = list(seen.keys())

            groups.append((gid, members))

            # Map element → group (only first occurrence retained)
            for m in members:
                if m not in elem2group:
                    elem2group[m] = gid

    return groups, elem2group


def main():
    # Load seed → family table
    seeds = load_seed_family_map(ALL_CSV)

    # Parse clustering file
    groups, elem2group = parse_gene_clusters(GENE_CLUSTER_TXT)
    group_map = {gid: members for gid, members in groups}

    rows = []
    missing = []

    """
    For each seed gene:
      - Identify its group via elem2group
      - For that group, label all members with the seed's family annotation
    This effectively propagates the seed's family label to all genes appearing
    in the same cluster.
    """
    for _, r in seeds.iterrows():
        seed_id = r["GeneID"]
        fam = r["Family"]
        gid = elem2group.get(seed_id)

        if gid is None:
            # Seed ID did not appear in any group
            missing.append(seed_id)
            continue

        members = group_map.get(gid, [])
        for mem in members:
            rows.append({
                "group": gid,
                "seed_gene": seed_id,   # Seed gene that triggered this group
                "family": fam,          # Seed's family annotation
                "member_gene": mem      # Member gene in the same group
            })

    out_df = pd.DataFrame(rows)

    # Remove duplicated lines caused by multiple seeds mapping to the same group
    if DEDUPLICATE and not out_df.empty:
        out_df = out_df.drop_duplicates(
            subset=["group", "family", "member_gene"],
            keep="first"
        )

    out_df.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")

    # Write seed genes that were not found in any cluster group
    if missing:
        with open(OUT_MISSING, "w", encoding="utf-8") as f:
            for x in missing:
                f.write(f"{x}\n")

    print("Saved output table:", OUT_CSV)
    if missing:
        print(f"{len(missing)} seed genes did not appear in any group. Written to:", OUT_MISSING)


if __name__ == "__main__":
    main()
