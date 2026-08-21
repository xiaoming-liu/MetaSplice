#!/usr/bin/env python3
"""
MetaSplice: Score splice-region variants for pathogenicity.

Input:  VCF or TSV/CSV file with chr, pos, ref, alt in hg19/GRCh37 coordinates
Output: TSV with splice pathogenicity scores and region classification

Three feature tiers:
  Minimal:  region + ESE/ESS + SMS (requires hg19 FASTA + GTF)
  Standard: + conservation + DANN + Eigen + gene-level + gnomAD + GPN-MSA
  Full:     + SpliceTransformer + Pangolin (requires GPU)
"""
import argparse
import csv
import gzip
import json
import os
import sys
import time
import bisect
import numpy as np
import pandas as pd
import pysam
import xgboost as xgb
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.join(SCRIPT_DIR, "model")

MODEL = xgb.XGBClassifier()
MODEL.load_model(os.path.join(MODEL_DIR, "metasplice.ubj"))
with open(os.path.join(MODEL_DIR, "features.json")) as f:
    FEAT_DEF = json.load(f)
ALL_FEATURES = FEAT_DEF["ALL_FEATURES"]
REGION_ONEHOT = FEAT_DEF["REGION_ONEHOT"]
ESE_FEATURES = FEAT_DEF["ESE_FEATURES"]
SMS_FEATURES = FEAT_DEF["SMS_FEATURES"]

with open(os.path.join(MODEL_DIR, "sms_7mer_scores.json")) as f:
    SMS_SCORES = json.load(f)

COMP = str.maketrans("ACGT", "TGCA")
BASES = {"A", "C", "G", "T"}

RESCUE_ESE = {
    "GAAGAA","AAGAAG","AGAAGA","GAAGAG","AAGAAA","GAAGAT","AAGAAC",
    "GAAGAC","AAGAAT","GAAGCA","AAGAGA","GAAGTA","AAGAGG","GAAGGA",
    "AAGAGT","GAAGGT","AAGAGC","GAAGGC","GAAGGG","GAAGTG","AAGATA",
    "GAAGTT","AAGATG","GAAGCT","AAGATC","GAAGCG","AAGACC","GAATGA",
    "AAGACG","AAGACT","AACAAG","AAGCAG","AATGAA","AAGCAA","AACGAA",
    "AAGGAG","AACAAC","AAGCAC","AAGGAA","AAGCGA","AAGGAC","AAGCGC",
    "AAGGTG","AAGGAT","AAGGCA","AAGGTC","AAGGCG","AAGGTA","AAGGCT",
    "AAGGGA","AAGGGC","AAGGGG","AAGGGT","AAGGTT","ACTGAA","ACAAAG",
    "ACAAGA","ACGAAG","ACAAGG","ACGAAC","ACTAAG","ACCAAG","AGAAAG",
    "AGAACG","AGAACC","AGAACT","AGAAGA","AGAAGC","AGAAGG","AGAAGT",
    "AGAATA","AGAATG","AGAATC","AGAATT","AGAACA","AGACAG","AGACGA",
    "AGACGC","AGAGAA","AGAGAC","AGAGAG","AGAGAT","AGAGCA","AGAGCG",
    "AGAGGA","AGAGGC","AGAGGG","AGAGGT","AGAGTG","AGATGA","AGATGC",
    "ATGAAG","CAAGAA","CAAGAC","CAAGAG","CAAGAT","CAAGCA","CAAGGA",
    "CAAGGC","CAAGGG","CAAGGT","CAAGTG","CAGAAG","CAGAAC","CAGAGA",
    "CAGAGC","CAGAGG","CTGAAG","CGAAGA","CGAAGC","CGAAGG","CGAAGT",
    "CGGAAG","CTGAAC","GAAACA","GAAACG","GAAAGA","GAAAGC","GAAAGG",
    "GAAAGT","GAAATA","GAAATC","GAAATG","GAAATT","GAACAA","GAACAC",
    "GAACAG","GAACAT","GAACCA","GAACCG","GAACGA","GAACGC","GAACGG",
    "GAACGT","GAACTA","GAACTG","GAACTT","GAAGCC","GAAGCG","GAAGCT",
    "GAAGGC","GAAGGG","GAAGGT","GAAGTA","GAAGTC","GAAGTG","GAAGTT",
    "GAATAA","GAATAC","GAATAG","GAATAT","GAATCA","GAATGA","GAATGC",
    "GAATGG","GAATGT","GAATTA","GAATTC","GAATTG","GAGAAA","GAGAAC",
    "GAGAAG","GAGAAT","GAGACA","GAGACC","GAGACG","GAGACT","GAGAGA",
    "GAGAGC","GAGAGG","GAGAGT","GAGATA","GAGATC","GAGATG","GAGATT",
    "GATGAA","GATGAC","GATGAG","GATGAT","GGATGA","GGAAGA","GGAAGC",
    "GGAAGG","GGAAGT","TGAAAG","TGAAGA","TGAAGC","TGAAGG","TGAAGT",
    "TTGAAG",
}

FAS_ESS = {
    "AAAAAA","AAAAAG","AAAAAT","AAAACA","AAAACC","AAAACG","AAAACT",
    "AAAAGA","AAAAGC","AAAAGG","AAAAGT","AAAATA","AAAATC","AAAATG",
    "AAAATT","AAACAA","AAACAC","AAACAG","AAACAT","AAACCA","AAACCC",
    "AAACCG","AAACCT","AAACGA","AAACGC","AAACTA","AAACTC","AAACTG",
    "AAACTT","AAAGTA","AAAGTT","AAATAA","AAATAC","AAATAG","AAATAT",
    "AAATCA","AAATCC","AAATCG","AAATGA","AAATGC","AAATGG","AAATGT",
    "AAATTA","AAATTC","AAATTG","AAATTT","AACAAA","AACAAT","AACCAA",
    "AACTAA","AAGAAA","AAGAAT","AAGTAA","AATAAA","AATAAC","AATAAT",
    "AATACA","AATACC","AATAGA","AATAGC","AATAGG","AATAGT","AATATA",
    "AATATC","AATATG","AATATT","AATCAA","AATGAA","AATGAT","AATTAA",
    "AATTAC","AATTAG","AATTAT","AATTCA","AATTGA","AATTTA","AATTTC",
    "AATTTG","AATTTT","ACAAAA","ACAAAT","ACATAA","ACTAAA","AGAAAA",
    "AGAAAT","AGATAA","AGTAAA","ATAAAA","ATAAAC","ATAAAG","ATAAAT",
    "ATAACA","ATAAGC","ATAAGG","ATAAGT","ATAATA","ATAATC","ATAATG",
    "ATAATT","ATACAA","ATAGAA","ATAGAT","ATATAT","ATCAAA","ATGAAA",
    "ATGAAT","ATGTAA","ATTAAA","ATTAAC","ATTAAT","ATTACA","ATTAGA",
    "ATTATA","ATTATC","ATTATG","ATTATT","ATTCAA","ATTGAA","ATTTAA",
    "ATTTAC","ATTTAG","ATTTAT","ATTTCA","ATTTGA","ATTTTA","ATTTTC",
    "ATTTTG","ATTTTT","CAAAAA","CAAAAT","CAATAA","GAAAAA","GAAAAT",
    "GAATAA","TAAAAA","TAAAAC","TAAAAG","TAAAAT","TAAACA","TAAAGA",
    "TAAATA","TAAATC","TAAATG","TAAATT","TAACAA","TAAGAA","TAATAA",
    "TAATAC","TAATAG","TAATAT","TAATCA","TAATGA","TAATTA","TACAAA",
    "TAGAAA","TATAAA","TATAAT","TATATA","TATATT","TATTAT","TATTTT",
    "TCAAAA","TGAAAA","TTAAAA","TTAAAC","TTAAAT","TTAATA","TTACAA",
    "TTAGAA","TTATAA","TTATAT","TTCAAA","TTGAAA","TTTAAA","TTTAAT",
    "TTTATA","TTTTAA","TTTTAT","TTTTTA","TTTTTT",
}

PRIORITY = {
    "Canonical Acceptor Site": 1, "Canonical Donor Site": 2,
    "Acceptor Region": 3, "Donor Region": 4,
    "PPT Region": 7, "Branch Point Region": 8, "Deep Intronic": 9,
}


def count_hexamers(seq, hexamer_set):
    return sum(1 for i in range(len(seq)-5) if seq[i:i+6] in hexamer_set)


# ============================================================
# Region classification
# ============================================================
def build_exon_index(gtf_path):
    """Parse GENCODE GTF and build boundary index."""
    exons_by_chr = defaultdict(list)
    gene_strand_map = {}
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split("\t")
            if fields[2] == "gene":
                gene = None
                for attr in fields[8].split(";"):
                    a = attr.strip()
                    if a.startswith("gene_name"):
                        gene = a.split('"')[1]; break
                if gene:
                    gene_strand_map[gene] = fields[6]
            if fields[2] != "exon": continue
            chrom = fields[0]
            if "_" in chrom: continue
            s, e, strand = int(fields[3]), int(fields[4]), fields[6]
            gene = None
            for attr in fields[8].split(";"):
                a = attr.strip()
                if a.startswith("gene_name"):
                    gene = a.split('"')[1]; break
            if gene:
                exons_by_chr[chrom].append((s, e, strand, gene))

    bnd_index = {}
    for chrom, exons in exons_by_chr.items():
        boundaries = []
        for s, e, strand, gene in exons:
            boundaries.append((s, s, e, strand, gene))
            boundaries.append((e, s, e, strand, gene))
        boundaries.sort()
        bnd_index[chrom] = (boundaries, [b[0] for b in boundaries])

    return bnd_index, gene_strand_map, exons_by_chr


def classify_variant(chrom, pos, bnd_index):
    if chrom not in bnd_index: return None, None
    boundaries, bnd_positions = bnd_index[chrom]
    idx = bisect.bisect_left(bnd_positions, pos)
    best_r, best_g, best_p = None, None, 999
    best_dist = float("inf")
    seen = set()
    for direction in (-1, 1):
        i = idx if direction == 1 else idx - 1
        while 0 <= i < len(boundaries):
            bp, s, e, strand, gene = boundaries[i]
            if abs(pos - bp) > best_dist and best_r is not None:
                break
            key = (s, e, strand, gene)
            if key not in seen:
                seen.add(key)
                if s <= pos <= e:
                    i += direction
                    continue
                if strand == "+":
                    if pos < s:
                        d = s - pos
                        if d <= 2: r = "Canonical Acceptor Site"
                        elif d <= 10: r = "Acceptor Region"
                        elif d <= 25: r = "PPT Region"
                        elif d <= 44: r = "Branch Point Region"
                        else: r = "Deep Intronic"
                    else:
                        d = pos - e
                        if d <= 2: r = "Canonical Donor Site"
                        elif d <= 8: r = "Donor Region"
                        else: r = "Deep Intronic"
                else:
                    if pos > e:
                        d = pos - e
                        if d <= 2: r = "Canonical Acceptor Site"
                        elif d <= 10: r = "Acceptor Region"
                        elif d <= 25: r = "PPT Region"
                        elif d <= 44: r = "Branch Point Region"
                        else: r = "Deep Intronic"
                    else:
                        d = s - pos
                        if d <= 2: r = "Canonical Donor Site"
                        elif d <= 8: r = "Donor Region"
                        else: r = "Deep Intronic"
                p = PRIORITY.get(r, 999)
                d_actual = min(abs(pos - s), abs(pos - e))
                if p < best_p or (p == best_p and d_actual < best_dist):
                    best_p, best_r, best_g = p, r, gene
                    best_dist = d_actual
            i += direction
    return best_r, best_g


# ============================================================
# Feature computation
# ============================================================
def compute_ese_sms(chrom, pos, ref, alt, strand, fasta):
    features = {}
    try:
        seq11 = fasta.fetch(chrom, pos-6, pos+5).upper()
    except: return {k: np.nan for k in ESE_FEATURES + SMS_FEATURES}
    if len(seq11) != 11:
        return {k: np.nan for k in ESE_FEATURES + SMS_FEATURES}
    mid5 = 5
    alt_seq11 = seq11[:mid5] + alt + seq11[mid5+1:]
    if strand == "-":
        seq11 = seq11.translate(COMP)[::-1]
        alt_seq11 = alt_seq11.translate(COMP)[::-1]
    re = count_hexamers(seq11, RESCUE_ESE); ae = count_hexamers(alt_seq11, RESCUE_ESE)
    rs = count_hexamers(seq11, FAS_ESS); as_ = count_hexamers(alt_seq11, FAS_ESS)
    total = re + rs
    features["delta_ESE"] = ae - re
    features["delta_ESS"] = as_ - rs
    features["ref_ESE_count"] = re
    features["ref_ESS_count"] = rs
    features["delta_ESE_ESS"] = (ae-re) - (as_-rs)
    features["ese_ess_ratio"] = round(re/total, 4) if total > 0 else 0.5

    try:
        seq13 = fasta.fetch(chrom, pos-7, pos+6).upper()
    except:
        features.update({k: np.nan for k in SMS_FEATURES}); return features
    if len(seq13) != 13:
        features.update({k: np.nan for k in SMS_FEATURES}); return features
    mid6 = 6
    alt_seq13 = seq13[:mid6] + alt + seq13[mid6+1:]
    if strand == "-":
        seq13 = seq13.translate(COMP)[::-1]
        alt_seq13 = alt_seq13.translate(COMP)[::-1]
    ref_s, alt_s = [], []
    for j in range(7):
        rk, ak = seq13[j:j+7], alt_seq13[j:j+7]
        if rk in SMS_SCORES: ref_s.append(SMS_SCORES[rk])
        if ak in SMS_SCORES: alt_s.append(SMS_SCORES[ak])
    if ref_s and alt_s and len(ref_s)==len(alt_s):
        deltas = [a-r for r,a in zip(ref_s, alt_s)]
        features["sms_abs_delta_mean"] = np.mean(np.abs(deltas))
        features["sms_abs_delta_max"] = max(np.abs(deltas))
        features["sms_ref_mean"] = np.mean(ref_s)
    else:
        features.update({k: np.nan for k in SMS_FEATURES})
    return features


def compute_gerp(variants, gerp_dir):
    """Compute GERP++ RS scores."""
    if not gerp_dir: return {}
    scores = {}
    by_chr = defaultdict(list)
    for i, (chrom, pos, _, _) in enumerate(variants):
        if not chrom: continue  # unmapped after reverse liftover
        by_chr[chrom.replace("chr","")].append((i, pos))

    for chrom_num, vlist in sorted(by_chr.items(), key=lambda x: (len(x[0]),x[0])):
        path = os.path.join(gerp_dir, f"chr{chrom_num}.maf.rates.gz")
        if not os.path.exists(path): continue
        positions = {pos for _, pos in vlist}
        min_pos, max_pos = min(positions), max(positions)
        gerp = {}
        with gzip.open(path, "rt") as f:
            ln = 0
            for line in f:
                ln += 1
                if ln < min_pos-10: continue
                if ln > max_pos+10: break
                if ln in positions:
                    parts = line.strip().split("\t")
                    if len(parts)>=2:
                        try: gerp[ln] = float(parts[1])
                        except: pass
        for idx, pos in vlist:
            if pos in gerp: scores[idx] = gerp[pos]
    return scores


def compute_dann_eigen(variants, dann_dir, eigen_dir):
    """Compute DANN and Eigen scores via sorted scan."""
    dann_scores, eigen_scores = {}, {}
    by_chr = defaultdict(list)
    for i, (chrom, pos, ref, alt) in enumerate(variants):
        if not chrom: continue  # unmapped after reverse liftover
        by_chr[chrom.replace("chr","")].append((i, pos, ref, alt))

    # Eigen files store coding and noncoding records as two separately sorted
    # blocks; intronic variants appear only in the noncoding block. Scanning must
    # not stop at the first position past max_pos or those records are never read.
    # DANN is a single sorted block, so its early break is kept.
    for db_name, db_dir, result_dict, score_col, multi_block in [
        ("DANN", dann_dir, dann_scores, 4, False),
        ("Eigen", eigen_dir, eigen_scores, 6, True),
    ]:
        if not db_dir: continue
        for chrom_num in sorted(by_chr.keys(), key=lambda x: (len(x),x)):
            vlist = by_chr[chrom_num]
            candidates = [f for f in os.listdir(db_dir)
                          if f"chr{chrom_num}." in f and f.endswith(".gz")]
            if not candidates: continue
            path = os.path.join(db_dir, candidates[0])
            lookup = {}
            for idx, pos, ref, alt in vlist:
                lookup[(pos, ref, alt)] = idx
                lookup[(pos, ref.translate(COMP), alt.translate(COMP))] = idx
            positions = {pos for _, pos, _, _ in vlist}
            min_pos, max_pos = min(positions), max(positions)
            remaining = set(idx for idx,_,_,_ in vlist)
            with gzip.open(path, "rt") as f:
                for line in f:
                    if line.startswith("#"): continue
                    parts = line.split("\t")
                    if len(parts) <= score_col: continue
                    try: pos = int(parts[1])
                    except: continue
                    if pos < min_pos: continue
                    if pos > max_pos:
                        if multi_block: continue   # a later block may restart lower
                        break
                    key = (pos, parts[2], parts[3])
                    if key in lookup:
                        idx = lookup[key]
                        if idx in remaining:
                            try:
                                val = parts[score_col].rstrip()
                                if val != ".": result_dict[idx] = float(val)
                            except: pass
                            remaining.discard(idx)
                            if not remaining: break
    return dann_scores, eigen_scores


def compute_conservation(variants, hg38_map, phylop_dir, phastcons_dir):
    """Compute phyloP and phastCons from wigFix files."""
    results = {k: {} for k in ["phyloP100way_vertebrate","phyloP470way_mammalian",
                                "phastCons100way_vertebrate","phastCons470way_mammalian"]}
    hg38_by_chr = defaultdict(set)
    for i, (hc, hp) in hg38_map.items():
        hg38_by_chr[f"chr{hc}"].add(hp)

    for feat_name, feat_dir, pattern in [
        ("phyloP100way_vertebrate", phylop_dir, "phyloP100way"),
        ("phyloP470way_mammalian", phylop_dir, "phyloP470way"),
        ("phastCons100way_vertebrate", phastcons_dir, "phastCons100way"),
        ("phastCons470way_mammalian", phastcons_dir, "phastCons470way"),
    ]:
        if not feat_dir: continue
        for chrom in hg38_by_chr:
            positions = hg38_by_chr[chrom]
            path = os.path.join(feat_dir, f"{chrom}.{pattern}.wigFix.gz")
            if not os.path.exists(path): continue
            min_pos, max_pos = min(positions), max(positions)
            scores = {}
            current_pos = 0
            with gzip.open(path, "rt") as f:
                for line in f:
                    if line.startswith("fixedStep"):
                        for part in line.strip().split():
                            if part.startswith("start="):
                                current_pos = int(part.split("=")[1])
                        continue
                    if current_pos > max_pos+1000: break
                    if current_pos in positions:
                        try: scores[current_pos] = float(line.strip())
                        except: pass
                    current_pos += 1
            for i, (hc, hp) in hg38_map.items():
                if f"chr{hc}" == chrom and hp in scores:
                    results[feat_name][i] = scores[hp]
    return results


def compute_gnomad(variants, hg38_map, gnomad_dir):
    """Compute gnomAD AF and POPMAX_AF."""
    af_scores, popmax_scores = {}, {}
    if not gnomad_dir: return af_scores, popmax_scores
    hg38_by_chr = defaultdict(list)
    for i, (hc, hp) in hg38_map.items():
        ref, alt = variants[i][2], variants[i][3]
        hg38_by_chr[hc].append((i, hp, ref, alt))

    for chrom_num, vlist in hg38_by_chr.items():
        path = os.path.join(gnomad_dir, f"gnomad.joint.v4.1.sites.chr{chrom_num}.vcf.bgz.snp.tsv.gz")
        if not os.path.exists(path): continue
        lookup = {}
        for idx, pos, ref, alt in vlist:
            lookup[(pos, ref, alt)] = idx
            lookup[(pos, ref.translate(COMP), alt.translate(COMP))] = idx
        positions = {pos for _,pos,_,_ in vlist}
        min_pos, max_pos = min(positions), max(positions)
        with gzip.open(path, "rt") as f:
            header = f.readline().strip().split("\t")
            af_idx = header.index("gnomAD_joint_AF")
            pm_idx = header.index("gnomAD_joint_POPMAX_AF")
            for line in f:
                parts = line.split("\t")
                try: pos = int(parts[1])
                except: continue
                if pos < min_pos: continue
                if pos > max_pos: break
                key = (pos, parts[2], parts[3])
                if key in lookup:
                    idx = lookup[key]
                    try: af_scores[idx] = float(parts[af_idx])
                    except: pass
                    try: popmax_scores[idx] = float(parts[pm_idx])
                    except: pass
                    del lookup[key]
                    if not lookup: break
    return af_scores, popmax_scores


def compute_gpn_msa(variants, hg38_map, gpn_msa_file):
    """Compute GPN-MSA scores via tabix."""
    scores = {}
    if not gpn_msa_file or not os.path.exists(gpn_msa_file): return scores
    tbx = pysam.TabixFile(gpn_msa_file)
    for i, (hc, hp) in hg38_map.items():
        ref, alt = variants[i][2], variants[i][3]
        rc, ac = ref.translate(COMP), alt.translate(COMP)
        try:
            for row in tbx.fetch(str(hc), hp-1, hp):
                parts = row.split("\t")
                if len(parts)>=5:
                    r, a = parts[2], parts[3]
                    if (r==ref and a==alt) or (r==rc and a==ac):
                        scores[i] = float(parts[4]); break
        except: continue
    tbx.close()
    return scores


def compute_gene_level(variants, gene_assignments, dbnsfp_path):
    """Compute gene-level features from dbNSFP gene file."""
    results = {}
    if not dbnsfp_path: return results
    gene_data = {}
    col_map = {"gnomAD_pLI":"gnomAD_pLI", "gnomAD_LOEUF":"gnomAD_LOEUF",
               "RVIS_EVS":"RVIS_EVS", "RVIS_percentile_EVS":"RVIS_percentile_EVS",
               "ClinGen_Haploinsufficiency_Score":"ClinGen_Haploinsufficiency_Score"}
    with gzip.open(dbnsfp_path, "rt") as f:
        header = f.readline().strip().split("\t")
        idx_map = {k: header.index(v) for k,v in col_map.items() if v in header}
        gname_idx = header.index("Gene_name")
        for line in f:
            fs = line.strip().split("\t")
            gn = fs[gname_idx]
            row = {}
            for k, ci in idx_map.items():
                v = fs[ci]
                try: row[k] = float(v) if v not in (".","") else np.nan
                except: row[k] = np.nan
            gene_data[gn] = row

    for i, gene in gene_assignments.items():
        if gene in gene_data:
            results[i] = gene_data[gene]
    return results


def compute_exon_context(variants, gene_assignments, exons_by_chr):
    """Compute exon context features."""
    results = {}
    flat = {}
    for c, exons in exons_by_chr.items():
        se = sorted(set((s,e) for s,e,_,_ in exons))
        flat[c] = se

    starts_by_chr = {c: [s for s,e in se] for c,se in flat.items()}

    for i, (chrom, pos, _, _) in enumerate(variants):
        exons = flat.get(chrom)
        starts = starts_by_chr.get(chrom)
        if not exons: continue
        idx = bisect.bisect_right(starts, pos)
        min_dist, nearest_len = float("inf"), np.nan
        for j in range(max(0,idx-50), min(len(exons),idx+50)):
            s, e = exons[j]
            if s <= pos <= e:
                min_dist, nearest_len = 0, e-s+1; break
            d = min(abs(pos-s), abs(pos-e))
            if d < min_dist: min_dist, nearest_len = d, e-s+1
        if min_dist < 10000:
            results[i] = {"dist_to_boundary": float(min_dist),
                          "exon_length": float(nearest_len),
                          "exon_length_mod3": float(nearest_len % 3)}
    return results


def load_liftover(chrom_num, liftover_dir):
    """Load hg19->hg38 position mapping for one chromosome."""
    path = os.path.join(liftover_dir, f"chr{chrom_num}.hg19tohg38.gz")
    if not os.path.exists(path): return {}
    mapping = {}
    with gzip.open(path, "rt") as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                try:
                    mapping[int(parts[1])] = (parts[2], int(parts[3]))
                except: continue
    return mapping


def load_reverse_liftover(chrom_num, rev_liftover_dir, target_positions):
    """Load hg38->hg19 position mapping for one chromosome, streaming and
    retaining only positions in target_positions (memory-bounded)."""
    path = os.path.join(rev_liftover_dir, f"chr{chrom_num}.hg38tohg19.gz")
    if not os.path.exists(path): return {}
    mapping = {}
    with gzip.open(path, "rt") as f:
        f.readline()  # header
        for line in f:
            parts = line.split("\t")
            if len(parts) < 4: continue
            try:
                hg38_pos = int(parts[1])
            except ValueError:
                continue
            if hg38_pos in target_positions:
                try:
                    mapping[hg38_pos] = (parts[2].strip(), int(parts[3]))
                except ValueError:
                    continue
    return mapping


# ============================================================
# Input parsing
# ============================================================
def parse_input(path):
    variants = []
    if path.endswith(".vcf"):
        with open(path) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                chrom = parts[0]
                if not chrom.startswith("chr"): chrom = f"chr{chrom}"
                pos = int(parts[1])
                ref, alt = parts[3], parts[4]
                if ref in BASES and alt in BASES:
                    variants.append((chrom, pos, ref, alt))
    else:
        sep = "\t" if path.endswith(".tsv") else ","
        with open(path) as f:
            reader = csv.DictReader(f, delimiter=sep)
            for row in reader:
                chrom = row.get("chr", row.get("CHROM", row.get("#CHROM", "")))
                if not chrom.startswith("chr"): chrom = f"chr{chrom}"
                pos = int(row.get("pos", row.get("POS", 0)))
                ref = row.get("ref", row.get("REF", ""))
                alt = row.get("alt", row.get("ALT", ""))
                if ref in BASES and alt in BASES:
                    variants.append((chrom, pos, ref, alt))
    return variants


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="MetaSplice: splice variant pathogenicity scorer")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    # --genome controls which coordinate system the input VCF is in.
    # hg19 path: hg19 FASTA + GENCODE v19; conservation/gnomAD/GPN-MSA looked up
    #           via hg19->hg38 liftover (--liftover-dir).
    # hg38 path: hg38 FASTA + GENCODE v49; conservation/gnomAD/GPN-MSA looked up
    #           directly. DANN/Eigen/GERP/sequence-context features looked up
    #           via hg38->hg19 reverse liftover (--reverse-liftover-dir + hg19
    #           FASTA/GTF for sequence/exon context).
    parser.add_argument("--genome", default="hg19", choices=["hg19","hg38"])
    parser.add_argument("--fasta", default=os.environ.get("METASPLICE_FASTA",""),
                        help="Reference FASTA matching --genome (hg19 or hg38). Used for "
                             "ESE/ESS hexamer + SMS 7-mer sequence context. Reference "
                             "nucleotides at any given biological position are identical "
                             "between hg19 and hg38, so either build is acceptable.")
    parser.add_argument("--gtf", default=os.environ.get("METASPLICE_GTF",""),
                        help="GENCODE GTF matching --genome. For --genome hg19 use v19 "
                             "(gencode.v19.annotation.gtf.gz); for --genome hg38 use v49 "
                             "(gencode.v49.annotation.gtf.gz).")
    parser.add_argument("--gerp-dir", default=os.environ.get("METASPLICE_GERP",""))
    parser.add_argument("--dann-dir", default=os.environ.get("METASPLICE_DANN",""))
    parser.add_argument("--eigen-dir", default=os.environ.get("METASPLICE_EIGEN",""))
    parser.add_argument("--phylop-dir", default=os.environ.get("METASPLICE_PHYLOP",""))
    parser.add_argument("--phastcons-dir", default=os.environ.get("METASPLICE_PHASTCONS",""))
    parser.add_argument("--gnomad-dir", default=os.environ.get("METASPLICE_GNOMAD",""))
    parser.add_argument("--gpn-msa", default=os.environ.get("METASPLICE_GPNMSA",""))
    parser.add_argument("--dbnsfp-gene", default=os.environ.get("METASPLICE_DBNSFP",""))
    parser.add_argument("--liftover-dir", default=os.environ.get("METASPLICE_LIFTOVER",""),
                        help="hg19->hg38 liftover files (used when --genome hg19)")
    parser.add_argument("--reverse-liftover-dir",
                        default=os.environ.get("METASPLICE_RLIFTOVER",""),
                        help="hg38->hg19 liftover files (used when --genome hg38). "
                             "Pre-compute with reproducibility/scripts/invert_liftover.py.")
    parser.add_argument("--spt-scores", default=os.environ.get("METASPLICE_SPT",""),
                        help="SpliceTransformer score TSV (columns: #CHROM, POS, REF, ALT, "
                             "score, plus 15 tissue columns). Variant order must match --input.")
    parser.add_argument("--pangolin-scores", default=os.environ.get("METASPLICE_PANGOLIN",""),
                        help="Pangolin score CSV (columns: CHROM, POS, REF, ALT, Pangolin). "
                             "Variant order must match --input.")
    parser.add_argument("--include-out-of-scope", action="store_true",
                        help="Also output variants that fall outside the 7 model regions "
                             "(exonic + unclassified). By default only in-scope variants "
                             "are written; the model is not trained on out-of-scope variants "
                             "and predictions for them are not meaningful.")
    args = parser.parse_args()

    TISSUES = ["Adipose Tissue","Blood","Blood Vessel","Brain","Colon","Heart",
               "Kidney","Liver","Lung","Muscle","Nerve","Small Intestine",
               "Skin","Spleen","Stomach"]

    variants = parse_input(args.input)
    n = len(variants)
    print(f"Loaded {n:,} SNVs")

    # Build exon index
    bnd_index, gene_strand_map, exons_by_chr = {}, {}, {}
    if args.gtf:
        print("Parsing GENCODE...")
        bnd_index, gene_strand_map, exons_by_chr = build_exon_index(args.gtf)

    # Classify regions
    print("Classifying regions...")
    regions, gene_assignments = {}, {}
    for i, (chrom, pos, _, _) in enumerate(variants):
        r, g = classify_variant(chrom, pos, bnd_index)
        if r: regions[i] = r
        if g: gene_assignments[i] = g

    in_scope = len(regions)
    print(f"  In scope (non-exonic splice regions): {in_scope:,}/{n:,}")

    # Initialize features
    feat_matrix = [{f: np.nan for f in ALL_FEATURES} for _ in range(n)]

    # Region one-hot
    for i in range(n):
        for r in REGION_ONEHOT:
            feat_matrix[i][f"r_{r}"] = 1 if regions.get(i) == r else 0

    # ESE/ESS + SMS
    fasta = pysam.FastaFile(args.fasta) if args.fasta else None
    if fasta:
        print("Computing ESE/ESS + SMS features...")
        for i, (chrom, pos, ref, alt) in enumerate(variants):
            gene = gene_assignments.get(i, "")
            strand = gene_strand_map.get(gene, "+")
            ese_sms = compute_ese_sms(chrom, pos, ref, alt, strand, fasta)
            feat_matrix[i].update(ese_sms)

    # Exon context
    if exons_by_chr:
        print("Computing exon context...")
        exon_ctx = compute_exon_context(variants, gene_assignments, exons_by_chr)
        for i, ctx in exon_ctx.items():
            feat_matrix[i].update(ctx)

    # ── Build coordinate maps based on input genome ──
    # hg19_variants: variants in hg19 coords (for GERP/DANN/Eigen lookups).
    # hg38_map: {idx: (hg38_chrom, hg38_pos)} for phyloP/phastCons/gnomAD/GPN-MSA.
    if args.genome == "hg19":
        hg19_variants = variants
        # Forward liftover for hg38-bridged features
        hg38_map = {}
        if args.liftover_dir and (args.phylop_dir or args.phastcons_dir
                                  or args.gnomad_dir or args.gpn_msa):
            print("Loading hg19->hg38 liftover...")
            by_chr = defaultdict(list)
            for i, (chrom, pos, _, _) in enumerate(variants):
                by_chr[chrom.replace("chr","")].append((i, pos))
            for chrom_num, vlist in by_chr.items():
                mapping = load_liftover(chrom_num, args.liftover_dir)
                for idx, pos in vlist:
                    if pos in mapping:
                        hg38_map[idx] = mapping[pos]
            print(f"  Mapped: {len(hg38_map):,}/{n:,}")
    else:  # args.genome == "hg38"
        # Input is hg38; conservation/gnomAD/GPN-MSA hit directly.
        hg38_map = {i: (chrom.replace("chr",""), pos)
                    for i, (chrom, pos, _, _) in enumerate(variants)}
        # Reverse liftover for hg19-direct features
        # Unmapped entries keep empty chrom -> skipped in compute_gerp/dann_eigen
        hg19_variants = [("", 0, ref, alt) for _, _, ref, alt in variants]
        if args.reverse_liftover_dir and (args.gerp_dir or args.dann_dir or args.eigen_dir):
            print("Loading hg38->hg19 reverse liftover...")
            by_chr = defaultdict(list)
            for i, (chrom, pos, _, _) in enumerate(variants):
                by_chr[chrom.replace("chr","")].append((i, pos))
            mapped_count = 0
            for chrom_num, vlist in by_chr.items():
                target = {pos for _, pos in vlist}
                mapping = load_reverse_liftover(chrom_num, args.reverse_liftover_dir, target)
                for idx, pos in vlist:
                    if pos in mapping:
                        hg19_chrom_raw, hg19_pos = mapping[pos]
                        hg19_chrom = (hg19_chrom_raw if hg19_chrom_raw.startswith("chr")
                                       else f"chr{hg19_chrom_raw}")
                        _, _, ref, alt = variants[idx]
                        hg19_variants[idx] = (hg19_chrom, hg19_pos, ref, alt)
                        mapped_count += 1
            print(f"  Mapped to hg19: {mapped_count:,}/{n:,}")
        elif (args.gerp_dir or args.dann_dir or args.eigen_dir):
            print("WARNING: --genome hg38 with hg19-direct feature dirs but no "
                  "--reverse-liftover-dir; DANN/Eigen/GERP will be skipped")

    # GERP++ (always queried against hg19 database via hg19_variants)
    if args.gerp_dir:
        print("Computing GERP++...")
        gerp = compute_gerp(hg19_variants, args.gerp_dir)
        for i, v in gerp.items(): feat_matrix[i]["GERP++_RS"] = v

    # DANN + Eigen (hg19 databases)
    if args.dann_dir or args.eigen_dir:
        print("Computing DANN/Eigen...")
        dann, eigen = compute_dann_eigen(hg19_variants, args.dann_dir, args.eigen_dir)
        for i, v in dann.items(): feat_matrix[i]["DANN_score"] = v
        for i, v in eigen.items(): feat_matrix[i]["Eigen-phred_coding"] = v

    # Conservation
    if args.phylop_dir or args.phastcons_dir:
        print("Computing conservation...")
        cons = compute_conservation(variants, hg38_map, args.phylop_dir, args.phastcons_dir)
        for feat_name, scores in cons.items():
            for i, v in scores.items(): feat_matrix[i][feat_name] = v

    # gnomAD
    if args.gnomad_dir:
        print("Computing gnomAD AF...")
        af, popmax = compute_gnomad(variants, hg38_map, args.gnomad_dir)
        for i, v in af.items(): feat_matrix[i]["gnomAD_AF"] = v
        for i, v in popmax.items(): feat_matrix[i]["gnomAD_POPMAX_AF"] = v

    # GPN-MSA
    if args.gpn_msa:
        print("Computing GPN-MSA...")
        gpn = compute_gpn_msa(variants, hg38_map, args.gpn_msa)
        for i, v in gpn.items(): feat_matrix[i]["gpn_msa_score"] = v

    # Gene-level
    if args.dbnsfp_gene:
        print("Computing gene-level features...")
        gene_feats = compute_gene_level(variants, gene_assignments, args.dbnsfp_gene)
        for i, feats in gene_feats.items():
            feat_matrix[i].update(feats)

    if fasta: fasta.close()

    # SpliceTransformer (optional, GPU-scored externally)
    if args.spt_scores:
        print(f"Loading SpliceTransformer scores: {args.spt_scores}")
        sep = "\t"
        with open(args.spt_scores) as fh:
            first = fh.readline()
            if "," in first and first.count(",") > first.count("\t"):
                sep = ","
        spt = pd.read_csv(args.spt_scores, sep=sep)
        if len(spt) != n:
            raise ValueError(
                f"SpT row count ({len(spt):,}) does not match input ({n:,}). "
                f"Score SpT on the SAME input VCF in the SAME order."
            )
        # Alignment check
        chrom_col = "#CHROM" if "#CHROM" in spt.columns else "CHROM"
        for i, (chrom, pos, ref, alt) in enumerate(variants):
            ec = str(spt[chrom_col].iloc[i]).lstrip("chr")
            bc = chrom.lstrip("chr")
            if ec != bc or int(spt["POS"].iloc[i]) != pos or \
               str(spt["REF"].iloc[i]) != ref or str(spt["ALT"].iloc[i]) != alt:
                raise ValueError(
                    f"SpT alignment mismatch at row {i}: "
                    f"input={chrom}:{pos}:{ref}:{alt}, "
                    f"SpT={spt[chrom_col].iloc[i]}:{spt['POS'].iloc[i]}:"
                    f"{spt['REF'].iloc[i]}:{spt['ALT'].iloc[i]}"
                )
        for i in range(n):
            feat_matrix[i]["spt_orig"] = float(spt["score"].iloc[i])
            for t in TISSUES:
                if t in spt.columns:
                    feat_matrix[i][f"spt_o_{t}"] = float(spt[t].iloc[i])

    # Pangolin (optional, GPU-scored externally)
    if args.pangolin_scores:
        print(f"Loading Pangolin scores: {args.pangolin_scores}")
        pg_inc = np.zeros(n)
        pg_dec = np.zeros(n)
        with open(args.pangolin_scores) as fh:
            fh.readline()  # header
            for i, line in enumerate(fh):
                if i >= n:
                    raise ValueError(
                        f"Pangolin has more rows than input ({n}). "
                        f"Score Pangolin on the SAME input in the SAME order."
                    )
                parts = line.strip().split(",", 4)
                # Alignment check on chr/pos/ref/alt
                if len(parts) >= 4:
                    ec = str(parts[0]).lstrip("chr")
                    bc = variants[i][0].lstrip("chr")
                    try:
                        ep = int(parts[1])
                    except ValueError:
                        ep = None
                    if ec != bc or ep != variants[i][1] or \
                       parts[2] != variants[i][2] or parts[3] != variants[i][3]:
                        raise ValueError(
                            f"Pangolin alignment mismatch at row {i}: "
                            f"input={variants[i]}, pangolin row={parts[:4]}"
                        )
                inc, dec = 0.0, 0.0
                if len(parts) >= 5:
                    for token in parts[4].replace(",", "|").split("|"):
                        if ":" not in token: continue
                        try:
                            s = float(token.rsplit(":", 1)[-1])
                            if s > inc: inc = s
                            if s < dec: dec = s
                        except ValueError:
                            pass
                pg_inc[i] = inc
                pg_dec[i] = dec
        for i in range(n):
            feat_matrix[i]["pangolin_inc"] = float(pg_inc[i])
            feat_matrix[i]["pangolin_dec"] = float(pg_dec[i])
            feat_matrix[i]["pangolin_max"] = float(max(pg_inc[i], abs(pg_dec[i])))

    # Predict
    print("Predicting...")
    X = np.array([[feat_matrix[i].get(f, np.nan) for f in ALL_FEATURES]
                   for i in range(n)], dtype=np.float32)
    preds = MODEL.predict_proba(X)[:, 1]

    # Feature coverage
    available = sum(1 for f in ALL_FEATURES
                    if sum(1 for i in range(n) if not np.isnan(feat_matrix[i].get(f, np.nan))) > n*0.5)
    print(f"Features with >50% coverage: {available}/{len(ALL_FEATURES)}")

    # Output
    rows = []
    for i, (chrom, pos, ref, alt) in enumerate(variants):
        row = {
            "chr": chrom, "pos": pos, "ref": ref, "alt": alt,
            "region": regions.get(i, "NA"),
            "gene": gene_assignments.get(i, "NA"),
            "splice_score": round(float(preds[i]), 6),
        }
        if args.include_out_of_scope:
            row["in_scope"] = i in regions
        rows.append(row)

    df = pd.DataFrame(rows)
    excluded = n - len(regions)
    if not args.include_out_of_scope:
        # Filter to in-scope variants (the 7 model regions); the model is not
        # trained on exonic / unclassified variants and predictions on them
        # are not meaningful. Pass --include-out-of-scope to keep them.
        in_scope_mask = [i in regions for i in range(n)]
        df = df[in_scope_mask].reset_index(drop=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved: {args.output}")
    in_scope_n = len(regions)
    if args.include_out_of_scope:
        print(f"  Reported: {len(df):,} variants ({in_scope_n:,} in-scope + "
              f"{excluded:,} out-of-scope with in_scope=False)")
        scope_mask = df["in_scope"].astype(bool)
    else:
        print(f"  Reported: {len(df):,} in-scope variants")
        if excluded:
            print(f"  Excluded: {excluded:,} out-of-scope (exonic / unclassified) variants")
        scope_mask = pd.Series([True] * len(df))
    score_col = df.loc[scope_mask, "splice_score"]
    print(f"  Score >= 0.1: {(score_col >= 0.1).sum():,}")
    print(f"  Score >= 0.2: {(score_col >= 0.2).sum():,}")
    print(f"  Score >= 0.5: {(score_col >= 0.5).sum():,}")


if __name__ == "__main__":
    main()
