# Brassica carinata Genome Annotation Pipeline

## Overview

**Goal:** Enrich your Brassica carinata GFF with Arabidopsis functional info (ortholog, description, GO IDs, GO names), then query fatty-acid/oil genes near GWAS SNPs.

---

### **Inputs**

- `genomic.gff` — B. carinata GFF
- `bcar.proteins.faa` — B. carinata protein FASTA (headers like `>KAL... hypothetical protein Bca###_###### [Brassica carinata]`)
- `ath.proteins.faa` — Arabidopsis thaliana protein FASTA

### **Outputs**

- `ath_functions_with_names.clean.tsv` — Arabidopsis functional table (AGI ID, description, GO IDs, GO names, etc.)
- `rbh.bcar_at.tsv` — Reciprocal best hits (B. carinata protein ↔ Arabidopsis AGI)
- `bcar_annotated.gff3` — Enriched GFF (Ortholog, Description, GO, GO_names, etc.)

### **Required Tools**

- `diamond` (≥ 2.0)
- `python3`
- `Rscript` with Bioconductor
- *(optional)* `bedtools`, `awk`

---

## Step 1: Create a Working Directory

```bash
mkdir -p proj_reannot && cd proj_reannot

# Copy your inputs here:
cp /path/to/genomic.gff .
cp /path/to/bcar.proteins.faa .
cp /path/to/ath.proteins.faa .
```
**Expected:** 3 files now in `proj_reannot`.

---

## Step 2: Build Arabidopsis Functional Table

### 2.1 `build_ath_functions.R`

Creates a table of Arabidopsis gene functions using Bioconductor.

<details>
<summary>Click to expand R script</summary>

```r
#!/usr/bin/env Rscript
# build_ath_functions.R

pkgs <- c("BiocManager","org.At.tair.db","AnnotationDbi","dplyr","readr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p == "BiocManager") install.packages("BiocManager", repos="https://cloud.r-project.org")
    else BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}
suppressPackageStartupMessages({
  library(org.At.tair.db); library(AnnotationDbi); library(dplyr); library(readr)
})
avail_cols <- columns(org.At.tair.db)
uni_col <- if ("UNIPROT"   %in% avail_cols) "UNIPROT" else
           if ("SWISSPROT" %in% avail_cols) "SWISSPROT" else
           if ("ACCNUM"    %in% avail_cols) "ACCNUM" else NA

tair_ids <- keys(org.At.tair.db, keytype = "TAIR")
sel_cols <- c("SYMBOL","GENENAME","GO","ONTOLOGY")
if (!is.na(uni_col)) sel_cols <- c(sel_cols, uni_col)

raw <- AnnotationDbi::select(org.At.tair.db, keys=tair_ids, columns=sel_cols, keytype="TAIR")

summ <- raw %>%
  group_by(TAIR) %>%
  summarise(
    description = dplyr::first(na.omit(GENENAME)),
    go_terms    = paste(sort(unique(na.omit(GO))), collapse="|"),
    uniprot     = if (!is.na(uni_col)) paste(sort(unique(na.omit(.data[[uni_col]]))), collapse="|") else "",
    .groups="drop"
  ) %>%
  mutate(
    description = ifelse(is.na(description), "", description),
    go_terms    = ifelse(go_terms=="", "", go_terms),
    uniprot     = ifelse(uniprot=="", "", uniprot)
  ) %>%
  filter(grepl("^AT[1-5MC]G[0-9]{5}(\\.[0-9]+)?$", TAIR)) %>%
  transmute(
    ath_gene = TAIR,
    description, go_terms,
    kegg = "", mapman = "", uniprot
  )

write_tsv(summ, "ath_functions.tsv")
cat("Wrote ath_functions.tsv with", nrow(summ), "rows\n")
```
</details>

**Run:**

```bash
Rscript build_ath_functions.R
head ath_functions.tsv
```

**Expected:**  
File `ath_functions.tsv` with header:  
`ath_gene description go_terms kegg mapman uniprot`

---

### 2.2 Add GO Names (human-readable)

#### `add_go_names.R`

<details>
<summary>Click to expand R script</summary>

```r
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GO.db)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
  library(stringr)
})
ath <- read_tsv("ath_functions.tsv", show_col_types = FALSE)

all_go <- unique(unlist(strsplit(paste(ath$go_terms, collapse="|"), "\\|")))
all_go <- all_go[!is.na(all_go) & all_go!=""]
go_map <- AnnotationDbi::select(GO.db, keys=all_go, columns=c("TERM"), keytype="GOID") %>%
  filter(!is.na(GOID)) %>% distinct(GOID, TERM)

lookup <- setNames(go_map$TERM, go_map$GOID)

ath$go_names <- vapply(ath$go_terms, function(x){
  if (is.na(x) || x=="") return("")
  ids <- unlist(strsplit(x, "\\|"))
  terms <- lookup[ids]
  paste(na.omit(terms), collapse="|")
}, FUN.VALUE = character(1))

write_tsv(ath, "ath_functions_with_names.tsv")
cat("Wrote ath_functions_with_names.tsv\n")
```
</details>

**Run & Clean:**

```bash
Rscript add_go_names.R
awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {for(i=1;i<=NF;i++) if($i=="NA") $i=""; print}' \
  ath_functions_with_names.tsv > ath_functions_with_names.clean.tsv
head ath_functions_with_names.clean.tsv
```

---

## Step 3: Compute Reciprocal Best Hits (RBH) with DIAMOND

#### `diamond_rbh.sh`

<details>
<summary>Click to expand Bash script</summary>

```bash
#!/usr/bin/env bash
set -euo pipefail
ATH=${1:-ath.proteins.faa}
BCAR=${2:-bcar.proteins.faa}
OUT=${3:-rbh}

command -v diamond >/dev/null || { echo "diamond not found"; exit 1; }
[[ -s "$ATH" ]] || { echo "Missing $ATH"; exit 1; }
[[ -s "$BCAR" ]] || { echo "Missing $BCAR"; exit 1; }

echo "==> Building DBs"
diamond makedb --in "$ATH" -d "${OUT}_ath"
diamond makedb --in "$BCAR" -d "${OUT}_bcar"

echo "==> Forward (Bcar -> Ath)"
diamond blastp -d "${OUT}_ath.dmnd" -q "$BCAR" -o "${OUT}.b2a.tsv" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --max-target-seqs 1 --evalue 1e-5 --query-cover 50 --subject-cover 50

echo "==> Reverse (Ath -> Bcar)"
diamond blastp -d "${OUT}_bcar.dmnd" -q "$ATH" -o "${OUT}.a2b.tsv" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --max-target-seqs 1 --evalue 1e-5 --query-cover 50 --subject-cover 50

echo "==> RBH"
awk 'NR==FNR {rev[$1]=$2; next} ($2 in rev) && rev[$2]==$1 {print $1"\t"$2}' \
  "${OUT}.a2b.tsv" "${OUT}.b2a.tsv" > "${OUT}.pairs_AT_to_KAL.tsv"

# Flip to KAL first (what our annotator prefers): KAL \t AT
awk 'BEGIN{FS=OFS="\t"} {print $2,$1}' "${OUT}.pairs_AT_to_KAL.tsv" > "${OUT}.bcar_at.tsv"

wc -l "${OUT}.b2a.tsv" "${OUT}.a2b.tsv" "${OUT}.pairs_AT_to_KAL.tsv" "${OUT}.bcar_at.tsv" | sed 's/^/   /'
echo "RBH written to ${OUT}.bcar_at.tsv (KAL first, AT second)"
```
</details>

**Run:**

```bash
chmod +x diamond_rbh.sh
./diamond_rbh.sh ath.proteins.faa bcar.proteins.faa rbh
head rbh.bcar_at.tsv
```

---

## Step 4: Annotate the B. carinata GFF

#### `annotate_gff_with_arabidopsis.py`

<details>
<summary>Click to expand Python script</summary>

```python
#!/usr/bin/env python3
import sys, re, urllib.parse
from collections import defaultdict

if len(sys.argv) < 6:
    sys.stderr.write(
        "Usage:\n  annotate_gff_with_arabidopsis.py bcar.gff3 rbh_pairs.tsv ath_functions.tsv bcar.proteins.faa out.gff3\n"
    ); sys.exit(1)

gff_in, rbh_tsv, ath_fn_tsv, prot_faa, gff_out = sys.argv[1:6]

def parse_attrs(attr):
    d={}
    for kv in attr.split(";"):
        if kv and "=" in kv:
            k,v = kv.split("=",1)
            d[k]=urllib.parse.unquote(v)
    return d

def enc(v): return urllib.parse.quote(v, safe="()|,:._- /")

def norm_prot(raw):
    if raw is None: return ""
    s = raw.strip()
    if not s: return ""
    s = s.split()[0]
    toks = s.split("|")
    cand = max(toks, key=lambda t: len(re.findall(r'[A-Za-z0-9._-]', t))) if len(toks)>1 else s
    m = re.search(r'[A-Za-z]{2,}\d{4,}\.\d+', cand)  # e.g., KAL0744661.1
    if m: return m.group(0)
    m = re.search(r'[A-Za-z0-9._-]+', cand)
    return m.group(0) if m else cand

def is_at(x): return re.match(r'^AT[1-5MC]G[0-9]+(\.[0-9]+)?$', x or "") is not None
def at_base(at): return at.split(".")[0] if at and "." in at else at

# 1) GFF: locus_tag -> gene_id
lines = []
gene_ids=set(); locus_to_gene={}
with open(gff_in) as fh:
    for line in fh:
        lines.append(line)
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9: continue
        if f[2].lower() in ("gene","pseudogene","mrna","transcript"):
            A=parse_attrs(f[8])
            gid = A.get("ID") or A.get("Name")
            if f[2].lower() in ("gene","pseudogene") and gid: gene_ids.add(gid)
            lt = A.get("locus_tag")
            par= A.get("Parent")
            if lt and gid: locus_to_gene.setdefault(lt, gid)
            if lt and par: locus_to_gene.setdefault(lt, par)

# 2) FASTA: protein -> gene via locus_tag in header
lt_pat = re.compile(r'(Bca\d{3,5}_\d{5,7})')
prot_to_gene={}
with open(prot_faa) as fh:
    for ln in fh:
        if not ln.startswith(">"): continue
        hdr = ln[1:].strip()
        pid = norm_prot(hdr.split()[0])       # e.g., KAL0744661.1
        m = lt_pat.search(hdr)
        if m:
            lt = m.group(1)                   # e.g., Bca101_101253
            gid = locus_to_gene.get(lt)
            if gid:
                prot_to_gene[pid]=gid
                pv=re.sub(r'\.\d+$','',pid)
                prot_to_gene.setdefault(pv,gid)

# 3) RBH (auto-detect columns)
prot_to_at={}
with open(rbh_tsv) as fh:
    for ln in fh:
        ln=ln.strip()
        if not ln: continue
        a,b = ln.split("\t")[:2]
        if is_at(a) and not is_at(b):
            at_gene, bcar_prot = a, norm_prot(b)
        elif is_at(b) and not is_at(a):
            at_gene, bcar_prot = b, norm_prot(a)
        else:
            at_gene, bcar_prot = b, norm_prot(a)  # fallback
        prot_to_at[bcar_prot]=at_gene
        if "." in bcar_prot:
            prot_to_at.setdefault(bcar_prot.split(".")[0], at_gene)

# 4) Arabidopsis function table (optionally has go_names)
ath_fn={}
with open(ath_fn_tsv) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    hidx = {h:i for i,h in enumerate(header)}
    req = ["ath_gene","description","go_terms","kegg","mapman","uniprot"]
    for r in req:
        if r not in hidx:
            sys.stderr.write(f"ERROR: column '{r}' missing in {ath_fn_tsv}\n"); sys.exit(2)
    has_go_names = ("go_names" in hidx)
    for ln in fh:
        parts = ln.rstrip("\n").split("\t")
        if len(parts) < len(header): parts += [""]*(len(header)-len(parts))
        agi = parts[hidx["ath_gene"]]
        if not agi: continue
        base = at_base(agi)
        rec = {
            "description": parts[hidx["description"]],
            "go_terms":    parts[hidx["go_terms"]],
            "kegg":        parts[hidx["kegg"]],
            "mapman":      parts[hidx["mapman"]],
            "uniprot":     parts[hidx["uniprot"]],
        }
        if has_go_names: rec["go_names"] = parts[hidx["go_names"]]
        ath_fn[base]=rec
        if agi!=base: ath_fn[agi]=rec

# 5) Build annotations per gene
from collections import defaultdict
gene_annot = defaultdict(lambda: {"Ortholog":"","Description":"","GO":"","KEGG":"","MapMan":"","UniProt":"","GO_names":""})
hit_prot=0
for prot, at in prot_to_at.items():
    gid = prot_to_gene.get(prot) or (prot_to_gene.get(prot.split(".")[0]) if "." in prot else None)
    if not gid: continue
    hit_prot+=1
    info = ath_fn.get(at) or ath_fn.get(at_base(at)) or {}
    gene_annot[gid]["Ortholog"]=at
    gene_annot[gid]["Description"]=info.get("description","")
    gene_annot[gid]["GO"]=info.get("go_terms","")
    gene_annot[gid]["KEGG"]=info.get("kegg","")
    gene_annot[gid]["MapMan"]=info.get("mapman","")
    gene_annot[gid]["UniProt"]=info.get("uniprot","")
    if "go_names" in info and info["go_names"]:
        gene_annot[gid]["GO_names"]=info["go_names"]

# 6) Write annotated GFF
ann_ct=0
with open(gff_out,"w") as out:
    out.write("##gff-version 3\n# + Arabidopsis RBH-based annotations (GO_names supported)\n")
    for line in lines:
        if not line or line.startswith("#"): out.write(line); continue
        f=line.rstrip("\n").split("\t")
        if len(f)<9: out.write(line); continue
        if f[2].lower() in ("gene","pseudogene"):
            A=parse_attrs(f[8]); gid=A.get("ID") or A.get("Name")
            ann=gene_annot.get(gid)
            if ann and (ann["Ortholog"] or ann["Description"] or ann["GO"] or ann["GO_names"]):
                add=[]
                if ann["Ortholog"]:   add.append(f"Ortholog={enc(ann['Ortholog'])}")
                if ann["Description"]:add.append(f"Description={enc(ann['Description'])}")
                if ann["GO"]:        add.append(f"GO={enc(ann['GO'])}")
                if ann["GO_names"]:  add.append(f"GO_names={enc(ann['GO_names'])}")
                if ann["KEGG"]:      add.append(f"KEGG={enc(ann['KEGG'])}")
                if ann["MapMan"]:    add.append(f"MapMan={enc(ann['MapMan'])}")
                if ann["UniProt"]:   add.append(f"UniProt={enc(ann['UniProt'])}")
                add.append("Source=RBH_Arabidopsis")
                sep = "" if f[8].endswith(";") or f[8]=="" else ";"
                f[8] = f[8] + sep + ";".join(add)
                out.write("\t".join(f) + "\n"); ann_ct+=1
            else:
                out.write(line)
        else:
            out.write(line)

sys.stderr.write(f"[INFO] GFF genes: {len(gene_ids)} | locus->gene: {len(locus_to_gene)} | prot->gene: {len(prot_to_gene)} | RBH prots: {len(prot_to_at)} | matched prots: {hit_prot} | genes annotated: {ann_ct}\n")
```
</details>

**Run:**

```bash
python3 annotate_gff_with_arabidopsis.py \
  genomic.gff rbh.bcar_at.tsv ath_functions_with_names.clean.tsv bcar.proteins.faa \
  bcar_annotated.gff3
```

---

## Step 5: Query Genes Near a SNP (With Fatty-Acid Keywords)

#### `view_region_genes.py`

<details>
<summary>Click to expand Python script</summary>

```python
#!/usr/bin/env python3
import sys, re, urllib.parse

def parse_attrs(s):
    d={}
    for kv in s.split(";"):
        if kv and "=" in kv:
            k,v=kv.split("=",1)
            d[k]=urllib.parse.unquote(v)
    return d

if len(sys.argv)<6:
    print("Usage: view_region_genes.py <gff> <chr> <start> <end> <out.tsv> [keyword_regex]", file=sys.stderr); sys.exit(1)

gff, chrom, start, end, out = sys.argv[1:6]
start, end = int(start), int(end)
kw_pat = re.compile(sys.argv[6], re.I) if len(sys.argv)>6 else None

rows=[]
with open(gff) as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<9 or f[0]!=chrom or f[2]!="gene": continue
        s,e = int(f[3]), int(f[4])
        if e<start or s>end: continue
        A = parse_attrs(f[8])
        text = " ".join([A.get("Description",""), A.get("GO",""), A.get("GO_names",""), A.get("Ortholog","")]).strip()
        if kw_pat and not kw_pat.search(text): continue
        rows.append([f[0],s,e,A.get("ID",""),A.get("locus_tag",""),A.get("Ortholog",""),
                     A.get("Description",""),A.get("GO",""),A.get("GO_names","")])

with open(out,"w") as fo:
    fo.write("chr\tstart\tend\tgene_id\tlocus_tag\tortholog\tdescription\tGO\tGO_names\n")
    for r in rows: fo.write("\t".join(map(str,r))+"\n")

print(f"Wrote {len(rows)} rows to {out}")
```
</details>

**Example:** Query ±100 kb window around a SNP at `CM081020.1:31630794`, fatty/lipid keywords

```bash
CHR="CM081020.1"; POS=31630794; WIN=100000
python3 view_region_genes.py \
  bcar_annotated.gff3 "$CHR" $((POS-WIN)) $((POS+WIN)) region_fatty.tsv \
  '(fatty|lipid|acyl|desaturase|elong|KCS|FAE1|FAD2|FAD3|DGAT|LPAAT|GPAT|OLEOSIN|ACC|KAS|FATA|FATB|triacylglycerol|acyltransferase)'

# Or export everything in the window
python3 view_region_genes.py bcar_annotated.gff3 "$CHR" $((POS-WIN)) $((POS+WIN)) region_all.tsv
```

---

## Step 6: Troubleshooting Cheatsheet

- **Annotated 0 gene features:**  
  Usually RBH file columns reversed or protein IDs not matching. Use `rbh.bcar_at.tsv` with KAL first, AT second.
- **No GO names appear:**  
  Use `ath_functions_with_names.clean.tsv` and latest annotator script.
- **Weird `%2B/%3B` characters:**  
  These are URL encodings required by GFF. `view_region_genes.py` decodes them.
- **Windows line endings:**  
  If you edited scripts on Windows, run `dos2unix *.sh *.py *.R`.

---

## Step 7: One-Shot Minimal Pipeline

```bash
# From empty folder with 3 inputs copied in:
# genomic.gff, bcar.proteins.faa, ath.proteins.faa

# A) Arabidopsis functions
Rscript build_ath_functions.R
Rscript add_go_names.R
awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {for(i=1;i<=NF;i++) if($i=="NA") $i=""; print}' \
  ath_functions_with_names.tsv > ath_functions_with_names.clean.tsv

# B) RBH
chmod +x diamond_rbh.sh
./diamond_rbh.sh ath.proteins.faa bcar.proteins.faa rbh

# C) Annotate GFF
python3 annotate_gff_with_arabidopsis.py \
  genomic.gff rbh.bcar_at.tsv ath_functions_with_names.clean.tsv bcar.proteins.faa \
  bcar_annotated.gff3

# D) Check
grep -c 'Source=RBH_Arabidopsis' bcar_annotated.gff3
grep -c 'GO_names=' bcar_annotated.gff3
grep -c 'Description=' bcar_annotated.gff3
```

---

## Contributions

Feel free to open issues or pull requests to improve this workflow or add new features.
