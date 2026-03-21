import csv
import time
import requests

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

input_file = "data/degsensembl.csv"
output_file = "data/degsensemblfinal.csv"


def normalize_gene(g):
    if not g:
        return None
    return g.replace("\xa0", "").replace("\r", "").strip().upper()


def robust_gene2ensembl(genes, batch_size=100):
    mapping = {}

    # --- batch lookup ---
    for i in range(0, len(genes), batch_size):
        batch = genes[i:i + batch_size]

        url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens"
        payload = {"symbols": batch}

        try:
            r = requests.post(url, headers=HEADERS, json=payload, timeout=20)
            data = r.json() if r.ok else {}

            for g in batch:
                entry = data.get(g)
                mapping[g] = entry.get("id") if entry else None

        except Exception:
            for g in batch:
                mapping[g] = None

        time.sleep(0.2)

    # --- fallback: xrefs ---
    missing = [g for g, v in mapping.items() if v is None]
    print("Missing after batch:", len(missing))

    for g in missing:
        try:
            url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{g}"
            r = requests.get(url, headers=HEADERS, timeout=10)

            if r.ok:
                xrefs = r.json()
                gene_ids = [x["id"] for x in xrefs if x.get("type") == "gene"]
                mapping[g] = gene_ids[0] if gene_ids else None

        except Exception:
            pass

        time.sleep(0.05)

    # --- final fallback ---
    still_missing = [g for g, v in mapping.items() if v is None]
    print("Still missing after xrefs:", len(still_missing))

    for g in still_missing:
        try:
            url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{g}"
            r = requests.get(url, headers=HEADERS, timeout=10)

            if r.ok:
                data = r.json()
                mapping[g] = data.get("id")

        except Exception:
            pass

        time.sleep(0.05)

    return mapping


def fill_ensembl_ids(input_path, output_path):
    with open(input_path, newline="", encoding="utf-8-sig") as infile:

        sample = infile.read(2048)
        infile.seek(0)

        dialect = csv.Sniffer().sniff(sample)
        reader = csv.DictReader(infile, dialect=dialect)

        fieldnames = [f.strip() for f in reader.fieldnames]
        if "ENSEMBLID" not in fieldnames:
            fieldnames.append("ENSEMBLID")

        rows = []
        genes = []

        for row in reader:
            row = {k.strip(): v for k, v in row.items()}
            gene_raw = row.get("GENEID")

            gene = normalize_gene(gene_raw) if gene_raw else None

            rows.append(row)
            if gene:
                genes.append(gene)

        genes = list(set(genes))

        print("Total rows:", len(rows))
        print("Unique genes:", len(genes))

        mapping = robust_gene2ensembl(genes)

        unmapped = [g for g, v in mapping.items() if v is None]
        print("Final unmapped:", len(unmapped))
        print("Examples:", unmapped[:10])

        with open(output_path, "w", newline="", encoding="utf-8") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=dialect.delimiter)
            writer.writeheader()

            for row in rows:
                gene_raw = row.get("GENEID")
                gene = normalize_gene(gene_raw) if gene_raw else None

                row["ENSEMBLID"] = mapping.get(gene, "") if gene else ""

                clean_row = {k: row.get(k, "") for k in fieldnames}
                writer.writerow(clean_row)


if __name__ == "__main__":
    fill_ensembl_ids(input_file, output_file)