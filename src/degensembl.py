import csv
import time
import requests

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

input_file = "data/degsensembl.csv"
output_file = "data/degsensemblfinal.csv"


def batch_gene2ensembl(genes, batch_size=100):
    """
    Maps gene symbols to Ensembl IDs using batch POST.
    Returns dict: {gene_symbol: ensembl_id or None}
    """
    mapping = {}

    for i in range(0, len(genes), batch_size):
        batch = genes[i:i + batch_size]

        url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens"
        payload = {"symbols": batch}

        try:
            r = requests.post(url, headers=HEADERS, json=payload, timeout=20)
            if not r.ok:
                for g in batch:
                    mapping[g] = None
                continue

            data = r.json()

            for g in batch:
                entry = data.get(g)
                mapping[g] = entry.get("id") if entry else None

        except Exception:
            for g in batch:
                mapping[g] = None

        time.sleep(0.2)  # mild rate control

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
            gene = row.get("GENEID")

            rows.append(row)
            if gene:
                genes.append(gene.strip())

        print("Total rows:", len(rows))
        print("Unique genes:", len(set(genes)))

        mapping = batch_gene2ensembl(list(set(genes)))

        with open(output_path, "w", newline="", encoding="utf-8") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=dialect.delimiter)
            writer.writeheader()

            for row in rows:
                gene = row.get("GENEID")
                row["ENSEMBLID"] = mapping.get(gene, "") if gene else ""

                clean_row = {k: row.get(k, "") for k in fieldnames}
                writer.writerow(clean_row)


if __name__ == "__main__":
    fill_ensembl_ids(input_file, output_file)