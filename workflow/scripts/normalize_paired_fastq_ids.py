#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def open_text(path, mode):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return path.open(mode)


def read_record(handle, path):
    record = [handle.readline() for _ in range(4)]
    if not record[0]:
        return None
    if any(line == "" for line in record):
        raise ValueError(f"Incomplete FASTQ record in {path}")
    if not record[0].startswith("@"):
        raise ValueError(f"Invalid FASTQ header in {path}: {record[0].rstrip()}")
    return record


def normalize_common_id(header):
    read_id = header[1:].strip().split()[0]
    if read_id.endswith("/1") or read_id.endswith("/2"):
        read_id = read_id[:-2]
    return read_id


def normalize_pair_ids(r1_in, r2_in, r1_out, r2_out):
    pair_count = 0
    with open_text(r1_in, "r") as r1_reader, open_text(r2_in, "r") as r2_reader, \
            open_text(r1_out, "w") as r1_writer, open_text(r2_out, "w") as r2_writer:
        while True:
            r1_record = read_record(r1_reader, r1_in)
            r2_record = read_record(r2_reader, r2_in)

            if r1_record is None and r2_record is None:
                break
            if r1_record is None or r2_record is None:
                raise ValueError("R1 and R2 FASTQ files contain different numbers of records")

            common_id = normalize_common_id(r1_record[0])
            r1_record[0] = f"@{common_id} 1/2\n"
            r2_record[0] = f"@{common_id} 2/2\n"
            r1_writer.writelines(r1_record)
            r2_writer.writelines(r2_record)
            pair_count += 1

    print(f"[normalize_paired_fastq_ids] normalized {pair_count} read pairs")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Normalize paired FASTQ mate headers to a shared read ID."
    )
    parser.add_argument("--r1-in", required=True)
    parser.add_argument("--r2-in", required=True)
    parser.add_argument("--r1-out", required=True)
    parser.add_argument("--r2-out", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    normalize_pair_ids(args.r1_in, args.r2_in, args.r1_out, args.r2_out)


if __name__ == "__main__":
    main()
