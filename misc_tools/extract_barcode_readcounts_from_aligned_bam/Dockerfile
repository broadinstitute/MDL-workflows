FROM python:3.10-slim

RUN pip install --no-cache-dir pysam numpy

COPY process_barcodes.py merge_barcode_counts.py /usr/local/bin/
RUN chmod +x /usr/local/bin/process_barcodes.py /usr/local/bin/merge_barcode_counts.py
