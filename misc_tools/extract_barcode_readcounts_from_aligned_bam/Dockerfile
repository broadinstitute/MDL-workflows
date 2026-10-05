FROM python:3.10-slim

RUN pip install --no-cache-dir pysam

COPY process_barcodes.py /usr/local/bin/process_barcodes.py
RUN chmod +x /usr/local/bin/process_barcodes.py
