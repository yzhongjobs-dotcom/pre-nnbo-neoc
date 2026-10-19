"""Verify the project's selected local papers; optionally restore missing files.

Downloading or extracting text does not change the recorded reading scope.
No research experiment, model training, or literature ranking is performed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def local_path(value):
    path = (ROOT / value).resolve()
    if not path.is_relative_to((ROOT / '文献').resolve()):
        raise ValueError('Material path outside 文献: ' + value)
    return path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true',
                        help='Restore missing selected PDFs, checking their recorded SHA256')
    parser.add_argument('--ids', nargs='+', help='Limit verification to material IDs')
    args = parser.parse_args()
    manifest = json.loads((ROOT / '文献/本地资料清单.json').read_text())
    records = manifest['records']
    if args.ids:
        missing = set(args.ids) - {r['id'] for r in records}
        if missing:
            raise ValueError('Unknown IDs: ' + ', '.join(sorted(missing)))
        records = [r for r in records if r['id'] in args.ids]
    verified = set()
    for record in records:
        pdf = local_path(record['local_path'])
        text = local_path(record['text_path'])
        if not pdf.exists() and args.download:
            url = record.get('recovery_url') or record.get('source_url')
            if not url:
                raise ValueError('No restoration URL: ' + record['id'])
            data = urllib.request.urlopen(url, timeout=45).read()
            if not data.startswith(b'%PDF') or hashlib.sha256(data).hexdigest() != record['sha256']:
                raise ValueError('Downloaded PDF identity/hash mismatch: ' + record['id'])
            pdf.parent.mkdir(parents=True, exist_ok=True)
            pdf.write_bytes(data)
        if sha256(pdf) != record['sha256']:
            raise ValueError('Local PDF hash mismatch: ' + record['id'])
        if not text.exists() and args.download:
            text.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(['pdftotext', '-layout', str(pdf), str(text)], check=True)
        if not text.is_file() or not text.read_text(errors='replace').strip():
            raise ValueError('Missing/empty text: ' + record['id'])
        if record.get('text_sha256') and sha256(text) != record['text_sha256']:
            raise ValueError('Text hash mismatch (pdftotext versions can differ): ' + record['id'])
        verified.add(record['local_path'])
    print(json.dumps({'material_records': len(records), 'unique_pdfs': len(verified),
                      'verified': True, 'new_scientific_reading': False}, ensure_ascii=False))


if __name__ == '__main__':
    main()
