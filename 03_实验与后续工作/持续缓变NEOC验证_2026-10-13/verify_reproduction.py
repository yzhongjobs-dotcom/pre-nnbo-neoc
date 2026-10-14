from pathlib import Path
import hashlib
import json
import subprocess
import sys

base=Path(__file__).resolve().parent
names=['results.json','summary.csv','trace.csv']
def hashes():
    return {name:hashlib.sha256((base/name).read_bytes()).hexdigest() for name in names}
first=hashes()
with (base/'run_log_reproduction.txt').open('w',encoding='utf-8') as log:
    subprocess.run([sys.executable,'-X','utf8',str(base/'run_continuous_check.py')],
                   stdout=log,stderr=subprocess.STDOUT,check=True)
second=hashes()
assert first==second
record=dict(status='TWO_RUNS_IDENTICAL',compared_files=names,sha256=second,
            limits='No physical closed loop, runtime advantage or floating-point interval certificate.')
(base/'reproducibility.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps(record,indent=2))
