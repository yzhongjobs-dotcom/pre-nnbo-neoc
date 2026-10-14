"""Verify retained originals and generated navigation after organization."""
from pathlib import Path
import csv
import hashlib
import json
import re
from urllib.parse import unquote

root = Path(__file__).resolve().parents[1]
with (root / '90_整理记录/文件迁移清单.csv').open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f))
for row in rows:
    p = root / row['新路径']
    assert p.is_file(), row
    assert hashlib.sha256(p.read_bytes()).hexdigest() == row['SHA256'], row
    assert p.stat().st_size == int(row['字节数']), row
    assert not (root / row['原路径']).exists(), row
names = {'README.md', '00_当前进度与阅读顺序.md', '00_研究现状阅读说明.md',
         '00_实验与后续工作清单.md', '文件索引.md'}
checked = 0
for p in root.rglob('*.md'):
    if p.name not in names:
        continue
    for target in re.findall(r'\]\(([^)]+)\)', p.read_text(encoding='utf-8')):
        if target.startswith(('https://', 'http://', 'mailto:')):
            continue  # External citations are not local filesystem paths.
        assert (p.parent / unquote(target)).exists(), (str(p), target)
        checked += 1
index = (root / '90_整理记录/文件索引.md').read_text(encoding='utf-8')
for row in rows:
    assert row['原路径'] in index and row['新路径'] in index
print(json.dumps({'original_files_verified': len(rows),
                  'navigation_links_verified': checked,
                  'index_rows_verified': len(rows)}, ensure_ascii=True, indent=2))
