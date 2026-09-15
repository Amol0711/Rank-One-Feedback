#!/usr/bin/env python3
"""Check the exact numerical-package inventory and content hashes."""
from __future__ import annotations
import hashlib,json
from pathlib import Path,PurePosixPath
ROOT=Path(__file__).resolve().parent

def verify(root=ROOT):
    root=Path(root)
    expected=json.loads((root/'MANIFEST.json').read_text())['files']
    actual={}
    for p in sorted(root.rglob('*')):
        if p.relative_to(root).parts[0]=='.git':continue
        if p.is_symlink():raise ValueError('Symlink in numerical package')
        if p.is_file() and p.relative_to(root).as_posix()!='MANIFEST.json':
            name=p.relative_to(root).as_posix()
            actual[name]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    if actual!=expected:raise ValueError('Numerical inventory or hash mismatch')
    for name in expected:
        q=PurePosixPath(name)
        if q.is_absolute() or '..' in q.parts:raise ValueError('Unsafe manifest path')
    return {'passed':True,'files':len(actual)+1,'method':'exact inventory and SHA-256 comparison'}

if __name__=='__main__':print(json.dumps(verify(),indent=2))
