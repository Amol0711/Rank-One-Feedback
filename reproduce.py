#!/usr/bin/env python3
"""Execute fixed-model numerical workflows in a fresh external output directory."""
from __future__ import annotations
import argparse,contextlib,hashlib,json,os,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent

def hashes(root):
    return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.rglob('*')) if p.is_file() and p.relative_to(root).parts[0]!='.git'}

def destination(path):
    raw=Path(path).absolute()
    if any(p.is_symlink() for p in (raw,*raw.parents)):raise ValueError('Output path contains a symlink')
    out=raw.resolve()
    if out.exists() or out==ROOT or ROOT in out.parents or out in ROOT.parents:
        raise ValueError('Use a new external directory that does not overlap the input package')
    return out

def install_guard(out):
    """Restrict local data access to this package and its external output tree."""
    dependency_roots={Path(sys.prefix).resolve(),Path(sys.base_prefix).resolve()}
    accesses=set()
    def category(path):
        for base,name in [(ROOT,'PACKAGE'),(out,'OUTPUT')]:
            if path==base or base in path.parents:return name+'/'+str(path.relative_to(base))
        for base in dependency_roots:
            if path==base or base in path.parents:return 'DEPENDENCY/'+str(path.relative_to(base))
        if path in (Path('/dev/null'),Path('/dev/urandom'),Path('/etc/ld.so.cache')):return 'SYSTEM/'+path.name
        raise PermissionError('File access outside numerical inputs, outputs and Python dependencies')
    def audit(event,args):
        if event.startswith('socket.') or event in ('subprocess.Popen','os.system','os.posix_spawn'):
            raise PermissionError('Network and subprocess execution are disabled')
        if event=='open' and isinstance(args[0],(str,bytes,os.PathLike)):
            p=Path(os.fsdecode(args[0])).resolve();label=category(p)
            mode=args[1];flags=args[2]
            write=(isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)))
            if write and not label.startswith('OUTPUT/'):
                raise PermissionError('Writing outside the output tree is disabled')
            accesses.add(('write' if write else 'read',label))
        if event in ('os.mkdir','os.remove','os.rmdir'):
            p=Path(os.fsdecode(args[0])).resolve()
            if not (p==out or out in p.parents):raise PermissionError('Mutation outside the output tree is disabled')
        if event in ('os.rename','os.symlink','os.link'):
            raise PermissionError('Renaming and linking are disabled during numerical execution')
    sys.addaudithook(audit)
    return accesses

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['check','export','responses','stationary','directions','comparison','finite-gain','intervals','integer-replay','all'],default='check')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=destination(args.output)
    before=hashes(ROOT);out.mkdir(parents=True);accesses=install_guard(out)
    from verify_release import verify
    from numerics.workflows import WORKFLOWS,check_reference,save
    verify();record={'mode':args.mode,'status':'running','workflows':{},'source_unchanged':False}
    save(out/'execution.json',record)
    try:
        t=time.perf_counter();record['reference_check']=check_reference()
        chosen=list(WORKFLOWS) if args.mode=='all' else ([] if args.mode=='check' else [args.mode])
        for name in chosen:
            dest=out/name;dest.mkdir();start=time.perf_counter()
            with (dest/'execution.log').open('w') as log,contextlib.redirect_stdout(log):
                result=WORKFLOWS[name](dest)
            record['workflows'][name]={'result':result,'seconds':time.perf_counter()-start}
            save(out/'execution.json',record);print(name+' PASS',flush=True)
        record.update(status='passed',elapsed_seconds=time.perf_counter()-t,source_unchanged=before==hashes(ROOT))
        if not record['source_unchanged']:raise ValueError('Input package changed')
    except Exception as exc:
        record.update(status='failed',error_type=type(exc).__name__,error=str(exc),source_unchanged=before==hashes(ROOT));raise
    finally:
        save(out/'execution.json',record)
        save(out/'file_access.json',{'network_and_subprocess_disabled':True,'accesses':[{'mode':m,'path':p} for m,p in sorted(accesses)],'path_tokens':'PACKAGE, OUTPUT, DEPENDENCY, SYSTEM'})
    print('All requested workflows passed.')

if __name__=='__main__':
    try:main()
    except (ValueError,ArithmeticError,OSError,KeyError,TypeError) as exc:
        print(type(exc).__name__+': '+str(exc),file=sys.stderr);sys.exit(1)
