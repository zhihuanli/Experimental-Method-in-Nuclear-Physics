from pathlib import Path
import subprocess,sys
from independent_pdf import portable
ROOT=Path(__file__).resolve().parent/'independent'
for k in map(int,sys.argv[1:] or ['2','3','4']):
    subprocess.run(['/opt/homebrew/bin/pandoc',f'L{k}.md','--from=markdown+raw_tex',f'--metadata-file=L{k}.yaml','--standalone','--number-sections','--top-level-division=chapter','-o',f'L{k}.tex'],cwd=ROOT,check=True)
    for iteration in range(3):
        proc=subprocess.run(['/Library/TeX/texbin/xelatex','-interaction=nonstopmode','-halt-on-error',f'L{k}.tex'],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        (ROOT/f'L{k}-compile-{iteration}.txt').write_text(proc.stdout)
        if proc.returncode:
            print(proc.stdout[-3500:]);raise SystemExit(proc.returncode)
    portable(ROOT/f'L{k}.pdf',ROOT/f'L{k}-portable.pdf')
