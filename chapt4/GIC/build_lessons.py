"""Build the GIC lessons, runnable notebooks, and their in-place ROOT outputs.

Run with the course PyROOT environment, from this directory:
    python3 build_lessons.py
Markdown sources contain paired python/cpp fences; outputs are generated, not typed.
"""
from pathlib import Path
import base64
import contextlib
import html
import io
import re
import subprocess
import tempfile

import nbformat
from nbconvert.filters.markdown import markdown2html
import ROOT

BASE = Path(__file__).resolve().parent
ROOT.gROOT.SetBatch(True)

STYLE = """
body{margin:0;background:#fff;color:#222;font:17px/1.7 system-ui,sans-serif}
main{max-width:980px;margin:auto;padding:30px 24px 70px}
h1{font-size:1.85rem}h2{font-size:1.4rem;border-bottom:1px solid #ddd;padding-top:22px}
h3{font-size:1.12rem}a{color:#145787}pre{background:#f6f7f8;padding:14px;overflow:auto;font:14px/1.55 monospace}
code{font-size:.9em}table{border-collapse:collapse;margin:18px 0;width:100%}td,th{border-bottom:1px solid #ddd;padding:8px;text-align:left}
figure{margin:18px auto 26px;max-width:780px}figure img{width:100%;height:auto}figcaption{font-size:15px;color:#555;text-align:center}
img{max-width:100%;height:auto}nav,.language-choice{padding:10px 0}button{font:inherit;padding:5px 14px;margin-right:8px;cursor:pointer}
button[aria-pressed=true]{background:#164a70;color:white;border-color:#164a70}.cpp{display:none}
body[data-language=cpp] .python{display:none}body[data-language=cpp] .cpp{display:block}
.output{background:#fff;border-left:3px solid #ccc}.downloads{font-size:15px}blockquote{margin-left:0;padding-left:16px;border-left:3px solid #aaa}
"""
SCRIPT = """
function setLanguage(lang){document.body.dataset.language=lang;
document.querySelectorAll('[data-lang]').forEach(b=>b.setAttribute('aria-pressed',b.dataset.lang===lang));
try{localStorage.setItem('gic-language',lang)}catch(e){}}
try{setLanguage(localStorage.getItem('gic-language')||'python')}catch(e){setLanguage('python')}
"""

def render_page(stem, title, body):
    page = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>{title}</title>
<style>{STYLE}</style><script>window.MathJax={{tex:{{inlineMath:[['$','$'],['\\\\(','\\\\)']]}}}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script></head>
<body data-language="python"><main><nav><a href="../../reference.html">学习资料</a> ·
<a href="GIC_simulation.html">GIC 信号形成</a> · <a href="Experimental_Data_Analysis.html">实验波形分析</a></nav>
{body}</main><script>{SCRIPT}</script></body></html>'''
    (BASE / f'{stem}.html').write_text(page)

def build(stem):
    (BASE/'figures').mkdir(exist_ok=True)
    source = (BASE / '_sources' / f'{stem}.md').read_text()
    # Each pair optionally names the figure that the displayed code generates.
    pattern = re.compile(r'```python\n(.*?)\n```\n\n```cpp\n(.*?)\n```(?:\n\n<!-- figure: ([\w-]+) -->)?', re.S)
    pieces, cells_py, cells_cpp, cpp_blocks = [], [], [], []
    namespace = {'__name__': '__lesson__'}
    start, execution = 0, 0
    for match in pattern.finditer(source):
        prose = source[start:match.start()].strip()
        if prose:
            pieces.append(markdown2html(prose))
            cells_py.append(nbformat.v4.new_markdown_cell(prose))
            cells_cpp.append(nbformat.v4.new_markdown_cell(prose))
        py, cpp, fig = match.groups()
        execution += 1
        output = io.StringIO()
        with tempfile.TemporaryDirectory(prefix='gic-root-output-') as temporary:
            root_log = Path(temporary)/'output.txt'
            ROOT.gSystem.RedirectOutput(str(root_log),'w')
            try:
                with contextlib.redirect_stdout(output):
                    exec(compile(py, f'{stem} cell {execution}', 'exec'), namespace)
            finally:
                ROOT.gInterpreter.ProcessLine('gSystem->RedirectOutput(nullptr);')
            output.write(root_log.read_text())
        outputs = []
        pieces.append(f'<div class="python"><pre><code>{html.escape(py)}</code></pre></div>')
        pieces.append(f'<div class="cpp"><pre><code>{html.escape(cpp)}</code></pre></div>')
        if output.getvalue():
            outputs.append(nbformat.v4.new_output('stream', name='stdout', text=output.getvalue()))
            pieces.append(f'<pre class="output">{html.escape(output.getvalue())}</pre>')
        if fig:
            namespace['c'].SaveAs(str(BASE / 'figures' / f'{fig}.png'))
            data = base64.b64encode((BASE / 'figures' / f'{fig}.png').read_bytes()).decode()
            outputs.append(nbformat.v4.new_output('display_data', data={'image/png':data}))
            pieces.append(f'<figure><img src="figures/{fig}.png" alt="{fig.replace("-", " ")}"></figure>')
        cells_py.append(nbformat.v4.new_code_cell(py, execution_count=execution, outputs=outputs))
        # C++ is run separately below; shared figure is a rendering of equivalent code.
        cells_cpp.append(nbformat.v4.new_code_cell(cpp))
        cpp_blocks.append(cpp)
        if fig:
            cpp_blocks.append(f'c->SaveAs("figures/{fig}-cpp.png");')
        start = match.end()
    prose = source[start:].strip()
    pieces.append(markdown2html(prose))
    cells_py.append(nbformat.v4.new_markdown_cell(prose))
    cells_cpp.append(nbformat.v4.new_markdown_cell(prose))
    title = source.splitlines()[0].removeprefix('# ')
    choice = f'''<div class="language-choice">代码语言：
<button data-lang="python" onclick="setLanguage('python')" aria-pressed="true">PyROOT</button>
<button data-lang="cpp" onclick="setLanguage('cpp')" aria-pressed="false">ROOT C++</button></div>'''
    for language, cells, kernel in [('PyROOT',cells_py,'python3'),('CPP',cells_cpp,'root')]:
        nb = nbformat.v4.new_notebook(cells=cells, metadata={'kernelspec':{'name':kernel,'display_name':language,'language':'python' if language=='PyROOT' else 'c++'}})
        nbformat.write(nb, BASE / f'{stem}_{language}.ipynb')
    includes = '#include <TFile.h>\n#include <TTree.h>\n#include <TGraph.h>\n#include <TH2D.h>\n#include <TCanvas.h>\n#include <TLegend.h>\n#include <TStyle.h>\n#include <TLine.h>\n#include <TBox.h>\n#include <TMath.h>\n#include <vector>\n#include <algorithm>\n#include <fstream>\n#include <iostream>\n#include <iomanip>\n#include <cmath>\n#include <numeric>\n'
    macro = includes + '#include <TSystem.h>\n' + f'void {stem}() {{\n' + 'gSystem->mkdir("figures",true);\n' + '\n\n'.join(cpp_blocks) + '\n}\n'
    (BASE / f'{stem}.C').write_text(macro)
    result = subprocess.run(['root','-l','-b','-q',f'{stem}.C'],cwd=BASE,text=True,capture_output=True)
    print(result.stdout)
    if result.returncode or 'error:' in result.stderr:
        raise RuntimeError(result.stderr)
    # Keep the plots made by the compiled C++ code in its notebook as well.
    nb_cpp=nbformat.read(BASE / f'{stem}_CPP.ipynb',as_version=4)
    code_cells=[cell for cell in nb_cpp.cells if cell.cell_type=='code']
    for cell,match in zip(code_cells,pattern.finditer(source)):
        fig=match.group(3)
        if fig:
            data=base64.b64encode((BASE/'figures'/f'{fig}-cpp.png').read_bytes()).decode()
            cell.outputs=[nbformat.v4.new_output('display_data',data={'image/png':data})]
    nbformat.write(nb_cpp,BASE / f'{stem}_CPP.ipynb')
    body='\n'.join(pieces)
    end_title=body.index('</h1>')+len('</h1>')
    render_page(stem, title, body[:end_title]+choice+body[end_title:])
    return namespace

if __name__ == '__main__':
    import os
    os.chdir(BASE)
    (BASE/'figures').mkdir(exist_ok=True)
    build('GIC_simulation')
    build('Experimental_Data_Analysis')
    readme=markdown2html((BASE/'ReadMe.md').read_text())
    render_page('data-guide','GIC 数据与运行说明',readme)
    for old, target, title in [
        ('Simulation_Pulses','GIC_simulation.html#persistence','模拟脉冲与角度依赖'),
        ('comparison_of_pulse_shapes','Experimental_Data_Analysis.html#comparison','实验与模拟脉冲比较')]:
        render_page(old,title,f'<h1>{title}</h1><p>本节已整合到<a href="{target}">{title}</a>，图、方法与实例代码在同一页面中连续阅读。</p>')
