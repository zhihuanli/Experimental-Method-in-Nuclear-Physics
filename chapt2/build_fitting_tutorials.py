"""Build two executed notebooks and their bilingual-code HTML pages.

Run: python3 build_fitting_tutorials.py
Requires an existing Python kernel with PyROOT, nbclient and nbconvert, plus root.
The paired code in _sources is the only source of displayed and executed code.
"""
from pathlib import Path
import html
import json
import re
import subprocess
import tempfile

import nbformat
from nbclient import NotebookClient
from nbconvert.filters.markdown import markdown2html
from bs4 import BeautifulSoup

from validate_fitting_tutorials import KERNEL_CHECKS, validate_outputs

BASE = Path(__file__).resolve().parent
PAIR = re.compile(r"```python\n(.*?)\n```\n\n```cpp\n(.*?)\n```(?:\n\n<!-- figure: ([\w-]+) -->)?", re.S)
STYLE = """
body{margin:0;color:#252525;background:white;font:17px/1.75 system-ui,sans-serif}
main{max-width:920px;margin:auto;padding:25px 24px 70px}h1{font-size:1.9rem;line-height:1.35}
h2{font-size:1.4rem;border-bottom:1px solid #ddd;padding-top:20px}h3{font-size:1.1rem}
a{color:#175c89}pre{overflow:auto;background:#f5f6f7;padding:14px;font:14px/1.55 monospace}
code{font-size:.9em}table{border-collapse:collapse;width:100%;font-size:15px}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}
figure{max-width:740px;margin:18px auto}img{max-width:100%;height:auto}figcaption{font-size:14px;color:#555}
button{padding:7px 16px;margin:5px 8px 5px 0;font:inherit;cursor:pointer;border:1px solid #aaa;border-radius:3px}
button[aria-pressed=true]{background:#175c89;color:white}.cpp{display:none}body[data-language=cpp] .python{display:none}
body[data-language=cpp] .cpp{display:block}.output{background:white;border-left:3px solid #bbb;white-space:pre-wrap}
nav{font-size:15px;margin-bottom:20px}.language-choice{position:sticky;top:0;background:#fffef5;z-index:5;padding:6px 10px}
blockquote{border-left:3px solid #bbb;margin-left:0;padding-left:18px}details{margin:10px 0}
"""
SCRIPT = """
function setLanguage(lang){document.body.dataset.language=lang;
document.querySelectorAll('[data-lang]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.lang===lang)));
try{localStorage.setItem('fit-language',lang)}catch(e){}}
try{setLanguage(localStorage.getItem('fit-language')||'python')}catch(e){setLanguage('python')}
"""
HEADERS = """#include <TROOT.h>
#include <TError.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Math/MinimizerOptions.h>
#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooWorkspace.h>
#include <RooFormulaVar.h>
#include <RooMinimizer.h>
#include <RooRandom.h>
#include <RooMsgService.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <stdexcept>
using namespace RooFit;
"""


def build(stem, macro_name):
    source = (BASE / '_sources' / f'{macro_name}.md').read_text()
    cells, cpp_blocks, pairs = [], [], []
    end = 0
    for match in PAIR.finditer(source):
        prose = source[end:match.start()].strip()
        if prose:
            cells.append(nbformat.v4.new_markdown_cell(prose))
        py, cpp, figure = match.groups()
        if figure:
            py += f'\nc.SaveAs("fit_figures/{figure}.png")\ndisplay(Image(filename="fit_figures/{figure}.png"))'
            cpp += f'\nc->SaveAs("fit_figures/{figure}.png");'
        cells.append(nbformat.v4.new_code_cell(py, metadata={'cpp': cpp, 'figure': figure or ''}))
        cells.append(nbformat.v4.new_markdown_cell(
            '<details><summary>ROOT C++ 对应代码</summary>\n\n```cpp\n' + cpp + '\n```\n\n</details>',
            metadata={'cpp_alternative': True}))
        cpp_blocks.append(cpp)
        pairs.append((py, cpp, figure))
        end = match.end()
    cells.append(nbformat.v4.new_markdown_cell(source[end:].strip()))
    for i, cell in enumerate(cells):
        cell.id = f'{macro_name.lower()}-{i:03d}'
    nb = nbformat.v4.new_notebook(cells=cells, metadata={
        'kernelspec': {'name': 'python3', 'display_name': 'Python 3 (PyROOT)', 'language': 'python'},
        'language_info': {'name': 'python'},
        'tutorial_source': 'L. Moneta, Fitting and Parameter Estimation, ROOT Tutorial at UERJ (2015)',
    })
    with tempfile.TemporaryDirectory(prefix='fitting-lessons-') as temporary:
        run_dir = Path(temporary)
        print(f'Executing {stem}: {len(pairs)} Python cells', flush=True)
        # Run numerical invariants in the same live kernel; omit this QA-only cell
        # from the teaching notebook, but retain its machine-readable result.
        nb.cells.append(nbformat.v4.new_code_cell(KERNEL_CHECKS[macro_name]))
        NotebookClient(nb, timeout=180, kernel_name='python3', resources={'metadata': {'path': temporary}}).execute()
        qa = nb.cells.pop()
        qa_text = ''.join(output.get('text', '') for output in qa.outputs)
        (BASE / 'fit_figures').mkdir(exist_ok=True)
        for figure in (run_dir / 'fit_figures').glob('*.png'):
            (BASE / 'fit_figures' / figure.name).write_bytes(figure.read_bytes())
        macro = HEADERS + f'\nvoid {macro_name}() {{\n' + '\n\n'.join(cpp_blocks) + '\n}\n'
        (BASE / f'{macro_name}.C').write_text(macro)
        (run_dir / f'{macro_name}.C').write_text(macro)
        result = subprocess.run(['root', '-l', '-b', '-q', f'{macro_name}.C'], cwd=temporary, text=True, capture_output=True)
        print(result.stdout[-5500:])
        if result.returncode or 'error:' in result.stderr:
            raise RuntimeError(result.stdout + result.stderr)
        # Displayed outputs are from the executed Python notebook, never typed results.
        validate_outputs(macro_name, nb, result.stdout + result.stderr, qa_text)
    nbformat.write(nb, BASE / f'{stem}.ipynb')
    render(stem, nb, source)


def render(stem, nb, source):
    fragments = []
    for cell in nb.cells:
        if cell.metadata.get('cpp_alternative'):
            continue
        if cell.cell_type == 'markdown':
            fragments.append(markdown2html(cell.source))
            continue
        fragments.append('<div class="python"><pre><code>' + html.escape(cell.source) + '</code></pre></div>')
        fragments.append('<div class="cpp"><pre><code>' + html.escape(cell.metadata.cpp) + '</code></pre></div>')
        for output in cell.outputs:
            if output.output_type == 'stream':
                rendered = '<pre class="output">' + html.escape(output.text) + '</pre>'
                if 'CurveNameSuffix is duplicated' in output.text:
                    rendered = '<details><summary>ROOT 误差带绘图提示（原始日志）</summary>' + rendered + '</details>'
                fragments.append(rendered)
            elif 'image/png' in output.get('data', {}):
                # PNG already belongs to this cell and sits immediately after its code/output.
                fragments.append(f'<figure><img src="fit_figures/{cell.metadata.figure}.png" alt="{cell.metadata.figure}"></figure>')
            elif 'text/plain' in output.get('data', {}):
                fragments.append('<pre class="output">' + html.escape(output.data['text/plain']) + '</pre>')
    body = '\n'.join(fragments)
    document = BeautifulSoup(body, 'html.parser')
    for anchor in document.select('a.anchor-link'):
        anchor.decompose()
    # nbconvert uses percent-encoded Chinese heading IDs. Use stable ASCII IDs
    # for the web TOC; explicit semantic anchors in the source remain intact.
    for i, heading in enumerate(document.select('h1,h2,h3')):
        heading['id'] = f'topic-{i}'
    toc = []
    for i, heading in enumerate(document.find_all('h2')):
        heading['id'] = heading.get('id', f'section-{i+1}')
        toc.append(f'<li><a href="#{heading["id"]}">{html.escape(heading.get_text())}</a></li>')
    body = str(document)
    choice = '<div class="language-choice">代码语言：<button data-lang="python" onclick="setLanguage(\'python\')" aria-pressed="true">PyROOT</button><button data-lang="cpp" onclick="setLanguage(\'cpp\')" aria-pressed="false">ROOT C++</button></div>'
    choice += '<details><summary>目录</summary><ol>' + ''.join(toc) + '</ol></details>'
    index = body.index('</h1>') + len('</h1>')
    body = body[:index] + choice + body[index:]
    title = source.splitlines()[0][2:]
    page = f'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title>
<style>{STYLE}</style><script>window.MathJax={{tex:{{inlineMath:[['$','$'],['\\\\(','\\\\)']]}}}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script></head>
<body data-language="python"><main><nav><a href="../reference.html">学习资料</a> · <a href="ROOT_Fit_Tutorial.html">ROOT Fit</a> · <a href="RooFit%20mini%20tutorial.html">RooFit</a></nav>
{body}</main><script>{SCRIPT}</script></body></html>'''
    (BASE / f'{stem}.html').write_text('\n'.join(line.rstrip() for line in page.splitlines())+'\n')
    print(f'Rendered {stem}.html', flush=True)


def update_related():
    """Add only an orientation link; preserve existing code, figures and outputs."""
    links = {
        'ML_fit': ('ROOT_Fit_Tutorial.html', 'ROOT Fit 教程', '这里保留少计数与高计数样本的专题比较。'),
        'linearfit_error band': ('ROOT_Fit_Tutorial.html#weighted-uncertainty', 'ROOT Fit：协方差与拟合值误差', '本页可作为误差传播的独立查阅实例。'),
        'likelihood_decay': ('RooFit%20mini%20tutorial.html#decay', 'RooFit：有限范围的衰变时间拟合', '完整的分步说明已纳入 RooFit 教程，本页保留原专题比较。'),
        'likelihood_signal_background': ('RooFit%20mini%20tutorial.html#signal-background', 'RooFit：signal + background', '完整的分步说明已纳入 RooFit 教程，本页保留原专题比较。'),
    }
    for stem, (target, title, text) in links.items():
        note = f'<p class="fitting-tutorial-guide">入门阅读：<a href="{target}">{title}</a>。{text}</p>'
        path = BASE / f'{stem}.html'
        page = path.read_text()
        page = re.sub(r'<p class="fitting-tutorial-guide">.*?</p>\s*', '', page)
        page = re.sub(r'(</h1>)', lambda m: m.group(1)+'\n'+note, page, count=1)
        path.write_text(page)
        path = BASE / f'{stem}.ipynb'
        notebook = nbformat.read(path, as_version=4)
        first = notebook.cells[0]
        first.source = re.sub(r'<p class="fitting-tutorial-guide">.*?</p>\s*', '', first.source)
        heading, content = first.source.split('\n', 1)
        first.source = heading + '\n\n' + note + '\n' + content
        nbformat.write(notebook, path)


if __name__ == '__main__':
    import sys
    lessons = [('ROOT_Fit_Tutorial', 'ROOT_Fit_Tutorial'), ('RooFit mini tutorial', 'RooFit_Tutorial')]
    for stem, macro_name in lessons:
        if '--render' in sys.argv:
            render(stem, nbformat.read(BASE / f'{stem}.ipynb', as_version=4), (BASE / '_sources' / f'{macro_name}.md').read_text())
        elif len(sys.argv) == 1 or macro_name in sys.argv[1:]:
            build(stem, macro_name)
    update_related()
