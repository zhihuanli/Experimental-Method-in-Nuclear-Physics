"""Check copied data, model invariants, histogram counts and generated lessons.

Run after build_lessons.py. Optional --original-data points to the GIC data/
directory containing experimental/raw/fall.root and experimental/processed/fana.root.
"""
import argparse
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

import nbformat
import numpy as np
import ROOT
from bs4 import BeautifulSoup

parser = argparse.ArgumentParser()
parser.add_argument('--original-data', type=Path)
args = parser.parse_args()
BASE = Path(__file__).resolve().parent
ROOT.gROOT.SetBatch(True)

if args.original_data:
    for local, original in [('gic_experiment.root','experimental/raw/fall.root'),
                            ('gic_filtered.root','experimental/processed/fana.root')]:
        source = ROOT.TFile.Open(str(args.original_data/original))
        copied = ROOT.TFile.Open(str(BASE/local))
        left, right = source.Get('tree'), copied.Get('tree')
        assert left.GetEntries() == right.GetEntries() == 5500
        names = [b.GetName() for b in left.GetListOfBranches()]
        assert names == [b.GetName() for b in right.GetListOfBranches()]
        for j in range(5500):
            left.GetEntry(j)
            right.GetEntry(j)
            for name in names:
                assert np.array_equal(np.asarray(getattr(left,name)),np.asarray(getattr(right,name))), (local,j,name)
        print(local, ': all branches and all 5500 entries match the original')
        source.Close()
        copied.Close()

f = ROOT.TFile.Open(str(BASE/'gic_grid.root'))
grid = f.Get('wave')
assert grid.GetEntries() == 171
angles, energies = set(), set()
templates = np.loadtxt(BASE/'gic_templates.txt')
for event in grid:
    angles.add(event.theta)
    energies.add(event.e)
    assert abs(event.qa-event.e) < 1e-8
    expected = event.e*(1-event.centroid/60.35*np.cos(np.deg2rad(event.theta)))
    assert abs(event.qc-expected) < 1e-8
    for current, charge in [('ic','qc'),('ia','qa')]:
        assert abs(np.sum(np.array(getattr(event,current)))*event.dt-getattr(event,charge)) < 1e-8
    for response, peak in [('uc','pc'),('ua','pa')]:
        assert abs(max(getattr(event,response))-getattr(event,peak)) < 1e-10
    if event.e == 5.805 and event.theta in (0,90):
        col = 1 if event.theta == 0 else 3
        assert np.max(np.abs(np.array(event.uc)-templates[:,col])) < 1e-10
        assert np.max(np.abs(np.array(event.ua)-templates[:,col+1])) < 1e-10
assert angles == set(range(0,91,5)) and len(energies) == 9
print('171 model records: charge integrals, angular dependence, peaks and 0/90-degree templates pass')
for expr, bins, low, high in [('ic:(Iteration$+.5)*dt',250,0,2.5),
                              ('uc:(Iteration$+1)*dt',5000,.005,50.005)]:
    h = ROOT.TH2D('bin_check','',bins,low,high,160,0,16)
    grid.Draw(expr+'>>bin_check',f'Iteration$<{bins}','goff')
    projection = h.ProjectionX('px',0,-1)
    assert h.GetEntries() == 171*bins
    assert all(projection.GetBinContent(j)==171 for j in range(1,bins+1))
    h.Delete()
    projection.Delete()
f.Close()

f = ROOT.TFile.Open(str(BASE/'gic_experiment.root'))
tree = f.Get('tree')
h = ROOT.TH2D('raw_bins','',1024,-.005,40.955,512,0,16384)
tree.Draw('cathod:Iteration$*.01>>raw_bins','','goff')
px = h.ProjectionX('raw_px',0,-1)
assert h.GetEntries() == 5500*4096
assert all(px.GetBinContent(j)==5500*4 for j in range(1,1025))
expected = np.loadtxt(BASE/'gic_amplitudes.txt')
actual = []
for event in tree:
    amplitudes = []
    for branch, sign in [(event.cathod,1),(event.anode,-1)]:
        y = np.array(branch,dtype=float)
        y = sign*(y-y[:1500].mean())
        q = y+.01/27*np.cumsum(y)
        mean = np.convolve(q,np.ones(100)/100,mode='full')[:4096]
        shaped = mean.copy()
        shaped[400:] -= mean[:-400]
        peak = 2200+np.argmax(shaped[2200:3000])
        amplitudes.append(shaped[peak-10:peak+11].mean())
    actual.append(amplitudes)
delta = np.max(np.abs(np.array(actual)-expected))
assert delta < 1e-7
print('Waveform time bins: exact expected occupancy; no missing or duplicated sample columns')
print(f'PyROOT vs C++: 5500 shaped amplitude pairs, max difference {delta:.3g} ADC')
f.Close()

for stem in ['GIC_simulation','Experimental_Data_Analysis']:
    source = (BASE/'_sources'/f'{stem}.md').read_text()
    page = BeautifulSoup((BASE/f'{stem}.html').read_text(),'html.parser')
    figures = re.findall(r'<!-- figure: ([\w-]+) -->', source)
    assert len(page.select('.python')) == len(page.select('.cpp'))
    assert not page.select('.downloads')
    assert 'scope8' not in page.get_text().lower()
    for name in figures:
        img = page.find('img',src=f'figures/{name}.png')
        assert img is not None
        previous = img.parent.find_previous_sibling()
        assert previous.name == 'pre' or 'cpp' in previous.get('class',[])
    for tag in page.select('[href], img[src]'):
        url = urlsplit(tag.get('href',tag.get('src','')))
        if url.scheme or not url.path:
            continue
        target = BASE/unquote(url.path)
        assert target.exists() or target.with_suffix('.md').exists(), target
    for language in ['PyROOT','CPP']:
        notebook = nbformat.read(BASE/f'{stem}_{language}.ipynb',as_version=4)
        nbformat.validate(notebook)
        outputs = [o for c in notebook.cells if c.cell_type=='code' for o in c.outputs]
        assert not any(o.output_type=='error' for o in outputs)
        assert sum('image/png' in o.get('data',{}) for o in outputs) == len(figures)
    print(stem, ': paired code,',len(figures),'generated plots, notebook outputs and links pass')

repo = BASE.parents[1]
layout = BeautifulSoup((repo/'_layouts/default.html').read_text(),'html.parser')
assert not layout.select('.course-sidebar')
assert 'grid-template-columns' not in (repo/'assets/css/style.scss').read_text()
resources = (repo/'reference.md').read_text()
assert 'coursework2.1/' not in resources and 'coursework3.1/' not in resources
assert 'chapt1/telescope/telescope.html' not in resources and 'chapt2/2.centerlimit.html' not in resources
print('Course sidebar removed; assignment-specific resource entries removed')
