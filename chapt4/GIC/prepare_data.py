"""Copy the original teaching inputs; never modifies the original files.

python3 prepare_data.py --experimental /path/to/fall.root --processed /path/to/fana.root
"""
import argparse
import hashlib
import json
from pathlib import Path
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument('--experimental', type=Path, required=True)
parser.add_argument('--processed', type=Path, required=True)
args = parser.parse_args()
base = Path(__file__).resolve().parent
rows = []
for line in (base/'alpha_0.922barAr+3.45%CO2.txt').read_text().splitlines():
    a = line.split()
    if len(a) != 10 or a[1] not in ('keV','MeV'):
        continue
    energy = float(a[0]) * (0.001 if a[1]=='keV' else 1)
    stopping = (float(a[2])+float(a[3]))*0.16499
    projected = float(a[4]) * (0.001 if a[5]=='um' else 1)
    rows.append(f'{energy:.8g} {stopping:.9g} {projected:.8g}\n')
(base/'srim_stopping.txt').write_text(''.join(rows))

original = ROOT.TFile.Open(str(args.experimental))
tree = original.Get('tree')
sample = ROOT.TFile(str(base/'gic_experiment.root'),'RECREATE')
entries = int(tree.GetEntries())
copy = tree.CloneTree()
copy.Write()
sample.Close()
original.Close()

original = ROOT.TFile.Open(str(args.processed))
tree = original.Get('tree')
if tree.GetEntries() != entries:
    raise ValueError('Raw and processed event counts differ')
sample = ROOT.TFile(str(base/'gic_filtered.root'),'RECREATE')
copy = tree.CloneTree()
copy.Write()
sample.Close()
original.Close()
manifest = {
 'experiment': {'source':args.experimental.name,'source_sha256':hashlib.sha256(args.experimental.read_bytes()).hexdigest(),
                'source_entries':entries,'entries':entries,'selection':'all entries, original order; no waveform transformation',
                'sample_period_ns':10,'branches':'cathod[4096]/F, anode[4096]/F'},
 'filtered': {'source':args.processed.name,'source_sha256':hashlib.sha256(args.processed.read_bytes()).hexdigest(),
              'entries':entries,'selection':'all original processed entries; no new filtering',
              'kind':'baseline subtraction, 101-point median, 21-point moving average; ea/ec are signed filtered peak amplitudes in ADC'},
 'simulation': {'generator':'generate_grid.C','energies_MeV':[1,2,3,4,4.775,5.155,5.499,5.805,6],
                'theta_degrees':'0 to 90 inclusive, step 5','entries':171,
                'kind':'deterministic energy-angle grid using the model developed in the lesson; no event-yield sampling'},
 'srim': {'source':'alpha_0.922barAr+3.45%CO2.txt','columns':['E_MeV','stopping_MeV_per_mm','projected_range_mm'],
          'conversion':'(electronic + nuclear) * 0.16499, from original SRIM unit conversion'}
}
(base/'data_provenance.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
