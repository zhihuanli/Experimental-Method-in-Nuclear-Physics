"""Prepare small, traceable teaching inputs; never modifies the original files.

python3 prepare_data.py --experimental /path/to/fall.root --scope8 /path/to/alpha.root
"""
import argparse
import hashlib
import json
from pathlib import Path
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument('--experimental', type=Path, required=True)
parser.add_argument('--scope8', type=Path, required=True)
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
copy = tree.CopyTree('', '', 1000, 0)
copy.Write()
sample.Close()
original.Close()

columns = ['event_id','energy_mev','theta_deg','track_range_mm',
           'anode_charge_mev','cathode_charge_mev']
frame = ROOT.RDataFrame('ChargeEvents',str(args.scope8))
frame.Snapshot('ChargeEvents',str(base/'scope8_gic_observables.root'),columns)
manifest = {
 'experiment': {'source':args.experimental.name,'source_sha256':hashlib.sha256(args.experimental.read_bytes()).hexdigest(),
                'source_entries':5500,'entries':1000,'selection':'first 1000 entries, original order; no waveform transformation',
                'sample_period_ns':10,'branches':'cathod[4096]/F, anode[4096]/F'},
 'scope8': {'source':args.scope8.name,'source_sha256':hashlib.sha256(args.scope8.read_bytes()).hexdigest(),
            'entries':int(frame.Count().GetValue()),'columns':columns,'selection':'all entries of the latest ChargeEvents cycle',
            'kind':'detector-level simulation; no electronics, ADC noise, or energy broadening added'},
 'srim': {'source':'alpha_0.922barAr+3.45%CO2.txt','columns':['E_MeV','stopping_MeV_per_mm','projected_range_mm'],
          'conversion':'(electronic + nuclear) * 0.16499, from original SRIM unit conversion'}
}
(base/'data_provenance.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(manifest,ensure_ascii=False,indent=2))
