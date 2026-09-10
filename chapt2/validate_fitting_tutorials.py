"""Numerical invariants and PyROOT/C++ agreement for the fitting lessons."""
from pathlib import Path
import json
import math
import re

KERNEL_CHECKS = {
    'ROOT_Fit_Tutorial': '''
import json
assert all(int(r) == 0 for r in (line_result, fit, r_poisson, r_ls, side_result, peak_result, minos_result))
assert fit.CovMatrixStatus() == 3 and peak_result.CovMatrixStatus() == 3
assert Cp[0][0] > 0 and Cp[2][2] > 0 and area_var > 0
assert abs(signal.Integral(40,60)/width - area) < 0.001
assert 0 < roi_area < area and 0 < roi_error < area
assert abs(deviance - 2*peak_result.MinFcnValue()) < 0.01
assert abs(C[0][1]) < math.sqrt(C[0][0]*C[1][1])
for i in range(ellipse.GetN()):
    da = ellipse.GetPointX(i)-line.GetParameter(0)
    db = ellipse.GetPointY(i)-line.GetParameter(1)
    q = (C[1][1]*da*da-2*C[0][1]*da*db+C[0][0]*db*db)/(C[0][0]*C[1][1]-C[0][1]**2)
    assert abs(q-2.30) < 1e-8
print(json.dumps({"area":area,"area_error":math.sqrt(area_var),"roi_area":roi_area,"roi_error":roi_error,"deviance":deviance,"fwhm":k_fwhm*sigma,"checks":"all passed"}))
''',
    'RooFit_Tutorial': '''
import json
assert all(r.status() == 0 for r in (gauss_result, shape_result, result, binned_result, time_result, time_binned_result, manual_result, sim_result))
assert result.covQual() == 3 and binned_result.covQual() == 3 and sim_result.covQual() == 3
assert tree.GetEntries() == from_tree.numEntries() == 1000
assert h_import.Integral() == from_hist.sumEntries() == 1000
assert data.numEntries() == binned.sumEntries()
assert abs(ns.getVal()+nb.getVal()-data.numEntries()) < 0.001*data.numEntries()
assert 0 < signal_fraction.getVal() < 1
assert 0 < signal_roi.getVal() < ns.getVal()
assert signal_roi.getPropagatedError(result) > 0
assert combined.numEntries() == data_a.numEntries()+data_b.numEntries() == 800
assert uv_data.numEntries() == 10000 and h2.Integral() == 10000
assert failed == 0 and toy_pull.GetEntries() == 200
assert accepted == 200 and math.isfinite(pull_sd)
assert migrad_status == hesse_status == minos_status == 0
print(json.dumps({"ns":ns.getVal(),"ns_error":ns.getError(),"roi_signal":signal_roi.getVal(),"roi_error":signal_roi.getPropagatedError(result),"events":data.numEntries(),"toy_successful":accepted,"toy_pull_mean":pull_mean,"toy_pull_sd":pull_sd,"checks":"all passed"}))
''',
}


def validate_outputs(name, notebook, cpp_log, qa_text):
    python_log = '\n'.join(o.get('text', '') for c in notebook.cells if c.cell_type == 'code' for o in c.outputs)
    # The only known ROOT 6.40 warning is a duplicate internal curve-name option.
    warnings = [line for line in cpp_log.splitlines() if 'WARNING' in line or 'Warning' in line]
    assert all('CurveNameSuffix is duplicated' in line for line in warnings), warnings
    patterns = {
        'ROOT_Fit_Tutorial': [r'Net peak area = ([\d.]+) \+/- ([\d.]+)', r'FWHM = ([\d.]+) \+/- ([\d.]+)', r'mean = ([\d.]+) \+/- ([\d.]+)', r'sigma = ([\d.]+) \+/- ([\d.]+)'],
        'RooFit_Tutorial': [r'Signal in ROI = ([\d.]+) \+/- ([\d.]+)', r'ns = ([\d.]+) \+/- ([\d.]+)', r'Unbinned: tau = ([\d.]+) \+/- ([\d.]+)', r'shared_mean = ([\d.]+) \+/- ([\d.]+)', r'Pull mean = (-?[\d.]+); sample SD = ([\d.]+)'],
    }
    for pattern in patterns[name]:
        py_match = re.search(pattern, python_log)
        cpp_match = re.search(pattern, cpp_log)
        assert py_match and cpp_match, pattern
        for p, c in zip(py_match.groups(), cpp_match.groups()):
            assert math.isclose(float(p), float(c), abs_tol=0.011), (pattern, p, c)
    summary = json.loads(qa_text.strip().splitlines()[-1])
    summary['language_agreement'] = 'passed at displayed precision'
    summary['python_cells'] = sum(c.cell_type == 'code' for c in notebook.cells)
    summary['known_plot_warnings'] = len(warnings)
    (Path(__file__).parent / '_sources' / f'{name}_validation.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(name, 'numerical invariants and PyROOT/C++ agreement passed', flush=True)
