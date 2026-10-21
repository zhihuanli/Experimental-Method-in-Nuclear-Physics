"""Numerical invariants and PyROOT/C++ agreement for the fitting lessons."""
from pathlib import Path
import json
import math
import re

KERNEL_CHECKS = {
    'ML_Fit_Comparison': '''
import json
low_counts = [int(h_low.GetBinContent(i)) for i in range(1,21)]
high_counts = [int(h_high.GetBinContent(i)) for i in range(1,21)]
assert len(low) == len(raw_low) == h_low.Integral() == 50
assert len(high) == len(raw_high) == h_high.Integral() == 10000
assert low_counts == [26,9,5,3,3,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0]
assert high_counts == [3933,2372,1462,935,487,332,187,116,60,50,20,21,9,7,4,1,1,0,1,2]
assert abs(sum(low)/50 - 0.8799174538946442) < 1e-12
assert abs(sum(high)/10000 - 0.997298694887933) < 1e-12
fits = [r_poisson_low,r_ls_low,r_poisson_high,r_ls_high]
assert all(int(r)==0 and r.CovMatrixStatus()==3 for r in fits)
assert r_ls_low.Ndf()+2 == sum(n>0 for n in low_counts) == 9
assert r_poisson_low.Ndf()+2 == 20
assert abs(tau_multinomial-poisson_low.GetParameter(1)) < 2e-5
assert abs(poisson_low.GetParameter(0)-50) < 1e-3
assert abs(poisson_high.GetParameter(0)-10000) < 0.02
for endpoint in [lower_low,upper_low]:
    assert abs(nll_low.Eval(endpoint)-nll_low.Eval(tau_low)-0.5) < 1e-6
for nll, tau, data, err in [(nll_low,tau_low,low,error_low),(nll_high,tau_high,high,error_high)]:
    # Independent analytic score and curvature at the optimum.
    n=len(data)
    score=n/tau-sum(data)/tau**2-n*T/(tau**2*math.expm1(T/tau))
    curvature=n/tau**2-n*T*T*math.exp(T/tau)/(tau**4*math.expm1(T/tau)**2)
    assert abs(score) < 0.002
    assert math.isclose(err,1/math.sqrt(curvature),rel_tol=1e-5)
for curve, model in zip(curves_low,[unbinned_low,poisson_low,ls_low]):
    assert abs(curve.Integral()-model.GetParameter(0)) < 1e-5
assert h_low.GetBinError(9)==0 and display_low.GetBinErrorUp(9)>0
for histogram,model,r in [(h_low,poisson_low,r_poisson_low),(h_high,poisson_high,r_poisson_high)]:
    counts=[histogram.GetBinContent(i) for i in range(1,21)]
    expected=[model.Integral((i-1)*width,i*width)/width for i in range(1,21)]
    half_deviance=sum(v-n+(n*math.log(n/v) if n else 0) for n,v in zip(counts,expected))
    assert abs(half_deviance-r.MinFcnValue())<1e-5
for histogram,model,r in [(h_low,ls_low,r_ls_low),(h_high,ls_high,r_ls_high)]:
    chi2=sum((histogram.GetBinContent(i)-model.Integral((i-1)*width,i*width)/width)**2/histogram.GetBinContent(i)
             for i in range(1,21) if histogram.GetBinContent(i)>0)
    assert abs(chi2-r.Chi2())<1e-5
print(json.dumps({"checks":"all passed","low_bin_counts":low_counts,"high_bin_counts":high_counts,
"unbinned_low":{"tau":tau_low,"error":error_low,"interval":[lower_low,upper_low]},
"poisson_low":{"tau":poisson_low.GetParameter(1),"error":poisson_low.GetParError(1)},
"multinomial_low":tau_multinomial,"ls_low":{"tau":ls_low.GetParameter(1),"error":ls_low.GetParError(1)},
"unbinned_high":{"tau":tau_high,"error":error_high},
"poisson_high":{"tau":poisson_high.GetParameter(1),"error":poisson_high.GetParError(1)},
"ls_high":{"tau":ls_high.GetParameter(1),"error":ls_high.GetParError(1)}}))
''',
    'Weighted_Fit_Uncertainty': '''
import json
assert int(result)==0 and result.CovMatrixStatus()==3 and result.Ndf()==4
assert list(x)==[0,2,4,6,8,10] and list(y)==[5,8,14,20,24,35] and list(ey)==[2,1.5,2.5,1,3,2]
# Independent closed-form weighted least-squares solution.
w=[1/e**2 for e in ey]
S=sum(w); Sx=sum(wi*xx for wi,xx in zip(w,x)); Sxx=sum(wi*xx**2 for wi,xx in zip(w,x))
Sy=sum(wi*yy for wi,yy in zip(w,y)); Sxy=sum(wi*xx*yy for wi,xx,yy in zip(w,x,y))
D=S*Sxx-Sx*Sx
analytic_a=(Sxx*Sy-Sx*Sxy)/D; analytic_b=(S*Sxy-Sx*Sy)/D
assert abs(line.GetParameter(0)-analytic_a)<1e-10 and abs(line.GetParameter(1)-analytic_b)<1e-10
assert abs(C[0][0]-Sxx/D)<1e-10 and abs(C[0][1]+Sx/D)<1e-10 and abs(C[1][1]-S/D)<1e-10
assert abs(sum(((yy-line.Eval(xx))/ee)**2 for xx,yy,ee in zip(x,y,ey))-result.Chi2())<1e-10
assert not result.NormalizedErrors()
assert abs(fx0-17.574505849132713)<1e-10
assert abs(sigma_fx0-0.6720369954664063)<1e-10
assert abs(ci[0]-sigma_fx0)<1e-8
assert abs(diagonal_only-1.7590045713859828)<1e-10
for i in range(band.GetN()):
    xx=band.GetPointX(i)
    assert abs(band.GetErrorY(i)**2-(C[0][0]+2*xx*C[0][1]+xx**2*C[1][1]))<1e-10
print(json.dumps({"checks":"all passed","parameters":[line.GetParameter(0),line.GetParameter(1)],
"covariance":[[C[i][j] for j in range(2)] for i in range(2)],"correlation":result.Correlation(0,1),
"chi2":result.Chi2(),"ndf":result.Ndf(),"f5":fx0,"sigma_f5":sigma_fx0,"root_interval":ci[0],"omitting_covariance":diagonal_only}))
''',
    'ROOT_Fit_Tutorial': '''
import json
assert all(int(r) == 0 for r in (line_result, fit, r_poisson, r_ls, peak_result, minos_result))
assert fit.CovMatrixStatus() == 3 and peak_result.CovMatrixStatus() == 3
assert Cp[0][0] > 0 and Cp[2][2] > 0 and area_var > 0
assert abs(signal.Integral(2.6,3.8)/width - area) < 0.001
assert 0 < roi_area < area and 0 < roi_error < area
assert abs(deviance - 2*peak_result.MinFcnValue()) < 0.01
assert abs(C[0][1]) < math.sqrt(C[0][0]*C[1][1])
for i in range(ellipse.GetN()):
    da = ellipse.GetPointX(i)-line.GetParameter(0)
    db = ellipse.GetPointY(i)-line.GetParameter(1)
    q = (C[1][1]*da*da-2*C[0][1]*da*db+C[0][0]*db*db)/(C[0][0]*C[1][1]-C[0][1]**2)
    assert abs(q-2.30) < 1e-8
from pathlib import Path
import hashlib
reference = json.loads(Path("data/lecture_mass_reference.json").read_text())
assert hashlib.sha256(Path("data/lecture_mass_events.csv").read_bytes()).hexdigest() == reference["sha256"]
assert len(mass_values) == hp.Integral() == 600
assert [hp.GetBinContent(i) for i in range(1,49)] == reference["bin_counts"]
assert int(peak_result.NPar()) == 5
lecture_parameters = {}
for j, name in enumerate(reference["parameters"]):
    expected = reference["binned"]["parameters"][name]
    value, error = model.GetParameter(j), math.sqrt(Cp[j][j])
    assert abs(value-expected["value"]) < 0.003*expected["error"], (name,value,expected)
    assert math.isclose(error, expected["error"], rel_tol=0.015), (name,error,expected)
    lecture_parameters[name] = {"value":value,"error":error}
assert abs(model.GetParameter(0)+model.GetParameter(1)-600) < 0.03
assert max(signal.Eval(2.6+1.2*i/600) for i in range(601)) > 15
print(json.dumps({"lecture_binned":lecture_parameters,"data_sha256":reference["sha256"],"area":area,"area_error":math.sqrt(area_var),"roi_area":roi_area,"roi_error":roi_error,"deviance":deviance,"fwhm":k_fwhm*sigma,"checks":"all passed"}))
''',
    'RooFit_Tutorial': '''
import json
assert all(r.status() == 0 for r in (gauss_result, shape_result, result, binned_shape_result, binned_result, time_result, time_binned_result, manual_result, sim_result))
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
from pathlib import Path
import hashlib
reference = json.loads(Path("data/lecture_mass_reference.json").read_text())
assert hashlib.sha256(Path("data/lecture_mass_events.csv").read_bytes()).hexdigest() == reference["sha256"]
assert len(mass_values) == data.numEntries() == binned.sumEntries() == 600
assert binned.numEntries() == 48
assert [binned.get(i) and binned.weight() for i in range(48)] == reference["bin_counts"]
assert result.floatParsFinal().getSize() == binned_result.floatParsFinal().getSize() == 5
assert shape_result.floatParsFinal().getSize() == binned_shape_result.floatParsFinal().getSize() == 4
lecture_results = {}
for kind, fit_result in [("unbinned",result),("binned",binned_result)]:
    parameters = {}
    for name in reference["parameters"]:
        par = fit_result.floatParsFinal().find(name.lower() if name in ("Ns","Nb") else name)
        expected = reference[kind]["parameters"][name]
        assert abs(par.getVal()-expected["value"]) < 0.003*expected["error"], (kind,name,par.getVal(),expected)
        assert math.isclose(par.getError(),expected["error"],rel_tol=0.015), (kind,name,par.getError(),expected)
        parameters[name] = {"value":par.getVal(),"error":par.getError()}
    lecture_results[kind] = parameters
for shape_fit, ext_fit in [(shape_result,result),(binned_shape_result,binned_result)]:
    for name in ("mu","sigma","beta"):
        sp, ep = shape_fit.floatParsFinal().find(name), ext_fit.floatParsFinal().find(name)
        assert abs(sp.getVal()-ep.getVal()) < 0.005*ep.getError()
    frac = shape_fit.floatParsFinal().find("frac").getVal()
    assert abs(frac-ext_fit.floatParsFinal().find("ns").getVal()/600) < 1e-4
print(json.dumps({"lecture":lecture_results,"data_sha256":reference["sha256"],"ns":ns.getVal(),"ns_error":ns.getError(),"roi_signal":signal_roi.getVal(),"roi_error":signal_roi.getPropagatedError(result),"events":data.numEntries(),"toy_successful":accepted,"toy_pull_mean":pull_mean,"toy_pull_sd":pull_sd,"checks":"all passed"}))
''',
}


def validate_outputs(name, notebook, cpp_log, qa_text):
    python_log = '\n'.join(o.get('text', '') for c in notebook.cells if c.cell_type == 'code' for o in c.outputs)
    # The only known ROOT 6.40 warning is a duplicate internal curve-name option.
    warnings = [line for line in cpp_log.splitlines() if 'WARNING' in line or 'Warning' in line]
    assert all('CurveNameSuffix is duplicated' in line for line in warnings), warnings
    patterns = {
        'ML_Fit_Comparison': [rf'{name}: tau = ([\d.]+) \+/- ([\d.]+)' for name in ('Unbinned low','Unbinned high','Poisson low','Least squares low','Poisson high','Least squares high')],
        'Weighted_Fit_Uncertainty': [r'a = ([\d.]+) \+/- ([\d.]+)', r'b = ([\d.]+) \+/- ([\d.]+)', r'f\(5\) = ([\d.]+) \+/- ([\d.]+)', r'ROOT interval half-width = ([\d.]+)'],
        'ROOT_Fit_Tutorial': [r'Net peak area = ([\d.]+) \+/- ([\d.]+)', r'FWHM = ([\d.]+) \+/- ([\d.]+)', r'mean = ([\d.]+) \+/- ([\d.]+)', r'sigma = ([\d.]+) \+/- ([\d.]+)', r'Ns\s*= ([\d.]+) \+/- ([\d.]+)', r'Nb\s*= ([\d.]+) \+/- ([\d.]+)'],
        'RooFit_Tutorial': [r'Signal in ROI = ([\d.]+) \+/- ([\d.]+)', r'ns\s*= ([\d.]+) \+/- ([\d.]+)', r'Unbinned: tau = ([\d.]+) \+/- ([\d.]+)', r'shared_mean = ([\d.]+) \+/- ([\d.]+)', r'Pull mean = (-?[\d.]+); sample SD = ([\d.]+)'],
    }
    for pattern in patterns[name]:
        py_match = re.search(pattern, python_log)
        cpp_match = re.search(pattern, cpp_log)
        assert py_match and cpp_match, pattern
        for p, c in zip(py_match.groups(), cpp_match.groups()):
            tolerance = 2e-6 if name in ('ML_Fit_Comparison','Weighted_Fit_Uncertainty') else 0.011
            assert math.isclose(float(p), float(c), abs_tol=tolerance), (pattern, p, c)
    summary = json.loads(qa_text.strip().splitlines()[-1])
    summary['language_agreement'] = 'passed at displayed precision'
    summary['python_cells'] = sum(c.cell_type == 'code' for c in notebook.cells)
    summary['known_plot_warnings'] = len(warnings)
    (Path(__file__).parent / '_sources' / f'{name}_validation.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(name, 'numerical invariants and PyROOT/C++ agreement passed', flush=True)
