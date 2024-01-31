# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 18:18:31 2023

@author: jhodges
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

from plotting import getJHcolors, getPlotLimits
from algorithms import getMaterials, processCaseData, sortCases
from algorithms import findFds, buildFdsFile, runModel, load_csv
from algorithms import calculateUncertainty, plotMaterialExtraction, calculateUncertaintyBounds
from algorithms import getTimeAveragedPeak
from algorithms import getTimeAveragedEnergy, getTimeAveraged_timeToEnergy

if __name__ == "__main__":
    
    # Initialize variables
    fs=16
    lw = 3
    windowSize = 60
    percentile = 90
    runSimulations = True
    colors = getJHcolors()
    
    # Read data
    spec_file_dict = getMaterials()
    materials = list(spec_file_dict.keys())
    
    # Prepare directories
    fdsdir, fdscmd = findFds()
    fileDir = os.path.dirname(os.path.abspath(__file__))
    
    output_statistics = dict()
    for material in materials:
        output_statistics[material] = dict()
        xlim, ylim = getPlotLimits(material)
        spec_file_dict[material] = processCaseData(spec_file_dict[material])
        
        properties = spec_file_dict[material]
        cases = properties['cases']
        
        totalEnergyMax = np.nanmax([cases[c]['totalEnergy'] for c in cases])
        
        if totalEnergyMax < 100:
            print("Total energy for %s is %0.1f < 100, skipping"%(material, totalEnergyMax))
            continue
            
        times_trimmed = [cases[c]['times_trimmed'] for c in cases]
        hrrs_trimmed = [cases[c]['hrrs_trimmed'] for c in cases]
        
        thicknesses = [cases[c]['delta'] for c in cases]
        fluxes = [cases[c]['cone'] for c in cases]
        
        chid = material.replace(' ','_')

        fluxes, deltas, tigns, cases_to_plot = sortCases(cases)
        
        workingDir = fileDir + os.sep +'..' + os.sep + 'input_files' + os.sep+ material + os.sep
        workingDir = workingDir.replace(' ', '_')
        
        if os.path.exists(workingDir) is False: os.mkdir(workingDir)
        # Calculate times to ignition
        Tign = 20
        front_h = 0
        txt = buildFdsFile(chid, cases, properties, Tign, front_h,
                           ignitionMode='Time', calculateDevcDt=False)
        
        with open("%s%s%s.fds"%(workingDir, os.sep, chid), 'w') as f:
            f.write(txt)
        
        if runSimulations:
            runModel(workingDir, chid+".fds", 1, fdsdir, fdscmd, printLiveOutput=False)
        
        
        
        data = load_csv(workingDir, chid)
        for i in range(0, len(cases)):
            c = cases_to_plot[i]
            namespace = ('CONE_%03.2f_%03d'%(cases[c]['delta']*1e3, cases[c]['cone'])).replace('.','p')
            label = r'%0.0f $\mathrm{kW/m^{2}}$'%(fluxes[i])
            delta0 = cases[c]['delta']
            times = data['Time'].values
            hrrpuas = data['"HRRPUA-'+namespace+'"'].values
            
            mod_peak = getTimeAveragedPeak(times, hrrpuas, windowSize)
            exp_peak = getTimeAveragedPeak(cases[c]['times'],cases[c]['HRRs'], windowSize) #, referenceTimes=times)
            
            output_statistics[material][c] = dict()
            
            energyThreshold, exp_t = getTimeAveragedEnergy(cases[c]['times']-cases[c]['tign'],cases[c]['HRRs'], windowSize, percentile)
            mod_t, timeAverage = getTimeAveraged_timeToEnergy(times-cases[c]['tign'], hrrpuas, windowSize, energyThreshold)
            
            output_statistics[material][c]['%0.0f_exp'%(percentile)] = exp_t
            output_statistics[material][c]['%0.0f_mod'%(percentile)] = mod_t
            print('%s\t\t%s\t\t%0.1f\t\t%0.1f\t\t%0.1f\t\t%0.1f'%(material, namespace, exp_peak, mod_peak, exp_t, mod_t))
            
            output_statistics[material][c]['delta'] = delta0
            output_statistics[material][c]['coneExposure'] = fluxes[i]
            output_statistics[material][c]['Qpeak_60s_exp'] = exp_peak
            output_statistics[material][c]['Qpeak_60s_mod'] = mod_peak
    
    metric_outputs = dict()
    metrics = ['Qpeak_60s', '90'] #, '10', '25','50','75','90']
    metricLabels = ['60s Avg', '90% Energy Consumed'] #'10% Energy Consumed', '25% Energy Consumed', '50% Energy Consumed', '75% Energy Consumed', '90% Energy Consumed']
    logs = [False, True, True, True, True, True]
    axmins = [0.0, 1e1, 1e1, 1e1, 1e1, 1e1]
    axmaxs = [2500, 1e4, 1e4, 1e4, 1e4, 1e4]
    
    for ii in range(0, len(metrics)): 
        metric = metrics[ii]
        label = metricLabels[ii]
        loglog = logs[ii]
        axmin = axmins[ii]
        axmax = axmaxs[ii]
        metric_outputs[metric] = dict()
        exp_metric = metric + '_exp'
        mod_metric = metric + '_mod'
        exp_points = []
        mod_points = []
        f = []
        ddd = []
        m = []
        mc = []
        c = []
        mask = []
        bases = []
        for material in materials:
            matClass = spec_file_dict[material]['materialClass']
            for case in list(output_statistics[material].keys()):
                if output_statistics[material][case]['Qpeak_60s_exp'] > 100:
                    mask.append(True)
                    bases.append(False) #bases.append(output_statistics[material][case]['basis'])
                    exp_points.append(output_statistics[material][case][exp_metric])
                    mod_points.append(output_statistics[material][case][mod_metric])
                    ddd.append(output_statistics[material][case]['delta'])
                    m.append(material)
                    mc.append(matClass)
                    c.append(output_statistics[material][case]['coneExposure'])
        mask = [x is False for x in bases]
        exp_points = np.array(exp_points)
        mod_points = np.array(mod_points)
        f = np.array(f)
        m = np.array(m)
        mc = np.array(mc)
        c = np.array(c)
        ddd = np.array(ddd)
        mask = np.array(mask)
        bias, sigma_m, sigma_e, points = calculateUncertainty(exp_points[mask], mod_points[mask])
        metric_outputs[metric]['bias'] = bias
        metric_outputs[metric]['sigma_m'] = sigma_m
        metric_outputs[metric]['sigma_e'] = sigma_e
        metric_outputs[metric]['points'] = points
        metric_outputs[metric]['exp_points'] = exp_points[mask]
        metric_outputs[metric]['mod_points'] = mod_points[mask]
        
        inds = np.argsort(abs(points))[::-1]
        worst_mats = m[inds]
        worst_cla = mc[inds]
        worst_c = c[inds]
        worst_ddd = ddd[inds]
        worst_pts = points[inds]
        for iiii in range(0, 10):
            print(worst_pts[iiii], worst_mats[iiii], worst_ddd[iiii], worst_c[iiii])
        print(metric)
        for uncertaintyStatistic in ['delta', 'exposure', 'materialClass']: # material
            if uncertaintyStatistic == 'delta':
                split = ddd
                for iii in range(0, len(m)):
                    ref_d = spec_file_dict[m[iii]]['case_basis']['case-1000']['delta']
                    case_d = ddd[iii]
                    
                    if ref_d < ddd[iii]:
                        split[iii] = 1
                    elif ref_d > ddd[iii]:
                        split[iii] = -1
                    else:
                        split[iii] = 0.1
                diff2 = split
                labelNames = {1 : r'$\mathrm{\delta > \delta_{ref}}$', 0.1 : r'$\mathrm{\delta = \delta_{ref}}$', -1 : r'$\mathrm{\delta < \delta_{ref}}$' }
                fname = 'uncertainty_statistics_delta_%s.png'%(metric)
            elif uncertaintyStatistic == 'material':
                split = m
                diff2 = m
                labelNames = {}
                for mat in materials: labelNames[mat] = mat.replace("2","")
                fname = 'uncertainty_statistics_materialClass_%s.png'%(metric)
            elif uncertaintyStatistic == 'materialClass':
                split = mc
                diff2 = mc
                labelNames = {'Others': 'Others', 'Polymers': 'Polymers', 'Wood-Based': 'Wood-Based', 'Mixtures': 'Mixtures'}
                fname = 'uncertainty_statistics_material_%s.png'%(metric)
            elif uncertaintyStatistic == 'exposure':
                split = c
                for iii in range(0, len(m)):
                    ref_q = spec_file_dict[m[iii]]['case_basis']['case-1000']['cone']
                    case_q = c[iii]
                    
                    if ref_q < c[iii]:
                        split[iii] = 1
                    elif ref_q > c[iii]:
                        split[iii] = -1
                    else:
                        split[iii] = 0.1
                diff2 = split
                #diff2 = [25 if x < 40 else x for x in c]
                #diff2 = [50 if (x > 40) and (x < 60) else x for x in diff2]
                #diff2 = [75 if (x > 60) else x for x in diff2]
                #labelNames = {25 : "$\mathrm{q_{cone}'' ~25 kW/m^{2}}$", 50 : "$\mathrm{q_{cone}'' ~50 kW/m^{2}}$", 75 : "$\mathrm{q_{cone}'' ~75 kW/m^{2}}$"}
                labelNames = {-1 : r"$\mathrm{q_{cone}''}$ < $\mathrm{q_{ref}''}$", 0.1 : r"$\mathrm{q_{cone}''}$ = $\mathrm{q_{ref}''}$", 1 : r"$\mathrm{q_{cone}''}$ > $\mathrm{q_{ref}''}$"}
                fname = 'uncertainty_statistics_exposure_%s.png'%(metric)
            
            fig, sigma_m, delta = plotMaterialExtraction(exp_points, mod_points, split, label, diff=diff2, axmin=axmin, axmax=axmax, loglog=loglog, labelName=labelNames, mask=mask)
            delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(exp_points, mod_points, split, split=True)
            
            metric_outputs[metric][uncertaintyStatistic] = dict()
            metric_outputs[metric][uncertaintyStatistic]['bias'] = delta
            metric_outputs[metric][uncertaintyStatistic]['sigma_m'] = sigma_m
            
            for key in list(delta.keys()):
                print('%s & %0.2f & %0.2f'%(key, delta[key], sigma_m[key]))
            plt.savefig('../figures/'+fname, dpi=300)
        detailed_output = pd.DataFrame(np.array([mc, m, ddd, c, exp_points, mod_points]).T, columns=['Class','Material','Thickness','Exposure','Exp','Mod'])
        detailed_output.to_csv('../output/output_%s.csv'%(metric))

    with pd.ExcelWriter('../figures/metrics_output.xlsx') as writer:
        tosave = dict()
        for key in list(metric_outputs.keys()):
            for v in ['bias', 'sigma_m']:
                tosave[key+'-'+v] = dict()
                tosave[key+'-'+v]['Mixtures'] = metric_outputs[key]['materialClass'][v]['Mixtures']
                tosave[key+'-'+v]['Others'] = metric_outputs[key]['materialClass'][v]['Others']
                tosave[key+'-'+v]['Polymers'] = metric_outputs[key]['materialClass'][v]['Polymers']
                tosave[key+'-'+v]['Wood-Based'] = metric_outputs[key]['materialClass'][v]['Wood-Based']
                tosave[key+'-'+v]['cone < ref'] = metric_outputs[key]['exposure'][v][-1]
                tosave[key+'-'+v]['cone = ref'] = metric_outputs[key]['exposure'][v][0.1]
                tosave[key+'-'+v]['cone > ref'] = metric_outputs[key]['exposure'][v][1]
                tosave[key+'-'+v]['delta < ref'] = metric_outputs[key]['delta'][v][-1]
                tosave[key+'-'+v]['delta = ref'] = metric_outputs[key]['delta'][v][0.1]
                tosave[key+'-'+v]['delta > ref'] = metric_outputs[key]['delta'][v][1]
                tosave[key+'-'+v]['total'] = metric_outputs[key][v]
        tosave_pd = pd.DataFrame(tosave)
        tosave_pd.to_excel(writer, sheet_name='FDS')
    
    