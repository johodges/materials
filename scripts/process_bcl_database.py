import os
import glob
import numpy as np
import pandas as pd
import subprocess
import matplotlib.pyplot as plt

from algorithms import findFds, load_csv

def getColors():
    colors = [[0,0.4470,0.7410],
              [0.8500,0.3250,0.0980],
              [0.9290,0.6940,0.1250],
              [0.4940,0.1840,0.5560],
              [0.4660,0.6740,0.1880],
              [0.3010,0.7450,0.9330],
              [0.6350,0.0780,0.1840]]
    colors = np.array(colors)
    return colors

def findLimits(times, HRRs, energyCutoff1, energyCutoff2):
    v = np.cumsum(HRRs)
    ind1 = 0 
    while ind1 == 0:
        try:
            ind1 = np.where(v < np.nanmax(v)*energyCutoff1)[0][-1]
        except:
            energyCutoff1 = energyCutoff1*2
            #print(energyCutoff1)
    ind2 = v.shape[0]
    while ind2 == v.shape[0]:
        try:
            ind2 = np.where(v > np.nanmax(v)*energyCutoff2)[0][0]
        except:
            energyCutoff2 = energyCutoff2*0.99
            #print(energyCutoff2)
    times_trimmed = times[ind1:ind2]
    hrrs_trimmed = HRRs[ind1:ind2]
    tign = times[ind1]
    return tign, times_trimmed, hrrs_trimmed


def buildFdsFile(chid, cases, properties, Tign, front_h,
                 ignitionMode='Temperature', outputTemperature=False,
                 calculateDevcDt=True, devc_dt=10.,
                 energyThreshold=0.0):
    ''' Generate a solid phase only FDS input file representing cone
    experimens at different exposures given a reference curve and
    material properties. The configuration can support up to 9
    thermal exposures, configured in a 3x3 grid.
    
    Notes:
        1. Zero heat transfer coefficient at the surface. This is
           intentional because the calculated flame heat flux which
           is included in the IGNITION-RAMPs includes convection.
        2. The estimated flame heat flux currently assumes a surface
           heat transfer coefficient of 15 W/m2-K and a gas phase
           radiative fraction of 0.35.
        3. If the ignition temperature is unknown, it can be calculated
           by comparing with experimental times to ignition. Changing
           the input of Tign to 'Calculated' will tell FDS to save
           out the WALL TEMPERATURE data needed to extract this
           information.
        4. All samples are assumed to have 0.5in / 12.7 mm of ceramic
           fiber insulation behind them.
    '''
    
    tend = np.nanmax([(cases[c]['times_trimmed'].max()+cases[c]['tign'])*2 for c in cases])
    
    density = properties['density']
    conductivity = properties['conductivity']
    emissivity = properties['emissivity']
    specific_heat = properties['specific_heat']
    soot_yield = max([properties['soot_yield'], 0.001])
    heat_of_combustion = properties['heat_of_combustion']
    
    tempOutput = '.TRUE.' if outputTemperature else '.FALSE.'
    DT_DEVC = devc_dt
    if calculateDevcDt:
        NFRAMES = 1200/1.
        DT_DEVC = tend/NFRAMES
    if ignitionMode == 'Time': Tign = -273.15
    txt = "&HEAD CHID='%s', /\n"%(chid)
    txt = txt+"&TIME DT=1., T_END=%0.1f /\n"%(tend)
    txt = txt+"&DUMP DT_CTRL=%0.1f, DT_DEVC=%0.1f, DT_HRR=%0.1f, SIG_FIGS=4, SIG_FIGS_EXP=2, /\n"%(DT_DEVC, DT_DEVC, DT_DEVC)
    txt = txt+"&MISC SOLID_PHASE_ONLY=.TRUE., TMPA=27., /\n"
    txt = txt+"&MESH ID='MESH', IJK=3,3,3, XB=0.,0.3,0.,0.3,0.,0.3, /\n"
    txt = txt+"&MESH ID='MESH', IJK=3,3,3, XB=0.6,0.9,0.,0.3,0.,0.3, /\n"
    txt = txt+"&REAC ID='PROPANE', FUEL='PROPANE', HEAT_OF_COMBUSTION=%0.0f, SOOT_YIELD=%0.3f /\n"%(heat_of_combustion*1e3, soot_yield)
    txt = txt+"&MATL ID='BACKING', CONDUCTIVITY=0.10, DENSITY=65., EMISSIVITY=0.9, SPECIFIC_HEAT=1.14, /\n"
    #txt = txt+"&MATL ID='BACKING', CONDUCTIVITY=0.2, DENSITY=585., EMISSIVITY=1., SPECIFIC_HEAT=0.8, /\n"
    txt = txt+"&MATL ID='SAMPLE', CONDUCTIVITY=%0.3f, DENSITY=%0.0f, EMISSIVITY=%0.3f, SPECIFIC_HEAT=%0.3f, /\n"%(conductivity, density, emissivity, specific_heat)
    
    totalEnergyMax = np.nanmax([cases[c]['totalEnergy'] for c in cases])
    filtered_cases = [c for c in cases if cases[c]['totalEnergy'] > totalEnergyMax*energyThreshold]
    all_names = []
    all_fluxes = [cases[c]['cone'] for c in filtered_cases]
    all_deltas = [cases[c]['delta'] for c in filtered_cases]
    all_tigns = [cases[c]['tign'] for c in filtered_cases]
    all_names = [('CONE_%03.2f_%03d'%(cases[c]['delta']*1e3, cases[c]['cone'])).replace('.','p') for c in filtered_cases]
    
    for i, c in enumerate(filtered_cases):
        time = cases[c]['times_trimmed']
        hrrpua = cases[c]['hrrs_trimmed']
        namespace = all_names[i]
        if cases[c]['totalEnergy'] < totalEnergyMax*energyThreshold: continue
        prevTime=-1e6
        for i in range(0, len(time)):
            if (time[i]-prevTime) < 0.0001:
                #txt = txt+"&RAMP ID='CONE-RAMP', T=%0.4f, F=%0.1f, /\n"%(time[i]-time[0]+0.0001, hrrpua[i])
                pass
            else:
                txt = txt+"&RAMP ID='%s', T=%0.4f, F=%0.1f, /\n"%(namespace, time[i]-time[0], hrrpua[i])
            prevTime = time[i]
        txt = txt+"&RAMP ID='%s', T=%0.4f, F=%0.1f, /\n"%(namespace, prevTime-time[0]+1, 0.0)
    
    ind = np.argsort(all_fluxes)
    all_names = [all_names[i] for i in ind]
    all_deltas = [all_deltas[i] for i in ind]
    all_tigns = [all_tigns[i] for i in ind]
    all_fluxes = [all_fluxes[i] for i in ind]
    
    ind = np.argsort(all_deltas)
    all_names = [all_names[i] for i in ind]
    all_fluxes = [all_fluxes[i] for i in ind]
    all_tigns = [all_tigns[i] for i in ind]
    all_deltas = [all_deltas[i] for i in ind]
    
    y = -0.05
    for i in range(0, len(all_names)):
        flux = all_fluxes[i]
        delta = all_deltas[i]
        namespace = all_names[i]
        
        filtered_flux = [c for c,n in zip(all_fluxes,all_names) if n != namespace]
        filtered_delta = [c for c,n in zip(all_deltas,all_names) if n != namespace]
        filtered_name = [c for c,n in zip(all_names,all_names) if n != namespace]
        
        if len(filtered_flux) == 0:
            filtered_flux = [c for c,n in zip(all_fluxes,all_names) ]
            filtered_delta = [c for c,n in zip(all_deltas,all_names) ]
            filtered_name = [c for c,n in zip(all_names,all_names) ]
        
        if i%3 == 0: y = y + 0.1
        XYZ = [((i % 3))*0.1+0.05, y, 0.0]
        XB = [XYZ[0]-0.05, XYZ[0]+0.05, XYZ[1]-0.05, XYZ[1]+0.05, 0.0,0.0]
        XB2 = [XYZ[0]-0.05+0.6, XYZ[0]+0.05+0.6, XYZ[1]-0.05, XYZ[1]+0.05, 0.0,0.0]
        
        # No ignition for surface temperature calculation
        txt = txt+"&SURF ID='SAMPLE-%s_noign', EXTERNAL_FLUX=%0.0f, "%(namespace, flux)
        txt = txt+"HEAT_TRANSFER_COEFFICIENT=%0.1f, HEAT_TRANSFER_COEFFICIENT_BACK=10., "%(10)
        txt = txt+"MATL_ID(1:2,1)='SAMPLE','BACKING', THICKNESS(1:2)=%0.2e,%0.2e, /\n"%(delta, 0.0254/2)
        
        # Ignition
        txt = txt+"&SURF ID='SAMPLE-%s', EXTERNAL_FLUX=%0.0f, "%(namespace, flux)
        txt = txt+"HEAT_TRANSFER_COEFFICIENT=%0.1f, HEAT_TRANSFER_COEFFICIENT_BACK=10., "%(front_h)
        txt = txt+"HRRPUA=1., IGNITION_TEMPERATURE=%0.0f, MATL_ID(1:2,1)='SAMPLE','BACKING', "%(Tign)
        txt = txt+"REFERENCE_HEAT_FLUX_TIME_INTERVAL=0.5, "
        txt = txt+"RAMP_Q="
        for j in range(0, len(filtered_name)):
            txt = txt + "'%s',"%(filtered_name[j])
        txt = txt+"REFERENCE_HEAT_FLUX="
        for j in range(0, len(filtered_name)):
            txt = txt + "%0.0f,"%(filtered_flux[j])
        txt = txt+"REFERENCE_THICKNESS="
        for j in range(0, len(filtered_name)):
            txt = txt + "%0.2e,"%(filtered_delta[j])
        txt = txt+'THICKNESS(1:2)=%0.2e,%0.2e, /\n'%(delta, 0.0254/2)
        
        # Add vent for before ignition
        txt = txt+"&VENT ID='SAMPLE-%s_noign', SURF_ID='SAMPLE-%s_noign', XB="%(namespace, namespace)
        for x in XB2:
            txt = txt+"%0.2f,"%(x)
        txt = txt+' /\n'
        
        # Add obst for after ignition
        txt = txt+"&OBST ID='SAMPLE-%s', SURF_ID='SAMPLE-%s', XB="%(namespace, namespace)
        for x in XB:
            txt = txt+"%0.2f,"%(x)
        if ignitionMode == 'Time':
            txt = txt+"DEVC_ID='TIGN-%s'"%(namespace)
        txt = txt+', /\n'
        
        txt = txt+"&DEVC ID='WALL TEMPERATURE-%s', INITIAL_STATE=.FALSE., IOR=3, OUTPUT=%s, "%(namespace, tempOutput)
        txt = txt+"QUANTITY='WALL TEMPERATURE', SETPOINT=%0.0f, XYZ=%0.2f,%0.2f,%0.2f, /\n"%(Tign, XYZ[0]+0.6, XYZ[1], XYZ[2])
        
        txt = txt+"&CTRL ID='IGNITION-CTRL-%s', FUNCTION_TYPE='ANY', INPUT_ID='WALL TEMPERATURE-%s', /\n"%(namespace, namespace)
        if ignitionMode == 'Time':
            txt = txt+"&DEVC ID='TIGN-%s', XYZ=0,0,0, SETPOINT=%0.2f, QUANTITY='TIME', INITIAL_STATE=.FALSE., /\n"%(namespace, all_tigns[i])
            
        txt = txt+"&DEVC ID='IGNITION_DEVC-%s', CTRL_ID='IGNITION-CTRL-%s', IOR=3, OUTPUT=.FALSE., QUANTITY='CONTROL', "%(namespace,namespace)
        txt = txt+"XYZ=%0.2f,%0.2f,%0.2f, /\n"%(XYZ[0], XYZ[1], XYZ[2])
        
        txt = txt+"&DEVC ID='HRRPUA-%s', IOR=3, QUANTITY='HRRPUA', SPEC_ID='PROPANE', "%(namespace)
        txt = txt+"XYZ=%0.2f,%0.2f,%0.2f, /\n"%(XYZ[0], XYZ[1], XYZ[2])
        
        txt = txt+"&DEVC ID='IGNITION-TIME-%s', NO_UPDATE_DEVC_ID='IGNITION_DEVC-%s', OUTPUT=.FALSE., "%(namespace,namespace)
        txt = txt+"QUANTITY='TIME', XYZ=%0.2f,%0.2f,%0.2f, /\n"%(XYZ[0], XYZ[1], XYZ[2])
                        
    return txt

def loadBclData(data_dir, plot=False):
    allData = dict()
    for thickness in [6, 10]:
        coneData = dict()
        for configuration in ['Horizontal','Vertical']:
            coneData[configuration] = dict()
            material = "Beskid_PMMA_CN_%s_%dmm.csv"%(configuration, thickness)
            d = pd.read_csv(data_dir + os.sep + material, sep=';')
            for flux in [25, 50, 75]:
                tmax = 0
                times = []
                hrrs = []
                tigns = []
                dt = 0.1
                for counter in [1, 2, 3]:
                    time = d['TIME_%d_%d'%(flux,counter)].values
                    HRR = d['HRR_%d_%d'%(flux,counter)].values
                    energyCutoff1 = 0.0001
                    energyCutoff2 = 0.99
                    ind = np.where(HRR > 30)[0][0]
                    ind = max(0, ind-2)
                    HRR[:ind] = 0
                    tmax = max(np.nanmax(time),tmax)
                    time2 = np.linspace(0, np.nanmax(time), int(np.nanmax(time)/0.1) + 1)
                    HRR2 = np.interp(time2, time, HRR)
                    tign, times_trimmed, hrrs_trimmed = findLimits(time2, HRR2, energyCutoff1, energyCutoff2)
                    times.append(times_trimmed)
                    hrrs.append(hrrs_trimmed)
                    tigns.append(tign)
                times_out = np.linspace(0, tmax, int(tmax/0.1) + 1)
                hrrs_out = np.zeros_like(times_out)
                tign_out = np.nanmean(tigns)
                for counter in [1, 2, 3]:
                    time2 = times[counter-1]-tigns[counter-1] + tign_out
                    tmp_hrr = np.interp(times_out, time2, hrrs[counter-1], left=0,right=0)
                    hrrs_out = hrrs_out + tmp_hrr
                    #plt.plot(time2, hrrs[counter-1])
                hrrs_out = hrrs_out/3
                hrrs_out[times_out < tign_out] = 0
                
                totalEnergy = np.trapz(hrrs_out, times_out)
                npoints = 10
                times_out2 = np.linspace(0, times_out[-1], npoints)
                times_out2 = np.append(times_out2, tign_out)
                times_out2 = np.append(times_out2, times_out[np.argmax(hrrs_out)])
                times_out2 = np.sort(times_out2)
                hrrs_out2 = np.interp(times_out2, times_out, hrrs_out)
                totalEnergy2 = np.trapz(hrrs_out2, times_out2)
                while abs(totalEnergy2-totalEnergy)/(totalEnergy+totalEnergy2)/2 > 0.001:
                    hrrs_out3 = np.interp(times_out, times_out2, hrrs_out2)
                    diff = abs(hrrs_out3-hrrs_out)
                    tind = np.argmax(diff)
                    times_out2 = np.append(times_out2,times_out[tind])
                    times_out2 = np.sort(times_out2)
                    hrrs_out2 = np.interp(times_out2, times_out, hrrs_out)
                    totalEnergy2 = np.trapz(hrrs_out2, times_out2)
                    
                
                coneData[configuration][flux] = {'time': times_out2, 'hrrs': hrrs_out2, 'tign': tign_out}
        if plot:
            plt.figure(figsize=(8,8))
            colors = getColors()
            ls = ['-','--']
            lw = 3
            fs = 12
            for j, configuration in enumerate(list(coneData.keys())):
                for i, flux in enumerate(list(coneData[configuration].keys())):
                    label = r'$\Delta$=%dmm flux=%0.0f kW/m$^{2}$'%(thickness, flux)
                    label = r'%s flux=%0.0f kW/m$^{2}$'%(configuration, flux)
                    plt.plot(coneData[configuration][flux]['time'], coneData[configuration][flux]['hrrs'], ls[j], label=label, color=colors[i], linewidth=lw)
            plt.legend(fontsize=fs, framealpha=1)
            plt.ylim(0, 1200)
            plt.xlim(0, 1500)
            plt.xlabel('Time (sec)', fontsize=fs)
            plt.ylabel('HRRPUA (kW/m$^{2}$)', fontsize=fs)
            plt.grid()
        allData[thickness] = coneData
    return allData

def runModel(outdir, outfile, mpiProcesses, fdsdir, fdscmd, printLiveOutput=False):
    ''' This function will run fds with an input file
    '''
    my_env = os.environ.copy()
    my_env['I_MPI_ROOT'] = fdsdir+"\\mpi"
    my_env['PATH'] = fdsdir + ';' + my_env['I_MPI_ROOT'] + ';' + my_env["PATH"]
    my_env['OMP_NUM_THREADS'] = '1'
    
    process = subprocess.Popen([fdscmd, outfile, ">&", "log.err"], cwd=r'%s'%(outdir), env=my_env, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    out, err = process.communicate()
    errcode = process.returncode   
    return out, err, errcode

if __name__ == "__main__":
    
    data_dir = "..\\data\\bcl_materials"
    allData = loadBclData(data_dir, False)
    properties = {'density': 1182, 'conductivity': 0.15, 'specific_heat': 1.48,
                  'emissivity': 1.0, 'heat_of_combustion': 29.5, 'soot_yield': 0.006}
    
    # Prepare directories
    os.environ['PATH'] = os.environ['PATH'] + ";C:\\Program Files\\firemodels\\FDS6_9_1\\FDS6\\bin"
    fdscmd = 'fds.exe'
    
    fdsdir, check = findFds()
    workingDir = 'working'
    runSimulations = True
    
    configuration_to_scale = 'Horizontal'
    configuration_to_compare = 'Vertical'
    
    chid = 'bcl_black_pmma'
    if configuration_to_scale == 'Horizontal': chid = chid + '_horizontal'
    if configuration_to_scale == 'Vertical': chid = chid + '_vertical'
    if configuration_to_compare == 'Horizontal': chid = chid + '_horizontal'
    if configuration_to_compare == 'Vertical': chid = chid + '_vertical'
    
    Tign = 436
    front_h = 10
    cases = dict()
    counter = 0
    for thickness in [6, 10]:
        for configuration in ['Horizontal']:
            for flux in [25, 50, 75]:
                tign = allData[thickness][configuration][flux]['tign']
                times_trimmed = allData[thickness][configuration][flux]['time']-tign
                hrrs_trimmed = allData[thickness][configuration][flux]['hrrs']
                inds = np.where(times_trimmed > -0.001)
                times_trimmed = times_trimmed[inds]
                hrrs_trimmed = hrrs_trimmed[inds]
                total_energy = np.trapz(hrrs_trimmed, times_trimmed)
                cases[counter] = {
                    'times_trimmed': times_trimmed,
                    'hrrs_trimmed': hrrs_trimmed,
                    'tign': tign,
                    'totalEnergy': total_energy,
                    'cone': flux,
                    'delta': thickness/1000}
                counter += 1

    txt = buildFdsFile(chid, cases, properties, Tign, front_h,
                     ignitionMode='Time', outputTemperature=False,
                     calculateDevcDt=True, devc_dt=10.,
                     energyThreshold=0.0)
    
    try: 
        os.mkdir(os.path.abspath(workingDir))
    except:
        pass
    with open("%s%s%s.fds"%(workingDir, os.sep, chid), 'w') as f:
        f.write(txt)

    if runSimulations:
        runModel(workingDir, chid+".fds", 1, fdsdir, fdscmd, printLiveOutput=False)
    
    data = load_csv(workingDir, chid)
    output_statistics = dict()
    plt.figure(figsize=(8,8))
    lw = 3
    fs = 12
    j=0
    ls = ['-','--']
    colors = getColors()
    cases_to_plot = list(cases.keys())
    for i in range(0, len(cases_to_plot)):
        c = cases_to_plot[i]
        thickness = cases[c]['delta']*1e3
        flux = cases[c]['cone']
        namespace = ('CONE_%03.2f_%03d'%(cases[c]['delta']*1e3, cases[c]['cone'])).replace('.','p')
        
        label = r'%0.0f $\mathrm{kW/m^{2}}$'%(cases[c]['cone'])
        delta0 = cases[c]['delta']
        times = data['Time'].values
        hrrpuas = data['"HRRPUA-'+namespace+'"'].values
        
        label = r'%s $\Delta$=%0dmm flux=%0.0f kW/m$^{2}$'%(configuration, thickness, flux)
        plt.plot(times, hrrpuas, ls[j], label=label, color=colors[i], linewidth=lw)
        
        t_raw = allData[thickness][configuration_to_compare][flux]['time']
        hrr_raw = allData[thickness][configuration_to_compare][flux]['hrrs']
        plt.plot(t_raw, hrr_raw, ls[1], color=colors[i], linewidth=lw)
    plt.legend(fontsize=fs, framealpha=1)
    plt.ylim(0, 1200)
    plt.xlim(0, 1500)
    plt.xlabel('Time (sec)', fontsize=fs)
    plt.ylabel('HRRPUA (kW/m$^{2}$)', fontsize=fs)
    plt.grid()
    plt.savefig(chid+'_comparison.png',dpi=300)