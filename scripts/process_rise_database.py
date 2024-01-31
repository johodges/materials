import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict

def extractConeAnalysisData(values, scaled_times, tigns, uniqueFluxes):
    coneAnalysisData = dict()
    for i, flux in enumerate(uniqueFluxes):
        
        v = values[i]
        t = scaled_times[i]
        tign = tigns[i]
        coneAnalysisData[flux] = dict()
        coneAnalysisData[flux]['peakHRRPUA'] = np.max(v)
        coneAnalysisData[flux]['timeToPeak'] = t[np.argmax(v)]-tign
        
        tMax = int(np.ceil(t.max()))
        tMax = max([tMax, 60])
        t2 = np.linspace(0, tMax, tMax+1)
        v2 = np.interp(t2, t, v)
        
        dt = t2[1] - t2[0]
        
        filters = [int(60/dt), int(180/dt), int(300/dt)]
        filteredValues = []
        for f in filters:
            fil = np.ones(f)/float(f)
            vfilt = np.convolve(v2, fil, mode='same')
            filteredValues.append(np.max(vfilt))
        
        coneAnalysisData[flux]['avg60s'] = filteredValues[0]
        coneAnalysisData[flux]['avg180s'] = filteredValues[1]
        coneAnalysisData[flux]['avg300s'] = filteredValues[2]
    return coneAnalysisData

def getConeData(material, directory):
    p = os.path.abspath(directory)
    files = glob.glob(p+os.sep+'cone_H*')
    files = [x for x in files if "Scalar" not in x]
    
    cad_prop = pd.read_csv(p+os.sep+material+'_Cone_Analysis_Data.csv', index_col=0)
    tmp_prop = pd.read_csv(p+os.sep+material+'_Ignition_Temp_Properties.csv', header=None, index_col=0).T
    density = tmp_prop['Density (kg/m3)'].values[0]
    coneData = dict()
    for j, file in enumerate(files):
        dc = pd.read_csv(file)
        rev = file.split('_')[-1].split('.csv')[0]
        HF = file.split('HF')[1].split('Scalar')[0].split('_')[0]
        f = material + '_Cone_HF' + HF + 'Scalar_*_' + rev + '.csv'
        d = glob.glob(os.path.join(directory, f))[0]
        ds = pd.read_csv(d, header=None, index_col=0)
        
        namespace = file.split('cone_')[-1].split('.csv')[0]
        
        coneData[j] = dict()
        coneData[j]['flux'] = float(ds.loc['HEAT FLUX', 1]) # kW/m2
        coneData[j]['mass'] = float(ds.loc['SPECIMEN MASS', 1]) / 1000 # kg
        coneData[j]['area'] = float(ds.loc['SURF AREA', 1]) # m2
        coneData[j]['thickness'] = coneData[j]['mass'] / (density*coneData[j]['area'])
        coneData[j]['time'] = dc['Time']
        coneData[j]['hrrpua'] = dc['HRRPUA']
        coneData[j]['timeMax'] = np.nanmax(dc['Time'])
        coneData[j]['timeInt'] = float(ds.loc['SCAN TIME', 1])
        coneData[j]['timeIgn'] = float(ds.loc['TIME TO IGN', 1]) # m2
        coneData[j]['peakHRRPUA'] = cad_prop.loc['Peak HRRPUA (kW/m2)', namespace]
        coneData[j]['timeToPeak'] = cad_prop.loc['Time to Peak HRRPUA (s)', namespace]
        coneData[j]['HRRPUA, 60s average'] = cad_prop.loc['Average HRRPUA over 60 seconds (kW/m2)', namespace]
        coneData[j]['HRRPUA, 180s average'] = cad_prop.loc['Average HRRPUA over 180 seconds (kW/m2)', namespace]
        coneData[j]['HRRPUA, 300s average'] = cad_prop.loc['Average HRRPUA over 300 seconds (kW/m2)', namespace]
        coneData[j]['HeatOfCombustion'] = cad_prop.loc['Avg. Effective Heat of Combustion (MJ/kg)', namespace]
    
    coneData = pd.DataFrame(coneData).T
    
    return coneData

def timeAverageFsriConeData(coneData, flux, outInt):
    times = coneData.loc[coneData['flux'] == flux, 'time'].values
    hrrpuas = coneData.loc[coneData['flux'] == flux, 'hrrpua'].values
    
    tMax = np.max(coneData.loc[coneData['flux'] == flux, 'timeMax'].values)
    dt = np.min(coneData.loc[coneData['flux'] == flux, 'timeInt'].values)
    

def importFsriMaterial(p, material, outInt, uniqueFluxes, filterWidth=11):
    material_dict = dict()
    tmp_prop = pd.read_csv(p+os.sep+material+'_Ignition_Temp_Properties.csv', header=None, index_col=0).T
    
    conductivity = tmp_prop['Thermal Conductivity (W/m-K)'].values[0]
    heatCapacity = tmp_prop['Heat Capacity (J/kg-K)'].values[0]
    density = tmp_prop['Density (kg/m3)'].values[0]
    
    material_dict['conductivity'] = conductivity
    material_dict['heatCapacity'] = heatCapacity
    material_dict['density'] = density
    material_dict['directory'] = p
    coneData = getConeData(material, p)
    
    #uniqueFluxes = np.unique(coneData['flux'].values)
    fil = np.ones(filterWidth)/filterWidth
    
    for flux in uniqueFluxes:
        outInt2 = outInt
        dt = np.min(coneData.loc[coneData['flux'] == flux, 'timeInt'].values)
        tMax = np.max(coneData.loc[coneData['flux'] == flux, 'timeMax'].values)
        
        thickness = np.median(coneData.loc[coneData['flux'] == flux, 'thickness'].values)
        mass = np.median(coneData.loc[coneData['flux'] == flux, 'mass'].values)
        area = np.median(coneData.loc[coneData['flux'] == flux, 'area'].values)
        tign = np.median(coneData.loc[coneData['flux'] == flux, 'timeIgn'].values)
        peakHRRPUA = np.median(coneData.loc[coneData['flux'] == flux, 'peakHRRPUA'].values)
        timeToPeak = np.median(coneData.loc[coneData['flux'] == flux, 'timeToPeak'].values)
        avg60s = np.nanmedian(coneData.loc[coneData['flux'] == flux, 'HRRPUA, 60s average'].values)
        avg180s = np.nanmedian(coneData.loc[coneData['flux'] == flux, 'HRRPUA, 180s average'].values)
        avg300s = np.nanmedian(coneData.loc[coneData['flux'] == flux, 'HRRPUA, 300s average'].values)
        avgHoC = np.nanmedian(coneData.loc[coneData['flux'] == flux, 'HeatOfCombustion'].values)
        
        times = coneData.loc[coneData['flux'] == flux, 'time'].values
        hrrpuas = coneData.loc[coneData['flux'] == flux, 'hrrpua'].values
        
        # Interpolate and filter raw data
        time = np.linspace(0, tMax, int(tMax/dt)+1)
        hrrpua = np.zeros_like(time)
        hrrpuas_interp = np.zeros((time.shape[0], len(times)))
        
        filterWidth=41
        fil = np.ones(filterWidth)/filterWidth
        
        for j in range(0, len(times)):
            hrrpuas_interp[:, j] = np.convolve(np.interp(time, times[j], hrrpuas[j]), fil, mode='same')
        
        # Find time to ignitions
        tIgns = []
        for j in range(0, len(times)):
            totalEnergy = np.trapz(hrrpuas_interp[:, j], time)
            totalEnergy2 = 0
            ind2 = 1
            while totalEnergy2 < 0.0001 * totalEnergy:
                ind2 = ind2 + 1
                totalEnergy2 = np.trapz(hrrpuas_interp[:ind2, j], time[:ind2])
            hrrpuas_interp[:, j] = np.interp(time, time, hrrpuas_interp[:, j])
            tIgns.append(time[ind2-1])
        tign = np.median(tIgns)
        
        # Average neglecting time to ignition
        hrrpuas_interp_notign = np.zeros_like(hrrpuas_interp)
        for j in range(0, len(times)):
            hrrpuas_interp_notign[:, j] = np.interp(time, time-tIgns[j], hrrpuas_interp[:, j])
        
        hrrpua = np.mean(hrrpuas_interp_notign, axis=1)
        
        #_, times_trim, hrrpuas_trim = findLimits(time+tign,hrrpua, 0.001, 0.99)
        
        #times_trim = np.append(np.array([0, tign*0.9]), times_trim)
        #hrrpuas_trim = np.append(np.array([0, 0]), hrrpuas_trim)
        
        #tMax = times_trim.max()+tign
        #targetTimes, HRRs_interp = interpolateExperimentalData(times, HR
        
        
        #hrrpua_mn = np.min(hrrpuas_interp_notign, axis=1)
        #hrrpua_mx = np.max(hrrpuas_interp_notign, axis=1)
        hrrpua_std = np.std(hrrpuas_interp_notign, axis=1)
        
        tMax = np.ceil((tMax-tign)/outInt2)*outInt2
        outTime = np.linspace(0, tMax, int(tMax/outInt2)+1)
        outHrrpua = np.interp(outTime, time, hrrpua)
        
        totalEnergy = np.trapz(hrrpua, time)
        totalEnergy2 = np.trapz(outHrrpua, outTime)
        
        peakHrrpua2 = np.nanmax(outHrrpua)
        peakHrrpua = np.nanmax(hrrpua)
        
        while abs(totalEnergy2 - totalEnergy)/totalEnergy > 0.00001 or abs(peakHrrpua2 - peakHrrpua)/peakHrrpua > 0.05:
            print(abs(totalEnergy2 - totalEnergy)/totalEnergy)
            outInt2 = outInt2*0.9
            outTime = np.linspace(0, tMax, int(tMax/outInt2)+1)
            outHrrpua = np.interp(outTime, time, hrrpua)
            totalEnergy2 = np.trapz(outHrrpua, outTime)
            peakHrrpua2 = np.nanmax(outHrrpua)
            if outInt2 <= 0.1:
                totalEnergy2 = totalEnergy
        
        outTime = np.append(np.array([0, tign*0.9]), outTime+tign)
        outHrrpua = np.append(np.array([0, 0]), outHrrpua)
        
        material_dict[flux] = dict()
        material_dict[flux]['time'] = outTime
        material_dict[flux]['hrrpua'] = outHrrpua
        material_dict[flux]['hrrpua_full_mean'] = hrrpua
        material_dict[flux]['hrrpua_full_std'] = hrrpua_std
        #material_dict[flux]['hrrpua_full_min'] = hrrpua_mn
        #material_dict[flux]['hrrpua_full_max'] = hrrpua_mx
        material_dict[flux]['time_full'] = time #-tign #-tStart
        material_dict[flux]['tIgn'] = tign
        material_dict[flux]['peakHRRPUA'] = peakHRRPUA
        material_dict[flux]['timeToPeak'] = timeToPeak
        material_dict[flux]['avg60s'] = avg60s
        material_dict[flux]['avg180s'] = avg180s
        material_dict[flux]['avg300s'] = avg300s
        material_dict[flux]['thickness'] = thickness
        material_dict[flux]['HeatOfCombustion'] = avgHoC
        material_dict[flux]['hrrpuas_interp'] = hrrpuas_interp
        material_dict[flux]['hrrpuas_interp_notign'] = hrrpuas_interp_notign
        material_dict[flux]['time_interp'] = time
        
    return material_dict

def checkMaterial(p, material, ignores):
    files = glob.glob(p+os.sep+'cone_H*')
    complete = True
    #if os.path.exists(p+os.sep+'ignition_temp.csv') is False: complete = False
    if os.path.exists(p+os.sep+material+'_Ignition_Temp_Properties.csv') is False: complete = False
    if len(files) == 0: complete = False
    if material in ignores: complete = False
    return complete

def importFsriDatabase(data_dir, outInt, Tinfty=300, ignores=['']):
    material_directories = glob.glob(os.path.join(os.path.abspath(data_dir),'*'))
    possible_materials = [d.split(os.sep)[-1] for d in material_directories]
    
    complete_materials = dict()
    materials = []
    for i in range(0, len(material_directories)):
        p = os.path.abspath(material_directories[i])
        check = checkMaterial(p, possible_materials[i], ignores)
        if check: materials.append(possible_materials[i])
    
    uniqueFluxes = [25, 50, 75]
    for i in range(0, len(materials)):
        material = materials[i]
        p = os.path.join(os.path.abspath(data_dir), material)
        material_dict = importFsriMaterial(p, material, outInt, uniqueFluxes)
        complete_materials[material] = material_dict
        
    return complete_materials

def interpolateExperimentalData(times, HRRs, targetDt=False, filterWidth=False):
    dt = np.nanmedian(times[1:]-times[:-1])
    if filterWidth is not False:
        filterWidth = int(filterWidth/dt)
        fil = np.ones(filterWidth)/filterWidth
        HRRs = np.convolve(HRRs, fil, mode='same')
    
    if targetDt is not False:
        dt = targetDt
    else:
        dt = np.nanmedian(times[1:]-times[:-1])
    tmax = np.round(times.max()/dt)*dt
    tmin = np.round(times.min()/dt)*dt
    targetTimes = np.linspace(tmin, tmax, int((tmax-tmin)/dt + 1))
    HRRs = np.interp(targetTimes, times, HRRs)
    
    return targetTimes, HRRs

def findLimits(times, HRRs, energyCutoff1, energyCutoff2):
    v = np.cumsum(HRRs)
    ind1 = 0 
    while ind1 == 0:
        try:
            ind1 = np.where(v < np.nanmax(v)*energyCutoff1)[0][-1]
        except:
            energyCutoff1 = energyCutoff1*2
            print(energyCutoff1)
    ind2 = v.shape[0]
    while ind2 == v.shape[0]:
        try:
            ind2 = np.where(v > np.nanmax(v)*energyCutoff2)[0][0]
        except:
            energyCutoff2 = energyCutoff2*0.99
            print(energyCutoff2)
    times_trimmed = times[ind1:ind2]
    hrrs_trimmed = HRRs[ind1:ind2]
    tign = times[ind1]
    return tign, times_trimmed, hrrs_trimmed

if __name__ == "__main__":
    
    systemPath = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(systemPath,'..','data','rise_materials')+os.sep
    out_dir = os.path.join(systemPath,'..','data','rise_materials_processed')+os.sep
    
    try:
        os.mkdir(out_dir)
    except:
        pass
    
    ignores = ['-29mm',
               'Carpet_Glue_Aluminum_plate_2_mm-05mm',
               'FR_EPS_Calcium_silicate_board-25mm',
               'FR_Particle_board-12mm',
               'Mineral_wool_faced_Combustible_face-30mm',
               'Plastic_faced_steel_sheet_Mineral_wool-25mm',
               'Plywood-12mm',
               'polyester_fabric-00mm',
               'RPPVC_EPR-16mm',
               'Solid_acrylic-12mm',
               ]
    
    material_database = defaultdict(bool)
    
    files = glob.glob(data_dir + "*.txt")
    
    for file in files:
        if True: #try:
            with open(file, 'r') as f:
                txt = f.read()
            
            tmp = txt.split('Keywords\n')[1].split('\nScalars')[0].split('\n')
            #txt[4].replace('\n','').split(';')
            
            #columns = txt[1].replace('\n','').split(';')
            #values = txt[2].replace('\n','').split(';')
            
            columns = tmp[0].split(';')
            values = tmp[1].split(';')
            
            material = '|'.join([x for x in values[:4] if x != ''])
            series = values[8]
    
            tmp = txt.split('Scalars\n')[1].split('VectorData\n')[0]
            
            tmp = tmp.replace('\n;',';').replace(';\n',';').split('\n')
            
            scalarColumns = tmp[0].split(';')
            tmp2 = '\n'.join(tmp[1:])
            scalarValues = tmp2.split(';')
    
            #scalarColumns = txt[4].replace('\n','').split(';')
            #scalarValues = txt[5].replace('\n','').split(';')
            
            scalars = dict()
            for c, v in zip(scalarColumns, scalarValues):
                try:
                    scalars[c] = float(v)
                except:
                    scalars[c] = v
            
            if 'Thickness (mm)' in scalarColumns:
                thickness = scalars['Thickness (mm)']
            elif 'Product thickness (mm)' in scalarColumns:
                thickness = scalars['Product thickness (mm)']
            else:
                print("Warning no thickness found for %s"%(file))
            
            if thickness == '':
                print("Warning no thickness found for %s"%(file))
                thickness = 0.0127
            elif type(thickness) == str:
                if '-' in thickness:
                    thickness = np.mean([float(x) for x in thickness.split('-')])
            
            HF = scalars['Flux (kW/m2)']
            tig = scalars['tig (s)']
            qpeak = scalars['qmax (kw/m2)']
            q180 = scalars['q180 (kw/m2)']
            q300 = scalars['q300 (kw/m2)']
            dHc = scalars['DHc (MJ/kg)']
            
            if 'Density (kg/m3)' in scalarColumns:
                rho = scalars['Density (kg/m3)']
            elif 'Product density (kg/m3)' in scalarColumns:
                rho = scalars['Product density (kg/m3)']
            
            if rho == '':
                if 'Sample mass (g)' in scalarColumns:
                    mass = scalars['Sample mass (g)']
                else:
                    print("Warning unable to find sample mass for %s"%(file))
                    if 'Mass loss (g)' in scalarColumns:
                        mass = scalars['Mass loss (g)']
                if 'Area (m2)' in scalarColumns:
                    area = scalars['Area (m2)']
                elif 'Area exposed (m2)' in scalarColumns:
                    area = scalars['Area exposed (m2)']
                try:
                    rho = mass / (area * thickness)
                except:
                    pass
                
                #assert False, "Stopped"
            #thickness_columns = [t for t in list(scalars.keys()) if 'Thickness' in t]
            #thicknesses = 
            
            #[name,idn,id2] = np.genfromtxt(file,delimiter=';',skip_header=2,max_rows=1,usecols=[0,4,9],dtype=str)
            #[HF,tig,qpeak,THR,q180,q300,dm,dHc,SPpeak,TSP,SEA,Area,t,rho] = np.genfromtxt(file,delimiter=';',skip_header=5,max_rows=1)
            
            material_dict = dict()
            material_dict['flux'] = HF
            material_dict['tIgn'] = tig
            material_dict['peakHRRPUA'] = qpeak
            material_dict['avg180s'] = q180
            material_dict['avg300s'] = q300
            
            material_dict['density'] = rho
            material_dict['thickness'] = thickness
            material_dict['HeatOfCombustion'] = dHc
            material_dict['FYI'] = series
            
            tmp = txt.split('VectorData\n')[1].split('\n')
            while tmp[-1] == '': tmp = tmp[:-1]
            '''
            for t in tmp[1:]:
                print(t)
                tt = t.split(';')
                v = [float(x) for x in tt]
            '''
            tmp2 = [[float(x) for x in y.split(';')] for y in tmp[1:] if (('Time' not in y) and (len(y.replace(';','')) > 0) and ('s' not in y)) ]
            hrrdata = np.array(tmp2)
            
            #hrrdata = np.genfromtxt(file,delimiter=';',skip_header=8,usecols=[0,1])
            
            times = hrrdata[:, 0]
            hrrpua = hrrdata[:, 1]
            
            material_dict['raw_times'] = times
            material_dict['raw_hrrpua'] = hrrpua
            
            if type(tig) == float:
                hrrpua[times < tig] = 0
            
            matind = material # + "_" + series
            if material_database[matind] is False:
                material_database[matind] = defaultdict(bool)
            
            if material_database[matind][thickness] is False:
                material_database[matind][thickness] = defaultdict(bool)
            
            if material_database[matind][thickness][HF] is False:
                material_database[matind][thickness][HF] = defaultdict(bool)
                material_database[matind][thickness][HF]['testCount'] = 0
            
            testCount = material_database[matind][thickness][HF]['testCount'] + 1
            material_database[matind][thickness][HF]['testCount'] = testCount
            
            material_database[matind][thickness][HF][testCount] = material_dict
        
        #except:
        #    print(file)
    
    material_database_filtered = defaultdict(bool)
    materials = list(material_database.keys())
    
    for material in materials:
        thicknesses = list(material_database[material].keys())
        for thickness in thicknesses:
            if type(thickness) is str: continue
            if thickness <= 0: continue
            fluxes = list(material_database[material][thickness].keys())
            if len(fluxes) < 2: continue
            densities = []
            for flux in fluxes:
                tigns = []
                testCount = material_database[material][thickness][flux]['testCount']
                for i in range(1, testCount+1):
                    density = material_database[material][thickness][flux][i]['density']
                    if type(density) == str:
                        pass
                    else:
                        densities.append(density)
                    tign = material_database[material][thickness][flux][i]['tIgn']
                    if type(tign) == str:
                        pass
                    else:
                        tigns.append(tign)
                if len(densities) == 0: continue
                material_database[material][thickness][flux]['tign'] = np.mean(tign)
            if material_database_filtered[material] is False:
                material_database_filtered[material] = defaultdict(bool)
                material_database_filtered[material]['Density'] = []
            material_database_filtered[material][thickness] = material_database[material][thickness]
            material_database_filtered[material]['Density'].append(np.mean(densities))
            #print(material, len(fluxes), fluxes)
            #print(densities)
    
    materials = list(material_database_filtered.keys())
    
    for material in materials:
        thicknesses = list(material_database_filtered[material].keys())
        thicknesses.remove('Density')
        for thickness in thicknesses:
            fluxes = list(material_database_filtered[material][thickness].keys())
            for flux in fluxes:
                testCount = material_database_filtered[material][thickness][flux]['testCount']
                tign = material_database_filtered[material][thickness][flux]['tign']
                outInt2 = 15
                if testCount > 1:
                    tMax = 0
                    for i in range(1, testCount+1):
                        raw_times = material_database_filtered[material][thickness][flux][i]['raw_times']
                        raw_hrrpuas = material_database_filtered[material][thickness][flux][i]['raw_hrrpua']
                        tign = material_database_filtered[material][thickness][flux][i]['tIgn']
                        tMax = max([tMax, np.nanmax(raw_times)-tign])
                    t_interp = np.linspace(0, np.round(tMax), int(np.round(tMax)*10+1))
                    hrrpua_interp = np.zeros_like(t_interp)
                    for i in range(1, testCount+1):
                        raw_times = material_database_filtered[material][thickness][flux][i]['raw_times']
                        raw_hrrpuas = material_database_filtered[material][thickness][flux][i]['raw_hrrpua']
                        tign = material_database_filtered[material][thickness][flux][i]['tIgn']
                        hrrpua_interp += np.interp(t_interp, raw_times-tign, raw_hrrpuas)
                    hrrpua_interp = hrrpua_interp / testCount
                    tign = material_database_filtered[material][thickness][flux]['tign']
                    times = np.append(np.array([0, tign-1]), t_interp+tign)
                    hrrpuas = np.append(np.array([0, 0]), hrrpua_interp)

                else:
                    tign = material_database_filtered[material][thickness][flux][1]['tIgn']
                    raw_times = material_database_filtered[material][thickness][flux][1]['raw_times']
                    raw_hrrpuas = material_database_filtered[material][thickness][flux][1]['raw_hrrpua']
                    tMax = np.nanmax(raw_times)-tign
                    t_interp = np.linspace(0, np.round(tMax), int(np.round(tMax)*10+1))
                    hrrpua_interp = np.interp(t_interp, raw_times-tign, raw_hrrpuas)
                    times = np.append(np.array([0, tign-1]), t_interp+tign)
                    hrrpuas = np.append(np.array([0, 0]), hrrpua_interp)
                
                hrrpuas[hrrpuas < 0] = 0
                
                hrrpuas[np.isnan(hrrpuas)] = 0.0
                
                outTime = np.linspace(0, tMax+tign, int(np.ceil((tMax+tign)/outInt2)+1))
                outHrrpua = np.interp(outTime, times, hrrpuas)
                
                totalEnergy = np.trapz(hrrpuas, times)
                totalEnergy2 = np.trapz(outHrrpua, outTime)
                
                peakHrrpua2 = np.nanmax(outHrrpua)
                peakHrrpua = np.nanmax(hrrpuas)
                
                while abs(totalEnergy2 - totalEnergy)/totalEnergy > 0.001 or abs(peakHrrpua2 - peakHrrpua)/peakHrrpua > 0.01:
                    outInt2 = outInt2*0.9
                    outTime = np.linspace(0, tMax+tign, int(np.round((tMax+tign)/outInt2)+1))
                    outHrrpua = np.interp(outTime, times, hrrpuas)
                    totalEnergy2 = np.trapz(outHrrpua, outTime)
                    peakHrrpua2 = np.nanmax(outHrrpua)
                    #print(outInt2, abs(totalEnergy2 - totalEnergy)/totalEnergy, abs(peakHrrpua2 - peakHrrpua)/peakHrrpua)
                    if outInt2 <= 0.1:
                        totalEnergy2 = totalEnergy
                        peakHrrpua2 = peakHrrpua
                
                    outHrrpua[outHrrpua < 0] = 0
                material_database_filtered[material][thickness][flux]['time'] = outTime
                material_database_filtered[material][thickness][flux]['hrrpua'] = outHrrpua
                
                d = pd.DataFrame(np.round(np.array([outTime, outHrrpua]).T, decimals=1), columns=['Time','HRRPUA'])
                dataFile = os.path.join(out_dir, '%s-%02dmm-%02d.csv'%(material.replace('|','_').replace(' ','_').replace(',','_'), thickness, flux))
                d.to_csv(dataFile, index=False)
    
    material_database = material_database_filtered
    resultDir = "../../../out/Scaling_Pyrolysis/"
    inputFileDir = "../../../fds/Validation/Scaling_Pyrolysis/"
    expFileDir = "../../../exp/RISE_Materials/"
    emissivity = 1
    txt = 'Code,Number,Series,Material,DataFile,ResultDir,InputFileDir,ExpFileDir,'
    txt = txt + 'ReferenceExposure,ReferenceTime,ReferenceHRRPUA,'
    txt = txt + 'ValidationTimes,ValidationHrrpuaColumns,ValidationFluxes,'
    txt = txt + 'Density,Conductivity,SpecificHeat,Emissivity,Thickness,'
    txt = txt + 'IgnitionTemperature,IgnitionTemperatureBasis,HeaderRows,FYI'
    
    for material in materials:
        thicknesses = list(material_database_filtered[material].keys())
        thicknesses.remove('Density')
        for thickness in thicknesses:
            conductivity = 0.4 #material_database[material]['conductivity']
            specific_heat = 1. #material_database[material]['heatCapacity']
            density = np.mean(material_database[material]['Density'])
            
            fluxes = [x for x in list(material_database[material][thickness].keys()) if type(x) is float]
            fluxes.sort()
            
            code ='d'
            number = 1
            mat = '%s-%02dmm'%(material, thickness)
            mat = mat.replace(',','_').replace(' ','_').replace('|','_')
            
            if mat in ignores:
                break
            
            dataFiles = ''
            for flux in fluxes:
                dataFile = os.path.join(expFileDir, '%s-%02d.csv'%(mat, flux))
                dataFiles = dataFiles + dataFile + '|'
            dataFiles = dataFiles[:-1]
            
            if (50 in fluxes):
                refFlux = 50
            else:
                ind = np.argmin([abs(x-50) for x in fluxes])
                refFlux = fluxes[ind]
            
            txt = txt + "\n" + "%s,%s,%s,%s,%s,%s,"%(code, number, "RISE_Materials", mat, dataFiles, resultDir)
            txt = txt + "%s,%s,%0.0f,%s-%0.0f.csv-Time,%s-%0.0f.csv-HRRPUA,"%(inputFileDir, expFileDir, refFlux, mat, refFlux, mat, refFlux)
            
            for flux in fluxes:
                txt = txt + '%s-%0.0f.csv-Time|'%(mat, flux)
            txt = txt[:-1] + ','
            for flux in fluxes:
                txt = txt + '%s-%0.0f.csv-HRRPUA|'%(mat, flux)
            txt = txt[:-1] + ','
            for flux in fluxes:
                txt = txt + '%0.0f|'%(flux)
            txt = txt[:-1] + ','
            txt = txt + '%0.1f,%0.4f,%0.4f,%0.4f,%0.8f,'%(density, conductivity, specific_heat, emissivity, thickness/1000)
            txt = txt + 'Calculate,'
            for flux in fluxes:
                txt = txt + '%0.0f|'%(flux)
            txt = txt[:-1] + ','
            for flux in fluxes:
                txt = txt + '1|'
            txt = txt[:-1] + ','
            txt = txt + 'RISE_materials'
            
    with open(os.path.join(systemPath,'..','data','rise_spec_file.csv'), 'w') as f:
        f.write(txt)

    
    