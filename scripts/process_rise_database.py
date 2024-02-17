import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict
from algorithms import getMaterialClass

if __name__ == "__main__":
    
    systemPath = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(systemPath,'..','data','rise_materials')+os.sep
    out_dir = os.path.join(systemPath,'..','data','rise_materials_processed')+os.sep
    out_dir_spec = 'rise_materials_processed' + os.sep
    
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
    
    warn_mass = False
    warn_thickness = False
    
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
                if warn_thickness:
                    print("Warning no thickness found for %s"%(file))
            
            if thickness == '':
                if warn_thickness:
                    print("Warning no thickness found for %s"%(file))
                thickness = 0.0127
            elif type(thickness) == str:
                if '-' in thickness:
                    thickness = np.mean([float(x) for x in thickness.split('-')])
            
            if 'Specific extinction area (avg) (m2/kg)' in scalarColumns:
                soot_yield = scalars['Specific extinction area (avg) (m2/kg)']
                if type(soot_yield) == str:
                    if (len(soot_yield) == 0) or (soot_yield in ['NA']):
                        #print("Warning no SEA found for %s"%(file))
                        soot_yield = -1
                    else:
                        print(file, soot_yield, len(soot_yield))
                    #assert False, "Stopped"
                else:
                    soot_yield = soot_yield/8700
                print(file, scalars['Specific extinction area (avg) (m2/kg)'], soot_yield)
            else:
                print("Warning no SEA found for %s"%(file))
                print(scalarColumns)
                soot_yield = -1
            
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
                    if warn_mass:
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
            material_dict['FYI'] = file.split('-')[-1].split('.txt')[0]
            material_dict['FYI2'] = series
            material_dict['SootYield'] = soot_yield
            
            tmp = txt.split('VectorData\n')[1].split('\n')
            while tmp[-1] == '': tmp = tmp[:-1]
            
            tmp2 = [[float(x) for x in y.split(';')] for y in tmp[1:] if (('Time' not in y) and (len(y.replace(';','')) > 0) and ('s' not in y)) ]
            hrrdata = np.array(tmp2)
            
            #hrrdata = np.genfromtxt(file,delimiter=';',skip_header=8,usecols=[0,1])
            
            times = hrrdata[:, 0]
            hrrpua = hrrdata[:, 1]
            
            material_dict['raw_times'] = times
            material_dict['raw_hrrpua'] = hrrpua
            
            if type(tig) == float:
                hrrpua[times < tig] = 0
            
            matind = material.replace('%','') # + "_" + series
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
        soot_yields = []
        for thickness in thicknesses:
            if type(thickness) is str: continue
            if thickness <= 0: continue
            fluxes = list(material_database[material][thickness].keys())
            if len(fluxes) < 2: continue
            densities = []
            HoCs = []
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
                    soot_yield = material_database[material][thickness][flux][i]['SootYield']
                    if soot_yield > 0:
                        soot_yields.append(soot_yield)
                    HoC = material_database[material][thickness][flux][i]['HeatOfCombustion']
                    if type(HoC) == str:
                        pass
                    else:
                        HoCs.append(HoC)
                if len(densities) == 0: continue
                material_database[material][thickness][flux]['tign'] = np.mean(tign)
            if material_database_filtered[material] is False:
                material_database_filtered[material] = defaultdict(bool)
                material_database_filtered[material]['Density'] = []
                material_database_filtered[material]['HeatOfCombustion'] = []
            material_database_filtered[material][thickness] = material_database[material][thickness]
            material_database_filtered[material]['Density'].append(np.mean(densities))
            material_database_filtered[material]['HeatOfCombustion'].append(np.mean(HoCs))
            #print(material, len(fluxes), fluxes)
            #print(densities)
        
        if material_database_filtered[material] is False:
            material_database_filtered.pop(material)
        else:
            pass
            if len(soot_yields) == 0:
                print("Warning no soot yields found for material %s"%(material))
                material_database_filtered[material]['SootYield'] = 0.05
            else:
                material_database_filtered[material]['SootYield'] = np.nanmean(soot_yields)
                print("Soot yield %0.8f found for material %s"%(np.nanmean(soot_yields), material))
        
    materials = list(material_database_filtered.keys())
    
    for material in materials:
        thicknesses = list(material_database_filtered[material].keys())
        thicknesses.remove('Density')
        thicknesses.remove('SootYield')
        thicknesses.remove('HeatOfCombustion')
        for thickness in thicknesses:
            fluxes = list(material_database_filtered[material][thickness].keys())
            for flux in fluxes:
                testCount = material_database_filtered[material][thickness][flux]['testCount']
                tign = material_database_filtered[material][thickness][flux]['tign']
                outInt2 = 15
                if testCount > 1:
                    tMax = 0
                    tests = ''
                    for i in range(1, testCount+1):
                        raw_times = material_database_filtered[material][thickness][flux][i]['raw_times']
                        raw_hrrpuas = material_database_filtered[material][thickness][flux][i]['raw_hrrpua']
                        tign = material_database_filtered[material][thickness][flux][i]['tIgn']
                        tMax = max([tMax, np.nanmax(raw_times)-tign])
                        tests=tests + '%s|'%(material_database_filtered[material][thickness][flux][i]['FYI'])
                    tests = tests[:-1]
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
                    tests = '%s'%(material_database_filtered[material][thickness][flux][1]['FYI'])
                
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
                material_database_filtered[material][thickness][flux]['tests'] = tests
                
                d = pd.DataFrame(np.round(np.array([outTime, outHrrpua]).T, decimals=1), columns=['Time','HRRPUA'])
                mat = '%s-%02dmm'%(material, thickness)
                mat = mat.replace(',','_').replace(' ','_').replace('|','_').replace('%','')
                mat_name = 'RISE_'+mat
                while len(mat_name) > 40:
                    tmp = mat_name.split('_')
                    long_ind = np.argmax([len(t) for t in tmp])
                    tmp[long_ind] = tmp[long_ind][:-1]
                    mat_name = '_'.join(tmp)
                material_database_filtered[material]['mat_name'] = mat_name
                if mat in ignores: continue
                dataFile = os.path.join(out_dir, mat_name+'-%0.0f.csv'%(flux))
                d.to_csv(dataFile, index=False)
    
    material_database = material_database_filtered
    resultDir = "../../../out/Scaling_Pyrolysis/"
    inputFileDir = "../../../fds/Validation/Scaling_Pyrolysis/"
    #expFileDir = "../../../exp/RISE_Materials/"
    emissivity = 1
    txt = 'Code,Number,Series,Material,MaterialClass,DataFile,ResultDir,InputFileDir,ExpFileDir,'
    txt = txt + 'ReferenceExposure,ReferenceThickness,ReferenceTime,ReferenceHRRPUA,'
    txt = txt + 'ValidationTimes,ValidationHrrpuaColumns,ValidationFluxes,'
    txt = txt + 'Density,Conductivity,SpecificHeat,Emissivity,Thickness,'
    txt = txt + 'CharFraction,HeatOfCombustion,SootYield,'
    txt = txt + 'IgnitionTemperature,IgnitionTemperatureBasis,HeaderRows,FYI'
    
    for material in materials:
        thicknesses = list(material_database_filtered[material].keys())
        thicknesses.remove('Density')
        thicknesses.remove('SootYield')
        thicknesses.remove('HeatOfCombustion')
        thicknesses.remove('mat_name')
        for thickness in thicknesses:
            conductivity = 0.4 #material_database[material]['conductivity']
            specific_heat = 1. #material_database[material]['heatCapacity']
            density = np.mean(material_database[material]['Density'])
            heat_of_combustion = np.mean(material_database[material]['HeatOfCombustion'])
            soot_yield = np.mean(material_database[material]['SootYield'])
            
            fluxes = [x for x in list(material_database[material][thickness].keys()) if type(x) is float]
            fluxes.sort()
            
            code ='d'
            number = 1
            mat = '%s-%02dmm'%(material, thickness)
            mat = mat.replace(',','_').replace(' ','_').replace('|','_')
            
            mat_name = material_database_filtered[material]['mat_name']
            
            #mat_name = 'RISE_'+mat
            #while len(mat_name) > 40:
            #    tmp = mat_name.split('_')
            #    long_ind = np.argmax([len(t) for t in tmp])
            #    tmp[long_ind] = tmp[long_ind][:-1]
            #    mat_name = '_'.join(tmp)
            
            if mat in ignores:
                continue
            
            dataFiles = ''
            for flux in fluxes:
                dataFile = out_dir_spec+ '%s-%02d.csv'%(mat_name, flux)
                dataFiles = dataFiles + dataFile + '|'
            dataFiles = dataFiles[:-1]
            
            if (50 in fluxes):
                refFlux = 50
            else:
                ind = np.argmin([abs(x-50) for x in fluxes])
                refFlux = fluxes[ind]
            matClass = getMaterialClass(material)
            if matClass == 'Unknown': code = 's'
            m = material.replace('%',' ').replace('-',' ').replace(',',' ').replace('.',' ')
            while ('  ' in m): m = m.replace('  ',' ')
            m = 'RISE_'+m + '-%0.0f'%(thickness)
            m = m.replace('|','-').replace('handlaminated','laminated').replace('Painted','Paint').replace('painted','painted')
            
            txt = txt + "\n" + "%s,%s,%s,%s,%s,%s,%s,"%(code, number, 'RISE_Materials', m, matClass, dataFiles, resultDir)
            txt = txt + "%s,%s,%0.0f,%0.8f,%s-%0.0f.csv-Time,%s-%0.0f.csv-HRRPUA,"%(inputFileDir, out_dir_spec, refFlux, thickness, mat_name, refFlux, mat_name, refFlux)
            
            for flux in fluxes:
                txt = txt + '%s-%0.0f.csv-Time|'%(mat_name, flux)
            txt = txt[:-1] + ','
            for flux in fluxes:
                txt = txt + '%s-%0.0f.csv-HRRPUA|'%(mat_name, flux)
            txt = txt[:-1] + ','
            for flux in fluxes:
                txt = txt + '%0.0f|'%(flux)
            txt = txt[:-1] + ','
            txt = txt + '%0.1f,%0.4f,%0.4f,%0.4f,%0.8f,'%(density, conductivity, specific_heat, emissivity, thickness/1000)
            
            txt = txt + '%0.4f,%0.4f,%0.8f,'%(0, heat_of_combustion, soot_yield)
            
            txt = txt + 'Calculate,'
            for flux in fluxes:
                txt = txt + '%0.0f|'%(flux)
            txt = txt[:-1] + ','
            for flux in fluxes:
                txt = txt + '1|'
            txt = txt[:-1] + ','
            for flux in fluxes:
                txt = txt + '_'.join([x for x in sorted([x for x in material_database_filtered[material][thickness][flux]['tests'].split('|')])]) + '|'
            txt = txt[:-1]
            
    with open(os.path.join(systemPath,'..','data','rise_spec_file.csv'), 'w') as f:
        f.write(txt)

    
    