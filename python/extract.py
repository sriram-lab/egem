"""
Extracting data from CCLE data obtained from Ghandi et al., 2019

The files that are needed as the input for the script:
  * CCLE_GCP.csv - obtained from the supplementary data under a different name, but it is the file that contains the global chromatin proteomics dataset

The files created from this script:
  * ccle_names.txt - the names of the CCLE cell lines corresponding to the H3 markers
  * h3marks.txt - the names of the H3 markers
  * h3_relval.txt - the values of the relative proteomics intensities obtained from Ghandi et al., 2019.
@author: Scott Campit
"""

import numpy as np
import pandas as pd
from openpyxl import load_workbook

def get_media(dictionary, key, lst):
    """
    This function requires a list of synonyms for a specific medium. It will construct a dictionary that maps the key to various synonyms of the medium.
    """
    for i in lst:
        try:
            dictionary[key].append(i)
        except KeyError:
            dictionary[key] = [i]
    return dictionary

def mapper(dict, series):
    """
    This function maps the values in a medium key to the dataframe elements to get the keys into the dataframe.
    """
    for k, v in dict.items():
        idx = series.isin(v)
        tmp = series[idx]
        tmp[idx] = k
        series[idx] = tmp
    return series

def extract():
    """
    Extract will get the medium conditions that were used for the CCLE histone proteomics paper to grow the cells and classify them to simpler medium conditions. The result will be outputted as a table.
    """

    # Basic data clean up
    df = pd.read_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/CCLE_GCP.csv')
    df = df.drop('BroadID', axis=1)

    media = pd.read_excel(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/summary.xlsx', sheet_name='Cell Line Annotations', usecols = ['CCLE_ID', 'Growth.Medium'])

    df = pd.merge(df, media, left_on='CellLineName', right_on='CCLE_ID')
    unqmed = df['Growth.Medium'].unique()

    to_remove = [
        ' +10%FBS',
        '+10%FBS',
        '+5%FBS',
        ' +5%FBS',
        ' with 10% fetal calf serum',
        ' with 10% fetal bovine serum',
        ' with 20% fetal bovine serum',
        ' (ATCC catalog # 30-2006) + 5%  FBS',
        '  + 5%  FBS',
        ' with 15%  fetal calf serum',
        ' + 10% h.i. FBS',
        ' + 10 % FBS',
        '+20% FBS',
        '; heat inactivated fetal bovine serum',
        ' + 5% FBS',
        ' + 10% FBS',
        ' +10 % FBS',
        ' +10 %FBS',
        '  +5% FBS',
        '  +10% FBS',
        ' + 10%FBS',
        ' (FBS),10%',
        ' +15%FBS',
        ' (FBs),10%',
        '+5% FBS',
        ' +10% FBS',
        '+10% FBS',
        '+ 10%FBS',
        '+20%FBS',
        '-20% FBS',
        ' + 20% FBS',
        '+ 10% FBS',
        ' +5-10 h.i. FBS',
        '+ 20% FBS',
        ' + 20% h.i. FBS',
        ' + 5 % FBS',
        '+15%FBS',
        ' +10 FBS',
        ' +20% h.i. FBS',
        ' +20 % FBS',
        ' +10% h.i. FBS',
        '+ 10-20% h. i. FBS',
        '+ 5% FBS',
        ' + 10% h.i.FBS',
        ' +15-20% h.i. FBS',
        ' +20% h.i FBS',
        ' +10% h. i. FBS',
        ' +10-20% h.i. FBS',
        ' +15% h. i. FBS',
        '+ 15% FBS',
        '+10% h.i FBS',
        ' + 10% h. i. FBS',
        ' +10-20% h. i. FBS',
        ' +20% h. i. FBS',
        ' +10%  h. i. FBS',
        ' + 10% h. i. FBS',
        ' +10% h.i.FBS',
        '+ 20% h. i. FBS',
        ' +10%h.i. FBS',
        ' +20% h.i. FBs',
        ' +15% h.i.FBS',
        ' +20% h.i.FBS',
        ' +10% h,i, FBS',
        '+ 10% h.i. FBS',
        '+10-20% h.i.FBS',
        ' 10% FBS',
        ' +15% h.i. FBS',
        ' +15% FBS',
        ' +0.5% human serum (+0.005 lu/ml TSH +5ug/ml human insulin-cells grow also without these latter supplements)',
        ' with fetal calf serum',
        ' heat inactivatedfetal bovine serum (FBS), 10%',
        ' +10%',
        '+10%fbs',
        ' with 15%  fetal bovine serum',
        ' with 15% fetal calf serum',
        ' with 20% heat inactivated fetal bovine serum',
        ' with 10% calf serum (FCS can be used)',
        ' with 5% fetal bovine serum',
        ' wiyh10% fetal bovine serum',
        ' with 10% heat inactivated fetal bovine serum',
        '. Refer cancer Res. 42:3858 (1982) for other medium with low serum concentration.',
        '+20 % FBS',
        '+15% FBS',
        ' + 20%FBS',
        '-10% FBS',
        "+20 % FBS"
    ]

    # Get synonyms for various medium components in dict with single key
    waymouth = {}
    waymouth_syn = [
        "Waymouth's",
        "Waymouth MB 7521 medium"
    ]
    waymouth = get_media(waymouth, 'Waymouth', waymouth_syn)

    l15 = {}
    l15_syn = [
        "L15",
        "Leibovitz's L-15 Medium"
    ]
    l15 = get_media(l15, "L15", l15_syn)

    rpmi = {}
    rpmi_syn = [
        "RPMI1640",
        "RPMI-1640",
        "RPMI 1640 medium",
        "RPMI 1640",
        "RPMI 1640(or DM-160AU)",
        "RPMI 1640 medium.",
        "RPMI ",
        "90-95% RPMI 1640",
        "80% RPMI 1640",
        "90% RPMI 1640",
        "85% RPMI 1640",
        "80-90% RPHI 1640",
        "80-90% RPMI 1640",
        "90% RPMI",
        "80-85% RPMI 1640",
        "80%RPMI-1640",
        "90% rPMI 1640",
        "90% RPMI1640",
        "80% RPMI 1640 "<
        "90% RPMI",
        "90 % Iscove's MDMD or RPMI 1640",
        "80-90% RPMI 1640",
        "90%RPMI 1640",
        "80-90% RPHI 1640",
        "RPI-1640",
        "80-90% RPHI 1640",
        "RPMI-1640 +1 mM Sodium pyruvate",
        "80% RPMI 1640 ",
        "RPMI-1640  ",
        "RPMI-1640",
        "RPMI-1640 ",
    ]
    rpmi = get_media(rpmi, "RPMI", rpmi_syn)

    rpmi_gln = {}
    rpmi_gln_syn = [
        "RPMI 1640  with L-glutamine (300mg/L), 90%",
        "RPMI 1640 with L-glutamine(300mg/L), 90%;",
        "RPMI w Gln"
    ]
    rpmi_gln = get_media(rpmi_gln, "RPMI w Gln", rpmi_gln_syn)

    ham_f12 = {}
    ham_f12_syn = [
        "Han's F12 medium 90%",
        "Ham's F12 medium",
        "HamF12",
        "Ham's F12",
        "HAMS F12 ",
        "F-12K  ATCC ",
        "F-12",
        "HAM F-10"
    ]
    ham_f12 = get_media(ham_f12, "HAM F-12", ham_f12_syn)

    ham_f10 = {}
    ham_f10_syn = [
        "HamF10",
        "HAMSF10",
        "Hams F10"
    ]
    ham_f10 = get_media(ham_f10, "HAM F-10", ham_f10_syn)

    aMEM = {}
    aMEM_syn = [
        "alpha-MEM",
        "80% alpha-MEM ",
        "70% alpha-MEM ",
        "alphaMEM"
    ]
    aMEM = get_media(aMEM, "alpha-MEM", aMEM_syn)

    dmem = {}
    dmem_syn = [
        "Dulbecco's modified Eagle's medium",
        "DMEM",
        "85% Dulbecco's MEM",
        "80% Dulbecco's MEM",
        "90% Dulbecco's MEM",
        "DEMEM",
        "DMEM ",
        "Eagle's minimal essential medium",
        "Eagle's minimal essential",
        "Eagle",
        "MEM",
        "EMEM",
        "Eagle MEM",
        "EMEM ",
        "80% Dulbecco's MEM(4.5g/L glucose)",
        "90% Dulbecco'sMEM(4.5g/l glucose)",
        "90% Dulbecco's MEM (4.5g/L glucose)",
        "90% Dulbecco's MEM(4.5 g/L glucose)",

    ]
    dmem = get_media(dmem, "DMEM", dmem_syn)

    mccoy = {}
    mccoy_syn = [
        "McCoy's 5A",
        "McCoy 5A"
    ]
    mccoy = get_media(mccoy, "McCoy 5A", mccoy_syn)

    dmem_f12_11 = {}
    dmem_f12_11_syn = [
        "DMEM:F12 (1:1)",
        "DMEM:HAM's F12  1:1",
        "DMEM:F12",
        "DMEM/F12(1:1)",
        "DMEM:F:12 (1:1)",
        "DMEM:F12(1:1)",
        "(DMEM: HamF12=1:1)",
        "(DMEM:HamF12=1:1)",
        "DMEM:HAM's F12  1:1 ",
        "80% mixture of Dulbecco's MEM + Ham's F 12 at 1:1",
        "80% mixture of Ham's F12 + Dulbecco's MEM (at 1:1)",
        "DMEM:HAM's F12  1:1  "
    ]
    dmem_f12_11 = get_media(dmem_f12_11, "DMEM:F12 (1:1)", dmem_f12_11_syn)

    dmem_rpmi_21 = {}
    dmem_rpmi_21_syn = [
        "2/3 DMEM:1/3 RPMI",
        "DMEM:RPMI (2:1)"
    ]
    dmem_rpmi_21 = get_media(dmem_rpmi_21, "DMEM:RPMI (2:1)", dmem_rpmi_21_syn)

    mcdb_m199 = {}
    mcdb_m199_syn = [
        "MCDB 105 (1:1)Medium 199",
        "MCDB 105:MEDIUM 199 (1:1)",
        "MCDB 105 1:1 Media 199",
    ]
    mcdb_m199 = get_media(mcdb_m199, "MCDB105:M199", mcdb_m199_syn)

    rpmi_dmem_11= {}
    rpmi_dmem_11_syn = [
        "RPMI:EMEM (1:1)",
        "RPMI:Eagles MEM",
        "RPM:MEM ",
        "RPMI:Eagles MEM"
    ]
    rpmi_dmem_11 = get_media(rpmi_dmem_11, "RPMI:Eagles MEM", rpmi_dmem_11_syn)

    dmem_iscove = {}
    dmem_iscove_syn = [
        "40% Dulbecco's MEM +40% Iscove's MDM"
    ]
    dmem_iscove = get_media(dmem_iscove, "DMEM:Iscove", dmem_iscove_syn)

    rpmi_iscove = {}
    rpmi_iscove_syn = [
        "90 % Iscove's MDMD or RPMI 1640",
        "80-90% Iscove's MDM or RPMI 1640"
    ]
    rpmi_iscove = get_media(rpmi_iscove, "RPMI:Iscove", rpmi_iscove_syn)

    iscove = {}
    iscove_syn = [
        "80% Iscove's MDM",
        "90% Iscove's MDM",
        "80-90% Iscove's MDM",
        "IMDM",
        "IMDM "
    ]
    iscove = get_media(iscove, "Iscove", iscove_syn)

    rpmi_f12 = {}
    rpmi_f12_syn = [
        "RPMI:F12 1:1",
        "RPMI:Ham F12"
    ]
    rpmi_f12 = get_media(rpmi_f12, "RPMI:F12", rpmi_f12_syn)

    list_of_medium = [
        waymouth,
        l15,
        rpmi,
        rpmi_gln,
        ham_f12,
        ham_f10,
        aMEM,
        dmem,
        mccoy,
        dmem_f12_11,
        dmem_rpmi_21,
        mcdb_m199,
        rpmi_dmem_11,
        dmem_iscove,
        rpmi_iscove,
        iscove,
        rpmi_f12
    ]

    import re
    unqmed = pd.Series(unqmed)
    unqmed = unqmed.str.replace('|'.join(map(re.escape, to_remove)), '')
    medium_map = df['Growth.Medium'].str.replace('|'.join(map(re.escape, to_remove)), '')

    for medium in list_of_medium:
        medium_map = mapper(medium, medium_map)

    # Create a map of all the unique medium in the CCLE study (media.txt)
    #unqmed = unqmed.unique()

    # Print out a file for all corresponding medium to common identifiers (ccle_media.txt)
    for medium in medium_map:
        print(medium)

    #df['CellLineName'] = df['CellLineName'].str.split('_').str[0]

    # output everything for matlab

    #df['CellLineName'].to_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/ccle_names.txt', header=False, index=False)

    #_ = df.pop('CellLineName')

    #h3marks = df.columns.tolist()

    #with open(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/ccle_h3marks2.txt', 'w') as f:
    #    for i in h3marks:
    #        print(i)
    #        f.write("%s\n" % i)

    #df.to_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/h3_relval.txt', header=False, index=False)

#extract()

def make_medium_conditions():
    """
    make_medium_conditions is an extension of what Shen et al., did with their analysis, but will now consider all medium components. To conduct this analysis, I will assume that RPMI is the default condition for current substrate uptake rates, and will scale the remaining substrate uptake rates for other medium conditions based on RPMI.
        For instance, if RPMI has glc = 2 and DMEM has glc = 4.5, I will scale the substrate uptake rate for glc in DMEM conditions to be 2.25x that of the RPMI condition.
    """

    # Get the default substrate uptake rates from RECON1. For those with uptake rates set to 0 initially, I changed it to a small number (-0.005 = default value)
    path = r'/mnt/c/Users/scampit/Desktop/MeGEM/data/'

    # For all medium, I will use RPMI as the base condition, and will scale the substrate uptake rates accordingly.

        # If the medium component is not present, the substrate uptake rate will be set to 0. Otherwise, it will be the default uptake rate * the scaling factor.

        # If the specific formulation was not provided, I simply divided the amount according to stoichiometric equivalents and added the medium concentrations together. If the amount of the medium component = infinity, I kept the infinity, as the amount will still be in excess. All infinity values will be converted into -100 for the substrate uptake rates.

    # Single medium conditions I already have that I need to make the mixed media that are not already measured...
    rpmi = pd.read_excel(path+'medium.xlsx', sheet_name='RPMI')
    dmem = pd.read_excel(path+'medium.xlsx', sheet_name='DMEM')
    mcdb105 = pd.read_excel(path+'medium.xlsx', sheet_name='MCDB105')
    m199 = pd.read_excel(path+'medium.xlsx', sheet_name='M199')
    f12 = pd.read_excel(path+'medium.xlsx', sheet_name='HAM F-12')
    iscove = pd.read_excel(path+'medium.xlsx', sheet_name='Iscove')

    # load the excel file so you don't overwrite the excel sheet (again...)
    book = load_workbook(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/medium.xlsx')
    writer = pd.ExcelWriter(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/medium.xlsx', engine='openpyxl')
    writer.book = book
    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

    # Media that I need to "make". Let's start with 2:1 DMEM: RPMI
    dmem_21 = dmem.copy(deep=True)
    dmem_21['mM'] = dmem_21['mM'].replace('Infinity', np.inf)
    dmem_21['mM'] = dmem_21['mM']*2/3
    rpmi_21 = rpmi.copy(deep=True)
    rpmi_21['mM'] = rpmi_21['mM']*1/3
    dmem_rpmi = pd.concat([dmem_21, rpmi_21]).groupby('Components')['mM'].sum().reset_index()
    #dmem_rpmi.to_excel(writer, sheet_name='DMEM-RPMI 2-1', index=False)

    # MCDB105-M199
    mcdb105_11 = mcdb105.copy(deep=True)
    mcdb105_11['mM'] = mcdb105_11['mM']*1/2
    m199_11 = m199.copy(deep=True)
    m199_11['mM'] = m199_11['mM'].replace('Infinity', np.inf)
    m199_11['mM'] = m199_11['mM']*1/2
    mcdb105_m199 = pd.concat([mcdb105_11, m199_11]).groupby('Components')['mM'].sum().reset_index()
    #mcdb105_m199.to_excel(writer, sheet_name='MCDB105-M199', index=False)

    # RPMI-F12
    rpmi_11 = rpmi.copy(deep=True)
    rpmi_11['mM'] = rpmi_11['mM']*1/2
    f12_11 = f12.copy(deep=True)
    f12_11['mM'] = f12_11['mM']*1/2
    rpmi_f12 = pd.concat([rpmi_11, f12_11]).groupby('Components')['mM'].sum().reset_index()
    #rpmi_f12.to_excel(writer, sheet_name='RPMI-F12', index=False)

    # DMEM-iscove
    dmem_11 = dmem.copy(deep=True)
    dmem_11['mM'] = dmem_11['mM'].replace('Infinity', np.inf)
    dmem_11['mM'] = dmem_11['mM']*1/2
    iscove_11 = iscove.copy(deep=True)
    iscove_11['mM'] = iscove_11['mM']*1/2
    dmem_iscove = pd.concat([dmem_11, iscove_11]).groupby('Components')['mM'].sum().reset_index()
    #dmem_iscove.to_excel(writer, sheet_name='DMEM-Iscove', index=False)

    # RPMI-Iscove
    rpmi_11 = rpmi.copy(deep=True)
    rpmi_11['mM'] = rpmi_11['mM']*1/2
    iscove_11 = iscove.copy(deep=True)
    iscove_11['mM'] = iscove_11['mM']*1/2
    rpmi_iscove = pd.concat([rpmi_11, iscove_11]).groupby('Components')['mM'].sum().reset_index()
    #rpmi_iscove.to_excel(writer, sheet_name='RPMI-Iscove', index=False)

    #writer.save()

    # Now that I made the medium conditions, I need to calculate those scaling coefficients that I will use for each medium condition. Divide all dataframes by RPMI.
    xlfil = pd.ExcelFile(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/medium.xlsx')
    media_conditions = xlfil.sheet_names

    for medium in media_conditions:
        # Hold this constant
        rpmi = pd.read_excel(path+'medium.xlsx', sheet_name='RPMI', index_col='Components')

        # New medium condition
        df = pd.read_excel(path+'medium.xlsx', sheet_name=medium, index_col='Components')
        df['mM'] = df['mM'].replace('Infinity', np.inf)

        # Create a column that will be written back into the excel file.
        df['Scaling factor'] = df['mM'].divide(rpmi['mM'], axis='index', fill_value=0)
        df['Scaling factor'] = df['Scaling factor'].replace(np.inf, 10)
        #df.to_excel(writer, sheet_name=medium, index=True)
        #writer.save()

    # I have scaling factors. Now I will add to the uptake.xlsx excel sheet to  map them back to specific substrate uptake rates in the metabolic model per medium condition.

    # load the excel file so you don't overwrite the excel sheet
    book = load_workbook(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/uptake.xlsx')
    writer = pd.ExcelWriter(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/uptake.xlsx', engine='openpyxl')
    writer.book = book
    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)

    for medium in media_conditions:
        # Keep this constant
        default_uptake = pd.read_excel(path+'uptake.xlsx', sheet_name='RPMI', index_col='Metabolite')

        # This will change based on medium component
        df = pd.read_excel(path+'medium.xlsx', sheet_name=medium, index_col='Components')
        print(df)

        # Create a column that will be written back into the excel file
        default_uptake['Adjusted uptake rate'] = default_uptake['Substrate Uptake Rate'].multiply(df['Scaling factor'], axis=0, fill_value=0)
        default_uptake.to_excel(writer, sheet_name=medium,index=True)
        writer.save()

make_medium_conditions()

def iterred(sheetnam=''):
    """
    iterred will iterate over the files outputted from the eGEM metabolic model and construct a list of dataframes that can be used for visualizations.
    """
    df = pd.read_excel('/mnt/c/Users/scampit/Desktop/MeGEM/matlab/tables/eGEMn.xlsx', sheet_name=sheetnam)
    df = df[df.columns.drop(list(df.filter(regex='.1')))]
    df = df.drop('Unnamed: 22', axis=1)
    df = df.head(50)
    return df
