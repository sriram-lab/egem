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

import pandas as pd

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
    extract will perform the actual extraction
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

extract()

def iterred(sheetnam=''):
    """
    iterred will iterate over the files outputted from the eGEM metabolic model and construct a list of dataframes that can be used for visualizations.
    """
    df = pd.read_excel('/mnt/c/Users/scampit/Desktop/MeGEM/matlab/tables/eGEMn.xlsx', sheet_name=sheetnam)
    df = df[df.columns.drop(list(df.filter(regex='.1')))]
    df = df.drop('Unnamed: 22', axis=1)
    df = df.head(50)
    return df
