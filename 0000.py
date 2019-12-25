import csv
import cobra.test
import cobra.io
from cobra.util.solver import linear_reaction_coefficients
from os.path import join
myfolder = 'C:/Users/15733/Desktop/模型构建/模型相关文件/'
myfile = 'iML1515-reverse-name.txt'
outfolder = 'C:/Users/15733/Desktop/模型构建/修改总酶量重新计算结果/单一碳源/1515计算结果/'
file = 'C:/Users/15733/Desktop/模型构建/修改总酶量重新计算结果/单一碳源/底物输入.txt'
d = {}
for one in open(file, 'r', encoding='utf-8'):
    tn = one.rstrip('\n').split('\t')
    if tn[0] in d:
        d[tn[0]].append(tn[2])
    else:
        d[tn[0]] = tn[2]
for i in d:
    model=cobra.io.read_sbml_model(myfolder + "iML1515.xml")
    model.reactions.get_by_id('EX_glc__D_e').bounds = (0.0,0.0)
    model.reactions.get_by_id(i).bounds = (-eval(d[i]),0.0)
    #model.reactions.get_by_id('BIOMASS_Ec_iML1515_core_75p37M').bounds = (0.2,0.2)
    model.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
    pfba_solution = cobra.flux_analysis.pfba(model)
    dict_fluxes = dict(pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-10])
    s1 = csv.reader(open(myfolder + myfile), delimiter='\t')  # 换行分隔符
    realist1 = []
    realist1.extend(s1)
    Dict = {}
    writefile = csv.writer(open(outfolder + i +'.txt', 'w', newline=''), dialect='excel',
                               delimiter='\t')
    for m in realist1:
        Dict[m[0]] = m
    for n in dict_fluxes:
        dict_fluxes[n] = round(dict_fluxes[n], 10)  # round（n，x）表示n值保留x位小数
        if n in Dict:
            Dict[n].insert(1, str(dict_fluxes[n]))  # insert(位置，值)
            writefile.writerow(Dict[n])
