# encoding: utf-8

#coding ='utf-8'

def excel_sbml(modelname,outputname):
    import datetime
    import cobra
    from cobra import Metabolite, Reaction, Model
    from cobra.io import write_sbml_model
    #from __future__ import print_function
    import pandas as pd
    #myfolder='/media/jupyter/zhang_x/study/'
    #mname="seed_1"

    #若需包括更多代谢物信息，则可先在模型中创建代谢物对象并添加到模型
    #A = Metabolite('A')
    #model.add_metabolites([A, B, C, D, E, P])
    starttime = datetime.datetime.now()
    model = Model('model')
    colnames=['id','reaction_eq','name','lower_bound','upper_bound','Object']
#将上述值改为excel表格中的相应列名称数据以便正确读出值，表格中不一定包括所有列,不需要改变列的顺序，通过列头识别
#objr='Objective'  #目标反应的名字
#rin=[['PGI',-10,10]] #可设置多个输入反应及其上下限
#若表格中不包括目标反应数据需人为给出，同样对输入反应需人为添加约束,若已有相应数据则这两个值不用改
#data = pd.read_csv(mname+'.csv', delimiter=",", na_values=['(none)']).fillna('')
    data = pd.read_excel(modelname, sheet_name='reactions', header=0, na_values=['(none)']).fillna('')  #'(none)'变NaN，NaN数据用''填充
#print(data)
#直接从excel文件读取不用再转换txt,header行为数据起始行和表头（因前面行可能有模型说明文件） 表中空值替换为空字符以便处理，不能用dropna，否则整行数据都会丢掉
    keys=data.keys()
    for index, row in data.iterrows():
        r = Reaction(row[colnames[0]].strip())  #r即每个反应对应的反应名字（name）
        model.add_reaction(r)
        r.build_reaction_from_string(row[colnames[1]],fwd_arrow='-->', rev_arrow='<--', reversible_arrow='<=>', term_split='+')  #excel中方程式colnames[1],转成SBML格式
        if colnames[2] in keys:
            r.subsystem=row[colnames[2]]
        if colnames[3]in keys:
            r.lower_bound=row[colnames[3]]
            r.upper_bound=row[colnames[4]]
        if colnames[4] in keys: #目标反应
            if row[colnames[5]]:
               r.objective_coefficient=row[colnames[5]]
    #if colnames[5] in keys:
    #    if row[colnames[5]]:
    #       r.objective_coefficient=1
    #s=row[colnames[6]] #对基因关系处理要看表格中是如何存储and/or关系的
    #if s:
    #    genes=s.split(", ")
    #   r.gene_reaction_rule='('+" or ".join(genes)+' )'
#if colnames[3] not in keys:  #需人为设定输入反应边界，注意根据反应方向确定是改下限还是上限及值的正负
#    for r in rin:
#        rea=model.reactions.get_by_id(r[0])
#        rea.lower_bound=r[1]
#        rea.upper_bound=r[2]
    write_sbml_model(model, outputname)
    #endtime = datetime.datetime.now()
    #interval=(endtime-starttime).seconds
    #print(str(interval)+'sec')
'''
def excel_sbml(modelname):
    import re
    import cobra
    import ast
    import collections
    from cobra import Metabolite, Reaction, Model, Gene
    import pandas as pd
    from cobra.io.dict import model_to_dict, model_from_dict,metabolite_from_dict,gene_from_dict,reaction_from_dict
    from cobra.io import read_sbml_model, write_sbml_model
    model = Model()
    excel = pd.ExcelFile(modelname)
    #处理代谢物
    mdictl=excel.parse('metabolites',index_col=None).T.to_dict().values()#将dataframe转化为list of dictionaries
    for m in mdictl:
        m['charge']=int(m['charge']) #字符串转化为整数
    model.add_metabolites([metabolite_from_dict(metabolite) for metabolite in mdictl]) #mdictl需要是一个代谢物对象的list而非词典list
    #处理反应，因涉及目标反应确定，逐一添加而未从反应list添加
    reactions = excel.parse('reactions',index_col=None).fillna('') #fillna('')很重要，否则会产生Nan值
    del reactions['reaction_eq'] #不需要反应方程列数据
    rdictl=reactions.T.to_dict().values() #将dataframe转化为listofdictionaries以便模型读取
    for r in rdictl:
        metstr=r['metabolites']
        values = re.search(r"OrderedDict\((.*)\)", metstr).group(1)
        metdict = collections.OrderedDict(ast.literal_eval(values))
        r['metabolites']=metdict
        annstr=r['notes']
        values = re.search(r"OrderedDict\((.*)\)", annstr).group(1)
        anndict = collections.OrderedDict(ast.literal_eval(values))
        r['notes']=anndict
        r['lower_bound']=float(r['lower_bound'])
        r['upper_bound']=float(r['upper_bound'])
        rea=reaction_from_dict(r, model) #由词典得到新反应并添加到模型
        model.add_reaction(rea)
        if r['objective_coefficient']: 
            rea.objective_coefficient=float(r['objective_coefficient']) #不将字符转化成浮点数也可以，但转化更安全
#           objr=r['objective_coefficient']=float(r['objective_coefficient'])
    model.add_reactions([reaction_from_dict(reaction, model) for reaction in rdictl])
    write_sbml_model(model, modelname[:-5]+"_excel_xml.xml")
'''