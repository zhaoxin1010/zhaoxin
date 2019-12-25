# encoding: utf-8

#input is sbml model and output-file' s name such as 'initialmdoel.xlsx'
def sbml_excel(modelname,output):
    from cobra.io.dict import model_to_dict, model_from_dict,metabolite_from_dict,gene_from_dict,reaction_from_dict
    from cobra.core import Model
    from cobra.io import read_sbml_model, write_sbml_model
    import pandas as pd
    model=read_sbml_model(modelname)
    a=model_to_dict(model,sort=False)
    #writer = pd.ExcelWriter(modelname[:-4]+'.xlsx')
    writer = pd.ExcelWriter(output)
    #pd.DataFrame(a['compartments']).to_excel(writer,'Sheet1',index=False)
    pd.DataFrame(a['metabolites']).to_excel(writer,'metabolites',index=False)
    pd.DataFrame(a['genes']).to_excel(writer,'genes',index=False)
    df_r=pd.DataFrame(a['reactions'])
    df_r['reaction_eq'] = df_r.id.apply(lambda x: model.reactions.get_by_id(x).build_reaction_string(use_metabolite_names=False))
    df_r['reaction_eq_name'] = df_r.id.apply(lambda x: model.reactions.get_by_id(x).build_reaction_string(use_metabolite_names=True))
    #del df_r['metabolites'] #生成反应方程列后将代谢物列删除，可保留以方便再次读取
    df_r.to_excel(writer,'reactions',index=False)
    #df_r.to_excel(writer,'reactions',columns=["id","name","reaction","lower_bound","upper_bound","gene_reaction_rule"],index=False)
    writer.save()


#input is cobrapy model and out
def model_excel(model,output):
    from cobra.io.dict import model_to_dict, model_from_dict,metabolite_from_dict,gene_from_dict,reaction_from_dict
    from cobra.core import Model
    from cobra.io import read_sbml_model, write_sbml_model
    import pandas as pd
    a=model_to_dict(model,sort=False)
    writer = pd.ExcelWriter(output)
    #pd.DataFrame(a['compartments']).to_excel(writer,'Sheet1',index=False)
    pd.DataFrame(a['metabolites']).to_excel(writer,'metabolites',index=False)
    pd.DataFrame(a['genes']).to_excel(writer,'genes',index=False)
    df_r=pd.DataFrame(a['reactions'])
    df_r['reaction_eq'] = df_r.id.apply(lambda x: model.reactions.get_by_id(x).build_reaction_string(use_metabolite_names=False))
    #del df_r['metabolites'] #生成反应方程列后将代谢物列删除，可保留以方便再次读取
    df_r.to_excel(writer,'reactions',index=False)
    #df_r.to_excel(writer,'reactions',columns=["id","name","reaction","lower_bound","upper_bound","gene_reaction_rule"],index=False)
    writer.save()
#if __name__ == "__main__":
    #sbml_excel(modelname)

def mkdir(path):
    # 引入模块
    import os

    # 去除首位空格
    path = path.strip()
    # 去除尾部 \ 符号
    path = path.rstrip("\\")

    # 判断路径是否存在
    # 存在     True
    # 不存在   False
    isExists = os.path.exists(path)

    # 判断结果
    if not isExists:
        # 如果不存在则创建目录
        # 创建目录操作函数
        os.makedirs(path)

        print(path + ' 创建成功')
        return True
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print(path + ' 目录已存在')
        return False

def pathway (model,fluxes,outputfile_name):
    import csv
    import pandas as pd
    flux = open(outputfile_name, 'w')
    for r,v in fluxes.iteritems():
        if abs(v)>1e-6:
            flux.write(r +'\t' + str(round(v,4)) + '\t' + model.reactions.get_by_id(r).build_reaction_string(use_metabolite_names=True) +'\n')
    flux.close()

def calculate (precursors,model,outpath,outfilename):
    import csv
    import output
    import cobra
    writefile = open(outfilename, 'w')
    model.reactions.get_by_id('rxn00062_c0').bounds=(0,1000)
    with model:
      for i in precursors:
        #print(i)
        i = str(i)
        product = model.metabolites.get_by_id(i)
        demand = model.add_boundary(product, type='demand')  # add demand reaction as the objective
        model.objective = demand
        outputfile_name = outpath + i + '.txt'
        fluxes = cobra.flux_analysis.pfba(model).fluxes
        a='DM_'+i
        obj=str(fluxes[a])
        reduce=reduce_degree(i)
        writefile.write(i+"\t"+model.metabolites.get_by_id(i).name+"\t"+obj+"\t"+str(reduce[0])+"\t"+str(reduce[1])+"\n")   
        output.pathway(model,fluxes,outputfile_name)
    writefile.close()
def reduce_degree (i):
    import csv
    import re
    #folder='./media/jupyter/yuan_qq/Quality control-new/metabolites_formula.txt'#存储代谢物化学式的文件
    s = csv.reader(open('metabolites_formula.txt'), delimiter="\t")
    realist = []
    realist.extend(s)
    for j in realist:
        if i==j[0]: 
            rf = re.findall (r'([A-Z][a-z]*)(\d*)', j[5])
            value_C=0
            value_H=0
            value_N=0   
            value_O=0
            value_P=0
            value_S=0
            for i in rf:
                if 'C' in i:
                    if i[1]=='':
                        value_C=1
                    else:
                        value_C=i[1]
                if 'H' in i:
                    if i[1]=='':
                        value_H=1
                    else:
                        value_H=i[1]
                if 'O' in i:
                    if i[1]=='':
                        value_O=1
                    else:
                        value_O=i[1]
                if 'N' in i:
                    if i[1]=='':
                        value_N=1
                    else:
                        value_N=i[1]
                if 'P' in i:
                    if i[1]=='':
                        value_P=1
                    else:
                        value_P=i[1]
                if 'S' in i:
                    if i[1]=='':
                        value_S=1
                    else:
                        value_S=i[1]
            Huanyuanli=4*int(value_C) + int(value_H)-int(j[6]) + (-2)*int(value_O) + (-3)*int(value_N) + 6*int(value_S) + 5*int(value_P)
            if not Huanyuanli==0:
                return (j[5],round(240/Huanyuanli,4))
            else:
                return (j[5],0.0)
