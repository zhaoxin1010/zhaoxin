# encoding: utf-8
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
    folder='/media/jupyter/yuan_qq/Quality control-new/metabolites_formula.txt'#存储代谢物化学式的文件
    s = csv.reader(open(folder), delimiter="\t")
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
'''
def calculate (precursors,model,outpath,outfilename):
    import csv
    import output
    import cobra
    writefile = open(outfilename, 'w')
    unsyn={}
    unsyn1=[]
    errorsyn={}
    errorsyn1=[]
    normal={}
    normal1=[]
    for i in precursors:
        #print(i)
        i = str(i)
        product = model.metabolites.get_by_id(i)
        demand = model.add_boundary(product, type='demand')  # add demand reaction as the objective
        model.objective = demand
        outputfile_name = outpath + i + '.txt'
        fluxes = cobra.flux_analysis.pfba(model).fluxes
        a='DM_'+i
        obj=str(round(fluxes[a],4))
        reduce=reduce_degree(i)
        if not reduce[1]=='' and reduce[1] < fluxes[a]:
            errorsyn1.append(i)
            errorsyn[i]={i,model.metabolites.get_by_id(i).name,obj,str(reduce)}
        elif not reduce[1]=='' and reduce[1] < 1e-6:
            unsyn1.append(i)
            unsyn[i]={i,model.metabolites.get_by_id(i).name,obj,str(reduce)}
        else:
            normal1.append(i)
            normal[i]={i,model.metabolites.get_by_id(i).name,obj,str(reduce)}         
        output.pathway(fluxes,'Seed1111.5.65232.xlsx',outputfile_name)
        model.reactions.get_by_id(a).knock_out()
    if not errorsyn1:
        writefile.write('Error synthetic biomass precursors'+"\n")
        for i in errorsyn1:
            writefile.write(errorsyn[i]+"\n")  
    if not unsyn1:
        writefile.write('Unsynthetic biomass precursors'+"\n")
        for i in unsyn:
            writefile.write(unsyn[i]+"\n")  
    if not normal1:
        writefile.write('Normal biomass precursors'+"\n")
        for i in normal:
            writefile.write(normal[i]+"\n")
    writefile.close()
'''
