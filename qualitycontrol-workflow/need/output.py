
def pathway (model,fluxes,outputfile_name):
    import csv
    import pandas as pd
    flux = open(outputfile_name, 'w')
    for r,v in fluxes.iteritems():
        if abs(v)>1e-6:
            flux.write(r +'\t' + str(round(v,4)) + '\t' + model.reactions.get_by_id(r).build_reaction_string(use_metabolite_names=True) +'\n')
    flux.close()
    

'''
def pathway (fluxes,outputfile_name):
    import csv
    import pandas as pd
    flux = open(outputfile_name, 'w')
    modelname = 'SEED-name.txt'
    s = csv.reader(open(modelname), delimiter="\t")
    realist = []
    realist.extend(s)
    for r,v in fluxes.iteritems():
        if abs(v)>1e-6:
            #print (r,v)
            for i in realist:
                if i[0]==r:
                    flux.write(i[0] +'\t' + str(v) + '\t' + i[1] +'\n')
    flux.close()
    
'''   
'''
def pathway (fluxes,inputfile_name,outputfile_name):
    import csv
    import pandas as pd
    flux = open(outputfile_name, 'w')
    modelname = pd.read_excel(inputfile_name, sheet_name=2)
    for r,v in fluxes.iteritems():
        if abs(v)>1e-6:
            #print (r,v)
            for i in modelname.id:
                if i==r:
                    flux.write(i[0] +'\t' + str(v) + '\t' + str(modelname[modelname['id']==i]['reaction_eq_name']) +'\n')
    flux.close()
'''
