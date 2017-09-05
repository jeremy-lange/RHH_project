'''
- to run this script from the command line, simply place this file, demographies.py, and ms-formatted simulations in the same directory and type "python run_dadi.py input_pattern func_str output"
- input_pattern is a common string in the file names of all ms-formatted simulations
- func_str is the string that identifies which demography you want to use in this maximization
- output is the name of the output file that this script will create
- note: the printing at the end of this script is less than ideal. It may be advantageous to edit to create better formatted output files
'''

import sys

import numpy

import dadi

import demographies

import glob

#get arguments from command line
args=sys.argv
input_pattern=args[1]
func_str=args[2]

if func_str=="trunk_2epoch_sizechange":
    params=([1,1])
    upper_bound = [ 100, 10]
    lower_bound = [ 1e-2,0]
elif func_str=="trunk_3epoch_sizechange":
    params=([1,1,1])
    upper_bound = [ 100,10, 10]
    lower_bound = [ 1e-2,0,0]
elif func_str=="bneck_3params":
    params=([0.047,0.042,0.0096])
    upper_bound = [ 100,10, 10]
    lower_bound = [ 1e-2,0,0]
elif func_str=="bneck_2params":
    params=([0.47,0.0096])
    upper_bound = [ 100, 10]
    lower_bound = [ 1e-2,0]
elif func_str=="trunk_2epoch_sizechange_bneck_3param":
    params=([1,1,0.042,0.047,0.0096])
    upper_bound = [ 100, 10,10,100,10]
    lower_bound = [ 1e-2,0,0,1e-2,0]
elif func_str=="trunk_2epoch_sizechange_bneck_2param":
    params=([1,1,0.047,0.0096])
    upper_bound = [ 100, 10,100,10]
    lower_bound = [ 1e-2,0,1e-2,0]
elif func_str=="trunk_3epoch_sizechange_bneck_3param":
    params=([1,1,1,0.042,0.047,0.0096])
    upper_bound = [ 100, 10,10,10,100,10]
    lower_bound = [ 1e-2,0,0,0,1e-2,0]
elif func_str=="trunk_3epoch_sizechange_bneck_2param":
    params=([1,1,1,0.047,0.0096])
    upper_bound = [ 100, 10,10,100,10]
    lower_bound = [ 1e-2,0,0,1e-2,0]
elif func_str=="IM_2params":
    params=([0.25, 1.25e-06])
    upper_bound = [ 10, 10]
    lower_bound = [ 0.0001,0]
elif func_str=="IM_2params_decpriors":
    func_str="IM_2params"
    params=([0.25, 1.25e-06])
    upper_bound = [ 10, 0.1]
    lower_bound = [ 0.0001,0]
elif func_str=="trunk_2epoch_sizechange_IM_2param":
    params=([1,1,0.25,1.25e-06])
    upper_bound = [ 100, 10,10,10]
    lower_bound = [ 1e-2,0,0.0001,0]
elif func_str=="trunk_2epoch_sizechange_IM_2param_decpriors":
    func_str="trunk_2epoch_sizechange_IM_2param"
    params=([1,1,0.25,1.25e-06])
    upper_bound = [ 100, 10,10,0.1]
    lower_bound = [ 1e-2,0,0.0001,0]
elif func_str=="trunk_3epoch_sizechange_IM_2param":
    params=([1,1,1,0.25,1.25e-6])
    upper_bound = [ 100, 10,10,10,10]
    lower_bound = [ 1e-2,0,0,0.0001,0]
elif func_str=="trunk_3epoch_sizechange_IM_2param_decpriors":
    func_str="trunk_3epoch_sizechange_IM_2param"
    params=([1,1,1,0.25,1.25e-6])
    upper_bound = [ 100, 10,10,10,10]
    lower_bound = [ 1e-2,0,0,0.0001,0]

func=getattr(demographies,func_str)
files=glob.glob(input_pattern +'*txt')
likelihoods=[]
optimized_params=[]

for f in range(0,len(files)):
    print f
    data=dadi.Spectrum.from_ms_file(files[f],average=False)
    ns = data.sample_sizes
    pts_l = [110,120,130]
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    model = func_ex(params, ns, pts_l)
    ll_model = dadi.Inference.ll_multinom(model, data)
    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound)
    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,lower_bound=lower_bound,upper_bound=upper_bound,maxiter=10)#,verbose=len(params))
    model=func_ex(popt, ns, pts_l)
    ll_opt=dadi.Inference.ll_multinom(model, data)
    likelihoods.append(ll_opt)
    optimized_params.append([popt])

output=args[3]
target=open(output,'w')
for i in range(0,len(likelihoods)):
    target.write(str(likelihoods[i]))
    target.write("\t")
    for j in range(0, len(optimized_params[i])):
        target.write(str(optimized_params[i][j]))
        target.write("\t")
    target.write("\n")
target.close()
