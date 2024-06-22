import matplotlib.pyplot as plt
from numpy import *           
from module.transient_model import *
import pickle

'''build pipeline model from the configuration file and run simulation'''
# specify the name of the configuration file
name='emoc_paper_case1'
# build transient model
model = PipelineSystem(name)
# save transient model
a_file = open(data_root_path+ name+"_class.pkl", "wb")
pickle.dump(model, a_file)
a_file.close()
# run simulation
model.run( mode=MODE.MOC)
# save results
model.write_recorder(filename=data_root_path+ name)

'''plot results'''
npipe=model.n_pipe 
t = model.pick_data('time')
fig1 = plt.figure(1, figsize=(8, 10))
# define nodes where the heads are plotted
plot_nodes=[0,1,2,3]
# plot
for i,node in enumerate(plot_nodes):
    num_subplots=len(plot_nodes)
    
    plt.subplots_adjust(left=0.135, wspace=1.5, hspace=0.4,
                        bottom=0.15, right=0.85, top=0.9)
    plt.subplot(num_subplots,1,i+1)
    plt.plot(t, model.nodesData[f'node{node}'],'k', linewidth=1)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Head (m)', fontsize=12)
    plt.xlim([0,model.T])
plt.savefig(figure_root_path + name +str(model.id)+ '_' + "pipe_" + '.png')
plt.show()

