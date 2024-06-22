


import matplotlib.pyplot as plt
from numpy import *           
from module.transient_model import *


'''build pipeline model from the configuration file and run simulation'''
# specify the name of the configuration file
name = 'network-112'
# build transient model
model = PipelineSystem(name)
# run simulation
model.run( mode=MODE.MOC)
# save results
model.write_recorder(filename=data_root_path+ name+'Model_1')
# plot results
npipe=model.n_pipe 
t = model.pick_data('time')
maxh=np.zeros(model.n_pipe)
index=range(model.n_pipe)
# plot 
index=[0,14,23,45]
fig1 = plt.figure(1, figsize=(6.4, 4))
plt.subplots_adjust(left=0.135, wspace=1.5, hspace=0.5,
                    bottom=0.15, right=0.85, top=0.9)
lenPlots=len(index)
for num,i  in enumerate( index):
    plt.subplot(lenPlots,1,num+1)
    y = model.pick_data(h=-1,pipenum=i)
    plt.plot(t, y, 'k', linewidth=1)
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Head (m)', fontsize=12)
    plt.xlim([0,model.T])
plt.show()