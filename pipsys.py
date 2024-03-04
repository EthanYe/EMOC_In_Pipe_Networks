
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import *           
from module.transient_model import *
import pandas as pd
import pickle


var = ['HNN', 'H0', 'QNN']

# name = 'j_tsnet3_emoc'
# name = '21dpipe_oneburst'
# name = 'sevenpipe_oneburst'
# name = 'j_experiment2_6pipe_deadend'
# name = 'j_experiment2_stepwave_6pipe'
# name = 'j_experiment2_stepwave_6pipe_10s'
# name='j_experiment2_continuous'
# name = 'j_kim_leak_ga'
# name = 'enkf_5_pipe'
# name = '4_pipe_2D'
# name='6_pipes_valve2'
# name='6_pipes_valve2_xld'
# name='GA_network'
# name = '4_pipes_valve_leak_embeded'
# name = '4_pipes_valve_leak_loc2_3'
# name = 'leak_experiment_5pipe'
# name = 'leak_experiment_outside'
# name = 'leak_experiment_outside_network'
# name="leak_experiment_close_sole_network"
# name = 'leak_experiment_branchpipe_demand'
# name = 'j_eem_7pipe'
# name = 'j_eem_8pipe1_emoc'
# name = 'j_eem_11pipe_emoc'
# name = 'j_eem_13pipe5_emoc'
# name = 'j_eem_7pipe_demand_emoc'
# name = 'j_eem_7pipe_1long_emoc'
# name = 'Tnet2_sub2_v2'
# name = 'Tnet6_like'
# name = 'Tnet23'
# name = 'Tnet2_sub3'
# name="Tnet_3inline_valves"
name = 'qt_7pipe'
# name = 'Tnet5_like'
# name = 'Tnet7'
# name = '10pipe'
# name = '10pipe_KV'
# name='3_pipes_valve_0plastics'
# # name = 'Tnet2_sub_removed2'
# name = 'RMOC_2pipes'
# name = 'Tnet2'
# x_EEM = np.load('seven_pipes_one_reach.npy')
# x_EEM2 = np.load('seven_pipes_two_reach.npy')
# x_EEM3 = np.load('seven_pipes_three_reach.npy')
# x_EEM5 = np.load('seven_pipes_five_reach.npy')
# X_EEM1_wei=pd.read_excel('seven_pipes.xls').to_numpy().T
# tlen = int(x_EEM.shape[1] / 2)
# model = PipelineSystem('single_pipe')
# model = PipelineSystem('reddy')
# model = PipelineSystem('21pipe_oneburst', isRecord=True)
# model = PipelineSystem(name, isRecord=True)
# name = '3_pipes_valve_1plastics_v'
model = PipelineSystem(name, isRecord=True,id=1)
# random_numbers = [random.randint(0, 65) for _ in range(20)]
# for i,id in enumerate(random_numbers):
#     if i<5:
#         model.pipes[id].a = model.pipes[i].a*1.1
#     elif i<10:
#         model.pipes[id].a = model.pipes[i].a*0.9
    # elif i<15:
    #     model.pipes[id].D = model.pipes[i].D*1.05
    # elif i<20:
    #     model.pipes[id].D = model.pipes[i].D*0.95
# model.update_pipes()
# a_file = open(data_root_path+ name+'id'+str(model.id) + "_class.pkl", "wb")
# pickle.dump(model, a_file)
# a_file.close()
model.run( mode=MODE.MOC)
# for node in [1,2,3,4,6]:  
#     for j in range(100):
#         cda= 0.0002*j/100
#         model.T=0
#         model.inv_state_var(model.Xs)
#         model.demands[0].update(model.nodes,node,cda)
#         for i in range(model.steps-1):
#             model.transient()
model.write_recorder(filename=data_root_path+ name)
        # print(node,j)
# model.demands[0].plot_demands()
# model.write_recorder(filename='GNN_traing_raw_data/'+ name)
npipe=model.n_pipe 
t = model.pick_data('time')
fignum = 1
maxh=np.zeros(model.n_pipe)
index=range(model.n_pipe)
index=[0,1,2]
fig1 = plt.figure(fignum, figsize=(6.4, 4))
plt.subplots_adjust(left=0.135, wspace=1.5, hspace=0.5,
                    bottom=0.15, right=0.85, top=0.9)
lenPlots=len(index)
for num,i  in enumerate( index):
    plt.subplot(lenPlots,1,num+1)
    y = model.pick_data(h=-1,pipenum=i)
    maxh[i]=max(y)
    plt.plot(t, y, 'k', linewidth=1, label='small')
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Head (m)', fontsize=12)
    plt.xlim([0,model.T])
    if num==0:
        plt.legend()
    
    fignum += 1

fig1.text(0.035, 0.92, '(a)', fontsize=12)
fig1.text(0.035, 0.45, '(b)', fontsize=12)
print(maxh)
plt.savefig(figure_root_path + name + '_' + "pipe_" + '.png')
plt.show()


a=1
