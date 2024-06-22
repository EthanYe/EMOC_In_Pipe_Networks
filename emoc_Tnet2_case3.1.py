from module.EMOC import *
def subnode2pipe(subnodes,sub_model):
    sp_dict={}
    ep_dict={}
    for node in range(subnodes):
        epa_node=sub_model.nodes[node].epa_node
        pip_sensor_node=node_epa2pip[epa_node]
        for pipe in data_model.pipes:
            if pipe.js==pip_sensor_node:
                sp_dict[node]=pipe.id
                break
            if pipe.je==pip_sensor_node:
                ep_dict[node]=pipe.id
                break
    return sp_dict, ep_dict

# name of the saved pipe network model
# name = 'leak_experiment_network'
mode='delta'
name = 'Tnet2_sub3id1'
sensor_data_name='Tnet23id1'
# node_epa2subpip = np.load(data_root_path+ 'Tnet2_node_epa2subpip.npy',allow_pickle=True).item() 
# node_epa2pip = np.load(data_root_path+ 'Tnet2_node_epa2pip.npy',allow_pickle=True).item()   
# pipe_epa2subpip = np.load(data_root_path+ 'Tnet2_pipe_epa2pip.npy',allow_pickle=True).item()  
# pipe_epa2pip = {}
node_epa2subpip = np.load(data_root_path+ 'Tnet2_sub3_node_epa2subpip.npy',allow_pickle=True).item() 
node_epa2pip = np.load(data_root_path+ 'Tnet23_node_epa2pip.npy',allow_pickle=True).item()   
pipe_epa2subpip = np.load(data_root_path+ 'Tnet23_pipe_epa2pip.npy',allow_pickle=True).item()  
# pipe_epa2subpip = {}
sensors_epa_nodes=[204,189,265,1015] #[204,1019,265,1015],[204,189,265,1015],[204,193,265,1015]
unknowns_epa_nodes=[197,113,159,1017]
sensors=[]
raw_sensors=[]
sensor_data=np.zeros((len(sensors_epa_nodes),500))
for sensor_epa in sensors_epa_nodes:
    sensors.append(node_epa2subpip[sensor_epa])
    raw_sensors.append(node_epa2pip[sensor_epa])
loaded_data=np.load(data_root_path + sensor_data_name + '_nodes_data.npz')
flag=0
nodeData=[]

# Valid path 1: [[[22, 16, 18, 21], [45, 46, 47, 48, 52, 14, 15, 13, 17, 24, 19, 1], [23, 4, 3, 2], [42, 44]]]
# Valid path 1: [[[22, 16, 18, 21], [45, 46, 47, 48, 52, 14, 15, 13, 17, 24, 19, 1], [23, 4, 3, 2], [42, 44]]]
for i in range(len(loaded_data)):
    
    array_name = f'node{i}'
    loaded_array = loaded_data[array_name]
    nodeData.append(loaded_array)
    if i in raw_sensors:
        sensor_data[raw_sensors.index(i)]=loaded_array
        flag+=1
# sensors=[3,19,18,39] # epa_node: 161, 193, 191, 1012 ,pipsys_node_entire:[33,49,48,75] pipsys_pipe_entire:[13,58,100,70]
# test_points=[54,29] #1028,1002 epa_node_entire:[90,65][92,81]
# extra_sensor=[3]
# define the location of unknown boundaries
unknowns=[]
for unknown_epa in unknowns_epa_nodes:
    unknowns.append(node_epa2subpip[unknown_epa])
# sensors=[1,2,21,44]
# unknowns=[1,2,21,44]

# # read sensor data
# a_file = open(data_root_path + name + "_class.pkl", "rb")

# data_model= pickle.load(a_file)
mode=''

mode='delta'

'''EMOC model'''
emoc = EMOC(name=name, sensors=sensors, sensor_data=sensor_data,unknowns=unknowns,mode=mode)
# for pipe in emoc.model.pipes:
#     if pipe.a==1100:
#         pass
#         # print(pipe.start_node.epa_node,pipe.end_node.epa_node)
#     elif pipe.a==900:
#         print(pipe.start_node.epa_node,pipe.end_node.epa_node)
emoc.emoc()
dataIndexs=[]
for i,node in enumerate(emoc.model.nodes):
    epa_index=node.epa_node
    dataIndexs.append(node_epa2pip[epa_index])
pipes = emoc.pipes
max_steps = emoc.plotSteps
# plot
plt_nodes=[]
fignum=1
errors=[]

# for node in emoc.nodes:
#     id=dataIndexs[node.id]
#     plt.figure(fignum)
    
#     if mode=='delta':
#         plt.plot(nodeData[id][:max_steps]-nodeData[id][ 0], label='Measured')
#         y1=nodeData[id][:max_steps]-nodeData[id][ 0] 
#     else:
#         plt.plot(nodeData[id][:max_steps], label='Measured')
#     for pipe in pipes:
#         if pipe.js==node.id:
#             plt.plot(pipe.hi[:max_steps, 0], '--', label='Predicted')
#             y2=pipe.hi[:max_steps, 0]
#         elif pipe.je==node.id:
#             plt.plot(pipe.hi[:max_steps, -1], '--', label='Predicted')
#             y2=pipe.hi[:max_steps, -1]
#     error= np.mean(np.square((y1-y2)))
#     errors.append(error)
#     fignum += 1
#     plt.xlabel('Time (s)')
#     plt.ylabel('Pressure head (m)')
#     # plt.title('Node '+str(js))
#     plt.legend(fontsize=16)
#     # plt.savefig('../figures/' + name + 'pressure_' + str(node) + '.png', dpi=300)
# print(errors)

# np.save('errors.npy',np.array(errors))

# plt.figure(fignum)
# # index = np.linspace(0,)
# # values=[0.54,0.5,0.68,4.30]
# # width=0.5
# # p2 = plt.bar(index, values, width, color="#87CEFA", label="hosflkksfjlsf")
# plt.legend()
# plt.show()
for pipe in pipes:
    node=pipe.js
    id=dataIndexs[node]
    plt.figure(fignum)
    if mode=='delta':
        plt.plot(nodeData[id][:max_steps]-nodeData[id][ 0], label='Measured') 
    else:
        plt.plot(nodeData[id][:max_steps], label='Measured')
    plt.plot(pipe.hi[:max_steps, 0], '--', label='Predicted')
    fignum += 1
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure head (m)')
    # plt.title('Node '+str(js))
    plt.legend(fontsize=16)
    plt.savefig('../figures/' + name + 'pressure_' + str(node) + '.png', dpi=300)
plt.legend()
plt.show()

plt_nodes=[]
fignum=1
data=emoc.data
for pipe in pipes:
    id=pipe.id
    plt.figure(fignum)
    if mode=='delta':
            plt.plot(data[id][:max_steps, 0]-data[id][0, 0], label='Measured') 
    else:
        plt.plot(data[id][:max_steps, 0], label='Measured')
    plt.plot(pipe.hi[:max_steps, 0], '--', label='Predicted')
    fignum += 1
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure head (m)')
    # plt.title('Node '+str(js))
    plt.legend(fontsize=16)
plt.legend()
plt.show()