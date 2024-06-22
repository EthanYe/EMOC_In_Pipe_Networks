from module.EMOC import *

name = 'subnetwork-65'
sensor_data_name='network-112Model_1'

sensors=[22,17,23,42] # sensor node (in code) + 1 = sensor (in paper) 
unknowns=[21,1,2,44]  # unknown node (in code) + 1 = unknown node (in paper) 
sensors_112=[52,47,57,78] # sensor node id in 112 network
sensor_data=np.zeros((len(sensors),500))
ini_data=[]
loaded_data=np.load(data_root_path + sensor_data_name + '_nodes_data.npz')
loaded_inidata=np.load(data_root_path + sensor_data_name + '_pipes_data.npz')
# read sensor data
nodeData=[]
for i in range(len(loaded_data)): 
    array_name = f'node{i}'
    loaded_array = loaded_data[array_name]
    nodeData.append(loaded_array)
    if i in sensors_112:
        sensor_data[sensors_112.index(i)]=loaded_array

pipeIDs=np.load('data/pipeIDs.npy') # load the 64 pipe IDs in 112 network
for i,pipeID in enumerate(pipeIDs):
    array_name = f'pipe{pipeID}'
    loaded_array = loaded_inidata[array_name]
    ini_data.append(loaded_array[0,1:])
# mode='delta'

'''EMOC model'''
emoc = EMOC(name=name, sensors=sensors, sensor_data=sensor_data,unknowns=unknowns,ini_data=ini_data)
emoc.emoc()
dataIndexs=np.load('data/nodeIDs.npy')
pipes = emoc.pipes
max_steps = emoc.plotSteps
# plot results
plt_nodes=[]
fignum=1
errors=[]
emoc.save_results()
for node in emoc.nodes:
    id=dataIndexs[node.id]
    plt.figure(fignum)
    plt.plot(nodeData[id][:max_steps], label='Measured')
    y1=nodeData[id][:max_steps]
    for pipe in pipes:
        if pipe.js==node.id:
            plt.plot(pipe.hi[:max_steps, 0], '--', label='Predicted')
            break
        elif pipe.je==node.id:
            plt.plot(pipe.hi[:max_steps, -1], '--', label='Predicted')
            break
    fignum += 1
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure head (m)')
    plt.title(f'Node{node.id}')
    plt.legend(fontsize=16)
    plt.savefig('figures/' + name + 'pressure_' + str(node.id) + '.png', dpi=300)
plt.show()
