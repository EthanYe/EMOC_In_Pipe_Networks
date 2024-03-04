import numpy as np
filename='Tnet_3inline_valves'
data_path='data/'
nodes = np.load('data/'+ filename + '_nodes_sorted.npy')
pipes = np.load('data/'+filename + '_pipes_sorted.npy')

node_epa2pip = {}
pipe_epa2pip = {}

for i in range(len(nodes)):
    node_epa2pip[nodes[i]]=i
for i in range(len(pipes)):
    pipe_epa2pip[pipes[i]]=i

np.save('data/'+filename+'_node_epa2pip.npy',node_epa2pip)
np.save('data/'+filename+'_pipe_epa2pip.npy',pipe_epa2pip)

a=1