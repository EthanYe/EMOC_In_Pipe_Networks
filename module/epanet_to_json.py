import pandas as pd
import numpy as np
filename='Tnet_3inline_valves_sub'
import json
data_path = 'data/'
config_path = 'config/'
with open(data_path+filename + '.inp', 'r+') as f:
    a=f.readlines()
    line_num=len(a)
    for i in range(line_num):
        j=i+2
        if a[i][:11]=='[JUNCTIONS]':
            junction_line_i=j
        elif a[i][:12]=='[RESERVOIRS]':
            reservoir_line_i=j
        elif a[i][:7]=='[TANKS]':
            tank_line_i=j
        elif a[i][:7]=='[PIPES]':
            pipe_line_i=j
        elif a[i][:8]=='[VALVES]':
            valve_line_i=j
        elif a[i][:7]=='[PUMPS]':
            pump_line_i=j
        elif a[i][:6]=='[TAGS]':
            tag_line_i=j
        elif a[i][:5]=='[COOR':
            coor_line_i=j
        elif a[i][:6] == '[VERTI':
            veti_line_i=j
    junction_num=reservoir_line_i-junction_line_i-3
    reservoir_num=tank_line_i-reservoir_line_i-3
    tank_num=pipe_line_i-tank_line_i-3
    pipe_num=pump_line_i-pipe_line_i-3
    pump_num=valve_line_i-pump_line_i-3
    valve_num=tag_line_i-valve_line_i-3
    coor_num = veti_line_i-coor_line_i-3
    # print(a[junction_line_i][:10],a[junction_line_i+junction_num-1][:10])
    # print(junction_line_i,reservoir_line_i,tank_line_i,pipe_line_i,valve_line_i,pump_line_i)
"""Create json data"""
json_dict={} # Create json dict
"""Create dict for control parameters"""
cp_index="control parameters"
cp_value={
        "time step": 0.1,
        "total time": 50,
        "gravity": 9.81,
        "friction factor": "f",
        "unit": "US"
    } # Any other parameters can be added


"""pipe info"""
pipes=[]
nodes_id=[]
pipe_info=[]
for i in range(pipe_line_i,pipe_line_i+pipe_num):
    pipe=a[i].split('\t')
    pipe_info+=[pipe]
    nodes_id+=[int(pipe[1])]
    nodes_id+=[int(pipe[2])]
"""Merge nodes of valves and pumps, and modify pipe info"""
nodes_vp=[]
modified_vp=[]
for i in range(valve_line_i,valve_line_i+valve_num):
    valve=a[i].split('\t')
    node_1=int( valve[1])
    node_2=int( valve[2])
    if node_1 in nodes_id and node_2 in nodes_id:
        # merge the two nodes to one node
        for i,info in enumerate(pipe_info):
            if int(info[1])==node_2:
                pipe_info[i][1]=node_1
            if int(info[2])==node_2:
                pipe_info[i][2]=node_1
        nodes_vp+=[node_2]
        modified_vp+=[node_1]
    # if node_1 in nodes_id and node_2 not in nodes_id:
    #     # merge the two nodes to one node
    #     for i,info in enumerate(pipe_info):
    #         if int(info[1])==node_2:
    #             pipe_info[i][1]=node_1
    #         if int(info[2])==node_2:
    #             pipe_info[i][2]=node_1
    #     nodes_vp+=[node_2]
    #     modified_vp+=[node_1]
for i in range(pump_line_i,pump_line_i+pump_num):
    pump=a[i].split('\t')
    node_1=int( pump[1])
    node_2=int( pump[2])
    if node_1 in nodes_id or node_2 in nodes_id:
        # merge the two nodes to one node
        for i,info in enumerate(pipe_info):
            if int(info[1])==node_2:
                pipe_info[i][1]=node_1
            if int(info[2])==node_2:
                pipe_info[i][2]=node_1
        nodes_vp+=[node_2]
        modified_vp+=[node_1]        
"""sort nodes"""
nodes_id=[]
pipes_id=[]
for pipe in pipe_info:
    pipes_id+=[int(pipe[0])]
    nodes_id+=[int(pipe[1])]
    nodes_id+=[int(pipe[2])]
nodes_id_sorted=np.sort(np.array(list(set(nodes_id)),dtype='int'))
pipes_id_sorted=np.array(pipes_id,dtype='int')
np.save(data_path+filename+'_pipes_sorted.npy',pipes_id_sorted)
np.save(data_path+filename+'_nodes_sorted.npy',nodes_id_sorted)
"""Create list for pipes"""
for i,pipe in enumerate(pipe_info):
    node_1,node_2=int(pipe[1]),int(pipe[2])
    n,m= np.where(nodes_id_sorted == node_1)[0],np.where(nodes_id_sorted == node_2)[0]
    pipes.append({"id":i,"js":int(n),"je":int(m),"wavespeed":1000,"diameter":float(pipe[4]),"length":float(pipe[3]),"f":0.022})

"""Create list for reservoirs"""
reservoirs=[]

for i in range(reservoir_line_i,reservoir_line_i+reservoir_num):
    reservoir=a[i].split('\t')
    r_node=int( reservoir[0])
    for node,modified_node in  zip(nodes_vp,modified_vp):
        if int( reservoir[0]) ==node:
            r_node=modified_node
    reservoirs.append({"id":len(reservoirs),"node":int(int( np.where(nodes_id_sorted ==r_node)[0])),"water level":float(reservoir[1])})


"""Create list for valves"""
inline_valves=[]
end_valves=[]
for i in range(valve_line_i,valve_line_i+valve_num):
    valve=a[i].split('\t')
    node_1=int( np.where(nodes_id_sorted == int(valve[1]))[0])
    # node_2=int(np.where(nodes_id_sorted == int(valve[2]))[0])
    # if the two nodes are connected to pipe, valve is inline valve. 
    # if only one node is connected to a pipe, than the valve is an end valve 
    if int(valve[2]) in nodes_vp:
        inline_valves.append({"id":len(inline_valves),"node":node_1,"motion":'sudden','Q0':0.1,"status":"open"})
    elif int(valve[2]) not in nodes_vp:
        end_valves.append({"id":len(end_valves),"node":node_1,"motion":'sudden','Q0':0.1})
    else:
        print("no pipe is connected to valve "+valve[0])
"""Create list for coordinate"""
coords=[]
for i in range(coor_line_i,coor_line_i+coor_num):
    coord=a[i].split('\t')
    if int(coord[0]) not in nodes_id_sorted:
        continue
    node=int( np.where(nodes_id_sorted == int(coord[0]))[0])
    coords.append({"node":node,"epa_node":int(coord[0]),'X-coord':float(coord[1]), 'Y-coord':float(coord[2])})



json_dict[cp_index]=cp_value # add dict for control parameters into json dict
json_dict["pipes"]=pipes # add dict for pipes into json dict
if len(reservoirs):
    json_dict["reservoirs"]=reservoirs # add dict for pipes into json dict
if len(inline_valves):
    json_dict["inline valve"]=inline_valves # add dict for pipes into json dict
if len(end_valves):
    json_dict["end valve"]=end_valves # add dict for pipes into json dict
if len(coords):
    json_dict["coordinates"] = coords  # add dict for pipes into json dict
with open(config_path+filename + '_config.json', 'w') as f:
    json.dump(json_dict, f)

# generate dict from epa nodes/pipes to pipsys nodes/pipes
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
# file.write({cp_index:cp_value})
# file.close()
