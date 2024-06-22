
import numpy as np



def shortest_paths(weighted_adjacent_matrix, sensors, unknowns):
    A = weighted_adjacent_matrix.copy()  # copy weighted_adjacent_matrix to store the distance of shortest paths
    N=weighted_adjacent_matrix.shape[0] # number of vertices
    shorest_path_list = [[[i,j] for j in range(N)] for i in range(N)] # list of shorest paths between vertices 
    # get shorest paths by Floyd
    # A represents the distance, shorest_path_list represents the paths
    for k in range(N):
        for i in range(N):
            for j in range(N):
                if A[i, j]>A[i,k]+A[k,j]:
                    A[i,j] =A[i,k]+A[k,j]
                    shorest_path_list[i][j] = shorest_path_list[i][k][:-1] + shorest_path_list[k][j]
    # get the full permutation of sensor arrays
    sensor_permutation = permute(sensors)
    # match the full permutation of sensors with unknown boundarys
    sus_combination_list=[]
    for i,sensors in enumerate (sensor_permutation):
        sus_combination=[]
        for j,sensor in enumerate( sensors):
            sus_combination.append([sensor,unknowns[j]])
        sus_combination_list.append(sus_combination)
    # get pahts that do not overlap
    valid_path=[]
    for sus in sus_combination_list:
        pipe =set()
        paths=[]
        valid=True
        for su in sus:
            path=shorest_path_list[su[0]][su[1]]
            paths.append(path)
            new_set = set(path)
            if pipe.isdisjoint(new_set):
                pipe=pipe|new_set
            else: 
                valid=False
                break
        if valid:
            valid_path.append(paths)
            print('Valid path '+ str(len(valid_path))+':',valid_path)
    if len(valid_path):
        return valid_path[0]
    else:
        print('No valid path!')
        exit(0)