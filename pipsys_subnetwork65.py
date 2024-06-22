import matplotlib.pyplot as plt
from numpy import *           
from module.transient_model import *
import pickle

'''build pipeline model from the configuration file and run simulation'''
# specify the name of the configuration file
name = 'subnetwork-65'
# build the transient model
model = PipelineSystem(name)
# save the model
a_file = open(data_root_path+ name+"_class.pkl", "wb")
pickle.dump(model, a_file)
a_file.close()