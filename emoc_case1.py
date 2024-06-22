from module.EMOC import *
import matplotlib.pyplot as plt
# name of the saved pipe network model

name = 'emoc_paper_case1'
# define the location of sensors
sensors=[2,3] # two sensors at node 2 and node 3
# define the location of unknown boundaries
unknowns=[0,4] # two unknown boundaries at node 0 and node 4
# define the name to save results
save_name='sensor_'+str(sensors)
# define the method to calculate initial condition

'''EMOC model'''
emoc = EMOC(name=name, sensors=sensors,unknowns=unknowns)
t=emoc.model.t
# MOC results, true values
data = emoc.data
# EMOC results are stored in pipe.hi and pipe.qi
# plot sensor data
emoc.emoc()
pipes = emoc.pipes
emoc.save_results(addname='')
# plot pressure head at unknown boundaries
max_steps = emoc.plotSteps
fig1=plt.figure(1,figsize=(8,4))
plt.subplots_adjust(left=0.135, wspace=1.5, hspace=0.35,
                    bottom=0.15, right=0.95, top=0.98)
plt.subplot(211)
plt.plot(t[:max_steps],data[0][:max_steps, -2], marker='o', markersize=5, markevery=np.arange(0, 985, 100),linewidth=1, label='True')
plt.plot(t[:max_steps],pipes[0].hi[:max_steps, -1],marker='D' , markersize=5, markevery=np.arange(30, 985, 100),linewidth=1,  label='Predicted')
plt.xlim([0,10])
plt.legend(ncol=2,loc=1, fontsize=11,frameon=False)
plt.xlabel('Time (s)',fontsize=12,)
plt.ylabel('Head (m)',fontsize=12,)
plt.subplot(212)
plt.plot(t[:max_steps],data[3][:max_steps, -2], marker='o', markersize=5, markevery=np.arange(0, 985, 100),linewidth=1, label='True')
plt.plot(t[:max_steps],pipes[3].hi[:max_steps, -1],marker='D' , markersize=5, markevery=np.arange(30, 985, 100),linewidth=1,  label='Predicted')
plt.xlim([0,10])
plt.xlabel('Time (s)',fontsize=12)
plt.ylabel('Head (m)',fontsize=12)
plt.legend(ncol=2,loc=1, fontsize=11,frameon=False)
fig1.text(0.03, 0.95, '(a)', fontsize=12)
fig1.text(0.03, 0.46, '(b)', fontsize=12)
plt.show()
