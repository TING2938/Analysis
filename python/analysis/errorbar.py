import matplotlib.pyplot as plt
from matplotlib import container
import pandas as pd
#%%
allcols = {'usecols':'C:AE', 'sheet_name':4, 'skiprows':0, 'nrows':63}

df = pd.read_excel(
        "file:///D:/Private/OneDrive - hust.edu.cn/Private/Gromacs/ILs/BmimBF4_SOL_LiBF4/BmimBF4.xlsx",
        header=None, **allcols)
#%%
x = [0, 0.5, 1, 2, 4]
data = ( {'x':x, 'y':[], 'err':[], 'title':'Li-Bmim'}, 
         {'x':x, 'y':[], 'err':[], 'title':'Li-BF4'}, 
         {'x':x, 'y':[], 'err':[], 'title':'SOL-Bmim'}, 
         {'x':x, 'y':[], 'err':[], 'title':'SOL-BF4'},
         {'x':x, 'y':[], 'err':[], 'title':'SOL-SOL'}, 
         {'x':x, 'y':[], 'err':[], 'title':'SOL-Li'}, 
         {'x':x, 'y':[], 'err':[], 'title':'Li-SOL'})

for i in range(len(df)):
    if df.iat[i, 0] == 'Average':
        j = 3
        for item in data:
            item['y'].append( (df.iloc[i-3:i, j] / df.iloc[i-3:i, j+1]).mean() )
            item['err'].append( (df.iloc[i-3:i, j] / df.iloc[i-3:i, j+1]).std() )
            #if j == 19: j -= 2
            j += 4
#%%  
def errbar(ax, data, legend=False, xlabel=False, ylabel=False, title=True, error=True):  
    
    if xlabel:
        ax.set_xlabel("Li : SOL")
    if ylabel:
        ax.set_ylabel("N (#)")
    if title:
        ax.set_title(data["title"])
    ax.set_xticks([0,1,2,3,4])

    if error:
        ax.errorbar(data['x'], data['y'][0:5], yerr=data['err'][0:5], fmt='-o', label='2489 ppm')
        ax.errorbar(data['x'], data['y'][5:10], yerr=data['err'][5:10], fmt='--o', label='5546 ppm')  
        ax.errorbar(data['x'], data['y'][10:15], yerr=data['err'][10:15], fmt=':o', label='8008 ppm')
    else:
        ax.plot(data['x'], data['y'][0:5], '-o', label='2489 ppm')
        ax.plot(data['x'], data['y'][5:10], '-.o', label='5546 ppm')  
        ax.plot(data['x'], data['y'][10:15], ':o', label='8008 ppm')
        
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
        #ax.legend(handles, labels, fontsize=16)
        return handles, labels


#%% plot
fig, ax = plt.subplots(ncols=4, nrows=2, sharex=True, figsize=(14, 7))
ax = ax.flatten()

errbar(ax[0], data[0], ylabel=True)
errbar(ax[1], data[1])
errbar(ax[2], data[2])
errbar(ax[3], data[3])
errbar(ax[4], data[4], xlabel=True, ylabel=True)
errbar(ax[5], data[5], xlabel=True)
leg = errbar(ax[6], data[6], xlabel=True, legend=True)
ax[7].legend(*leg, frameon=False, fontsize=16, loc=6)
ax[7].axis("off")
#fig.text(0.806, 0.2, "$Energy(VDM+COU)$")

#fig.savefig("Nr.png", dpi=600)
fig.show()