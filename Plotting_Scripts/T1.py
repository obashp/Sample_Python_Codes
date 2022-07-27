import matplotlib.pyplot as plt
import matplotlib as mpl

fig, ax = plt.subplots(figsize=(6, 1))
fig.subplots_adjust(right=0.2)

cmap = mpl.cm.cool
norm = mpl.colors.Normalize(vmin=5, vmax=10)


fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax, orientation='vertical', label='Some Units')
fig.add_subplot(figsize = (5,1))
