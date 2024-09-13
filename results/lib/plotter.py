import matplotlib.pyplot as plt

class DataEntry:
    def __init__(self, data, label):
        self.data = data
        self.label = label
        
class Plotter:
    def __init__(self, title, x_label, y_label, plot_type='line', legend_position='lower center', ncol=4, save_as='plot.png'):
        self.title = title
        self.x_label = x_label
        self.y_label = y_label
        self.plot_type = plot_type
        self.legend_position = legend_position
        self.ncol = ncol
        self.save_as = save_as
        
    def _plot_group(self, ax, group, x_axis_data_type, y_axis_data_type, linestyle_colors):
        for entry, (linestyle, color) in zip(group, linestyle_colors):
            if self.plot_type == 'line':
                ax.plot(entry.data[x_axis_data_type], entry.data[y_axis_data_type], label=entry.label, linestyle=linestyle, color=color)
            elif self.plot_type == 'scatter':
                ax.scatter(entry.data[x_axis_data_type], entry.data[y_axis_data_type], label=entry.label, color=color)
            elif self.plot_type == 'bar':
                ax.bar(entry.data[x_axis_data_type], entry.data[y_axis_data_type], label=entry.label, color=color)
        
    def create_plot(self, groups, x_axis_data_type, y_axis_data_type):
        fig, axs = plt.subplots(1, len(groups), figsize=(25, 12), sharex=True)
        if len(groups) == 1:
            axs = [axs]
            
        linestyle_colors = [('-', 'red'), ('--', 'green'), ('-.', 'black'), (':', 'blue')]
        
        for i, group in enumerate(groups):
            self._plot_group(axs[i], group, x_axis_data_type, y_axis_data_type, linestyle_colors)
            axs[i].set_xlabel(self.x_label)
            axs[i].set_ylabel(self.y_label)
            
        fig.suptitle(self.title, fontsize=20)
        
        # single legend for both subplots
        handles = [axs[0].lines[i] for i in range(len(groups[0]))] 
        labels = [handle.get_label() for handle in handles]
        fig.legend(handles, labels, loc=self.legend_position, bbox_to_anchor=(0.5, -0.05), ncol=self.ncol, fontsize='17', borderaxespad=3)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.1) 
        plt.savefig(self.save_as)