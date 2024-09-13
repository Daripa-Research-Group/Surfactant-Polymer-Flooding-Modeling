import warnings
from data_reader import DataReader
from plotter import DataEntry, Plotter

# Suppress all FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

if __name__ == '__main__':
    sh_off_file = 'shear-thinning-off.xlsx'
    sh_on_file = 'shear-thinning-on.xlsx'
    
    sh_off_reader = DataReader(sh_off_file)
    sh_on_reader = DataReader(sh_on_file)
    
    # retrieving data from homogenous permeability fields
    group_data_xanthane_homogenous = [
        DataEntry(sh_off_reader.get_sheet_data('no-surfactant', 'homogenous', 'xanthane'), 'No Shear-thinning + No Surfactant'),
        DataEntry(sh_off_reader.get_sheet_data('surfactant-0.9', 'homogenous', 'xanthane'), 'No Shear-thinning + Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('no-surfactant', 'homogenous', 'xanthane'), 'Shear-thinning + No Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('surfactant-0.9', 'homogenous', 'xanthane'), 'Shear-thinning + Surfactant')
    ]
    group_data_schizophyllan_homogenous = [
        DataEntry(sh_off_reader.get_sheet_data('no-surfactant', 'homogenous', 'schizophyllan'), 'No Shear-thinning + No Surfactant'),
        DataEntry(sh_off_reader.get_sheet_data('surfactant-0.9', 'homogenous', 'schizophyllan'), 'No Shear-thinning + Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('no-surfactant', 'homogenous', 'schizophyllan'), 'Shear-thinning + No Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('surfactant-0.9', 'homogenous', 'schizophyllan'), 'Shear-thinning + Surfactant')
    ]
    homogenous_plot = Plotter(
        title='Cumulative Oil Recovery for Homogenous Permeability Field',
        x_label='Time (s)',
        y_label='Cumulative Oil Recovered',
        plot_type='line',
        save_as='figures/homogenous_cumulative_oil_recovery_vs_time.png'
    )
    homogenous_plot.create_plot([group_data_xanthane_homogenous, group_data_schizophyllan_homogenous], 'timestamp', 'coc')
    
    # retrieving data from heterogenous permeability fields
    group_data_xanthane_heterogenous = [
        DataEntry(sh_off_reader.get_sheet_data('no-surfactant', 'heterogenous', 'xanthane'), 'No Shear-thinning + No Surfactant'),
        DataEntry(sh_off_reader.get_sheet_data('surfactant-0.9', 'heterogenous', 'xanthane'), 'No Shear-thinning + Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('no-surfactant', 'heterogenous', 'xanthane'), 'Shear-thinning + No Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('surfactant-0.9', 'heterogenous', 'xanthane'), 'Shear-thinning + Surfactant')
    ]
    group_data_schizophyllan_heterogenous = [
        DataEntry(sh_off_reader.get_sheet_data('no-surfactant', 'heterogenous', 'schizophyllan'), 'No Shear-thinning + No Surfactant'),
        DataEntry(sh_off_reader.get_sheet_data('surfactant-0.9', 'heterogenous', 'schizophyllan'), 'No Shear-thinning + Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('no-surfactant', 'heterogenous', 'schizophyllan'), 'Shear-thinning + No Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('surfactant-0.9', 'heterogenous', 'schizophyllan'), 'Shear-thinning + Surfactant')
    ]
    heterogenous_plot = Plotter(
        title='Cumulative Oil Recovery for Heterogenous Permeability Field',
        x_label='Time (s)',
        y_label='Cumulative Oil Recovered',
        plot_type='line',
        save_as='figures/heterogenous_cumulative_oil_recovery_vs_time.png'
    )
    heterogenous_plot.create_plot([group_data_xanthane_heterogenous, group_data_schizophyllan_heterogenous], 'timestamp', 'coc')

    # retrieving data from quarter-five-spot permeability fields
    group_data_xanthane_quarter_five = [
        DataEntry(sh_off_reader.get_sheet_data('no-surfactant', 'quarter-five-spot', 'xanthane'), 'No Shear-thinning + No Surfactant'),
        DataEntry(sh_off_reader.get_sheet_data('surfactant-0.9', 'quarter-five-spot', 'xanthane'), 'No Shear-thinning + Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('no-surfactant', 'quarter-five-spot', 'xanthane'), 'Shear-thinning + No Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('surfactant-0.9', 'quarter-five-spot', 'xanthane'), 'Shear-thinning + Surfactant')
    ]
    group_data_schizophyllan_quarter_five = [
        DataEntry(sh_off_reader.get_sheet_data('no-surfactant', 'quarter-five-spot', 'schizophyllan'), 'No Shear-thinning + No Surfactant'),
        DataEntry(sh_off_reader.get_sheet_data('surfactant-0.9', 'quarter-five-spot', 'schizophyllan'), 'No Shear-thinning + Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('no-surfactant', 'quarter-five-spot', 'schizophyllan'), 'Shear-thinning + No Surfactant'),
        DataEntry(sh_on_reader.get_sheet_data('surfactant-0.9', 'quarter-five-spot', 'schizophyllan'), 'Shear-thinning + Surfactant')
    ]
    heterogenous_plot = Plotter(
        title='Cumulative Oil Recovery for Quarter-Five Spot Permeability Field',
        x_label='Time (s)',
        y_label='Cumulative Oil Recovered',
        plot_type='line',
        save_as='figures/quarter-five_cumulative_oil_recovery_vs_time.png'
    )
    heterogenous_plot.create_plot([group_data_xanthane_quarter_five, group_data_schizophyllan_quarter_five], 'timestamp', 'coc')