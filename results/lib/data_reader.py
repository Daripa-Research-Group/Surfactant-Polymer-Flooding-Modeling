import pandas as pd

class DataReader:
    def __init__(self, file_name):
        self.file_name = file_name
        
    def get_sheet_data(self, sheet, field_name, polymer_type):
        file = pd.read_excel(self.file_name, sheet)
        timestamp = []
        coc = []
        mfw = []
        
        # finding start columns
        timestamp_col = 0
        coc_col = 1
        mfw_col = 2
        if field_name == 'heterogenous':
            timestamp_col += 7
            coc_col += 7
            mfw_col += 7
        elif field_name == 'quarter-five-spot':
            timestamp_col += 14
            coc_col += 14
            mfw_col += 14
            
        if polymer_type != 'xanthane':
            timestamp_col += 3
            coc_col += 3
            mfw_col += 3
            
        cur_row = 2
        
        # Ensure columns are within bounds
        max_columns = file.shape[1]
        if timestamp_col >= max_columns or coc_col >= max_columns or mfw_col >= max_columns:
            raise IndexError("Column index exceeds available columns in the sheet")
        
        
        max_rows = file.shape[0]
        while cur_row < max_rows:
            if (pd.isna(file.iloc[cur_row, timestamp_col]) or 
            pd.isna(file.iloc[cur_row, coc_col]) or 
            pd.isna(file.iloc[cur_row, mfw_col])):
                break 
            
            timestamp.append(file.iloc[cur_row][timestamp_col])
            coc.append(file.iloc[cur_row][coc_col])
            mfw.append(file.iloc[cur_row][mfw_col])
            cur_row += 1

        return dict(timestamp=timestamp, coc=coc, mfw=mfw)