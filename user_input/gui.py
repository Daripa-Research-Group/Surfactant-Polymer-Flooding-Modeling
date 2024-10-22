import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd

class UserInputGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("User Input GUI")
        self.master.geometry("800x600")

        self.main_frame = tk.Frame(self.master)
        self.main_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(self.main_frame)
        self.scrollbar = tk.Scrollbar(self.main_frame, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = tk.Frame(self.canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.inputs = []
        self.row_count = 0

        self.create_table_headers()
        self.add_simulation()
        self.create_buttons()
        
        self.all_inputs = None

    def create_table_headers(self):
        headers = ["Simulation", "Figure Number", "Model Type", "C0iter", "G0iter", "Plot Type", 
                   "Polymer", "Permeability Flag", "Alpha Value", 
                   "Nsim", "Todd Longstaff Mixing Parameter", "Shear Flag"]
        tk.Label(self.scrollable_frame, text="Simulation Parameters", font=("Arial", 14, "bold")).grid(row=0, column=0, padx=10, pady=10, columnspan=12)

        for col, header in enumerate(headers):
            label = tk.Label(self.scrollable_frame, text=header, font=("Arial", 10, "bold"))
            label.grid(row=1, column=col, padx=5, pady=5, sticky="nsew")
            self.scrollable_frame.columnconfigure(col, weight=1)

    def add_simulation(self):
        row = self.row_count + 2
        input_vars = {}

        simulation_label = tk.Label(self.scrollable_frame, text=str(self.row_count + 1))
        simulation_label.grid(row=row, column=0, padx=5, pady=5, sticky="nsew")

        for col, header in enumerate([
            "figure_number", "model_type", "c0iter", "g0iter", "plot_type",
            "polymer", "permeability_flag", "alpha_value",
            "nsim", "Todd_Longstaff_mixing_parameter", "shear_flag"]):
            var = tk.DoubleVar() if "iter" in header or "parameter" in header else tk.IntVar()
            input_vars[header] = var
            
            entry = tk.Entry(self.scrollable_frame, textvariable=var, width=10)
            entry.grid(row=row, column=col + 1, padx=5, pady=5, sticky="ew")
            self.scrollable_frame.columnconfigure(col + 1, weight=1)

        self.inputs.append(input_vars)
        self.row_count += 1
        self.scrollable_frame.rowconfigure(row, weight=1)

    def create_buttons(self):
        button_frame = tk.Frame(self.master)
        button_frame.pack(side=tk.BOTTOM, fill=tk.X)

        tk.Button(button_frame, text="Submit All", command=self.submit_all).pack(side=tk.LEFT, padx=10, pady=10)
        tk.Button(button_frame, text="Add Simulation", command=self.add_simulation).pack(side=tk.LEFT, padx=10)
        tk.Button(button_frame, text="Upload CSV", command=self.load_csv).pack(side=tk.LEFT, padx=10)

        for button in button_frame.winfo_children():
            button.pack(side=tk.LEFT, padx=(10, 0), pady=10)

    def load_csv(self):
        file_path = filedialog.askopenfilename(title="Select a CSV file", filetypes=[("CSV files", "*.csv")])
        if file_path:
            try:
                data = pd.read_csv(file_path)
                messagebox.showinfo("CSV Loaded", data.head().to_string())
                print(data)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load CSV file: {e}")

    def submit_all(self):
        all_inputs = []
        for input_vars in self.inputs:
            input_values = {key: var.get() for key, var in input_vars.items()}
            all_inputs.append(input_values)

        # print("All Simulations:")
        # print(all_inputs)

        if all(all_inputs):
            print("All inputs received:")
            self.all_inputs = all_inputs
            self.master.destroy()
        else:
            messagebox.showwarning("Input Error", "Please fill in all fields before submitting.")
            
    def get_input(self):
        return self.all_inputs

