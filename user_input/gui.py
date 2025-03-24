import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd


class UserInputGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("User Input GUI")
        self.master.geometry("1050x600")

        self.main_frame = tk.Frame(self.master)
        self.main_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(self.main_frame)
        self.scrollbar = tk.Scrollbar(
            self.main_frame, orient="vertical", command=self.canvas.yview
        )
        self.scrollable_frame = tk.Frame(self.canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")),
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

    def create_table_headers(self) -> None:
        headers = [
            "Simulation",
            "Model",
            "Geometry",
            "Permeability",
            "Polymer Type",
            "Polymer Concentration",
            "Surfactant Type",
            "Surfactant Concentration",
        ]
        tk.Label(
            self.scrollable_frame,
            text="Simulation Parameters",
            font=("Arial", 14, "bold"),
        ).grid(row=0, column=0, padx=10, pady=10, columnspan=len(headers))

        for col, header in enumerate(headers):
            label = tk.Label(
                self.scrollable_frame, text=header, font=("Arial", 10, "bold")
            )
            label.grid(row=1, column=col, padx=5, pady=5, sticky="nsew")
            self.scrollable_frame.columnconfigure(col, weight=1)

    def add_simulation(self) -> None:
        row = self.row_count + 2
        input_vars = {}

        simulation_label = tk.Label(self.scrollable_frame, text=str(self.row_count + 1))
        simulation_label.grid(row=row, column=0, padx=5, pady=5, sticky="nsew")

        # Model dropdown
        model_options = {
            "No Shear Thinning": 1,
            "Sourav Implementation": 2,
            "Shear Thinning": 3,
        }
        model_var = tk.StringVar(value=list(model_options.keys())[0])
        input_vars["model_type"] = model_var

        model_menu = tk.OptionMenu(
            self.scrollable_frame, model_var, *model_options.keys()
        )
        model_menu.config(width=18)
        model_menu.grid(row=row, column=1, padx=5, pady=5, sticky="ew")
        self.scrollable_frame.columnconfigure(1, weight=1)
        
        # Geometry dropdown
        reservoir_geometry_options = {
            "Rectilinear": 1,
            "Quarter Five Spot": 2,
        }
        reservoir_geometry_var = tk.StringVar(value=list(reservoir_geometry_options.keys())[0])
        input_vars["reservoir_geometry"] = reservoir_geometry_var
         
        geometry_menu = tk.OptionMenu(
            self.scrollable_frame, reservoir_geometry_var, *reservoir_geometry_options.keys()
        )
        geometry_menu.config(width=18)
        geometry_menu.grid(row=row, column=2, padx=5, pady=5, sticky="ew")
        self.scrollable_frame.columnconfigure(2, weight=1)
         

        # Permeability dropdown
        permeability_options = {
            "Homogeneous": 1,
            "Heterogeneous": 2,
        }
        permeability_var = tk.StringVar(value=list(permeability_options.keys())[0])
        input_vars["permeability"] = permeability_var

        perm_menu = tk.OptionMenu(
            self.scrollable_frame, permeability_var, *permeability_options.keys()
        )
        perm_menu.config(width=15)
        perm_menu.grid(row=row, column=3, padx=5, pady=5, sticky="ew")
        self.scrollable_frame.columnconfigure(3, weight=1)

        # Polymer Type dropdown
        polymer_type_options = {
            "Xanthane": 1,
            "Schizophyllan": 2,
            "No Polymer": 3,
        }
        polymer_type_var = tk.StringVar(value=list(polymer_type_options.keys())[0])
        input_vars["polymer_type"] = polymer_type_var

        polymer_menu = tk.OptionMenu(
            self.scrollable_frame, polymer_type_var, *polymer_type_options.keys()
        )
        polymer_menu.config(width=15)
        polymer_menu.grid(row=row, column=4, padx=5, pady=5, sticky="ew")
        self.scrollable_frame.columnconfigure(4, weight=1)

        # Polymer Concentration Entry
        poly_conc_var = tk.DoubleVar()
        input_vars["polymer_concentration"] = poly_conc_var
        poly_entry = self.create_numeric_input(
            self.scrollable_frame, row, 5, poly_conc_var
        )
        self.scrollable_frame.columnconfigure(5, weight=1)

        # Surfactant Type dropdown
        surfactant_type_options = {
            "Alkyl Ether Sulfate": 1,
            "No Surfactant": 2,
        }
        surfactant_type_var = tk.StringVar(
            value=list(surfactant_type_options.keys())[0]
        )
        input_vars["surfactant_type"] = surfactant_type_var

        surfactant_menu = tk.OptionMenu(
            self.scrollable_frame, surfactant_type_var, *surfactant_type_options.keys()
        )
        surfactant_menu.config(width=15)
        surfactant_menu.grid(row=row, column=6, padx=5, pady=5, sticky="ew")
        self.scrollable_frame.columnconfigure(6, weight=1)

        # Surfactant Concentration Entry
        surf_conc_var = tk.DoubleVar()
        input_vars["surfactant_concentration"] = surf_conc_var
        surf_entry = self.create_numeric_input(
            self.scrollable_frame, row, 7, surf_conc_var
        )
        self.scrollable_frame.columnconfigure(7, weight=1)

        self.inputs.append(input_vars)
        self.row_count += 1
        self.scrollable_frame.rowconfigure(row, weight=1)

    def create_buttons(self) -> None:
        button_frame = tk.Frame(self.master)
        button_frame.pack(side=tk.BOTTOM, fill=tk.X)

        tk.Button(button_frame, text="Submit All", command=self.submit_all).pack(
            side=tk.LEFT, padx=10, pady=10
        )
        tk.Button(
            button_frame, text="Add Simulation", command=self.add_simulation
        ).pack(side=tk.LEFT, padx=10)
        tk.Button(button_frame, text="Upload CSV", command=self.load_csv).pack(
            side=tk.LEFT, padx=10
        )

        for button in button_frame.winfo_children():
            button.pack(side=tk.LEFT, padx=(10, 0), pady=10)

    def load_csv(self) -> None:
        file_path = filedialog.askopenfilename(
            title="Select a CSV file", filetypes=[("CSV files", "*.csv")]
        )
        if file_path:
            try:
                data = pd.read_csv(file_path)
                messagebox.showinfo("CSV Loaded", data.head().to_string())
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load CSV file: {e}")

    def submit_all(self) -> None:
        model_options = {
            "No Shear Thinning": 1,
            "Sourav Implementation": 2,
            "Shear Thinning": 3,
        }

        permeability_options = {
            "Homogeneous": 1,
            "Heterogeneous": 2,
        }

        polymer_type_options = {
            "Xanthane": 1,
            "Schizophyllan": 2,
            "No Polymer": 3,
        }

        surfactant_type_options = {
            "Alkyl Ether Sulfate": 1,
            "No Surfactant": 2,
        }

        all_inputs = []
        for input_vars in self.inputs:
            input_values = {}
            for key, var in input_vars.items():
                if key == "model_type":
                    input_values[key] = model_options[var.get()]
                elif key == "permeability":
                    input_values[key] = permeability_options[var.get()]
                elif key == "polymer_type":
                    input_values[key] = polymer_type_options[var.get()]
                elif key == "surfactant_type":
                    input_values[key] = surfactant_type_options[var.get()]
                else:
                    input_values[key] = var.get()

            # Validate polymer and surfactant concentrations
            polymer_conc = input_values.get("polymer_concentration")
            surfactant_conc = input_values.get("surfactant_concentration")

            if not self._is_valid_concentration(polymer_conc):
                messagebox.showwarning(
                    "Input Error",
                    "Polymer Concentration must be a numeric value between 0 and 1 inclusive.",
                )
                return

            if not self._is_valid_concentration(surfactant_conc):
                messagebox.showwarning(
                    "Input Error",
                    "Surfactant Concentration must be a numeric value between 0 and 1 inclusive.",
                )
                return

            all_inputs.append(input_values)

        if all(all_inputs):
            print("All inputs received:")
            for sim in all_inputs:
                print(sim)
            self.all_inputs = all_inputs
            self.master.destroy()
        else:
            messagebox.showwarning(
                "Input Error", "Please fill in all fields before submitting."
            )

    def get_input(self) -> list[dict]:
        return self.all_inputs

    def _is_valid_concentration(self, concentration):
        """Helper method to check if the concentration is numeric and between 0 and 1."""
        try:
            concentration_value = float(concentration)
            return 0 <= concentration_value <= 1
        except ValueError:
            return False

    def create_numeric_input(self, parent, row, column, input_var, width=10):
        """Create an input field that only accepts numeric values."""

        # Create the validation command
        validate_cmd = (parent.register(self.validate_numeric_input), "%P")

        entry = tk.Entry(
            parent,
            textvariable=input_var,
            width=width,
            validate="key",
            validatecommand=validate_cmd,
        )
        entry.grid(row=row, column=column, padx=5, pady=5, sticky="ew")
        parent.columnconfigure(column, weight=1)
        return entry

    def validate_numeric_input(self, value):
        """Validate that the input is a number (integer or float)."""
        try:
            # Check if the value can be converted to a float
            float(value)
            return True
        except ValueError:
            return False
